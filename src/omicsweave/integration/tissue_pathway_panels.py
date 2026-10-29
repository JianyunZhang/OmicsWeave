"""配置驱动的组织通路证据复核与绘图数据导出。

复用已有差异分析；只在动物层重算脂质类别。MFA的效应量入选标记
与FDR判断分开保存。原始拟合迭代不作为箱线图的生物学重复。
"""
from __future__ import annotations

import argparse
from datetime import datetime
import hashlib
import json
from pathlib import Path
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import yaml

from omicsweave.core.project import project_dir


CLASS_PATTERNS = {
    "TG": r"^TG\(", "DG": r"^DG\(", "Ceramide": r"^Cer",
    "Acylcarnitine": "Carnitine", "PC/LPC": r"^(?:PC|LPC)\(",
    "PE/LPE": r"^(?:PE|LPE)[\(O ]", "PA/LPA": r"^(?:PA|LPA)[\(\s]",
}


def save_csv(frame, path):
    frame.to_csv(path, index=False, encoding="utf-8-sig")


def descriptive_loo(x, y):
    differences = [np.delete(x, i).mean() - y.mean() for i in range(len(x))]
    differences += [x.mean() - np.delete(y, i).mean() for i in range(len(y))]
    delta = float(x.mean() - y.mean())
    return {"visual_mean_difference": delta, "loo_difference_min": min(differences),
            "loo_difference_max": max(differences),
            "loo_direction_retained": bool(all(np.sign(d) == np.sign(delta) for d in differences))}


def extract_brain_reference(reference, project_root, output, read, check_bh):
    """脑HFD/LFD参照独立导出，绝不追加到治疗合并表。"""
    paths = {key: project_root / value for key, value in reference["paths"].items()}
    out = output / reference["output_subdir"]
    out.mkdir(parents=True, exist_ok=True)
    comparison = f"{reference['case']}_vs_{reference['control']}"
    group_map = {reference["case"]: reference["case_label"], reference["control"]: reference["control_label"]}
    rna_path = paths["transcriptomics"] / "04_de_results/brain" / f"de_{comparison}.csv"
    protein_path = paths["proteinomics"] / "01_dp_results/brain" / f"dp_brain_{comparison}.csv"
    met_path = paths["metabolomics"] / "00_dm_results/brain" / f"dm_brain_{reference['metabolite_comparison']}.csv"
    rna = read(rna_path).set_index("gene")
    protein = read(protein_path).set_index("Protein")
    met = read(met_path).set_index("Index")
    assert rna.index.is_unique and protein.index.is_unique and met.index.is_unique
    check_bh(protein, "pvalue", "padj", "brain_reference_protein_BH")
    check_bh(met, "pvalue", "padj", "brain_reference_metabolite_BH")
    rna_meta = read(paths["transcriptomics"] / "sample_metadata.csv")
    rna_meta = rna_meta[rna_meta.tissue.eq("brain") & rna_meta.condition.isin(group_map)].copy()
    counts = read(paths["transcriptomics"] / "04_de_results/gene_counts_matrix.csv").set_index("gene")
    libraries = counts[rna_meta["sample"].tolist()].sum()
    cpm = counts[rna_meta["sample"].tolist()].div(libraries, axis=1) * 1e6
    protein_meta = read(paths["proteinomics"] / "00_qc/brain/sample_metadata.csv")
    protein_meta = protein_meta[protein_meta.condition.isin(group_map)].copy()
    protein_values = read(paths["proteinomics"] / "01_dp_results/brain/normalized_matrix.csv").set_index("Feature")
    filtered = read(paths["metabolomics"] / "01_data_preprocessing/brain/filtered_data.csv").set_index("Index")
    frames, evidence_rows = [], []
    for order, node in enumerate(reference["nodes"]):
        layer, node_id = node["layer"], node["id"]
        rows = []
        if layer == "RNA":
            ev, source_file = rna.loc[node_id], rna_path
            effect = ev.log2FoldChange_unshrunk
            family = "reference_whole_transcriptome_source"
            for _, sample in rna_meta.iterrows():
                rows.append({"sample_id": sample["sample"], "mouse_id": str(sample.mouse_id),
                             "group": group_map[sample.condition], "cpm": cpm.loc[node_id, sample["sample"]],
                             "value": np.log2(cpm.loc[node_id, sample["sample"]] + 1), "value_scale": "log2_CPM_plus_1"})
        elif layer == "protein":
            assert ";" not in node_id, "共享蛋白组不能自动作为单基因测量"
            ev, source_file = protein.loc[node_id], protein_path
            effect = ev.log2FoldChange
            family = "reference_all_testable_protein_groups"
            for _, sample in protein_meta.iterrows():
                value = protein_values.loc[node_id, sample["sample"]]
                if pd.notna(value):
                    rows.append({"sample_id": sample["sample"], "mouse_id": str(sample.mouse_id),
                                 "group": group_map[sample.condition], "value": value,
                                 "value_scale": "log2_normalized_unimputed_protein_signal"})
        else:
            ev, source_file = met.loc[node_id], met_path
            effect = ev.log2FoldChange
            family = f"reference_all_{len(met)}_metabolite_features"
            for _, sample in rna_meta.iterrows():
                column = f"{sample.mouse_id}-brain"
                rows.append({"sample_id": column, "mouse_id": str(sample.mouse_id),
                             "group": group_map[sample.condition], "value": filtered.loc[node_id, column],
                             "value_scale": "log2_PQN_signal"})
        frame = pd.DataFrame(rows)
        expected = {reference["case_label"]: reference["expected_counts"]["case"], reference["control_label"]: reference["expected_counts"]["control"]}
        assert frame.groupby("group").size().to_dict() == expected
        if layer == "protein":
            assert len(frame[frame.group.eq(reference["case_label"])]) == ev.n_numerator
            assert len(frame[frame.group.eq(reference["control_label"])]) == ev.n_denominator
            delta = frame[frame.group.eq(reference["case_label"])].value.mean() - frame[frame.group.eq(reference["control_label"])].value.mean()
            assert np.isclose(delta, effect, atol=1e-8), "蛋白逐动物值与源对比效应不符"
        for key, value in {"tissue": "brain_reference", "layer": layer, "node_id": node_id, "node_label": node["label"],
                           "node_order": order, "comparison": comparison, "source_pvalue": ev.pvalue, "source_padj": ev.padj,
                           "source_effect": effect, "source_effect_scale": "unshrunk_log2FC" if layer == "RNA" else "log2FC",
                           "fdr_family": family, "source_file": str(source_file.relative_to(project_root)),
                           "data_role": "HFD_LFD_obesity_reference", "is_treatment_comparison": False,
                           "experiment_id": reference["experiment_id"]}.items():
            frame[key] = value
        frames.append(frame)
        evidence_rows.append({"tissue": "brain_reference", "layer": layer, "node_id": node_id, "label": node["label"],
                              "annotation": ev.get("Compounds", node_id), "effect": effect, "effect_scale": frame.source_effect_scale.iloc[0],
                              "pvalue": ev.pvalue, "padj": ev.padj, "n_case": expected[reference["case_label"]],
                              "n_control": expected[reference["control_label"]], "fdr_family": family,
                              "source_file": str(source_file.relative_to(project_root)), "comparison": comparison,
                              "data_role": "HFD_LFD_obesity_reference", "is_treatment_comparison": False})
    data = pd.concat(frames, ignore_index=True)
    evidence = pd.DataFrame(evidence_rows)
    evidence["fdr_significant"] = evidence.padj.lt(.05)
    evidence["display_direction"] = np.where(evidence.padj.lt(.05), np.where(evidence.effect.gt(0), "up", "down"), "NS")
    assert not data.duplicated(["layer", "node_id", "sample_id"]).any()
    save_csv(data, out / "brain_reference_molecular_plot_data.csv")
    save_csv(evidence, out / "brain_reference_pathway_evidence.csv")
    summary = {"experiment_id": reference["experiment_id"], "comparison": comparison, "is_treatment_comparison": False,
               "brain_region": reference["metadata_brain_region"], "n_nodes": len(reference["nodes"]), "n_plot_records": len(data),
               "rna_fdr_total": int(rna.padj.lt(.05).sum()), "protein_fdr_total": int(protein.padj.lt(.05).sum()),
               "metabolite_fdr_total": int(met.padj.lt(.05).sum()), "source_dir": reference["output_subdir"]}
    (out / "brain_reference_scope.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"脑HFD/LFD参照：独立导出{len(reference['nodes'])}个节点、{len(data)}条动物-节点记录。")
    return summary


def extract(config_path: Path, project_root: Path):
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    paths = {key: project_root / value for key, value in config["paths"].items()}
    output = paths["output"]
    output.mkdir(parents=True, exist_ok=True)
    case, control = config["case"], config["control"]
    comparison = f"{case}_vs_{control}"
    case_label, control_label = config["case_label"], config["control_label"]
    expected_molecular = {case_label: config["expected_counts"]["molecular"]["case"], control_label: config["expected_counts"]["molecular"]["control"]}
    expected_mfa = {case_label: config["expected_counts"]["mfa"]["case"], control_label: config["expected_counts"]["mfa"]["control"]}
    input_paths = [config_path, Path(__file__).resolve(), Path(__file__).with_name("plot_tissue_pathway_panels.py").resolve()]
    validations = []

    def read(path):
        input_paths.append(path)
        return pd.read_csv(path, low_memory=False)

    def check_bh(frame, p_column, q_column, name):
        finite = np.isfinite(frame[p_column])
        expected = multipletests(frame.loc[finite, p_column], method="fdr_bh")[1]
        actual = frame.loc[finite, q_column].to_numpy(float)
        error = float(np.max(np.abs(expected - actual)))
        if not np.allclose(expected, actual, atol=1e-10, rtol=1e-8):
            raise ValueError(f"{name}: BH校正来源不一致，最大差={error}")
        validations.append({"check": name, "status": "pass", "n_tests": int(finite.sum()), "max_error": error})

    meta = read(paths["transcriptomics"] / "sample_metadata.csv")
    meta = meta[meta.condition.isin([case, control]) & meta.tissue.isin(config["tissues"])].copy()
    counts = read(paths["transcriptomics"] / "04_de_results/gene_counts_matrix.csv").set_index("gene")
    assert counts.index.is_unique
    libraries = counts[meta["sample"].tolist()].sum()
    cpm = counts[meta["sample"].tolist()].div(libraries, axis=1) * 1e6
    gene_modules = {}
    for module, genes in config["gene_modules"].items():
        for gene in genes:
            gene_modules.setdefault(gene, []).append(module)
    save_csv(pd.DataFrame([{"gene": g, "modules": ";".join(v)} for g, v in gene_modules.items()]), output / "candidate_gene_panel.csv")

    mfa_values = read(paths["mfa"] / "01_df_results/all_tissues/flux_matrix/sample_net_flux.csv")
    mfa_all = []
    # 现有MFA先在同一比较的所有组织×反应上BH，再拆分文件，须重建整个family核查。
    for tissue in sorted(mfa_values.tissue.unique()):
        path = paths["mfa"] / "01_df_results" / tissue / "differential_flux" / f"differential_flux_{comparison.upper()}.csv"
        frame = read(path)
        frame["source_file"] = str(path.relative_to(project_root))
        mfa_all.append(frame)
    mfa_all = pd.concat(mfa_all, ignore_index=True)
    assert mfa_all.flux_id.is_unique
    check_bh(mfa_all, "pvalue_ttest", "fdr_ttest", "MFA_Welch_BH_all_tissue_reactions")
    check_bh(mfa_all, "pvalue_mannwhitney", "fdr_mannwhitney", "MFA_MW_BH_all_tissue_reactions")
    save_csv(mfa_all[["flux_id", "tissue", "reaction", "pvalue_ttest", "fdr_ttest", "pvalue_mannwhitney",
                      "fdr_mannwhitney", "primary_statistic", "evidence_level", "is_significant"]], output / "mfa_fdr_family_audit.csv")
    mid = read(paths["mfa"] / "01_df_results/all_tissues/mid/mid_fit_residuals.csv")
    legacy = read(paths["legacy_lipid_classes"])
    legacy = legacy[(legacy.cohort == config["legacy_cohort"]) & legacy.comparison.eq(case_label)].copy()
    legacy["lipid_class"] = legacy.lipid_class.replace({"PA": "PA/LPA"})
    all_molecular, all_flux, all_evidence, summaries = [], [], [], []
    all_rna_candidates, all_lipid_tests = [], []

    for tissue, settings in config["tissues"].items():
        tissue_out = output / tissue
        tissue_out.mkdir(exist_ok=True)
        molecular_available = settings.get("molecular_available", True)
        if molecular_available:
            rna_path = paths["transcriptomics"] / "04_de_results" / tissue / f"de_{comparison}.csv"
            rna = read(rna_path).set_index("gene")
            assert rna.index.is_unique
            candidates = rna.loc[rna.index.intersection(gene_modules)].copy().reset_index()
            candidates.insert(0, "tissue", tissue)
            candidates["module"] = candidates.gene.map(lambda gene: ";".join(gene_modules[gene]))
            candidates["source_file"] = str(rna_path.relative_to(project_root))
            candidates["fdr_significant"] = candidates.padj.lt(.05)
            save_csv(candidates, tissue_out / f"{tissue}_RNA_candidate_evidence.csv")
            all_rna_candidates.append(candidates)

            met_path = paths["metabolomics"] / "00_dm_results" / tissue / f"dm_{tissue}_{comparison}.csv"
            met = read(met_path)
            assert met.Index.is_unique
            check_bh(met, "pvalue", "padj", f"{tissue}_metabolite_BH_all_features")
            met = met.set_index("Index")
            filtered_path = paths["metabolomics"] / "01_data_preprocessing" / tissue / "filtered_data.csv"
            filtered_all = read(filtered_path)
            assert filtered_all.Index.is_unique
            sample_columns = [c for c in filtered_all if re.match(config["metabolite_sample_regex"], c.lower())]
            case_cols = [c for c in sample_columns if f"-{config['metabolite_case_token']}-" in c]
            control_cols = [c for c in sample_columns if f"-{config['metabolite_control_token']}-" in c]
            assert len(case_cols) == expected_molecular[case_label] and len(control_cols) == expected_molecular[control_label]
            filtered = filtered_all.drop_duplicates(subset=sample_columns)
            class_scores, class_tests, class_members = {}, [], []
            for name, pattern in CLASS_PATTERNS.items():
                subset = filtered[filtered.Compounds.str.contains(pattern, regex=True, na=False)]
                means = subset[sample_columns].mean()
                class_scores[name] = means
                x, y = means[case_cols].to_numpy(float), means[control_cols].to_numpy(float)
                p = stats.ttest_ind(x, y, equal_var=False).pvalue
                se = np.sqrt(x.var(ddof=1) / len(x) + y.var(ddof=1) / len(y))
                v1, v2 = x.var(ddof=1) / len(x), y.var(ddof=1) / len(y)
                dof = (v1 + v2) ** 2 / (v1 ** 2 / (len(x) - 1) + v2 ** 2 / (len(y) - 1))
                interval = stats.t.ppf(.975, dof) * se
                class_tests.append({"tissue": tissue, "node_id": name, "n_features": len(subset),
                                    "n_case": len(x), "n_control": len(y), "delta": x.mean() - y.mean(),
                                    "pvalue": p, "p_mw_exact_sensitivity": stats.mannwhitneyu(x, y, alternative="two-sided", method="exact").pvalue,
                                    "ci95_low": x.mean() - y.mean() - interval, "ci95_high": x.mean() - y.mean() + interval,
                                    **descriptive_loo(x, y)})
                members = subset[["Index", "Compounds", "Class I", "Class II"]].copy()
                members.insert(0, "node_id", name)
                class_members.append(members)
            class_tests = pd.DataFrame(class_tests).set_index("node_id")
            class_tests["padj"] = multipletests(class_tests.pvalue, method="fdr_bh")[1]
            class_tests["fdr_significant"] = class_tests.padj.lt(.05)
            for name, row in class_tests.iterrows():
                old = legacy[legacy.tissue.eq(tissue) & legacy.lipid_class.eq(name)]
                if len(old) == 1:
                    same = np.isclose(old.delta.iloc[0], row.delta, atol=1e-8) and np.isclose(old.qvalue_bh.iloc[0], row.padj, atol=1e-8)
                    validations.append({"check": f"legacy_lipid_{tissue}_{name}", "status": "pass" if same else "updated_from_current_data",
                                        "old_q": old.qvalue_bh.iloc[0], "current_q": row.padj})
            save_csv(class_tests.reset_index(), tissue_out / f"{tissue}_lipid_class_tests.csv")
            save_csv(pd.concat(class_members, ignore_index=True), tissue_out / f"{tissue}_lipid_class_members.csv")
            all_lipid_tests.append(class_tests.reset_index())
            if settings.get("ceramide_definition_sensitivity", False):
                strict = filtered[filtered.Compounds.str.match(r"^Cer\(", na=False)]
                means = strict[sample_columns].mean()
                x, y = means[case_cols].to_numpy(float), means[control_cols].to_numpy(float)
                p_strict = stats.ttest_ind(x, y, equal_var=False).pvalue
                family = class_tests.pvalue.copy()
                family.loc["Ceramide"] = p_strict
                q_strict = multipletests(family, method="fdr_bh")[1][family.index.get_loc("Ceramide")]
                sensitivity = {"definition": "^Cer( excludes CerP, retains dihydro/phytoceramides", "n_features": len(strict),
                               "delta": float(x.mean() - y.mean()), "pvalue": p_strict, "q_BH_seven_classes": q_strict,
                               **descriptive_loo(x, y)}
                save_csv(pd.DataFrame([sensitivity]), tissue_out / f"{tissue}_ceramide_definition_sensitivity.csv")

            evidence = []
            for gene, row in rna.iterrows():
                if gene in gene_modules:
                    evidence.append({"tissue": tissue, "layer": "RNA", "node_id": gene, "label": gene,
                                     "effect": row.log2FoldChange_unshrunk, "effect_scale": "unshrunk_log2FC",
                                     "pvalue": row.pvalue, "padj": row.padj, "n_case": row.n_case, "n_control": row.n_control,
                                     "fdr_family": "whole_transcriptome_source", "source_file": str(rna_path.relative_to(project_root))})
            for node in settings["molecular_nodes"]:
                if node["layer"] == "metabolite":
                    row = met.loc[node["id"]]
                    evidence.append({"tissue": tissue, "layer": "metabolite", "node_id": node["id"], "label": node["label"],
                                     "annotation": row.Compounds, "effect": row.log2FoldChange, "effect_scale": "log2FC_relative_signal",
                                     "pvalue": row.pvalue, "padj": row.padj, "n_case": len(case_cols), "n_control": len(control_cols),
                                     "fdr_family": f"all_{len(met)}_tissue_metabolite_features", "source_file": str(met_path.relative_to(project_root))})
            for name, row in class_tests.iterrows():
                evidence.append({"tissue": tissue, "layer": "lipid_class", "node_id": name, "label": name,
                                 "effect": row.delta, "effect_scale": "difference_mean_log2_class_signal",
                                 "pvalue": row.pvalue, "padj": row.padj, "n_case": row.n_case, "n_control": row.n_control,
                                 "n_features": row.n_features, "fdr_family": "seven_lipid_classes_per_tissue_comparison",
                                 "source_file": str(filtered_path.relative_to(project_root))})
            evidence = pd.DataFrame(evidence)
            evidence["comparison"] = comparison
            evidence["data_role"] = "cohort_B_treatment"
            evidence["is_treatment_comparison"] = True
            evidence["fdr_significant"] = evidence.padj.lt(.05)
            evidence["display_direction"] = np.where(evidence.padj.isna(), "NA", np.where(evidence.padj.ge(.05), "NS", np.where(evidence.effect.gt(0), "up", "down")))
            save_csv(evidence, tissue_out / f"{tissue}_pathway_evidence.csv")
            all_evidence.append(evidence)
            lookup = evidence.set_index(["layer", "node_id"])
            molecular = []
            selected_meta = meta[meta.tissue.eq(tissue)]
            for order, node in enumerate(settings["molecular_nodes"]):
                layer, node_id = node["layer"], node["id"]
                ev = lookup.loc[(layer, node_id)]
                rows = []
                if layer == "RNA":
                    for _, sample in selected_meta.iterrows():
                        row = {"sample_id": sample["sample"], "mouse_id": str(sample.mouse_id),
                               "group": case_label if sample.condition == case else control_label,
                               "value": np.log2(cpm.loc[node_id, sample["sample"]] + 1),
                               "cpm": cpm.loc[node_id, sample["sample"]], "gene_count": counts.loc[node_id, sample["sample"]],
                               "library_gene_count_sum": libraries.loc[sample["sample"]], "value_scale": "log2_CPM_plus_1"}
                        rows.append(row)
                else:
                    values = class_scores[node_id] if layer == "lipid_class" else filtered_all.set_index("Index").loc[node_id, sample_columns]
                    for column in control_cols + case_cols:
                        rows.append({"sample_id": column, "mouse_id": column.split("-")[0],
                                     "group": case_label if column in case_cols else control_label,
                                     "value": float(values[column]), "value_scale": "mean_log2_PQN_signal" if layer == "lipid_class" else "log2_PQN_signal"})
                frame = pd.DataFrame(rows)
                assert frame.groupby("group").size().to_dict() == expected_molecular
                for key, value in {"tissue": tissue, "layer": layer, "node_id": node_id, "node_label": node["label"],
                                   "node_order": order, "comparison": comparison, "source_pvalue": ev.pvalue,
                                   "source_padj": ev.padj, "source_effect": ev.effect, "source_effect_scale": ev.effect_scale,
                                   "fdr_family": ev.fdr_family, "source_file": ev.source_file,
                                   "data_role": "cohort_B_treatment", "is_treatment_comparison": True,
                                   "experiment_id": config["experiment_id"]}.items():
                    frame[key] = value
                molecular.append(frame)
            molecular = pd.concat(molecular, ignore_index=True)
            assert not molecular.duplicated(["layer", "node_id", "sample_id"]).any()
            save_csv(molecular, tissue_out / f"{tissue}_molecular_plot_data.csv")
            all_molecular.append(molecular)
            loo_rows = []
            for (layer, node_id), frame in molecular.groupby(["layer", "node_id"]):
                x = frame[frame.group.eq(case_label)].value.to_numpy(float)
                y = frame[frame.group.eq(control_label)].value.to_numpy(float)
                loo_rows.append({"tissue": tissue, "layer": layer, "node_id": node_id, **descriptive_loo(x, y)})
            save_csv(pd.DataFrame(loo_rows), tissue_out / f"{tissue}_molecular_LOO_direction.csv")

        else:
            # 不以独立HFD/LFD队列替代当前治疗比较的缺失分子层。
            rna = pd.DataFrame(columns=["padj"])
            candidates = pd.DataFrame(columns=["fdr_significant"])
            met = pd.DataFrame(columns=["padj"])
            class_tests = pd.DataFrame(columns=["padj"])
            molecular = pd.DataFrame()
            availability = pd.DataFrame([
                {"tissue": tissue, "layer": layer, "comparison": comparison,
                 "status": "not_available_for_treatment",
                 "task_status": "not_applicable_without_new_source_data",
                 "reason": settings["missing_molecular_reason"],
                 "separate_reference": config.get("brain_reference", {}).get("output_subdir", "")}
                for layer in ["RNA", "protein", "metabolite"]
            ])
            save_csv(availability, tissue_out / f"{tissue}_molecular_data_availability.csv")

        native = settings["mfa_tissue"]
        native_stats = mfa_all[mfa_all.tissue.eq(native)].copy()
        save_csv(native_stats, tissue_out / f"{tissue}_MFA_all_reaction_evidence.csv")
        if not molecular_available:
            selected_mfa = native_stats[native_stats.reaction.isin(settings["mfa_nodes"])].copy()
            evidence = pd.DataFrame({
                "tissue": tissue, "layer": "MFA", "node_id": selected_mfa.reaction,
                "label": selected_mfa.reaction, "effect": selected_mfa["diff"],
                "effect_scale": "difference_of_animal_model_net_flux_medians",
                "pvalue": selected_mfa.pvalue_mannwhitney, "padj": selected_mfa.fdr_mannwhitney,
                "fdr_ttest": selected_mfa.fdr_ttest, "n_case": expected_mfa[case_label],
                "n_control": expected_mfa[control_label], "fdr_family": f"all_{len(mfa_all)}_tissue_reactions_in_comparison",
                "source_file": selected_mfa.source_file, "comparison": comparison,
                "evidence_level": selected_mfa.evidence_level, "identifiability": selected_mfa.identifiability,
                "data_role": "cohort_B_treatment", "is_treatment_comparison": True,
            })
            evidence["fdr_significant"] = evidence.padj.lt(.05)
            evidence["display_direction"] = np.where(evidence.padj.lt(.05), "FDR_model", "NS")
            save_csv(evidence, tissue_out / f"{tissue}_pathway_evidence.csv")
            all_evidence.append(evidence)
        mfa_subset = mfa_values[mfa_values.tissue.eq(native) & mfa_values.condition.isin([case.lower(), control.lower()])
                                & mfa_values.reaction.isin(settings["mfa_nodes"])].copy()
        assert mfa_subset.reaction.nunique() == len(settings["mfa_nodes"])
        mfa_subset["native_tissue"] = native
        mfa_subset["tissue"] = tissue
        mfa_subset["group"] = mfa_subset.condition.map({case.lower(): case_label, control.lower(): control_label})
        mfa_subset["node_id"] = mfa_subset.reaction
        mfa_subset["node_label"] = mfa_subset.reaction
        mfa_subset["node_order"] = mfa_subset.reaction.map({r: i for i, r in enumerate(settings["mfa_nodes"])})
        mfa_subset["layer"] = "MFA"
        mfa_subset["value"] = mfa_subset.net_flux
        mfa_subset["value_scale"] = "animal_median_filtered_model_net_flux"
        mfa_subset["comparison"] = comparison
        mfa_subset["data_role"] = "cohort_B_treatment"
        mfa_subset["is_treatment_comparison"] = True
        mfa_subset["experiment_id"] = config["experiment_id"]
        mfa_subset["near_one_flag"] = np.isclose(mfa_subset.value, 1, atol=1e-3, rtol=0)
        mfa_subset["near_zero_flag"] = np.isclose(mfa_subset.value, 0, atol=1e-6, rtol=0)
        mfa_subset["source_values_file"] = str((paths["mfa"] / "01_df_results/all_tissues/flux_matrix/sample_net_flux.csv").relative_to(project_root))
        native_lookup = native_stats.set_index("reaction")
        for column in ["pvalue_ttest", "fdr_ttest", "pvalue_mannwhitney", "fdr_mannwhitney", "evidence_level",
                       "primary_statistic", "identifiability", "description", "is_significant", "source_file"]:
            mfa_subset[column] = mfa_subset.reaction.map(native_lookup[column])
        mfa_subset.rename(columns={"is_significant": "source_effect_size_selected"}, inplace=True)
        mfa_subset["fdr_family"] = f"all_{len(mfa_all)}_tissue_reactions_in_comparison"
        for reaction, frame in mfa_subset.groupby("reaction"):
            assert frame.groupby("group").size().to_dict() == expected_mfa
            ev = native_lookup.loc[reaction]
            for group, condition in [(case_label, case.lower()), (control_label, control.lower())]:
                assert np.isclose(frame[frame.group.eq(group)].value.mean(), ev[f"mean_{condition}"], rtol=1e-8, atol=1e-7)
        save_csv(mfa_subset.sort_values(["node_order", "group", "sample_id"]), tissue_out / f"{tissue}_MFA_plot_data.csv")
        save_csv(mfa_subset.pivot(index=["group", "sample_id"], columns="node_id", values="value").reset_index(), tissue_out / f"{tissue}_MFA_plot_data_wide.csv")
        all_flux.append(mfa_subset)
        tissue_mid = mid[mid.condition.isin([case.lower(), control.lower()]) & mid.metabolite.str.startswith(native + "--")]
        save_csv(tissue_mid, tissue_out / f"{tissue}_MID_fit_residuals.csv")
        summary = {"tissue": tissue, "display": settings["display"], "native_mfa_tissue": native,
                   "molecular_available": molecular_available,
                   "rna_total_tested_rows": len(rna) if molecular_available else None,
                   "rna_fdr_total": int(rna.padj.lt(.05).sum()) if molecular_available else None,
                   "rna_candidate_rows": len(candidates) if molecular_available else None,
                   "rna_candidate_fdr": int(candidates.fdr_significant.sum()) if molecular_available else None,
                   "metabolite_tested_features": len(met) if molecular_available else None,
                   "metabolite_fdr": int(met.padj.lt(.05).sum()) if molecular_available else None,
                   "lipid_class_fdr": int(class_tests.padj.lt(.05).sum()) if molecular_available else None,
                   "mfa_total_reactions": len(native_stats), "mfa_source_effect_size_selected": int(native_stats.is_significant.sum()),
                   "mfa_mw_fdr_significant": int(native_stats.fdr_mannwhitney.lt(.05).sum()),
                   "mfa_welch_fdr_significant": int(native_stats.fdr_ttest.lt(.05).sum()),
                   "mfa_min_q_mw": float(native_stats.fdr_mannwhitney.min()), "mfa_min_q_welch": float(native_stats.fdr_ttest.min()),
                   "molecular_plot_rows": len(molecular), "mfa_plot_rows": len(mfa_subset)}
        summaries.append(summary)
        molecular_note = f"RNA候选FDR={summary['rna_candidate_fdr']}" if molecular_available else "同期治疗分子数据不适用"
        output_note = "分子/MFA逐动物CSV" if molecular_available else "MFA逐动物CSV和分子层数据覆盖说明"
        print(f"{tissue}: {molecular_note}，MFA效应量标记={summary['mfa_source_effect_size_selected']}，MFA FDR通过={summary['mfa_mw_fdr_significant']}；已导出{output_note}。")

    reference_summary = None
    if config.get("brain_reference"):
        reference_summary = extract_brain_reference(config["brain_reference"], project_root, output, read, check_bh)
    save_csv(pd.concat(all_molecular, ignore_index=True), output / "all_tissues_molecular_plot_data.csv")
    save_csv(pd.concat(all_flux, ignore_index=True), output / "all_tissues_MFA_plot_data.csv")
    save_csv(pd.concat(all_evidence, ignore_index=True), output / "all_tissues_pathway_evidence.csv")
    save_csv(pd.concat(all_rna_candidates, ignore_index=True), output / "all_tissues_RNA_candidate_evidence.csv")
    save_csv(pd.concat(all_lipid_tests, ignore_index=True), output / "all_tissues_lipid_class_tests.csv")
    save_csv(pd.DataFrame(summaries), output / "evidence_summary.csv")
    save_csv(pd.DataFrame(validations), output / "source_validation.csv")
    manifest = {"analysis_date": config["analysis_date"], "generated_at_local": datetime.now().isoformat(),
                "experiment_id": config["experiment_id"], "comparison": comparison, "config": str(config_path.relative_to(project_root)),
                "rna_metabolite_n": expected_molecular, "mfa_n": expected_mfa,
                "mfa_bh_family_size": len(mfa_all), "mfa_bh_family": "all tissue-reaction pairs in the same contrast, before tissue splitting",
                "boxplot_unit": "one biological animal", "mfa_value_definition": "median of existing per-animal filtered fitted net fluxes",
                "RNA_plot_scale": "log2(CPM+1) from gene count matrix for visualization; p-values remain source DESeq2 model results",
                "boxplot_definition": "animal values: median, Q1/Q3, Tukey 1.5 IQR whiskers; all animals overlaid",
                "source_mfa_effect_size_flag_is_not_fdr": True, "full_pipeline_refit": False,
                "candidate_panel_is_exploratory": True,
                "available_treatment_layers": {tissue: (["RNA", "metabolite", "lipid_class", "MFA"] if settings.get("molecular_available", True) else ["MFA"]) for tissue, settings in config["tissues"].items()},
                "brain_reference": reference_summary,
                "input_files": [{"path": str(p.relative_to(project_root)), "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
                                for p in dict.fromkeys(input_paths)]}
    (output / "source_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    (output / "resolved_config.json").write_text(json.dumps(config, ensure_ascii=False, indent=2), encoding="utf-8")
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--project-root", type=Path, default=None,
                        help="项目目录；默认按 omicsweave.core.project.project_dir() 查找")
    args = parser.parse_args()
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    root = (args.project_root or project_dir()).resolve()
    config_path = args.config.resolve()
    output = extract(config_path, root)
    print("证据与CSV导出完成：", output)


if __name__ == "__main__":
    main()
