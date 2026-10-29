"""核查正式通路图/箱线图/逐动物CSV，并可打包发布本地结果。"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import sys
import xml.etree.ElementTree as ET
import zipfile

import numpy as np
import pandas as pd

from omicsweave.common.project import project_dir


def check_plot(folder, prefix, mode, expected, role):
    data = pd.read_csv(folder / f"{prefix}_{mode}_plot_data.csv")
    box = pd.read_csv(folder / f"{prefix}_{mode}_box_statistics.csv")
    keys = ["layer", "node_id"]
    assert not data.duplicated(keys + ["sample_id"]).any()
    assert np.isfinite(data.value).all()
    assert data.is_treatment_comparison.eq(role).all()
    for (layer, node_id), group in data.groupby(keys):
        assert group.groupby("group").size().to_dict() == expected
        for label, sample in group.groupby("group"):
            rows = box[box.layer.eq(layer) & box.node_id.eq(node_id) & box.group.eq(label)]
            assert len(rows) == 1
            stats = rows.iloc[0]
            assert stats.n_animals == len(sample)
            assert np.allclose(np.quantile(sample.value, [.25, .5, .75]), [stats.q1, stats.med, stats.q3], atol=1e-10, rtol=1e-11)
            assert stats.axis_ymin <= sample.value.min() and stats.axis_ymax >= sample.value.max()
    for extension in ["svg", "pdf", "png"]:
        assert (folder / f"{prefix}_{mode}_key_nodes_boxplots.{extension}").is_file()
    return {"prefix": prefix, "mode": mode, "n_nodes": len(data[keys].drop_duplicates()),
            "n_animal_node_records": len(data), "n_box_groups": len(box), "is_treatment_comparison": role, "status": "pass"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence-dir", type=Path, required=True)
    parser.add_argument("--project-root", type=Path, default=None,
                        help="项目目录；默认按 omicsweave.common.project.project_dir() 查找")
    parser.add_argument("--package", action="store_true")
    args = parser.parse_args()
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    root, out = (args.project_root or project_dir()).resolve(), args.evidence_dir.resolve()
    config = json.loads((out / "resolved_config.json").read_text(encoding="utf-8"))
    manifest = json.loads((out / "source_manifest.json").read_text(encoding="utf-8"))
    for source in manifest["input_files"]:
        assert hashlib.sha256((root / source["path"]).read_bytes()).hexdigest() == source["sha256"], source["path"]
    plot_checks, requirements = [], []
    for tissue, settings in config["tissues"].items():
        folder = out / tissue
        expected_mol = {config["case_label"]: config["expected_counts"]["molecular"]["case"], config["control_label"]: config["expected_counts"]["molecular"]["control"]}
        expected_mfa = {config["case_label"]: config["expected_counts"]["mfa"]["case"], config["control_label"]: config["expected_counts"]["mfa"]["control"]}
        for component, name in [("evidence_review", f"{tissue}_pathway_evidence.csv"), ("pathway", f"panelB_{tissue}_pathway.svg"),
                                ("MFA_boxplot", f"{tissue}_MFA_key_nodes_boxplots.svg"), ("MFA_fit_QA", f"{tissue}_MFA_fit_QA.svg"),
                                ("MFA_animal_CSV", f"{tissue}_MFA_plot_data.csv"), ("MFA_box_statistics", f"{tissue}_MFA_box_statistics.csv")]:
            assert (folder / name).is_file()
            requirements.append({"tissue": tissue, "component": component, "status": "complete", "formal_file": str((folder / name).relative_to(out)), "role": "treatment"})
        plot_checks.append(check_plot(folder, tissue, "MFA", expected_mfa, True))
        if settings.get("molecular_available", True):
            plot_checks.append(check_plot(folder, tissue, "molecular", expected_mol, True))
            for component, name in [("molecular_boxplot", f"{tissue}_molecular_key_nodes_boxplots.svg"), ("molecular_animal_CSV", f"{tissue}_molecular_plot_data.csv"), ("molecular_box_statistics", f"{tissue}_molecular_box_statistics.csv")]:
                requirements.append({"tissue": tissue, "component": component, "status": "complete", "formal_file": str((folder / name).relative_to(out)), "role": "treatment"})
        else:
            assert not (folder / f"{tissue}_molecular_plot_data.csv").exists()
            assert not (folder / f"{tissue}_molecular_key_nodes_boxplots.svg").exists()
            coverage = pd.read_csv(folder / f"{tissue}_molecular_data_availability.csv")
            assert coverage.status.eq("not_available_for_treatment").all()
            requirements.append({"tissue": tissue, "component": "treatment_molecular_boxplot_and_CSV", "status": "not_applicable_no_treatment_source", "formal_file": str((folder / f"{tissue}_molecular_data_availability.csv").relative_to(out)), "role": "treatment"})
        annotations = pd.read_csv(folder / f"{tissue}_pathway_annotations.csv")
        for _, row in annotations.iterrows():
            mode = "MFA" if row.layer == "MFA" else "molecular"
            points = pd.read_csv(folder / f"{tissue}_{mode}_plot_data.csv")
            selected = points[points.layer.eq(row.layer) & points.node_id.eq(row.node_id)]
            assert len(selected) == sum(expected_mfa.values() if mode == "MFA" else expected_mol.values())
            if row.display_direction in ["↑", "↓"]:
                assert row.padj < .05 and ((row.effect > 0) == (row.display_direction == "↑"))
        if tissue == "brain":
            assert annotations.layer.eq("MFA").all()
            assert not annotations.display_direction.isin(["↑", "↓"]).any()

    primary_mol = pd.read_csv(out / "all_tissues_molecular_plot_data.csv")
    primary_mfa = pd.read_csv(out / "all_tissues_MFA_plot_data.csv")
    assert primary_mol.is_treatment_comparison.eq(True).all() and primary_mfa.is_treatment_comparison.eq(True).all()
    assert primary_mol.experiment_id.eq(config["experiment_id"]).all()
    assert primary_mfa.experiment_id.eq(config["experiment_id"]).all()
    assert "brain" not in primary_mol.tissue.unique()
    if config.get("brain_reference"):
        ref = config["brain_reference"]
        folder = out / ref["output_subdir"]
        expected = {ref["case_label"]: ref["expected_counts"]["case"], ref["control_label"]: ref["expected_counts"]["control"]}
        plot_checks.append(check_plot(folder, "brain_reference", "molecular", expected, False))
        data = pd.read_csv(folder / "brain_reference_molecular_plot_data.csv")
        assert data.experiment_id.eq(ref["experiment_id"]).all()
        assert data.data_role.eq("HFD_LFD_obesity_reference").all()
        assert not data.group.isin([config["case_label"], config["control_label"]]).any()
        for name in ["brain_reference_pathway.svg", "brain_reference_molecular_key_nodes_boxplots.svg", "brain_reference_molecular_plot_data.csv", "brain_reference_molecular_box_statistics.csv", "brain_reference_pathway_evidence.csv"]:
            assert (folder / name).is_file()
            requirements.append({"tissue": "brain_reference", "component": name, "status": "complete", "formal_file": str((folder / name).relative_to(out)), "role": "HFD_LFD_reference_not_treatment"})

    svg_files = list(out.rglob("*.svg"))
    ns = {"s": "http://www.w3.org/2000/svg"}
    for file in svg_files:
        node = ET.parse(file).getroot()
        assert len(node.findall(".//s:text", ns)) > 0 and not node.findall(".//s:image", ns)
        if "grayscale" in file.name:
            colors = re.findall(r"(?:fill|stroke):\s*#([0-9a-fA-F]{6})", file.read_text(encoding="utf-8"))
            assert all(c[:2] == c[2:4] == c[4:] for c in colors)
    checks = pd.read_csv(out / "source_validation.csv")
    assert checks.status.eq("pass").all()
    consistency = []
    audit = out / "completion_audit"
    audit.mkdir(exist_ok=True)
    for mode in ["molecular", "MFA"]:
        old_path = audit / f"pre_extension_four_tissue_{mode}_data.csv"
        if not old_path.exists():
            continue
        old = pd.read_csv(old_path)
        new = pd.read_csv(out / f"all_tissues_{mode}_plot_data.csv")
        new = new[new.tissue.isin(["heart", "liver", "sq_adipose", "gast_muscle"])]
        keys = ["tissue", "layer", "node_id", "group", "sample_id"]
        values = ["value", "source_padj", "source_effect"] if mode == "molecular" else ["value", "fdr_mannwhitney", "fdr_ttest", "fit_p25", "fit_p75", "n_fit_points"]
        merged = old[keys + values].merge(new[keys + values], on=keys, how="outer", suffixes=("_before", "_after"), indicator=True, validate="one_to_one")
        assert merged["_merge"].eq("both").all()
        for column in values:
            a, b = merged[column + "_before"].to_numpy(float), merged[column + "_after"].to_numpy(float)
            assert np.allclose(a, b, atol=1e-10, rtol=1e-10, equal_nan=True)
            consistency.append({"mode": mode, "field": column, "n_records": len(merged), "max_abs_change": float(np.nanmax(abs(a - b))), "status": "pass"})
    pd.DataFrame(consistency).to_csv(audit / "original_four_tissue_consistency.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(requirements).to_csv(audit / "requirements_after.csv", index=False, encoding="utf-8-sig")
    report = {"status": "pass", "n_treatment_tissues": len(config["tissues"]),
              "treatment_molecular_records": len(primary_mol), "treatment_MFA_records": len(primary_mfa),
              "reference_records": sum(r["n_animal_node_records"] for r in plot_checks if not r["is_treatment_comparison"]),
              "source_hashes_verified": len(manifest["input_files"]), "source_statistical_checks": len(checks),
              "editable_vector_svg_files": len(svg_files), "plot_checks": plot_checks,
              "original_four_tissue_numeric_equivalence": bool(consistency), "cohort_separation_verified": True,
              "brain_treatment_molecular_status": "not_applicable_no_source_data", "all_available_data_deliverables_complete": True,
              "validation_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (out / "final_validation.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    if args.package:
        archive = out / "panelB_six_tissues_figures_and_csv.zip"
        assert archive.resolve().parent == out
        files = [p for p in out.rglob("*") if p.is_file() and p.suffix.lower() in [".svg", ".pdf", ".png", ".csv", ".json", ".md"]
                 and not p.name.startswith("pre_extension_")]
        with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as z:
            for file in files:
                z.write(file, file.relative_to(out))
            for file in [Path(__file__).with_name("tissue_pathway_panels.py"), Path(__file__).with_name("plot_tissue_pathway_panels.py"), Path(__file__), root / manifest["config"]]:
                z.write(file, Path("reproduce") / file.name)
        with zipfile.ZipFile(archive) as z:
            assert z.testzip() is None
        print(f"已打包：{archive.name}，{archive.stat().st_size / 1024**2:.2f} MB")
    print(json.dumps({key: value for key, value in report.items() if key != "plot_checks"}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
