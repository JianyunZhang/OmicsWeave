#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reusable Metabolic Flux Analysis (MFA) Module

Analyses flux_raw_data.xlsx and mid_raw_data.xlsx produced by isotope-tracing
MFA fitting pipelines (e.g. INCA, Metran, 13C-MFA).

Designed to be project-agnostic. All inputs are provided via constructor / config.

分析计算委托给 df_analysis 模块，CSV 输出到 01_df_results/ 子目录。
绘图委托给 plotting 子包。

Typical usage:
    analyser = MFAAnalyser(output_dir="path/to/output")
    analyser.load_flux_data("path/to/flux_raw_data.xlsx")
    analyser.load_mid_data("path/to/mid_raw_data.xlsx")
    analyser.run_all(comparisons=[("hfd", "lfd"), ("hfd", "hsd"), ("hsd", "lfd")])
"""

import pandas as pd
import numpy as np
from pathlib import Path
from itertools import combinations
import warnings

warnings.filterwarnings("ignore")

TISSUE_COLORS = {
    "serum": "#FF4500", # Red
    "heart": "#A83C09", # Dark Mandarin
    "liver": "#8B008B", # Purple
    "kidney": "#663399", # RebeccaPurple
    "brain": "#FFDAB9", # PeachPuff
    "muscle": "#FA8072",
    "white_adipose": "#FFEFD5",

}

CONDITION_COLORS = {
    "hfd": "#FFA500",  # Orange
    "lfd": "#008000",  # Green (Office Green)
    "hsd": "#BA55D3",  # Purple
}

PATHWAY_GROUPS = {
    "Glycolysis": [
        "HEX_c", "PGI_c", "PFK_c", "FBA_c", "TPI_c",
        "GAPD_c", "PGK_c", "PGM_c", "ENO_c", "PYK_c",
    ],
    "TCA_Cycle": [
        "CS_m", "ACONT_m", "ICDH_m", "AKGD_m", "SUCOAS_m",
        "SUCD_m", "FUMH_m", "MDH_m",
    ],
    "Pentose_Phosphate": [
        "G6PDH2R_PGL_GND_c", "RPI_c", "RPE_c",
        "TKT1_c", "TKT2_c", "TALA_c",
    ],
    "Pyruvate_Metabolism": [
        "PDH_m", "LDH_c", "PC_m", "PEPCK_c", "ME2_c",
    ],
    "Amino_Acid": [
        "GLUD_m", "GLND_m", "GLNS_c", "ASPTA_m", "ASPTA_c", "GPT_c",
    ],
    "Lipid_Synthesis": [
        "ACITL_c", "LIPID_c", "ACSS2_m", "ACOT_c", "FA_usage",
    ],
    "Gluconeogenesis": [
        "G6Pase_c", "FBPase_c", "PEPCK_c",
    ],
    "Ketogenesis": [
        "ACAC_production_m", "BHBDH_m", "BHB_trans", "BHB_output",
    ],
    "Glycerol_Metabolism": [
        "GK_c", "G3PP_c", "GPS_c",
    ],
    "Transport": [
        "PYR_trans", "ASPGLU_m", "AKGMAL_m", "CIT_trans",
        "GLN_trans", "GLU_trans",
    ],
}


def _parse_flux_column(col_name: str):
    """Parse 'tissue--reaction' or 'MIX:tissue--metabolite_idx' into components."""
    if col_name.startswith("MIX:"):
        rest = col_name[4:]
        parts = rest.split("--", 1)
        if len(parts) == 2:
            return "MIX", parts[0], parts[1]
        return "MIX", "unknown", rest
    parts = col_name.split("--", 1)
    if len(parts) == 2:
        return "flux", parts[0], parts[1]
    return "flux", "unknown", col_name


class MFAAnalyser:
    """
    MFA 分析引擎：数据加载 + 流程编排。

    分析计算委托给 df_analysis 模块（输出 CSV 到 01_df_results/），
    可视化委托给 plotting 子包（输出 PNG 到 02_result_plot/）。

    输出目录结构（组织为中心）:
      01_df_results/
        all_tissues/
          pca/  qc/  overview/  mid/   跨组织分析
        <tissue>/differential_flux/      差异通量
        <tissue>/pathway_flux/           通路通量
        <tissue>/heatmap/                热图数据
      02_result_plot/
        <tissue>/1_volcano/   ...        组织特异性图表
        all_tissues/0_qc/  ...           跨组织图表
    """

    # 跨组织分析类别（始终放在 01_df_results/all_tissues/ 下）
    _CROSS_TISSUE_DF_CATEGORIES = {"pca", "qc", "overview", "mid"}
    # 组织特异性分析类别（放在 <tissue>/ 子目录下）
    _TISSUE_DF_CATEGORIES = {"differential_flux", "pathway_flux", "heatmap"}

    # 跨组织图表类别（始终放在 02_result_plot/all_tissues/ 下）
    _CROSS_TISSUE_PLOT_CATEGORIES = {
        "0_qc",
        "3_pca_clustering",
        "6_tissue_overview",
        "7_mid_analysis",
    }
    # 组织特异性图表类别（放在 <tissue>/ 子目录下）
    _TISSUE_PLOT_CATEGORIES = {"1_volcano", "2_heatmap", "4_pathway_flux"}

    def __init__(self, output_dir: str, condition_colors: dict = None,
                 tissue_colors: dict = None, pathway_groups: dict = None,
                 analysis_config: dict = None, run_context: dict = None):
        self.output_dir = Path(output_dir)
        self.output_root = self._infer_output_root(self.output_dir)

        # 跨组织分析目录（始终存在）
        for sub in self._CROSS_TISSUE_DF_CATEGORIES:
            (self.output_dir / "01_df_results" / "all_tissues" / sub).mkdir(parents=True, exist_ok=True)
        for sub in self._CROSS_TISSUE_PLOT_CATEGORIES:
            (self.output_dir / "02_result_plot" / "all_tissues" / sub).mkdir(parents=True, exist_ok=True)

        # 多同位素联合分析目录骨架（输出到 mfa_results/integration/ 下）
        (self.output_root / "integration" / "tracer_harmonized").mkdir(parents=True, exist_ok=True)
        (self.output_root / "integration" / "cross_tracer_scoring").mkdir(parents=True, exist_ok=True)

        # 组织特异性目录在 tissues 已知后延迟创建（见 _ensure_output_dirs）
        self._output_dirs_initialized = False

        # 创建解读报告目录
        (self.output_dir / "03_interpretation_report").mkdir(parents=True, exist_ok=True)

        # 统一小写：load_flux_data 对 condition 做 .lower()，
        # 因此颜色映射的 key 也必须小写才能匹配
        self.condition_colors = {k.lower(): v for k, v in (condition_colors or CONDITION_COLORS).items()}
        self.tissue_colors = {k.lower(): v for k, v in (tissue_colors or TISSUE_COLORS).items()}
        self.pathway_groups = pathway_groups or PATHWAY_GROUPS
        self.analysis_config = analysis_config or {}
        # 批次上下文：血糖记录、灌注 QC 阈值、示踪剂与检测元数据。
        # 供稳态 QC 与运行溯源使用；缺省时对应步骤自动跳过。
        self.run_context = run_context or {}

        self.flux_df = None          # fluxes x samples (best-fit)
        self.flux_all_fits = None    # all fits: dict {sample_id: DataFrame}
        self.flux_info = None        # sheet_information DataFrame
        self.sample_conditions = {}  # sample_id -> condition

        self.mid_df = None           # MIDs best-fit: multi-index (metabolite, mid_index) x samples
        self.mid_measured_df = None  # 实测 MID（experiments 列），同索引同列
        self.mid_info = None

        self.tissues = []
        self.conditions = []

    @staticmethod
    def _infer_output_root(output_dir: Path) -> Path:
        """从 runs/<protocol>/<batch> 结构中回溯 mfa_results 根目录。"""
        parts = list(output_dir.parts)
        for i, part in enumerate(parts):
            if part.lower() == "runs":
                if i == 0:
                    return output_dir
                return Path(*parts[:i])
        return output_dir

    def _ensure_output_dirs(self):
        """为已知的组织列表创建组织特异性输出目录。幂等操作。"""
        if self._output_dirs_initialized:
            return
        for tissue in self.tissues:
            for cat in self._TISSUE_DF_CATEGORIES:
                (self.output_dir / "01_df_results" / tissue / cat).mkdir(parents=True, exist_ok=True)
            for cat in self._TISSUE_PLOT_CATEGORIES:
                (self.output_dir / "02_result_plot" / tissue / cat).mkdir(parents=True, exist_ok=True)
        self._output_dirs_initialized = True

    # ── Data Loading ─────────────────────────────────────────────────────

    def load_flux_data(self, file_path: str, sheet_info_name="sheet_information"):
        """Load flux_raw_data.xlsx: extract best-fit row per sample."""
        fp = Path(file_path)
        if not fp.exists():
            print(f"[WARN] Flux file not found: {fp}")
            return

        # 只取工作表名并立即关闭文件：pd.ExcelFile 不关闭会一直占用文件句柄（Windows 上无法移动或删除）
        with pd.ExcelFile(fp) as xl:
            sheet_names = set(xl.sheet_names)
        if sheet_info_name not in sheet_names:
            print(f"[WARN] Sheet '{sheet_info_name}' not in {fp.name}")
            return

        info = pd.read_excel(fp, sheet_name=sheet_info_name)
        self.flux_info = info.copy()

        skip_cols = {"Index", "index", "loss", "Loss"}
        samples, vectors, all_fits = [], [], {}
        cols_ref = None

        for _, row in info.iterrows():
            sn = row["sheet_name"]
            if sn not in sheet_names:
                continue
            sid = str(row.get("complete_result_name", sn))
            cond = str(row.get("condition", "unknown")).lower()

            df = pd.read_excel(fp, sheet_name=sn)
            if df.empty:
                continue

            # MIX:* 是样本特异的模拟 MID 观测量，不是真实通量。不同样本可测
            # 代谢物不同，因此这些列不能进入跨样本固定的通量列模板。
            fcols = [c for c in df.columns
                     if c not in skip_cols
                     and not str(c).startswith("MIX:")
                     and pd.api.types.is_numeric_dtype(df[c])]
            if not fcols:
                continue

            if cols_ref is None:
                cols_ref = list(fcols)

            all_fits[sid] = df

            if "loss" in df.columns or "Loss" in df.columns:
                loss_col = "loss" if "loss" in df.columns else "Loss"
                ls = pd.to_numeric(df[loss_col], errors="coerce")
                # 找到该样本 150 次拟合中 Loss 最小（最优）的那一行
                best_idx = ls.idxmin() if ls.notna().any() else 0
                best = df.loc[best_idx]
            else:
                best = df[fcols].mean()

            # 最终 flux_df 只存储每个样本的这一个"代表值"
            vec = pd.to_numeric(best[cols_ref], errors="coerce").values.astype(float)
            samples.append(sid)
            vectors.append(vec)
            self.sample_conditions[sid] = cond

        if not samples:
            print("[WARN] No flux samples parsed")
            return

        self.flux_df = pd.DataFrame(
            np.vstack(vectors), index=samples, columns=cols_ref
        )
        self.flux_all_fits = all_fits

        self._extract_tissues()
        self.conditions = sorted(set(self.sample_conditions.values()))

        print(f"[OK] Flux data: {len(samples)} samples x {len(cols_ref)} fluxes")
        print(f"     Tissues: {', '.join(self.tissues)}")
        print(f"     Conditions: {', '.join(self.conditions)}")

    def load_mid_data(self, file_path: str, sheet_info_name="sheet_information"):
        """Load mid_raw_data.xlsx: extract best-fit MID per sample."""
        fp = Path(file_path)
        if not fp.exists():
            print(f"[WARN] MID file not found: {fp}")
            return

        # 只取工作表名并立即关闭文件：pd.ExcelFile 不关闭会一直占用文件句柄（Windows 上无法移动或删除）
        with pd.ExcelFile(fp) as xl:
            sheet_names = set(xl.sheet_names)
        if sheet_info_name not in sheet_names:
            print(f"[WARN] Sheet '{sheet_info_name}' not in {fp.name}")
            return

        info = pd.read_excel(fp, sheet_name=sheet_info_name)
        self.mid_info = info.copy()

        all_samples = {}
        for _, row in info.iterrows():
            sn = row["sheet_name"]
            if sn not in sheet_names:
                continue
            sid = str(row.get("complete_result_name", sn))
            cond = str(row.get("condition", "unknown")).lower()

            df = pd.read_excel(fp, sheet_name=sn)
            if df.empty:
                continue

            loss_row = df[df["metabolite_name"] == "loss"]
            mid_rows = df[df["metabolite_name"] != "loss"].copy()

            iter_cols = [c for c in df.columns
                         if c not in ("metabolite_name", "mid_index", "experiments")
                         and pd.api.types.is_numeric_dtype(df[c])]

            if loss_row.empty or not iter_cols:
                best_col = "experiments"
            else:
                losses = pd.to_numeric(loss_row[iter_cols].iloc[0], errors="coerce")
                best_col = losses.idxmin() if losses.notna().any() else iter_cols[0]

            mid_rows = mid_rows.copy()
            mid_rows["best_fit"] = pd.to_numeric(mid_rows[best_col], errors="coerce")
            # experiments 列是**实测** MID。此前它被读入后直接丢弃，
            # 于是"模型拟合得好不好"这个最基本的问题无从回答。
            if "experiments" in mid_rows.columns:
                mid_rows["measured"] = pd.to_numeric(
                    mid_rows["experiments"], errors="coerce"
                )
            else:
                mid_rows["measured"] = np.nan
            mid_rows["sample_id"] = sid
            mid_rows["condition"] = cond
            all_samples[sid] = mid_rows[
                ["metabolite_name", "mid_index", "best_fit", "measured"]
            ].copy()

            if sid not in self.sample_conditions:
                self.sample_conditions[sid] = cond

        if not all_samples:
            print("[WARN] No MID samples parsed")
            return

        # Build a unified MID matrix by merging on (metabolite_name, mid_index),
        # since samples may have different sets of MID entries.
        def _build_matrix(value_col: str) -> pd.DataFrame:
            first_sid = list(all_samples.keys())[0]
            matrix = all_samples[first_sid][
                ["metabolite_name", "mid_index", value_col]
            ].rename(columns={value_col: first_sid})
            for sid, sdf in list(all_samples.items())[1:]:
                matrix = matrix.merge(
                    sdf[["metabolite_name", "mid_index", value_col]].rename(
                        columns={value_col: sid}
                    ),
                    on=["metabolite_name", "mid_index"],
                    how="outer",
                )
            return matrix.set_index(["metabolite_name", "mid_index"])

        mid_matrix = _build_matrix("best_fit")
        self.mid_df = mid_matrix
        self.mid_measured_df = _build_matrix("measured")

        self.conditions = sorted(set(self.sample_conditions.values()))
        print(f"[OK] MID data: {len(all_samples)} samples, {len(mid_matrix)} MID entries")

    def _extract_tissues(self):
        """Extract tissue list from flux column names."""
        tissues = set()
        if self.flux_df is not None:
            for col in self.flux_df.columns:
                kind, tissue, _ = _parse_flux_column(col)
                if kind == "flux" and tissue != "unknown":
                    tissues.add(tissue)
        self.tissues = sorted(tissues)
        self._ensure_output_dirs()

    # ── 分析方法（委托给 df_analysis）───────────────────────────────────

    def quality_control(self):
        """拟合质量评估。委托给 df_analysis.quality_control()。"""
        from .df_analysis import quality_control
        return quality_control(self.flux_all_fits, self.sample_conditions,
                               self.output_dir)

    def flux_pca(self, n_components=None):
        """通量 PCA 分析。委托给 df_analysis.flux_pca()。"""
        from .df_analysis import flux_pca
        return flux_pca(self.flux_df, self.sample_conditions,
                        self.output_dir, self.tissues, n_components)

    def differential_flux(self, comparisons=None, tissue_filter=None):
        """差异通量分析。委托给 df_analysis.differential_flux()。"""
        from .df_analysis import differential_flux
        return differential_flux(
            self.flux_df, self.flux_all_fits, self.sample_conditions,
            self.conditions, self.tissues, self.output_dir,
            comparisons=comparisons, tissue_filter=tissue_filter,
            analysis_config=self.analysis_config,
        )

    def pathway_flux_comparison(self, comparisons=None):
        """通路通量比较。委托给 df_analysis.pathway_flux_comparison()。"""
        from .df_analysis import pathway_flux_comparison
        return pathway_flux_comparison(
            self.flux_df, self.flux_all_fits, self.sample_conditions,
            self.conditions, self.tissues, self.pathway_groups,
            self.output_dir, comparisons=comparisons,
            analysis_config=self.analysis_config,
        )

    def flux_heatmap(self, comparisons=None, top_n=None, p_thresh=None):
        """热图数据计算。委托给 df_analysis.flux_heatmap()。"""
        from .df_analysis import flux_heatmap
        heatmap_cfg = self.analysis_config.get("heatmap", {})
        top_n = top_n if top_n is not None else int(heatmap_cfg.get("top_n", 50))
        p_thresh = p_thresh if p_thresh is not None else float(heatmap_cfg.get("pvalue_threshold", 0.05))
        return flux_heatmap(
            self.flux_df, self.flux_all_fits, self.sample_conditions,
            self.conditions, self.tissues, self.output_dir,
            comparisons=comparisons, top_n=top_n, p_thresh=p_thresh,
            analysis_config=self.analysis_config,
        )

    def tca_cycle_analysis(self):
        """TCA 循环分析（兼容 API，默认不再输出独立 tca_cycle 文件夹）。"""
        from .df_analysis import tca_cycle_analysis
        return tca_cycle_analysis(
            self.flux_df, self.sample_conditions, self.conditions,
            self.tissues, self.pathway_groups, self.output_dir,
        )

    def mid_fit_residuals(self):
        """实测 vs 模拟 MID 残差 QC。委托给 df_analysis.mid_fit_residuals()。"""
        from .df_analysis import mid_fit_residuals
        return mid_fit_residuals(
            self.mid_df, self.mid_measured_df, self.sample_conditions,
            self.output_dir,
        )

    def mid_analysis(self, top_metabolites=20):
        """MID 分析。委托给 df_analysis.mid_analysis()。"""
        from .df_analysis import mid_analysis
        return mid_analysis(
            self.mid_df, self.sample_conditions, self.conditions,
            self.output_dir, top_metabolites=top_metabolites,
        )

    def infusion_steady_state_qc(self):
        """灌注期血糖稳态 QC。委托给 infusion_qc 模块。

        MFA 的全部推断都以采样时已达代谢/同位素稳态为前提，这一步是对
        该前提最直接的检验。

        血糖表按动物编号记录，所以优先用 ``groups.yaml: mouse_groups``；
        没有时退回 ``sample_groups``（20260312 批次以动物编号为键）。
        """
        from .infusion_qc import infusion_steady_state_qc, select_animal_groups
        animal_groups = select_animal_groups(
            self.run_context.get("mouse_groups"),
            self.run_context.get("sample_groups"),
        )
        return infusion_steady_state_qc(
            self.run_context.get("glycemia_file"),
            self.output_dir,
            self.run_context.get("infusion_qc"),
            animal_groups,
        )

    def write_provenance(self):
        """写出本次运行的技术链指纹，供整合阶段判断可合并性。"""
        from .provenance import write_run_provenance

        reactions = set()
        if self.flux_df is not None:
            for col in self.flux_df.columns:
                kind, _, reaction = _parse_flux_column(col)
                if kind == "flux" and not reaction.endswith("__R"):
                    reactions.add(reaction)

        n_fit_iterations = None
        if self.flux_all_fits:
            counts = [len(df) for df in self.flux_all_fits.values()]
            if counts:
                n_fit_iterations = int(np.median(counts))

        return write_run_provenance(
            self.output_dir,
            experiment_id=self.run_context.get("experiment_id", ""),
            protocol_id=self.run_context.get("protocol_id"),
            batch_id=self.run_context.get("batch_id"),
            batch_metadata=self.run_context.get("batch_metadata"),
            tracer=self.run_context.get("tracer"),
            reactions=reactions,
            tissues=self.tissues,
            conditions=self.conditions,
            n_samples=len(self.flux_df.index) if self.flux_df is not None else 0,
            n_fit_iterations=n_fit_iterations,
            analysis_config=self.analysis_config,
        )

    def sample_flux_matrix(self):
        """样本级净通量矩阵。委托给 df_analysis.sample_flux_matrix()。

        跨批次桥接样本一致性检查（bridging_qc）依赖该表。
        """
        from .df_analysis import sample_flux_matrix
        return sample_flux_matrix(
            self.flux_df, self.flux_all_fits, self.sample_conditions,
            self.tissues, self.output_dir,
            analysis_config=self.analysis_config,
        )

    def tissue_flux_overview(self):
        """组织通量概览。委托给 df_analysis.tissue_flux_overview()。"""
        from .df_analysis import tissue_flux_overview
        return tissue_flux_overview(
            self.flux_df, self.sample_conditions, self.conditions,
            self.tissues, self.output_dir,
        )

    # ── Report Generation ────────────────────────────────────────────────

    def generate_report(self):
        """
        生成 MFA 解读报告（Step 9）。

        调用 generate_mfa_interpretation_report 模块，
        按组织（tissue）生成独立解读报告到 03_interpretation_report/ 目录，
        替代旧的单一 mfa_report.md 文件清单。
        """
        print("\n--- Generating Interpretation Report (Step 9) ---")
        from .generate_mfa_interpretation_report import generate_report_from_output_dir

        # 构建条件显示名映射（小写 key → 显示名）
        cond_display = {c: c.upper() for c in self.conditions}

        generate_report_from_output_dir(
            results_dir=self.output_dir,
            tissues=self.tissues,
            comp_tuples=None,
            cond_display=cond_display,
        )

    # ── Convenience Runner ───────────────────────────────────────────────

    def _generate_all_plots(self, comparisons=None):
        """Call the plotting subpackage to generate all visualizations."""
        from .plotting import generate_all_plots

        if comparisons is None:
            comparisons = list(combinations(self.conditions, 2)) if len(self.conditions) >= 2 else []
        generate_all_plots(
            output_dir=self.output_dir,
            condition_colors=self.condition_colors,
            tissue_colors=self.tissue_colors,
            comparisons=comparisons,
            tissues=self.tissues,
            conditions=self.conditions,
            sample_conditions=self.sample_conditions,
            analysis_config=self.analysis_config,
        )

    def run_all(self, comparisons=None):
        """Run the full MFA analysis pipeline."""
        if comparisons is None and len(self.conditions) >= 2:
            comparisons = list(combinations(self.conditions, 2))

        # 统一小写：确保 comparisons 的条件名与 sample_conditions 一致
        comparisons = [(a.lower(), b.lower()) for a, b in comparisons]

        # 确保组织特异性输出目录已创建（幂等）
        self._ensure_output_dirs()

        print("=" * 70)
        print("MFA Analysis Pipeline")
        print("=" * 70)

        self.quality_control()
        self.infusion_steady_state_qc()
        self.write_provenance()
        self.sample_flux_matrix()
        self.flux_pca()
        self.differential_flux(comparisons=comparisons)
        self.flux_heatmap(comparisons=comparisons)
        self.pathway_flux_comparison(comparisons=comparisons)
        self.tissue_flux_overview()
        self.mid_analysis()
        self.mid_fit_residuals()
        self._generate_all_plots(comparisons=comparisons)
        self.generate_report()

        print("\n" + "=" * 70)
        print("MFA analysis complete!")
        print(f"Results: {self.output_dir}")
        print("=" * 70)
