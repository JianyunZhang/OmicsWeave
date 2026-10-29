"""
转录组学数据加载器 —— DE 差异表达 & GSEA / Enrichr 富集结果。

从 transcriptomics_results 标准目录结构加载数据，统一列名为:
  - gene
  - log2fc  (原始列 log2FoldChange)
  - pvalue
  - padj
  - significant
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import pandas as pd


# ── 已知组织名（长名优先，避免 "sq_adipose" 误切为 "sq"）──────────
KNOWN_TISSUES = [
    "sq_adipose",
    "gast_muscle",
    "go_adipose",
    "renal_cortex",
    "cerebral_cortex",
    "brain",
    "heart",
    "kidney",
    "liver",
]


def _split_tissue_comp(stem: str) -> tuple[Optional[str], str]:
    """
    将 'liver_DIO_vs_WT' 拆成 ('liver', 'DIO_vs_WT')。
    多下划线组织名（sq_adipose / gast_muscle）优先匹配。
    """
    for t in sorted(KNOWN_TISSUES, key=len, reverse=True):
        if stem.startswith(t + "_"):
            return t, stem[len(t) + 1 :]
    return None, stem


# ── 列名映射表 ────────────────────────────────────────────────────
_DE_COL_RENAME = {
    "log2FoldChange": "log2fc",
    "pooled_log2FoldChange": "log2fc",
    "pooled_se": "lfcSE",
    "pooled_pvalue": "pvalue",
    "pooled_fdr": "padj",
    "i2": "heterogeneity_i2",
}


class TranscriptomicsLoader:
    """
    Parameters
    ----------
    config : object
        需具备以下属性:
        - transcriptomics_results : Path | None  转录组分析结果根目录
        - tissues : list[str]                     组织列表
        - comparison_map : dict[str, str]         多组学比较名 → 转录组比较名
        - mfa_tissue_map : dict[str, str]         (可选) MFA 组织名映射
    """

    def __init__(self, config) -> None:
        self.tx_dir: Optional[Path] = (
            Path(config.transcriptomics_results)
            if config.transcriptomics_results is not None
            else None
        )
        self.tissues: list[str] = list(config.tissues)
        # Prefer transcriptomics_comparison_map (explicit name); fall back to comparison_map
        tx_cmap = getattr(config, "transcriptomics_comparison_map", None)
        if tx_cmap:
            self.comparison_map: dict[str, str] = dict(tx_cmap)
        else:
            self.comparison_map: dict[str, str] = dict(config.comparison_map)
        self.mfa_tissue_map: dict[str, str] = getattr(config, "mfa_tissue_map", {})
        integration = getattr(config, "integration", {}) or {}
        self._integration_input = (
            ((integration.get("inputs", {}) or {}).get("transcriptomics", {}) or {})
        )
        self._canonical_comparison_map = {
            str(source): str(canonical)
            for canonical, source in self.comparison_map.items()
            if source is not None
        }
        self._configured_comparisons = {str(value) for value in self.comparison_map}
        self._configured_source_comparisons = {
            str(value) for value in self.comparison_map.values() if value is not None
        }
        self._sample_metadata: Optional[pd.DataFrame] = None
        self._sample_metadata_loaded = False

        # 存放加载结果
        self.de: dict[str, dict[str, pd.DataFrame]] = {}   # {comparison: {tissue: df}}
        self.gsea: dict[str, dict[str, pd.DataFrame]] = {}  # {tissue: {comparison: df}}
        self.enrichment: dict[str, dict[str, pd.DataFrame]] = {}  # {tissue: {comp_dir: df}}

    # ── 公开接口 ──────────────────────────────────────────────────

    def load(self) -> dict:
        """
        加载全部转录组学数据并返回标准化结果字典。

        Returns
        -------
        dict
            {"de": {comparison: {tissue: DataFrame}},
             "gsea": {tissue: {comparison: DataFrame}},
             "enrichment": {tissue: {comparison_dir: DataFrame}}}
        """
        if self.tx_dir is None or not self.tx_dir.exists():
            print(f"[WARN] 转录组结果目录不存在或未配置: {self.tx_dir}")
            return {"de": self.de, "gsea": self.gsea, "enrichment": self.enrichment}

        if str(self._integration_input.get("mode", "single_run")).lower() in {
            "integrated", "auto",
        }:
            loaded_integrated = self._load_integrated_de_results()
            if not loaded_integrated and str(self._integration_input.get("mode", "")).lower() == "integrated":
                raise FileNotFoundError(
                    "Configured integrated transcriptomics results were not found under "
                    f"{self.tx_dir / 'integration' / 'within_protocol'}"
                )
            if not loaded_integrated:
                self._load_de_results()
        else:
            self._load_de_results()
        self._load_gsea_results()
        self._load_enrichment_results()

        n_de = sum(len(v) for v in self.de.values())
        n_gsea = sum(len(v) for v in self.gsea.values())
        n_enr = sum(len(v) for v in self.enrichment.values())
        print(
            f"[OK] 转录组加载完成: "
            f"{n_de} 个 DE 文件, "
            f"{n_gsea} 个 GSEA 文件, "
            f"{n_enr} 个 Enrichment 文件"
        )

        return {
            "de": self.de,
            "gsea": self.gsea,
            "enrichment": self.enrichment,
        }

    # ── DE 差异表达 ───────────────────────────────────────────────

    def _load_de_results(self) -> None:
        """
        加载 {tx_dir}/04_de_results/{tissue}/de_{tx_comparison}.csv。
        优先使用 de_for_multiomics_* 标准化文件（固定列名），
        回退到原始 de_*.csv 并做列名标准化。
        使用 comparison_map 将多组学比较名映射到转录组文件名中的比较名。
        """
        de_dir = self.tx_dir / "04_de_results"
        if not de_dir.exists():
            print(f"[WARN] DE 结果目录不存在: {de_dir}")
            return

        missing_by_comparison: dict[str, list[str]] = {}
        for multi_comp, tx_comp in self.comparison_map.items():
            self.de.setdefault(multi_comp, {})
            for tissue in self.tissues:
                # 将 MFA 组织名映射到转录组组织名
                tx_tissue = self.mfa_tissue_map.get(tissue, tissue)

                # Priority 1: multi-omics standardized export. 兼容当前 pipeline
                # 的 de_{tissue}_for_multiomics_* 与早期 de_for_multiomics_* 命名。
                standardized_candidates = [
                    de_dir / tx_tissue / f"de_{tx_tissue}_for_multiomics_{tx_comp}.csv",
                    de_dir / tx_tissue / f"de_{tissue}_for_multiomics_{tx_comp}.csv",
                    de_dir / tx_tissue / f"de_for_multiomics_{tx_comp}.csv",
                ]
                fp_mo = next(
                    (candidate for candidate in standardized_candidates if candidate.exists()),
                    None,
                )
                if fp_mo is not None:
                    df = pd.read_csv(fp_mo)
                    # Standardized file already has: gene_symbol, log2fc, padj, etc.
                    if "gene_symbol" in df.columns:
                        df = df.rename(columns={"gene_symbol": "gene"})
                    df = self._annotate_sample_counts(df, tx_tissue, tx_comp)
                    self.de[multi_comp][tissue] = df
                    continue

                # Priority 2: primary DE file — try multiple naming patterns
                fp = None
                candidates = [
                    de_dir / tx_tissue / f"de_{tx_comp}.csv",
                    de_dir / tx_tissue / f"de_{tx_tissue}_{tx_comp}.csv",
                    de_dir / tx_tissue / f"de_{tissue}_{tx_comp}.csv",
                    de_dir / f"de_{tx_tissue}_{tx_comp}.csv",
                ]
                for c in candidates:
                    if c.exists():
                        fp = c
                        break

                if fp is None:
                    # Try wilcoxon variants
                    for c in [
                        de_dir / tx_tissue / f"de_{tx_comp}_wilcoxon.csv",
                        de_dir / tx_tissue / f"de_{tx_tissue}_{tx_comp}_wilcoxon.csv",
                        de_dir / tx_tissue / f"de_{tissue}_{tx_comp}_wilcoxon.csv",
                    ]:
                        if c.exists():
                            fp = c
                            break

                if fp is None:
                    missing_by_comparison.setdefault(str(multi_comp), []).append(
                        str(tissue)
                    )
                    continue

                df = pd.read_csv(fp)
                df = self._normalize_de_columns(df)
                df = self._annotate_sample_counts(df, tx_tissue, tx_comp)
                self.de[multi_comp][tissue] = df
        for comparison, tissues in missing_by_comparison.items():
            print(
                f"[WARN] 转录组比较 {comparison} 缺少 {len(tissues)} 个组织的 DE 数据: "
                + ", ".join(tissues)
            )

    def _load_integrated_de_results(self) -> bool:
        """读取同一 protocol 内的转录组批次 Meta 结果。"""
        protocol_id = str(self._integration_input.get("protocol_id", "")).strip()
        base = self.tx_dir / "integration" / "within_protocol"
        protocol_dirs: list[Path] = []
        if protocol_id:
            protocol_dirs.append(base / protocol_id)
        if base.exists():
            protocol_dirs.extend(sorted(p for p in base.iterdir() if p.is_dir()))
        root = next((p for p in protocol_dirs if p.exists()), None)
        if root is None:
            return False

        loaded = False
        for multi_comp, tx_comp in self.comparison_map.items():
            for tissue in self.tissues:
                tx_tissue = self.mfa_tissue_map.get(tissue, tissue)
                candidates = [
                    root / tx_tissue / f"de_{tx_comp}_integrated.csv",
                    root / tx_tissue / f"de_{multi_comp}_integrated.csv",
                ]
                fp = next((path for path in candidates if path.exists()), None)
                if fp is None:
                    continue
                df = self._normalize_de_columns(pd.read_csv(fp))
                df["source_mode"] = "within_protocol_meta"
                df["protocol_id"] = root.name
                if "meta_note" not in df.columns:
                    df["meta_note"] = "random_effects_meta"
                self.de.setdefault(multi_comp, {})[tissue] = df
                loaded = True
        if loaded:
            print(f"[INFO] 使用转录组批次整合结果: {root}")
        return loaded

    # ── GSEA ─────────────────────────────────────────────────────

    def _load_gsea_results(self) -> None:
        """
        加载 {tx_dir}/05_result_plot/5_gsea/gsea_{tissue}_{comparison}_results.csv。
        通过扫描目录自动发现组织与比较名。
        """
        gsea_dir = self.tx_dir / "05_result_plot" / "5_gsea"
        if not gsea_dir.exists():
            print(f"[WARN] GSEA 结果目录不存在: {gsea_dir}")
            return

        for f in sorted(gsea_dir.glob("gsea_*_results.csv")):
            # 文件名形如 gsea_liver_DIO_vs_WT_results.csv
            stem_part = f.stem.replace("gsea_", "").replace("_results", "")
            tissue, comp = _split_tissue_comp(stem_part)
            if tissue is None:
                print(f"[WARN] 无法解析 GSEA 文件名: {f.name}")
                continue

            if self.comparison_map and not self._is_configured_comparison(comp):
                continue
            canonical_comp = self._canonical_comparison_map.get(comp, comp)
            self.gsea.setdefault(tissue, {})[canonical_comp] = pd.read_csv(f)

            if tissue not in self.tissues:
                self.tissues.append(tissue)

    # ── Enrichr 富集 ─────────────────────────────────────────────

    def _load_enrichment_results(self) -> None:
        """
        加载 {tx_dir}/05_result_plot/4_enrichment/enrichr_*.csv （可选）。
        """
        enr_dir = self.tx_dir / "05_result_plot" / "4_enrichment"
        if not enr_dir.exists():
            # enrichment 目录缺失不视为错误
            return

        for f in sorted(enr_dir.glob("enrichr_*.csv")):
            # 文件名形如 enrichr_liver_DIO_vs_WT_up.csv
            stem_part = f.stem.replace("enrichr_", "")
            tissue, comp_dir = _split_tissue_comp(stem_part)
            if tissue is None:
                print(f"[WARN] 无法解析 Enrichr 文件名: {f.name}")
                continue

            canonical_comp = comp_dir
            matched_configured = False
            for source, target in sorted(
                self._canonical_comparison_map.items(), key=lambda item: len(item[0]), reverse=True
            ):
                if comp_dir == source or comp_dir.startswith(source + "_"):
                    canonical_comp = target + comp_dir[len(source):]
                    matched_configured = True
                    break
            if not matched_configured:
                for canonical in sorted(self._configured_comparisons, key=len, reverse=True):
                    if comp_dir == canonical or comp_dir.startswith(canonical + "_"):
                        matched_configured = True
                        break
            if self.comparison_map and not matched_configured:
                continue
            self.enrichment.setdefault(tissue, {})[canonical_comp] = pd.read_csv(f)

            if tissue not in self.tissues:
                self.tissues.append(tissue)

    # ── 列名标准化 ───────────────────────────────────────────────

    @staticmethod
    def _normalize_de_columns(df: pd.DataFrame) -> pd.DataFrame:
        """
        将 DE DataFrame 的列名统一为:
          gene, log2fc, pvalue, padj, significant
        保留其它列不变。
        """
        df = df.rename(columns=_DE_COL_RENAME)
        return df

    def _is_configured_comparison(self, comparison: str) -> bool:
        """仅接受配置中显式声明的 canonical 或上游比较名。"""
        return (
            comparison in self._configured_comparisons
            or comparison in self._configured_source_comparisons
        )

    def _get_sample_metadata(self) -> Optional[pd.DataFrame]:
        """延迟读取 RNA 样本元数据；缺少显式字段时不猜测样本量。"""
        if self._sample_metadata_loaded:
            return self._sample_metadata
        self._sample_metadata_loaded = True
        if self.tx_dir is None:
            return None
        metadata_path = self.tx_dir / "sample_metadata.csv"
        if not metadata_path.exists():
            return None
        metadata = pd.read_csv(metadata_path)
        required = {"tissue", "condition"}
        if not required.issubset(metadata.columns):
            print(
                f"[WARN] RNA 样本元数据缺少字段 {sorted(required - set(metadata.columns))}: "
                f"{metadata_path}"
            )
            return None
        if "mouse_id" not in metadata.columns and "sample" not in metadata.columns:
            print(f"[WARN] RNA 样本元数据缺少 mouse_id/sample，不能审计 n_eff: {metadata_path}")
            return None
        self._sample_metadata = metadata
        return metadata

    def _annotate_sample_counts(
        self,
        df: pd.DataFrame,
        tissue: str,
        source_comparison: str,
    ) -> pd.DataFrame:
        """从显式样本表按组织与条件回填生物学重复数。"""
        if "_vs_" not in source_comparison:
            return df
        case_condition, control_condition = source_comparison.split("_vs_", 1)
        metadata = self._get_sample_metadata()
        if metadata is None:
            return df

        tissue_rows = metadata.loc[
            metadata["tissue"].astype(str).eq(str(tissue))
        ].copy()
        id_column = "mouse_id" if "mouse_id" in tissue_rows.columns else "sample"
        count_rows = tissue_rows.dropna(subset=[id_column]).copy()
        counts = count_rows.groupby(
            count_rows["condition"].astype(str)
        )[id_column].nunique()
        n_case = counts.get(case_condition)
        n_control = counts.get(control_condition)
        if pd.isna(n_case) or pd.isna(n_control):
            return df

        result = df.copy()
        filled_from_metadata = pd.Series(False, index=result.index)
        for column, value in (("n_case", int(n_case)), ("n_control", int(n_control))):
            if column not in result.columns:
                result[column] = value
                filled_from_metadata[:] = True
            else:
                existing = pd.to_numeric(result[column], errors="coerce")
                missing = existing.isna()
                result[column] = existing.mask(missing, value)
                filled_from_metadata |= missing

        calculated_n = (
            pd.to_numeric(result["n_case"], errors="coerce")
            + pd.to_numeric(result["n_control"], errors="coerce")
        )
        if "n_eff" not in result.columns:
            result["n_eff"] = calculated_n
            filled_from_metadata[:] = True
        else:
            existing_n = pd.to_numeric(result["n_eff"], errors="coerce")
            missing_n = existing_n.isna()
            result["n_eff"] = existing_n.mask(missing_n, calculated_n)
            filled_from_metadata |= missing_n

        if "n_eff_source" not in result.columns:
            result["n_eff_source"] = "de_result"
        result.loc[filled_from_metadata, "n_eff_source"] = "sample_metadata"
        return result
