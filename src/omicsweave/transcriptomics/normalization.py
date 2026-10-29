"""
基因层聚合与文库归一化
========================

本模块提供两类纯函数（不依赖全局配置，便于单元测试）：

1. **转录本 → 基因层聚合**：``build_gene_level_counts()``
   默认使用 tximport 的 ``lengthScaledTPM`` 策略，把转录本有效长度
   （kallisto ``abundance.tsv`` 的 ``eff_length`` 列）纳入计算。

   直接对 ``est_counts`` 求和会在两组之间发生**异构体切换**
   （isoform switching）时产生假阳性差异表达：基因总 count 变化了，
   但真实表达量没变。tximport 存在的全部理由就是消除这一偏差。

2. **文库大小归一化**：``median_of_ratios_size_factors()`` / ``normalize_counts()``
   DESeq2 的 median-of-ratios 方法。非参数检验（Wilcoxon）必须先做
   归一化，否则测序深度差异会在**所有基因**上造成系统性 log2FC 偏移。
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd

# 关闭 pandas 的 numexpr 加速路径。本模块的 lengthScaledTPM 聚合会做
# `DataFrame.mul(Series, axis=0)`；在 35,983 基因 × 95 样本（341 万元素）
# 下它越过 pandas 的 numexpr 阈值，实测约 93% 的输出被静默填成 0，
# 且每次运行错的位置都不同。19 样本（68 万元素）时完全正常，所以小规模
# 单元测试测不到。详见 omicsweave/common/pandas_safety.py。
from omicsweave.common.pandas_safety import harden_pandas

harden_pandas()

# 聚合方法标识，会写入 DE 结果的溯源列
AGGREGATION_LENGTH_SCALED_TPM = "lengthScaledTPM"
AGGREGATION_SIMPLE_SUM = "simple_sum"
VALID_AGGREGATION_METHODS = (AGGREGATION_LENGTH_SCALED_TPM, AGGREGATION_SIMPLE_SUM)


def map_transcripts_to_genes(
    transcript_ids: pd.Index,
    tx2gene: dict[str, dict[str, str]],
    key: str = "gene_symbol",
) -> pd.Series:
    """
    把转录本 ID 映射到基因标识。

    优先精确匹配；未命中时去掉 Ensembl 版本号后缀（``ENSMUST0000000001.5``
    → ``ENSMUST0000000001``）再匹配。仍未命中的转录本保留原 ID，
    使其自成一个「基因」而不是被静默丢弃。

    这里用一次性构建的去版本号索引做 O(1) 查找；旧实现对每个未命中的
    转录本都要全表扫描 ``tx2gene``，在 ~12 万转录本规模下是 O(N×M)。

    Parameters
    ----------
    transcript_ids : pd.Index
        计数矩阵的转录本 ID。
    tx2gene : dict
        ``{transcript_id: {"gene_id": ..., "gene_symbol": ...}}``。
    key : str
        取 ``gene_symbol`` 还是 ``gene_id``。

    Returns
    -------
    pd.Series
        以 transcript_id 为索引、基因标识为值。
    """
    # 去版本号索引只建一次
    versionless: dict[str, str] = {}
    for tx_id, info in tx2gene.items():
        base = tx_id.split(".")[0]
        # 首次写入即可；同一 base 对应多版本时取任意一个（基因归属一致）
        if base not in versionless:
            versionless[base] = info.get(key) or info.get("gene_id") or tx_id

    mapped = []
    for tx_id in transcript_ids:
        info = tx2gene.get(tx_id)
        if info is not None:
            mapped.append(info.get(key) or info.get("gene_id") or tx_id)
            continue
        mapped.append(versionless.get(str(tx_id).split(".")[0], tx_id))

    return pd.Series(mapped, index=transcript_ids, name="gene")


def _gene_effective_lengths(
    tpm_df: pd.DataFrame,
    efflen_df: pd.DataFrame,
    gene_key: pd.Series,
) -> pd.DataFrame:
    """
    计算基因层有效长度：以转录本丰度为权重的加权平均。

    对丰度全为 0 的基因（权重和为 0），退回该基因所有转录本有效长度的
    简单平均——否则会得到 0/0。
    """
    gene_abundance = tpm_df.groupby(gene_key, sort=False).sum()
    weighted = (tpm_df * efflen_df).groupby(gene_key, sort=False).sum()

    with np.errstate(divide="ignore", invalid="ignore"):
        gene_len = weighted / gene_abundance

    unweighted = efflen_df.groupby(gene_key, sort=False).mean()
    unweighted = unweighted.reindex(index=gene_len.index, columns=gene_len.columns)

    gene_len = gene_len.replace([np.inf, -np.inf], np.nan)
    gene_len = gene_len.where(gene_abundance > 0, unweighted)
    gene_len = gene_len.fillna(unweighted)

    # 仍为 NaN/0 的（转录本长度缺失）用全局中位数兜底，避免下游除零
    finite_vals = gene_len.to_numpy(dtype=float)
    finite_vals = finite_vals[np.isfinite(finite_vals) & (finite_vals > 0)]
    fallback = float(np.median(finite_vals)) if finite_vals.size else 1.0
    gene_len = gene_len.fillna(fallback)
    gene_len = gene_len.mask(gene_len <= 0, fallback)

    return gene_len


# 长度加权只改变每个基因**分到多少** count，不会让基因整体消失。实测正常
# 运行该比值恒在 0.95 以上；2026-09-01 的 numexpr 事故中它塌到 0.11。
# 0.5 因此是一条很宽的安全线，不会误报。
MIN_NONZERO_RATIO = 0.5


def _nonzero_ratio(
    gene_counts: pd.DataFrame,
    raw_gene_counts: pd.DataFrame,
    samples: list[str],
) -> pd.Series:
    """每样本：长度加权聚合的非零基因数 ÷ 直接求和的非零基因数。"""
    raw_nonzero = (raw_gene_counts[samples] > 0).sum(axis=0)
    out_nonzero = (gene_counts[samples] > 0).sum(axis=0)
    return out_nonzero / raw_nonzero.where(raw_nonzero > 0)


def _assert_aggregation_sane(
    gene_counts: pd.DataFrame,
    raw_gene_counts: pd.DataFrame,
    samples: list[str],
    min_ratio: float = MIN_NONZERO_RATIO,
) -> None:
    """
    对聚合**输出**做不变量校验：非零基因数不应远低于直接求和的结果。

    这道闸门刻意不去判断成因。它抓的是症状，因此对 numexpr 静默算错、
    未来的 pandas 广播 bug、索引对齐写错等等一视同仁。之所以必须校验
    「非零个数」而不是「总和」：``factors = old_sum / new_sum`` 这一步会
    把列和重标定回正确的文库大小，总和自检在事故中是**通过**的。
    """
    ratio = _nonzero_ratio(gene_counts, raw_gene_counts, samples).dropna()
    bad = ratio[ratio < min_ratio]
    if bad.empty:
        return

    worst = bad.idxmin()
    raise RuntimeError(
        "基因层聚合结果不可信：长度加权后非零基因数相对直接求和大幅塌缩，"
        "拒绝返回该矩阵。\n"
        f"  受影响样本 : {len(bad)} / {len(ratio)}（阈值 {min_ratio:.0%}）\n"
        f"  最差样本   : {worst} 保留 {bad.min():.1%} 的非零基因\n"
        "  常见成因   : pandas 的 numexpr 路径在 >100 万元素时静默算错；\n"
        "               请确认 omicsweave/common/pandas_safety.py 的 "
        "harden_pandas() 已生效，\n"
        "               或用 assert_broadcast_ops_sane() 探测当前环境。"
    )


def build_gene_level_counts(
    counts_df: pd.DataFrame,
    tpm_df: pd.DataFrame,
    efflen_df: pd.DataFrame,
    tx2gene: dict[str, dict[str, str]],
    method: str = AGGREGATION_LENGTH_SCALED_TPM,
    gene_key_field: str = "gene_symbol",
) -> tuple[pd.DataFrame, dict[str, object]]:
    """
    把转录本层矩阵聚合到基因层。

    ``lengthScaledTPM``（默认，等价于 tximport ``countsFromAbundance``）：

    1. 基因 TPM = 该基因所有转录本 TPM 之和
    2. 基因有效长度 = 以转录本 TPM 为权重的有效长度加权平均
    3. counts = 基因 TPM × 基因有效长度在样本间的均值
    4. 按列重标定，使每个样本的总 count 与原始求和结果一致
       （保留文库大小信息，供 DESeq2 估计离散度）

    ``simple_sum`` 保留旧行为（直接对 est_counts 求和），仅用于回溯比对。

    Returns
    -------
    (gene_counts, info)
        ``gene_counts`` 为整数矩阵（基因 × 样本）；``info`` 记录方法与
        规模，用于写入溯源信息。
    """
    if method not in VALID_AGGREGATION_METHODS:
        raise ValueError(
            f"Unknown gene aggregation method {method!r}; "
            f"expected one of {VALID_AGGREGATION_METHODS}"
        )

    gene_key = map_transcripts_to_genes(counts_df.index, tx2gene, key=gene_key_field)

    raw_gene_counts = counts_df.groupby(gene_key, sort=False).sum()

    if method == AGGREGATION_SIMPLE_SUM:
        out = raw_gene_counts.round(0).astype(int)
        return out, {
            "gene_aggregation": AGGREGATION_SIMPLE_SUM,
            "n_transcripts": int(counts_df.shape[0]),
            "n_genes": int(out.shape[0]),
            "length_corrected": False,
        }

    # 列对齐后再做逐元素运算，防止样本顺序不一致导致错位
    common_samples = [c for c in counts_df.columns if c in tpm_df.columns and c in efflen_df.columns]
    if not common_samples:
        raise ValueError("counts/tpm/eff_length matrices share no common samples")
    tpm_sub = tpm_df.loc[counts_df.index, common_samples]
    efflen_sub = efflen_df.loc[counts_df.index, common_samples]

    gene_len = _gene_effective_lengths(tpm_sub, efflen_sub, gene_key)
    gene_tpm = tpm_sub.groupby(gene_key, sort=False).sum()

    # counts ∝ 丰度 × 长度
    mean_len = gene_len.mean(axis=1)
    scaled = gene_tpm.mul(mean_len, axis=0)

    # 按列重标定回原始文库大小
    old_sum = raw_gene_counts[common_samples].sum(axis=0)
    new_sum = scaled.sum(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        factors = old_sum / new_sum
    factors = factors.replace([np.inf, -np.inf], np.nan).fillna(1.0)
    scaled = scaled.mul(factors, axis=1)

    gene_counts = scaled.round(0).fillna(0).astype(int)
    gene_counts[gene_counts < 0] = 0

    _assert_aggregation_sane(gene_counts, raw_gene_counts, common_samples)

    info = {
        "gene_aggregation": AGGREGATION_LENGTH_SCALED_TPM,
        "n_transcripts": int(counts_df.shape[0]),
        "n_genes": int(gene_counts.shape[0]),
        "length_corrected": True,
        "mean_gene_eff_length": float(mean_len.mean()),
        "min_nonzero_ratio": float(
            _nonzero_ratio(gene_counts, raw_gene_counts, common_samples).min()
        ),
    }

    nonzero_ratio = (gene_counts > 0).sum(axis=0) / (raw_gene_counts[common_samples] > 0).sum(axis=0)
    if (nonzero_ratio < 0.5).any():
        raise RuntimeError(
            f"聚合后非零基因数不足原始的 50%（最低 {nonzero_ratio.min():.1%}）——"
            "疑似 pandas/numexpr 广播损坏，请检查 compute.use_numexpr"
        )

    return gene_counts, info


# ---------------------------------------------------------------------------
# 文库大小归一化
# ---------------------------------------------------------------------------

def median_of_ratios_size_factors(counts_df: pd.DataFrame) -> pd.Series:
    """
    DESeq2 median-of-ratios 大小因子（基因 × 样本 → 每样本一个因子）。

    仅使用在所有样本中都为正的基因计算几何均值参考。若不存在这样的
    基因（低深度或稀疏矩阵），退回按文库大小归一化，并保证返回值为正。
    """
    counts = counts_df.astype(float)
    if counts.shape[1] == 0:
        return pd.Series(dtype=float)

    positive_everywhere = (counts > 0).all(axis=1)
    if positive_everywhere.any():
        log_counts = np.log(counts[positive_everywhere])
        log_geo_mean = log_counts.mean(axis=1)
        ratios = log_counts.sub(log_geo_mean, axis=0)
        size_factors = np.exp(ratios.median(axis=0))
    else:
        lib = counts.sum(axis=0)
        mean_lib = lib[lib > 0].mean() if (lib > 0).any() else 1.0
        size_factors = lib / mean_lib if mean_lib else pd.Series(1.0, index=counts.columns)

    size_factors = size_factors.replace([np.inf, -np.inf], np.nan).fillna(1.0)
    size_factors[size_factors <= 0] = 1.0
    size_factors.name = "size_factor"
    return size_factors


def normalize_counts(counts_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.Series]:
    """
    用 median-of-ratios 大小因子归一化计数矩阵。

    Returns
    -------
    (normalized_counts, size_factors)
    """
    size_factors = median_of_ratios_size_factors(counts_df)
    normalized = counts_df.astype(float).div(size_factors, axis=1)
    return normalized, size_factors
