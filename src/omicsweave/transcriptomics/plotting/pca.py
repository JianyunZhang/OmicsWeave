"""
PCA and Sample Clustering Plots
================================
Generates global PCA, per-tissue PCA, sample-sample correlation heatmap,
批次着色 PCA，以及主成分与协变量的关联检验。

合并多批次数据前必须先看到批次效应有多大。此前本模块完全没有批次相关
逻辑，等于把最关键的诊断留空。
"""

import os

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from . import _config as cfg

# 参与「主成分 ~ 协变量」关联检验的候选协变量
_COVARIATES = ("condition", "tissue", "batch")
_N_PCS_FOR_ASSOCIATION = 5


def generate_all():
    """Generate PCA and hierarchical clustering plots."""
    out_dir = os.path.join(cfg.PLOT_DIR, '3_pca_clustering')
    os.makedirs(out_dir, exist_ok=True)
    cfg.flush_print("\n[3/5] Generating PCA & Clustering Plots...")

    counts_path = os.path.join(cfg.DE_DIR, 'gene_counts_matrix.csv')
    counts_df = pd.read_csv(counts_path, index_col=0, low_memory=False)
    counts_df = counts_df.apply(pd.to_numeric, errors='coerce').fillna(0)
    # Drop samples (columns) that are entirely non-numeric
    counts_df = counts_df.loc[:, counts_df.sum() > 0]

    sample_batches = cfg.get_sample_batches()

    # ---- Global PCA (all samples, top 5000 variable genes) ----
    _plot_global_pca(counts_df, out_dir, sample_batches)

    # ---- Per-tissue PCA with mouse_id annotations ----
    _plot_per_tissue_pca(counts_df, out_dir)

    # ---- Sample-sample correlation heatmap (Pearson) ----
    _plot_correlation_heatmap(counts_df, out_dir)

    # ---- PC ~ 协变量关联检验（批次效应诊断） ----
    _pc_covariate_association(counts_df, out_dir, sample_batches)

    cfg.flush_print("  PCA & Clustering done.")


def _build_sample_annotation(sample_ids, sample_batches):
    """从样本名与配置构造注释表：mouse_id / tissue / condition / batch。"""
    ann = pd.DataFrame(index=pd.Index(sample_ids, name='sample'))
    ann['mouse_id'] = [str(s).split('_')[0] for s in sample_ids]
    ann['tissue'] = [
        str(s).split('_', 1)[1] if '_' in str(s) else 'unknown' for s in sample_ids
    ]
    ann['condition'] = ann['mouse_id'].map(cfg.MOUSE_GROUPS)
    if sample_batches:
        ann['batch'] = [
            sample_batches.get(str(s)) or sample_batches.get(str(s).split('_')[0])
            for s in sample_ids
        ]
    return ann


def _log_cpm_top_variable(counts_df, n_top=5000):
    """
    log2-CPM 变换后取高变基因。

    先做 CPM 再取对数，可以避免文库大小差异主导第一主成分——
    否则 PC1 往往只是测序深度。
    """
    lib = counts_df.sum(axis=0)
    lib = lib.replace(0, np.nan)
    cpm = counts_df.div(lib, axis=1) * 1e6
    log_cpm = np.log2(cpm.fillna(0) + 1).T  # samples x genes
    gene_var = log_cpm.var()
    top = gene_var.nlargest(min(n_top, len(gene_var))).index
    return log_cpm[top]


def _plot_global_pca(counts_df, out_dir, sample_batches=None):
    """Global PCA of all samples, colored by group / tissue / batch."""
    cfg.flush_print(f"  Global PCA (all {counts_df.shape[1]} samples)...")
    log_sub = _log_cpm_top_variable(counts_df, n_top=5000)

    scaler = StandardScaler()
    scaled = scaler.fit_transform(log_sub)

    pca = PCA(n_components=min(10, len(log_sub)))
    pcs = pca.fit_transform(scaled)
    pc_df = pd.DataFrame(pcs[:, :3], index=log_sub.index, columns=['PC1', 'PC2', 'PC3'])

    # Add metadata
    ann = _build_sample_annotation(log_sub.index, sample_batches or {})
    pc_df['mouse_id'] = ann['mouse_id']
    pc_df['tissue'] = ann['tissue']
    pc_df['group'] = ann['condition']
    has_batch = 'batch' in ann.columns and ann['batch'].notna().any()
    if has_batch:
        pc_df['batch'] = ann['batch']

    var_explained = pca.explained_variance_ratio_ * 100

    n_panels = 3 if has_batch else 2
    fig, axes = plt.subplots(1, n_panels, figsize=(7 * n_panels, 6))

    # PCA colored by group
    for grp, color in cfg.GROUP_COLORS.items():
        mask = pc_df['group'] == grp
        axes[0].scatter(
            pc_df.loc[mask, 'PC1'], pc_df.loc[mask, 'PC2'],
            c=color, label=grp, s=40, alpha=0.8,
            edgecolors='white', linewidths=0.5,
        )
    axes[0].set_xlabel(f'PC1 ({var_explained[0]:.1f}%)')
    axes[0].set_ylabel(f'PC2 ({var_explained[1]:.1f}%)')
    axes[0].set_title('PCA \u2014 Colored by Group')
    axes[0].legend(fontsize=9)
    axes[0].axhline(0, ls='--', lw=0.3, color='grey')
    axes[0].axvline(0, ls='--', lw=0.3, color='grey')

    # PCA colored by tissue
    tissue_cmap = plt.get_cmap(cfg.PLOT_COLORS.get('pca_tissue_cmap', 'Set2'))(
        np.linspace(0, 1, len(cfg.TISSUES))
    )
    tissue_color_map = dict(zip(cfg.TISSUES, tissue_cmap))
    for tis in cfg.TISSUES:
        mask = pc_df['tissue'] == tis
        axes[1].scatter(
            pc_df.loc[mask, 'PC1'], pc_df.loc[mask, 'PC2'],
            c=[tissue_color_map[tis]], label=cfg.TISSUE_LABELS.get(tis, tis),
            s=40, alpha=0.8, edgecolors='white', linewidths=0.5,
        )
    axes[1].set_xlabel(f'PC1 ({var_explained[0]:.1f}%)')
    axes[1].set_ylabel(f'PC2 ({var_explained[1]:.1f}%)')
    axes[1].set_title('PCA \u2014 Colored by Tissue')
    axes[1].legend(fontsize=8, loc='upper right')
    axes[1].axhline(0, ls='--', lw=0.3, color='grey')
    axes[1].axvline(0, ls='--', lw=0.3, color='grey')

    # PCA colored by batch —— 与按分组着色并排，便于直接判断批次是否主导分离
    if has_batch:
        batches = sorted(pc_df['batch'].dropna().astype(str).unique())
        batch_cmap = plt.get_cmap('tab10')(np.linspace(0, 1, max(len(batches), 2)))
        batch_color_map = dict(zip(batches, batch_cmap))
        for bat in batches:
            mask = pc_df['batch'].astype(str) == bat
            axes[2].scatter(
                pc_df.loc[mask, 'PC1'], pc_df.loc[mask, 'PC2'],
                c=[batch_color_map[bat]], label=bat,
                s=40, alpha=0.8, edgecolors='white', linewidths=0.5,
            )
        axes[2].set_xlabel(f'PC1 ({var_explained[0]:.1f}%)')
        axes[2].set_ylabel(f'PC2 ({var_explained[1]:.1f}%)')
        axes[2].set_title('PCA \u2014 Colored by Batch')
        axes[2].legend(fontsize=8, loc='upper right')
        axes[2].axhline(0, ls='--', lw=0.3, color='grey')
        axes[2].axvline(0, ls='--', lw=0.3, color='grey')

    fig.suptitle(
        f'Global PCA \u2014 All {len(log_sub)} Samples (log2-CPM, Top 5000 Variable Genes)',
        fontsize=13,
    )
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, 'pca_global.png'))
    plt.close(fig)


def _plot_per_tissue_pca(counts_df, out_dir):
    """Per-tissue PCA with mouse_id annotations."""
    for tissue in cfg.TISSUES:
        tissue_samples = [s for s in counts_df.columns if s.split('_', 1)[1] == tissue]
        if len(tissue_samples) < 4:
            continue

        sub = np.log2(counts_df[tissue_samples] + 1).T
        gene_var = sub.var()
        top_genes = gene_var.nlargest(min(3000, len(gene_var))).index
        sub = sub[top_genes]

        scaled = StandardScaler().fit_transform(sub)
        pca_t = PCA(n_components=min(3, len(sub)))
        pcs_t = pca_t.fit_transform(scaled)

        fig, ax = plt.subplots(figsize=(7, 5))
        for s_idx, sample in enumerate(tissue_samples):
            mouse_id = sample.split('_')[0]
            grp = cfg.MOUSE_GROUPS.get(mouse_id, 'Unknown')
            color = cfg.GROUP_COLORS.get(grp, '#999')
            ax.scatter(
                pcs_t[s_idx, 0], pcs_t[s_idx, 1],
                c=color, s=80, edgecolors='white', linewidths=0.8, zorder=3,
            )
            ax.annotate(
                mouse_id, (pcs_t[s_idx, 0], pcs_t[s_idx, 1]),
                fontsize=8, ha='left', va='bottom',
                xytext=(4, 4), textcoords='offset points',
            )

        var_exp = pca_t.explained_variance_ratio_ * 100
        ax.set_xlabel(f'PC1 ({var_exp[0]:.1f}%)')
        ax.set_ylabel(f'PC2 ({var_exp[1]:.1f}%)')
        ax.set_title(f'PCA \u2014 {cfg.TISSUE_LABELS.get(tissue, tissue)}')
        ax.axhline(0, ls='--', lw=0.3, color='grey')
        ax.axvline(0, ls='--', lw=0.3, color='grey')

        # Legend
        legend_elements = [
            Patch(facecolor=cfg.GROUP_COLORS[g_], label=g_)
            for g_ in sorted(cfg.GROUP_COLORS.keys())
        ]
        ax.legend(handles=legend_elements, fontsize=9)

        fig.tight_layout()
        fig.savefig(os.path.join(out_dir, f'pca_{tissue}.png'))
        plt.close(fig)
        cfg.flush_print(f"  PCA {tissue}: done")


def _plot_correlation_heatmap(counts_df, out_dir):
    """Sample-sample Pearson correlation heatmap with group color annotation."""
    cfg.flush_print("  Sample correlation heatmap...")
    log_all = np.log2(counts_df + 1)
    corr = log_all.corr(method='pearson')

    # Color annotation
    row_colors = []
    for s in corr.index:
        mouse_id = s.split('_')[0]
        grp = cfg.MOUSE_GROUPS.get(mouse_id, 'Unknown')
        row_colors.append(cfg.GROUP_COLORS.get(grp, '#999'))

    g = sns.clustermap(
        corr,
        cmap='YlOrRd',
        vmin=0.5, vmax=1.0,
        figsize=(16, 14),
        row_colors=row_colors,
        col_colors=row_colors,
        xticklabels=True,
        yticklabels=True,
        linewidths=0,
    )
    g.ax_heatmap.tick_params(labelsize=5)
    g.fig.suptitle('Sample-Sample Pearson Correlation', fontsize=14, y=1.01)
    g.savefig(os.path.join(out_dir, 'sample_correlation_heatmap.png'))
    plt.close(g.fig)


def _pc_covariate_association(counts_df, out_dir, sample_batches=None):
    """
    主成分与协变量的关联检验（批次效应的定量诊断）。

    对前若干个主成分，逐一检验其得分在各协变量水平间是否存在系统差异
    （Kruskal-Wallis 非参数检验，对小样本和非正态更稳健）。

    读法：
      * ``batch`` 在 PC1/PC2 上显著 → 批次效应主导全局变异，
        合并分析必须在模型中加入 batch 项。
      * ``batch`` 与 ``condition`` 同时在同一主成分上显著 → 两者可能混杂，
        需要回到设计层面确认是否可分离。

    输出 ``pc_covariate_association.csv`` 与 ``pc_covariate_association.png``。
    """
    cfg.flush_print("  PC ~ covariate association...")

    log_sub = _log_cpm_top_variable(counts_df, n_top=5000)
    if log_sub.shape[0] < 3:
        cfg.flush_print("    Skipped: fewer than 3 samples.")
        return None

    scaled = StandardScaler().fit_transform(log_sub)
    n_pcs = int(min(_N_PCS_FOR_ASSOCIATION, log_sub.shape[0] - 1, log_sub.shape[1]))
    if n_pcs < 1:
        cfg.flush_print("    Skipped: not enough components.")
        return None

    pca = PCA(n_components=n_pcs)
    pcs = pca.fit_transform(scaled)
    var_explained = pca.explained_variance_ratio_ * 100

    ann = _build_sample_annotation(log_sub.index, sample_batches or {})

    rows = []
    for cov in _COVARIATES:
        if cov not in ann.columns:
            continue
        labels = ann[cov]
        if labels.isna().all():
            continue
        for pc_idx in range(n_pcs):
            scores = pd.Series(pcs[:, pc_idx], index=log_sub.index)
            groups = [
                scores[labels == lvl].to_numpy()
                for lvl in labels.dropna().unique()
            ]
            groups = [g for g in groups if len(g) > 0]
            # Kruskal-Wallis 至少需要 2 组，且不能所有组都只有 1 个样本
            if len(groups) < 2 or sum(len(g) for g in groups) < 3:
                pvalue = np.nan
            else:
                try:
                    _, pvalue = stats.kruskal(*groups)
                except ValueError:
                    # 所有取值完全相同时 scipy 会抛错
                    pvalue = np.nan
            rows.append({
                'covariate': cov,
                'PC': f'PC{pc_idx + 1}',
                'variance_explained_pct': round(float(var_explained[pc_idx]), 3),
                'n_levels': len(groups),
                'pvalue': pvalue,
                'neg_log10_p': (
                    float(-np.log10(pvalue)) if pvalue and np.isfinite(pvalue) and pvalue > 0
                    else np.nan
                ),
                'test': 'kruskal_wallis',
            })

    if not rows:
        cfg.flush_print("    Skipped: no usable covariates.")
        return None

    assoc = pd.DataFrame(rows)
    csv_path = os.path.join(out_dir, 'pc_covariate_association.csv')
    assoc.to_csv(csv_path, index=False)
    cfg.flush_print(f"    Saved: {csv_path}")

    # 热图：行=协变量，列=主成分，值=-log10(p)
    pivot = assoc.pivot(index='covariate', columns='PC', values='neg_log10_p')
    pivot = pivot.reindex(columns=[f'PC{i + 1}' for i in range(n_pcs)])

    fig, ax = plt.subplots(figsize=(1.6 * n_pcs + 3, 1.1 * len(pivot) + 2.2))
    sns.heatmap(
        pivot, annot=True, fmt='.2f', cmap='YlOrRd', vmin=0,
        cbar_kws={'label': '-log10(p)'}, linewidths=0.5, linecolor='white', ax=ax,
    )
    # p=0.05 对应 -log10(p)=1.30
    ax.set_title(
        'PC ~ Covariate Association (Kruskal-Wallis)\n'
        'values > 1.30 correspond to p < 0.05',
        fontsize=11,
    )
    ax.set_xlabel('')
    ax.set_ylabel('')
    xlabels = [
        f'{pc}\n({assoc.loc[assoc["PC"] == pc, "variance_explained_pct"].iloc[0]:.1f}%)'
        for pc in pivot.columns
    ]
    ax.set_xticklabels(xlabels, rotation=0)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, 'pc_covariate_association.png'))
    plt.close(fig)

    # 控制台提示最需要注意的信号
    batch_hits = assoc[(assoc['covariate'] == 'batch') & (assoc['pvalue'] < 0.05)]
    if not batch_hits.empty:
        pcs_hit = ', '.join(
            f"{r.PC}(p={r.pvalue:.3g}, {r.variance_explained_pct:.1f}% var)"
            for r in batch_hits.itertuples()
        )
        cfg.flush_print(
            f"    [WARN] 批次在以下主成分上显著关联: {pcs_hit}\n"
            f"           合并分析需在 DE 模型中加入 batch 项；"
            f"ComBat 校正矩阵只可用于可视化，不可作为 DESeq2 输入。"
        )

    return assoc
