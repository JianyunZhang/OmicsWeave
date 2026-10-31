"""
Bulk RNA-Seq Analysis Pipeline (Windows Native)
=================================================
Complete pipeline using:
  Step 1: FASTQ Quality Control (Pure Python)
  Step 2: Adapter Trimming (cutadapt)
  Step 3: Pseudo-alignment & Quantification (kallisto - Windows native)
  Step 4: Differential Expression Analysis (PyDESeq2)

Environment: conda activate tutorial-mfa-2025 (Python 3.12)
Author: Auto-generated for jianyunzhang's obesity research project
Date: 2026-02-07
"""

import os
import sys
import gzip
import json
import subprocess
import collections
import argparse
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd

from .experiment_config import (
    load_transcriptomics_config,
    config_arg_help,
    resolve_comparisons,
    discover_experiment_protocols,
)
from .run_wilcoxon_de import run_wilcoxon_analysis
from .gene_annotate import annotate_csv_file
from .transcriptomics_integration import run_integration as run_integration_analysis
from .normalization import (
    build_gene_level_counts,
    AGGREGATION_LENGTH_SCALED_TPM,
    AGGREGATION_SIMPLE_SUM,
)
from . import qc_gate
from . import provenance as provenance_mod
from omicsweave.core.project import PROJECT_DIR_ENV, project_dir

# ============================================================================
# CONFIGURATION — all values are loaded from experiment YAML config.
# No hardcoded experiment-specific defaults.
# ============================================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Placeholders — will be overwritten by _apply_config()
RAW_DATA_ROOT = ""
SAMPLE_LIST_PATH = ""
OUTPUT_BASE = ""
KALLISTO_EXE = ""
TRANSCRIPTOME_FASTA = ""
KALLISTO_INDEX = ""
QC_SAMPLE_SIZE = None
THREADS = 8
MOUSE_GROUPS: dict[str, str] = {}
SAMPLE_BATCHES: dict[str, str] = {}
TISSUES: list[str] = []
COMPARISONS: list[tuple[str, str, str]] = []
PROTOCOL_ID: str | None = None
BATCH_ID: str | None = None
QC_THRESHOLDS: dict = {}
ON_QC_FAIL: str = "flag"
TECHNICAL_PARAMS: dict = {}
STRANDEDNESS: str = "unstranded"
GENE_AGGREGATION: str = "lengthScaledTPM"


def _apply_config(experiment: "str | None", config_path: "str | None", outdir: "str | None") -> dict:
    """
    Apply experiment YAML config to module-level globals.
    Returns a small dict of resolved config fields for logging.
    """
    cfg = load_transcriptomics_config(
        experiment_id=experiment,
        config_path=config_path,
        outdir_override=outdir,
    )

    global RAW_DATA_ROOT, SAMPLE_LIST_PATH, OUTPUT_BASE
    global KALLISTO_EXE, TRANSCRIPTOME_FASTA, KALLISTO_INDEX
    global QC_SAMPLE_SIZE, THREADS, TISSUES, MOUSE_GROUPS, SAMPLE_BATCHES, COMPARISONS
    global PROTOCOL_ID, BATCH_ID
    global QC_THRESHOLDS, ON_QC_FAIL, TECHNICAL_PARAMS, STRANDEDNESS, GENE_AGGREGATION

    RAW_DATA_ROOT = cfg.raw_data_root
    SAMPLE_LIST_PATH = cfg.sample_list_tsv
    OUTPUT_BASE = os.path.abspath(str(cfg.output_dir))
    os.environ["TRANSCRIPTOMICS_OUTDIR"] = OUTPUT_BASE

    KALLISTO_EXE = os.path.abspath(str(cfg.kallisto_exe))
    TRANSCRIPTOME_FASTA = os.path.abspath(str(cfg.transcriptome_fasta))
    KALLISTO_INDEX = os.path.abspath(str(cfg.kallisto_index))

    QC_SAMPLE_SIZE = cfg.qc_sample_size
    THREADS = int(cfg.threads)
    TISSUES = list(cfg.tissues)
    MOUSE_GROUPS = dict(cfg.mouse_groups)
    SAMPLE_BATCHES = dict(cfg.sample_batches)
    COMPARISONS = resolve_comparisons(cfg)
    PROTOCOL_ID = cfg.protocol_id
    BATCH_ID = cfg.batch_id
    QC_THRESHOLDS = dict(cfg.qc_thresholds or {})
    ON_QC_FAIL = str(cfg.on_qc_fail or "flag")
    TECHNICAL_PARAMS = dict(cfg.technical_params or {})
    STRANDEDNESS = str(cfg.strandedness or "unstranded")
    GENE_AGGREGATION = str(cfg.gene_aggregation or AGGREGATION_LENGTH_SCALED_TPM)

    return {
        "experiment_id": cfg.experiment_id,
        "config_path": str(cfg.config_path),
        "groups_path": str(cfg.groups_path),
        "protocol_id": cfg.protocol_id,
        "batch_id": cfg.batch_id,
        "output_root": str(cfg.output_root),
        "output_dir": OUTPUT_BASE,
        "strandedness": STRANDEDNESS,
        "gene_aggregation": GENE_AGGREGATION,
        "qc_thresholds": QC_THRESHOLDS,
        "on_qc_fail": ON_QC_FAIL,
    }


def flush_print(*args, **kwargs):
    """Print with immediate flush for real-time output over network."""
    print(*args, **kwargs, flush=True)


# ============================================================================
# STEP 1: PURE PYTHON FASTQ QUALITY CONTROL
# ============================================================================
def parse_fastq_gz(filepath, max_reads=None):
    """Parse a gzipped FASTQ file and yield (name, sequence, quality) tuples."""
    count = 0
    with gzip.open(filepath, 'rt', encoding='ascii') as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            _ = f.readline()  # +
            qual = f.readline().strip()
            yield header, seq, qual
            count += 1
            if max_reads and count >= max_reads:
                break


def phred_to_quality(qual_string):
    """Convert Phred+33 quality string to quality scores."""
    return [ord(c) - 33 for c in qual_string]


def fastq_qc(filepath, max_reads=None):
    """Perform quality control on a FASTQ.GZ file. Returns dict of QC metrics."""
    flush_print(f"    QC: {os.path.basename(filepath)}...", end=" ")

    total_reads = 0
    total_bases = 0
    gc_count = 0
    read_lengths = []
    quality_per_position = collections.defaultdict(list)
    quality_scores_all = []
    n_count = 0

    for header, seq, qual in parse_fastq_gz(filepath, max_reads):
        total_reads += 1
        seq_len = len(seq)
        total_bases += seq_len
        read_lengths.append(seq_len)
        gc_count += seq.count('G') + seq.count('C')
        n_count += seq.count('N')
        scores = phred_to_quality(qual)
        quality_scores_all.extend(scores)
        for i, s in enumerate(scores):
            quality_per_position[i].append(s)

        if total_reads % 50000 == 0:
            flush_print(f"{total_reads}...", end=" ")

    read_lengths = np.array(read_lengths)
    quality_scores_all = np.array(quality_scores_all)

    max_pos = max(quality_per_position.keys()) + 1
    mean_qual_per_pos = []
    for i in range(max_pos):
        if i in quality_per_position:
            mean_qual_per_pos.append(float(np.mean(quality_per_position[i])))
        else:
            mean_qual_per_pos.append(0)

    qc_result = {
        'file': os.path.basename(filepath),
        'total_reads': total_reads,
        'total_bases': total_bases,
        'mean_read_length': float(np.mean(read_lengths)),
        'min_read_length': int(np.min(read_lengths)),
        'max_read_length': int(np.max(read_lengths)),
        'gc_content_pct': round(gc_count / total_bases * 100, 2) if total_bases > 0 else 0,
        'n_content_pct': round(n_count / total_bases * 100, 4) if total_bases > 0 else 0,
        'mean_quality': round(float(np.mean(quality_scores_all)), 2),
        'q20_pct': round(float(np.sum(quality_scores_all >= 20) / len(quality_scores_all) * 100), 2),
        'q30_pct': round(float(np.sum(quality_scores_all >= 30) / len(quality_scores_all) * 100), 2),
        'mean_quality_per_position': [round(x, 2) for x in mean_qual_per_pos],
    }

    flush_print(f"done ({total_reads} reads, Q={qc_result['mean_quality']})")
    return qc_result


def run_qc_for_sample(sample_id, r1_path, r2_path, output_dir, max_reads=None):
    """Run QC for a pair of FASTQ files and save results."""
    os.makedirs(output_dir, exist_ok=True)
    report_path = os.path.join(output_dir, f"{sample_id}_qc_report.json")

    # Skip if report already exists
    if os.path.exists(report_path):
        flush_print(f"  [{sample_id}] QC report exists, loading...")
        with open(report_path) as f:
            return json.load(f)

    flush_print(f"  [{sample_id}]")
    results = {}
    for label, fpath in [('R1', r1_path), ('R2', r2_path)]:
        qc = fastq_qc(fpath, max_reads=max_reads)
        results[label] = qc

    with open(report_path, 'w') as f:
        json.dump(results, f, indent=2)

    return results


# ============================================================================
# STEP 2: TRIMMING (cutadapt)
# ============================================================================
def run_trimming(sample_id, r1_path, r2_path, output_dir):
    """Trim adapters and low-quality bases using cutadapt."""
    os.makedirs(output_dir, exist_ok=True)

    r1_out = os.path.join(output_dir, f"{sample_id}_R1.trimmed.fq.gz")
    r2_out = os.path.join(output_dir, f"{sample_id}_R2.trimmed.fq.gz")
    report = os.path.join(output_dir, f"{sample_id}_cutadapt.json")

    if os.path.exists(r1_out) and os.path.exists(r2_out):
        flush_print(f"  [{sample_id}] Trimming: skipped (output exists)")
        return r1_out, r2_out

    cmd = [
        sys.executable, "-m", "cutadapt",
        "-j", str(THREADS),
        "-a", "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",  # Illumina TruSeq R1
        "-A", "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",  # Illumina TruSeq R2
        "-q", "20",    # Quality cutoff
        "-m", "20",    # Minimum length
        "-o", r1_out,
        "-p", r2_out,
        "--json", report,
        r1_path, r2_path
    ]

    flush_print(f"  [{sample_id}] Trimming...")
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        flush_print(f"  [{sample_id}] Trimming: done")
    except FileNotFoundError:
        flush_print(f"  [{sample_id}] ERROR: cutadapt not found (pip install cutadapt)")
        return r1_path, r2_path
    except subprocess.CalledProcessError as e:
        flush_print(f"  [{sample_id}] ERROR: {e.stderr[:200] if e.stderr else 'Unknown'}")
        return r1_path, r2_path

    return r1_out, r2_out


# ============================================================================
# STEP 3: KALLISTO QUANTIFICATION
# ============================================================================
def build_kallisto_index():
    """Build kallisto index from transcriptome FASTA."""
    if os.path.exists(KALLISTO_INDEX):
        flush_print(f"  kallisto index exists: {KALLISTO_INDEX}")
        return True

    os.makedirs(os.path.dirname(KALLISTO_INDEX), exist_ok=True)

    if not os.path.exists(TRANSCRIPTOME_FASTA):
        flush_print(f"  ERROR: Transcriptome not found: {TRANSCRIPTOME_FASTA}")
        flush_print(f"  Run download_transcriptome.py first.")
        return False

    flush_print(f"  Building kallisto index (may take a few minutes)...")
    cmd = [
        KALLISTO_EXE, "index",
        "-i", KALLISTO_INDEX,
        TRANSCRIPTOME_FASTA
    ]

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        flush_print(f"  Index built: {KALLISTO_INDEX}")
        return True
    except subprocess.CalledProcessError as e:
        flush_print(f"  ERROR building index: {e.stderr[:300] if e.stderr else ''}")
        return False


STRANDEDNESS_FLAGS = {
    "unstranded": [],
    "fr-stranded": ["--fr-stranded"],
    "rf-stranded": ["--rf-stranded"],
}


def resolve_strandedness_flags(strandedness: str) -> list[str]:
    """
    把配置里的链特异性设置翻译成 kallisto 命令行开关。

    默认 ``unstranded`` 保持 kallisto 原生行为。若文库实为 dUTP 反向链
    特异（目前主流 mRNA 建库方式），按非链特异处理会让反义链与重叠基因
    的计数被错误归属。
    """
    key = str(strandedness or "unstranded").strip().lower()
    if key not in STRANDEDNESS_FLAGS:
        flush_print(f"  [WARN] Unknown strandedness {strandedness!r}; "
                    f"expected one of {sorted(STRANDEDNESS_FLAGS)}. Using 'unstranded'.")
        return []
    return list(STRANDEDNESS_FLAGS[key])


def run_kallisto_quant(sample_id, r1_path, r2_path, output_dir, strandedness="unstranded"):
    """Run kallisto quantification for a sample."""
    quant_dir = os.path.join(output_dir, sample_id)
    abundance_file = os.path.join(quant_dir, "abundance.tsv")

    if os.path.exists(abundance_file):
        flush_print(f"  [{sample_id}] Quant: skipped (output exists)")
        return abundance_file

    os.makedirs(quant_dir, exist_ok=True)

    cmd = [
        KALLISTO_EXE, "quant",
        "-i", KALLISTO_INDEX,
        "-o", quant_dir,
        "-t", str(THREADS),
        *resolve_strandedness_flags(strandedness),
        r1_path, r2_path
    ]

    flush_print(f"  [{sample_id}] Quantifying...")
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        # Parse mapping rate from stderr
        for line in result.stderr.split('\n'):
            if 'pseudoaligned' in line.lower():
                flush_print(f"  [{sample_id}] {line.strip()}")
                break
        flush_print(f"  [{sample_id}] Quant: done")
        return abundance_file
    except subprocess.CalledProcessError as e:
        flush_print(f"  [{sample_id}] ERROR: {e.stderr[:300] if e.stderr else ''}")
        return None


# ============================================================================
# STEP 4: MERGE COUNTS & DIFFERENTIAL EXPRESSION (PyDESeq2)
# ============================================================================
def merge_kallisto_counts(quant_dir, sample_ids):
    """
    Merge kallisto abundance.tsv files into count, TPM and effective-length matrices.

    Returns
    -------
    (counts_df, tpm_df, efflen_df)
        Transcripts as rows, samples as columns.  ``counts_df`` keeps kallisto's
        fractional ``est_counts`` (rounding is deferred to gene-level aggregation,
        which needs the un-rounded values for length scaling).
    """
    flush_print("\nMerging kallisto quantification results...")
    counts_dict = {}
    tpm_dict = {}
    efflen_dict = {}

    for sid in sample_ids:
        af = os.path.join(quant_dir, sid, "abundance.tsv")
        if os.path.exists(af):
            df = pd.read_csv(af, sep='\t').set_index('target_id')
            counts_dict[sid] = df['est_counts']
            tpm_dict[sid] = df['tpm']
            if 'eff_length' in df.columns:
                efflen_dict[sid] = df['eff_length']

    if not counts_dict:
        flush_print("  ERROR: No abundance files found.")
        return None, None, None

    counts_df = pd.DataFrame(counts_dict).fillna(0)
    tpm_df = pd.DataFrame(tpm_dict).fillna(0)
    efflen_df = pd.DataFrame(efflen_dict).fillna(0) if efflen_dict else None

    flush_print(f"  Merged {len(counts_dict)} samples, {len(counts_df)} transcripts.")
    if efflen_df is None:
        flush_print("  [WARN] abundance.tsv has no eff_length column; "
                    "length-scaled aggregation will be unavailable.")
    return counts_df, tpm_df, efflen_df


def build_tx2gene_from_fasta(fasta_path):
    """
    Build transcript-to-gene mapping from Ensembl cDNA FASTA headers.
    Headers look like: >ENSMUST00000103301.3 cdna ... gene:ENSMUSG00000076500.3 ...
    Returns dict {transcript_id: gene_id}
    """
    flush_print("  Building transcript-to-gene mapping from FASTA...")
    tx2gene = {}
    
    open_func = gzip.open if fasta_path.endswith('.gz') else open
    with open_func(fasta_path, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                parts = line.strip().split()
                tx_id = parts[0][1:]  # Remove >
                gene_id = None
                gene_symbol = None
                for part in parts:
                    if part.startswith('gene:'):
                        gene_id = part.split(':')[1]
                    if part.startswith('gene_symbol:'):
                        gene_symbol = part.split(':')[1]
                if gene_id:
                    tx2gene[tx_id] = {
                        'gene_id': gene_id,
                        'gene_symbol': gene_symbol or gene_id
                    }
    
    flush_print(f"  Mapped {len(tx2gene)} transcripts to genes.")
    return tx2gene


def aggregate_counts_to_gene(counts_df, tx2gene, tpm_df=None, efflen_df=None,
                             method=AGGREGATION_LENGTH_SCALED_TPM):
    """
    Aggregate transcript-level counts to gene level.

    Defaults to tximport's ``lengthScaledTPM``, which folds transcript effective
    length into the gene-level counts.  Plain summation of ``est_counts`` is
    biased whenever isoform usage shifts between conditions: the gene total
    changes while true expression does not.

    Falls back to ``simple_sum`` (with a warning) when TPM/effective-length
    matrices are unavailable, so the pipeline still runs on legacy outputs.
    """
    flush_print("  Aggregating transcript counts to gene level...")

    resolved = method
    if resolved == AGGREGATION_LENGTH_SCALED_TPM and (tpm_df is None or efflen_df is None):
        flush_print("  [WARN] TPM/eff_length unavailable — falling back to simple_sum "
                    "(gene counts will NOT be length-corrected).")
        resolved = AGGREGATION_SIMPLE_SUM

    if resolved == AGGREGATION_SIMPLE_SUM:
        gene_counts, info = build_gene_level_counts(
            counts_df, counts_df, counts_df, tx2gene, method=AGGREGATION_SIMPLE_SUM,
        )
    else:
        gene_counts, info = build_gene_level_counts(
            counts_df, tpm_df, efflen_df, tx2gene, method=AGGREGATION_LENGTH_SCALED_TPM,
        )

    flush_print(f"  Gene-level counts: {info['n_genes']} genes "
                f"(method={info['gene_aggregation']}, "
                f"length_corrected={info['length_corrected']}).")
    return gene_counts, info


def run_pydeseq2_analysis(counts_df, metadata_df, design_formula='condition',
                          contrast=None, output_prefix="de", shrink_lfc=True):
    """
    Run differential expression analysis using PyDESeq2.

    Args:
        counts_df: Gene count matrix (genes x samples)
        metadata_df: Sample metadata with 'condition' column
        design_formula: Design formula column name
        contrast: Tuple of (factor, numerator, denominator) for comparison
        output_prefix: Prefix for output files
        shrink_lfc: Apply apeGLM LFC shrinkage (recommended; p-values unchanged)

    Returns:
        DataFrame with DE results including test_method, n_case, n_control,
        evidence_level columns.
    """
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    flush_print(f"\n  Running PyDESeq2 analysis...")
    flush_print(f"  Design: ~{design_formula}")
    if contrast:
        flush_print(f"  Contrast: {contrast[1]} vs {contrast[2]}")

    # Ensure sample order matches
    common_samples = list(set(counts_df.columns) & set(metadata_df.index))
    if len(common_samples) == 0:
        flush_print("  ERROR: No common samples between counts and metadata.")
        return None

    counts_sub = counts_df[common_samples].T  # PyDESeq2 expects samples x genes
    metadata_sub = metadata_df.loc[common_samples]

    design_terms = _validate_de_design(metadata_sub, design_formula)

    flush_print(f"  Samples: {len(common_samples)}, Genes: {counts_sub.shape[1]}")

    # Filter low-count genes (at least 10 total counts)
    gene_sums = counts_sub.sum(axis=0)
    keep = gene_sums >= 10
    counts_filtered = counts_sub.loc[:, keep]
    flush_print(f"  After filtering (sum >= 10): {counts_filtered.shape[1]} genes")

    # Run DESeq2
    dds = DeseqDataSet(
        counts=counts_filtered,
        metadata=metadata_sub,
        design=f"~{design_formula}",
    )

    flush_print("  Fitting DESeq2 model...")
    dds.deseq2()

    # Get results
    if contrast:
        stat_res = DeseqStats(dds, contrast=list(contrast))
    else:
        stat_res = DeseqStats(dds)

    stat_res.summary()
    results_df = stat_res.results_df

    # LFC 收缩（apeGLM）——不改变 p 值，只稳定效应量估计
    shrink_info = {"lfc_shrinkage": "none", "lfc_shrink_coeff": "", "lfc_prior": ""}
    if shrink_lfc and contrast:
        unshrunk_lfc = results_df['log2FoldChange'].copy()
        results_df, shrink_info = _apply_lfc_shrinkage(stat_res, contrast, unshrunk_lfc)
        results_df['log2FoldChange_unshrunk'] = unshrunk_lfc.reindex(results_df.index)

    # Add significance labels
    results_df['significant'] = (
        (results_df['padj'] < 0.05) &
        (results_df['log2FoldChange'].abs() > 1)
    )

    # Add metadata columns
    numerator_name = contrast[1] if contrast else "unknown"
    denominator_name = contrast[2] if contrast else "unknown"
    n_num = int((metadata_sub['condition'] == numerator_name).sum())
    n_den = int((metadata_sub['condition'] == denominator_name).sum())

    results_df['test_method'] = 'pydeseq2'
    results_df['n_case'] = n_num
    results_df['n_control'] = n_den
    results_df['evidence_level'] = 'statistical_test'
    results_df['design_formula'] = f"~{design_formula}"
    results_df['batch_adjusted'] = 'batch' in design_terms
    results_df['n_batches'] = (
        int(metadata_sub['batch'].nunique()) if 'batch' in design_terms else 1
    )
    results_df['batch_ids'] = (
        ";".join(sorted(metadata_sub['batch'].astype(str).unique()))
        if 'batch' in design_terms else ""
    )
    results_df['source_mode'] = (
        'joint_count_model_batch_adjusted'
        if 'batch' in design_terms else 'single_results_directory'
    )
    results_df['lfc_shrinkage'] = shrink_info.get('lfc_shrinkage', 'none')
    results_df['lfc_shrink_coeff'] = shrink_info.get('lfc_shrink_coeff', '')

    sig = results_df[results_df['significant']]
    up = sig[sig['log2FoldChange'] > 0]
    down = sig[sig['log2FoldChange'] < 0]

    flush_print(f"  Results: {len(results_df)} genes tested")
    flush_print(f"  Significant DEGs (padj<0.05, |log2FC|>1): {len(sig)}")
    flush_print(f"    Up-regulated: {len(up)}")
    flush_print(f"    Down-regulated: {len(down)}")

    return results_df


def _resolve_shrinkage_coeff(lfc_columns, factor, numerator, denominator):
    """
    找出与 contrast 对应的 LFC 系数名，并判断方向是否需要取反。

    PyDESeq2 0.5.4 的 ``lfc_shrink()`` 只能收缩**设计矩阵里实际存在的
    系数**。参考水平由 formulaic 按字母序决定，``ref_level`` 参数已废弃
    且无效。因此当 contrast 的分母不是参考水平时，可用系数
    ``condition[T.<分母>]`` 表示的方向与 contrast **恰好相反**。

    实测确认：contrast=(hfd, lfd) 而系数为 ``condition[T.lfd_14w]`` 时，
    收缩后 log2FoldChange 符号 100% 翻转。若不处理，所有上/下调结论
    都会被静默颠倒。

    Returns
    -------
    (coeff_name, needs_negation) 或 (None, False) —— 找不到匹配系数时。
    """
    prefix = f"{factor}[T."
    candidates = [c for c in lfc_columns if str(c).startswith(prefix) and str(c).endswith("]")]
    for coeff in candidates:
        level = str(coeff)[len(prefix):-1]
        if level == numerator:
            return coeff, False
        if level == denominator:
            return coeff, True
    return None, False


def _apply_lfc_shrinkage(stat_res, contrast, unshrunk_lfc):
    """
    对 DE 结果应用 apeGLM LFC 收缩，并校验方向未被翻转。

    未收缩的 log2FC 对低表达基因会剧烈膨胀，而 ``|log2FC| > 1`` 的筛选
    阈值恰好优先选中这些噪声基因。收缩不改变 p 值，只改变效应量估计
    与其标准误（变为后验标准差）。

    Returns
    -------
    (results_df, shrinkage_info)
    """
    factor, numerator, denominator = contrast
    info = {"lfc_shrinkage": "none", "lfc_shrink_coeff": "", "lfc_prior": ""}

    coeff, needs_negation = _resolve_shrinkage_coeff(
        list(stat_res.LFC.columns), factor, numerator, denominator
    )
    if coeff is None:
        flush_print("    [WARN] No matching LFC coefficient found; skipping shrinkage.")
        info["lfc_shrinkage"] = "skipped_no_coeff"
        return stat_res.results_df, info

    try:
        stat_res.lfc_shrink(coeff=coeff)
    except Exception as exc:  # noqa: BLE001 — 收缩失败不应中断整个比较
        flush_print(f"    [WARN] lfc_shrink failed ({type(exc).__name__}: {exc}); "
                    "keeping unshrunken LFC.")
        info["lfc_shrinkage"] = "failed"
        return stat_res.results_df, info

    results_df = stat_res.results_df
    if needs_negation:
        results_df["log2FoldChange"] = -results_df["log2FoldChange"]

    # 方向校验：与收缩前比较，强效应基因的符号必须保持一致。
    # 这道断言让 PyDESeq2 未来改变参考水平约定时立刻暴露，而不是静默出错。
    aligned = unshrunk_lfc.reindex(results_df.index)
    strong = aligned.abs() > 0.5
    if strong.sum() >= 10:
        agreement = float(
            (np.sign(aligned[strong]) == np.sign(results_df.loc[strong, "log2FoldChange"])).mean()
        )
        if agreement < 0.8:
            flush_print(
                f"    [ERROR] LFC shrinkage direction check failed "
                f"(sign agreement {agreement:.2f} < 0.80). Reverting to unshrunken LFC."
            )
            results_df["log2FoldChange"] = aligned
            info["lfc_shrinkage"] = "reverted_direction_check"
            return results_df, info
        info["lfc_direction_agreement"] = round(agreement, 4)

    flush_print(f"    LFC shrinkage: apeGLM on '{coeff}'"
                f"{' (negated to match contrast)' if needs_negation else ''}")
    info["lfc_shrinkage"] = "apeglm"
    info["lfc_shrink_coeff"] = str(coeff)
    info["lfc_prior"] = "apeglm_cauchy"
    return results_df, info


def _validate_de_design(metadata_df: pd.DataFrame, design_formula: str) -> list[str]:
    """验证 DE 设计列完整且满秩，返回解析后的设计项。"""
    design_terms = [term.strip() for term in str(design_formula).split("+") if term.strip()]
    if not design_terms:
        raise ValueError("DE design requires at least one metadata term")
    missing_design = [term for term in design_terms if term not in metadata_df.columns]
    if missing_design:
        raise ValueError("Metadata missing design columns: " + ", ".join(missing_design))
    missing_values = metadata_df[design_terms].isna().any(axis=1)
    if missing_values.any():
        raise ValueError(
            "DE design metadata contains missing values for samples: "
            + ", ".join(metadata_df.index[missing_values].astype(str))
        )

    design_parts = [pd.DataFrame({"intercept": np.ones(len(metadata_df))}, index=metadata_df.index)]
    for term in design_terms:
        encoded = pd.get_dummies(
            metadata_df[term].astype(str), prefix=term, drop_first=True, dtype=float
        )
        design_parts.append(encoded)
    audit_design = pd.concat(design_parts, axis=1)
    if np.linalg.matrix_rank(audit_design.to_numpy(dtype=float)) < audit_design.shape[1]:
        raise ValueError(
            f"DE design '~{design_formula}' is rank deficient/confounded; "
            "use within-batch effects and Meta-analysis"
        )
    return design_terms


def export_multiomics_csv(de_df, output_dir, comp_name):
    """
    Export a standardized DE result file for multi-omics integration.

    Creates ``de_for_multiomics_<comp>.csv`` with fixed column names so
    downstream consumers (FluxOmicsScorer, TranscriptomicsLoader) never
    need column-name adapters.

    Columns: gene_symbol, ensembl_gene_id, ncbi_gene_id, log2fc, lfcSE,
             pvalue, padj, significant, direction, test_method, evidence_level,
             n_case, n_control, n_eff, n_eff_source, design_formula,
             batch_adjusted, n_batches, batch_ids, source_mode

    Parameters
    ----------
    de_df : pd.DataFrame
        DE result (index = gene symbol, or has 'gene' column).
    output_dir : str
        Directory to write the file.
    comp_name : str
        Comparison name (e.g. 'hfd_14w_vs_lfd_14w').
    """
    from .gene_annotate import annotate_gene_dataframe

    df = de_df.copy()
    # Normalize gene column
    if 'gene' in df.columns:
        df['gene_symbol'] = df['gene'].astype(str)
    else:
        df['gene_symbol'] = df.index.astype(str)

    # Normalize log2fc column
    if 'log2fc' not in df.columns and 'log2FoldChange' in df.columns:
        df['log2fc'] = df['log2FoldChange']

    # Direction column
    def _direction(row):
        if row.get('significant', False):
            fc = row.get('log2fc', 0)
            if pd.notna(fc) and fc > 0:
                return 'up'
            if pd.notna(fc) and fc < 0:
                return 'down'
            return 'ns'
        return 'ns'
    df['direction'] = df.apply(_direction, axis=1)

    # Try to add ID annotations
    try:
        annotated = annotate_gene_dataframe(df, gene_col='gene_symbol')
        if 'ensembl_gene_id' in annotated.columns:
            df['ensembl_gene_id'] = annotated['ensembl_gene_id'].values
        else:
            df['ensembl_gene_id'] = None
        if 'ncbi_gene_id' in annotated.columns:
            df['ncbi_gene_id'] = annotated['ncbi_gene_id'].values
        else:
            df['ncbi_gene_id'] = None
    except Exception:
        df['ensembl_gene_id'] = None
        df['ncbi_gene_id'] = None

    # Ensure all metadata columns exist
    for col, default in [
        ('test_method', 'unknown'),
        ('evidence_level', 'unknown'),
        ('n_case', np.nan),
        ('n_control', np.nan),
        ('design_formula', '~condition'),
        ('batch_adjusted', False),
        ('n_batches', 1),
        ('batch_ids', ''),
        ('source_mode', 'single_results_directory'),
    ]:
        if col not in df.columns:
            df[col] = default

    if 'n_eff' not in df.columns:
        df['n_eff'] = (
            pd.to_numeric(df['n_case'], errors='coerce')
            + pd.to_numeric(df['n_control'], errors='coerce')
        )
    if 'n_eff_source' not in df.columns:
        df['n_eff_source'] = np.where(
            pd.to_numeric(df['n_eff'], errors='coerce').notna(),
            'de_design_metadata',
            'not_available',
        )

    # Fixed output columns
    out_cols = [
        'gene_symbol', 'ensembl_gene_id', 'ncbi_gene_id',
        'log2fc', 'lfcSE', 'pvalue', 'padj', 'significant',
        'direction', 'test_method', 'evidence_level',
        'n_case', 'n_control', 'n_eff', 'n_eff_source', 'design_formula',
        'batch_adjusted', 'n_batches', 'batch_ids', 'source_mode',
    ]
    export_df = df[[c for c in out_cols if c in df.columns]].copy()

    out_path = os.path.join(output_dir, f"de_for_multiomics_{comp_name}.csv")
    export_df.to_csv(out_path, index=False)
    flush_print(f"  Multi-omics export: {out_path}")
    return out_path


# ============================================================================
# MAIN PIPELINE
# ============================================================================
def check_environment():
    """Check Python version and available tools."""
    flush_print("=" * 60)
    flush_print("Environment Check")
    flush_print("=" * 60)
    flush_print(f"Python: {sys.version}")
    flush_print(f"Platform: {sys.platform}")

    # Check kallisto
    try:
        result = subprocess.run([KALLISTO_EXE, "version"],
                                capture_output=True, text=True)
        flush_print(f"kallisto: {result.stdout.strip()}")
    except FileNotFoundError:
        flush_print(f"kallisto: NOT FOUND at {KALLISTO_EXE}")
        flush_print("  Run setup_kallisto.py first.")

    # Check cutadapt
    try:
        result = subprocess.run([sys.executable, "-m", "cutadapt", "--version"],
                                capture_output=True, text=True)
        flush_print(f"cutadapt: {result.stdout.strip()}")
    except Exception:
        flush_print("cutadapt: NOT FOUND (pip install cutadapt)")

    # Check pydeseq2
    try:
        import pydeseq2
        flush_print(f"pydeseq2: available")
    except ImportError:
        flush_print("pydeseq2: NOT FOUND (pip install pydeseq2)")

    # Check scipy/statsmodels
    try:
        import scipy
        flush_print(f"scipy: {scipy.__version__}")
    except ImportError:
        flush_print("scipy: NOT FOUND")

    try:
        import statsmodels
        flush_print(f"statsmodels: {statsmodels.__version__}")
    except ImportError:
        flush_print("statsmodels: NOT FOUND")

    flush_print("=" * 60)


def run_all_batches(experiment_id: str, outdir: str | None = None) -> list[tuple[str | None, str | None, str]]:
    """
    Run all enabled protocol/batch entries discovered from experiment.yaml.

    Each batch is executed as an independent subprocess in single-batch mode.
    """
    protocols = discover_experiment_protocols(experiment_id)
    results: list[tuple[str | None, str | None, str]] = []

    n_run = 0
    n_skip = 0
    for proto in protocols:
        protocol_id = proto.get("protocol_id")
        if not proto.get("enabled", False):
            for batch in proto.get("batches", []):
                results.append((protocol_id, batch.get("batch_id"), "disabled_protocol"))
                n_skip += 1
            continue

        for batch in proto.get("batches", []):
            batch_id = batch.get("batch_id")
            config_path = batch.get("config_path")
            has_data = bool(batch.get("has_data", False))

            if not has_data or not config_path:
                results.append((protocol_id, batch_id, "no_data"))
                n_skip += 1
                flush_print(f"  [skip] {protocol_id}/{batch_id}: 无可用配置或数据")
                continue

            cmd = [sys.executable, "-m", "omicsweave.transcriptomics.rnaseq_pipeline",
                   "--config", str(config_path)]
            if outdir:
                cmd.extend(["--outdir", outdir])
            # 子进程不继承上下文中显式指定的项目目录（--project-dir），经环境变量传递
            child_env = {**os.environ, PROJECT_DIR_ENV: str(project_dir())}

            flush_print(f"\n{'=' * 70}")
            flush_print(f"  Running batch: {protocol_id}/{batch_id}")
            flush_print(f"  Config: {config_path}")
            flush_print(f"{'=' * 70}")

            ret = subprocess.run(cmd, env=child_env).returncode
            if ret == 0:
                results.append((protocol_id, batch_id, "success"))
                n_run += 1
            else:
                results.append((protocol_id, batch_id, f"failed({ret})"))

    flush_print(f"\n{'=' * 70}")
    flush_print(f"  Batch Runner Summary: {n_run} ran, {n_skip} skipped")
    for pid, bid, status in results:
        flush_print(f"    {pid}/{bid}: {status}")
    flush_print(f"{'=' * 70}")

    return results


def _run_qc_gate(all_qc, quant_summary, qc_dir):
    """
    按配置阈值执行 QC 门控，落盘判定结果，并按 ``on_qc_fail`` 策略处置。

    在此之前 ``analysis_defaults.qc`` / ``alignment`` / ``batch_qc`` 里的
    阈值从未被任何代码读取——QC 是纯描述性的，不合格样本不会被拦下。

    Returns
    -------
    list[str]
        需要排除的样本（仅 ``exclude`` 策略下非空）。
    """
    if not QC_THRESHOLDS:
        flush_print("  [QC gate] 未配置阈值，跳过门控（仅生成描述性 QC 报告）。")
        return []

    qc_summary = qc_gate.summarize_qc_reports(all_qc)
    if qc_summary.empty:
        flush_print("  [QC gate] 无可用 QC 报告，跳过门控。")
        return []

    gate_df = qc_gate.evaluate_qc_gate(qc_summary, quant_summary, QC_THRESHOLDS)
    gate_path = os.path.join(qc_dir, "qc_gate_report.csv")
    gate_df.to_csv(gate_path, index=False)

    n_fail = int((~gate_df["qc_pass"]).sum())
    flush_print(f"\n  [QC gate] 阈值: "
                + ", ".join(f"{k}={v}" for k, v in sorted(QC_THRESHOLDS.items())))
    flush_print(f"  [QC gate] {len(gate_df) - n_fail}/{len(gate_df)} 样本通过 → {gate_path}")

    if n_fail == 0:
        return []

    for _, row in gate_df[~gate_df["qc_pass"]].iterrows():
        flush_print(f"    [FAIL] {row['sample']}: {row['qc_failures']}")

    excluded, policy = qc_gate.apply_qc_policy(gate_df, ON_QC_FAIL)
    if policy == qc_gate.POLICY_FAIL:
        raise RuntimeError(
            f"QC gate failed for {n_fail} sample(s) and on_qc_fail='fail'. "
            f"See {gate_path}."
        )
    if policy == qc_gate.POLICY_EXCLUDE:
        flush_print(f"  [QC gate] on_qc_fail='exclude' → 从下游分析中移除 {len(excluded)} 个样本。")
        return excluded

    flush_print("  [QC gate] on_qc_fail='flag' → 仅告警，样本仍参与分析。")
    return []


def get_sample_metadata(samples_df):
    """
    Build metadata DataFrame from sample list.
    Extracts mouse ID and tissue from sample names like '52_brain'.
    """
    meta_rows = []
    for _, row in samples_df.iterrows():
        sid = row['sample']
        parts = sid.split('_', 1)
        mouse_id = parts[0]
        tissue = parts[1] if len(parts) > 1 else 'unknown'
        group = MOUSE_GROUPS.get(sid) or MOUSE_GROUPS.get(mouse_id) or 'Unknown'
        batch = SAMPLE_BATCHES.get(sid) or SAMPLE_BATCHES.get(mouse_id)
        meta_rows.append({
            'sample': sid,
            'mouse_id': mouse_id,
            'tissue': tissue,
            'condition': group,
            'batch': batch,
        })

    meta_df = pd.DataFrame(meta_rows).set_index('sample')
    if SAMPLE_BATCHES:
        missing_batch = meta_df['batch'].isna()
        if missing_batch.any():
            raise ValueError(
                "sample_batches is configured but missing samples: "
                + ", ".join(meta_df.index[missing_batch].astype(str))
            )
    elif 'batch' in meta_df.columns:
        meta_df = meta_df.drop(columns=['batch'])
    return meta_df


def normalize_sample_list_columns(samples_df: pd.DataFrame) -> pd.DataFrame:
    """
    Normalize sample list columns to: sample, file1, file2.
    """
    sample_col_candidates = ["sample", "sample_id", "complete_result_name"]
    r1_col_candidates = ["file1", "fastq_1", "r1", "read1"]
    r2_col_candidates = ["file2", "fastq_2", "r2", "read2"]

    sample_col = next((c for c in sample_col_candidates if c in samples_df.columns), None)
    r1_col = next((c for c in r1_col_candidates if c in samples_df.columns), None)
    r2_col = next((c for c in r2_col_candidates if c in samples_df.columns), None)

    missing = []
    if not sample_col:
        missing.append("sample")
    if not r1_col:
        missing.append("file1/fastq_1")
    if not r2_col:
        missing.append("file2/fastq_2")
    if missing:
        raise ValueError(f"sample list is missing required columns: {', '.join(missing)}")

    normalized = samples_df.rename(columns={
        sample_col: "sample",
        r1_col: "file1",
        r2_col: "file2",
    }).copy()
    return normalized


def main():
    """Main pipeline entry point."""
    parser = argparse.ArgumentParser(
        description="Bulk RNA-Seq pipeline (kallisto + PyDESeq2) driven by experiment config.",
        epilog=config_arg_help(),
    )
    parser.add_argument("--experiment", default=None, help="Experiment id under experiments/<id>/transcriptomics_data/")
    parser.add_argument("--config", default=None, help="Path to transcriptomics run-entry config.yaml or batch.yaml")
    parser.add_argument("--outdir", default=None, help="Override output_dir (transcriptomics_results root)")
    parser.add_argument("--all-batches", action="store_true",
                        help="Run all enabled protocol/batch entries from experiment.yaml")
    parser.add_argument("--integrate", action="store_true",
                        help="Run multi-batch/multi-protocol integration after analysis "
                             "(joint ~batch + condition model + meta-analysis)")
    parser.add_argument("--joint-model-only", action="store_true",
                        help="Run only the cross-batch joint model on existing run outputs "
                             "(no re-quantification, no meta-analysis)")
    args = parser.parse_args()

    if args.joint_model_only:
        if not args.experiment:
            parser.error("--joint-model-only requires --experiment <id>.")
        from .joint_model import run_joint_model
        summary = run_joint_model(args.experiment, outdir=args.outdir)
        flush_print("\n" + json.dumps(summary, indent=2, ensure_ascii=False, default=str))
        return

    if args.all_batches:
        if not args.experiment:
            parser.error("--all-batches requires --experiment <id>.")
        run_all_batches(args.experiment, outdir=args.outdir)
        if args.integrate:
            run_integration_analysis(args.experiment)
        return

    cfg_info = _apply_config(args.experiment, args.config, args.outdir)

    start_time = datetime.now()

    flush_print("\n" + "=" * 60)
    flush_print("  Bulk RNA-Seq Analysis Pipeline")
    flush_print("  kallisto + PyDESeq2 (Windows Native)")
    flush_print("=" * 60)
    flush_print(f"Start time: {start_time}")
    flush_print(f"Experiment: {cfg_info.get('experiment_id')}")
    if cfg_info.get("protocol_id") and cfg_info.get("batch_id"):
        flush_print(f"Run:        {cfg_info.get('protocol_id')}/{cfg_info.get('batch_id')}")
    flush_print(f"Config:     {cfg_info.get('config_path')}")
    flush_print(f"Raw data:   {RAW_DATA_ROOT}")
    flush_print(f"Output root:{cfg_info.get('output_root')}")
    flush_print(f"Output:     {OUTPUT_BASE}")

    check_environment()

    # Read sample list
    try:
        samples_df = pd.read_csv(SAMPLE_LIST_PATH, sep='\t')
        samples_df = normalize_sample_list_columns(samples_df)
        flush_print(f"\nLoaded {len(samples_df)} samples from sample list.")
    except Exception as e:
        flush_print(f"ERROR: Cannot read sample list: {e}")
        return

    # Build metadata
    metadata_df = get_sample_metadata(samples_df)
    flush_print(f"\nExperimental groups:")
    for grp, cnt in metadata_df['condition'].value_counts().items():
        mice = sorted(metadata_df[metadata_df['condition'] == grp]['mouse_id'].unique())
        flush_print(f"  {grp}: {cnt} samples (mice: {', '.join(mice)})")

    # Create output directories
    qc_dir = os.path.join(OUTPUT_BASE, "01_qc")
    trim_dir = os.path.join(OUTPUT_BASE, "02_trimmed")
    quant_dir = os.path.join(OUTPUT_BASE, "03_kallisto_quant")
    de_dir = os.path.join(OUTPUT_BASE, "04_de_results")
    for d in [qc_dir, trim_dir, quant_dir, de_dir]:
        os.makedirs(d, exist_ok=True)

    # Save metadata
    metadata_path = os.path.join(OUTPUT_BASE, "sample_metadata.csv")
    metadata_df.to_csv(metadata_path)
    flush_print(f"Metadata saved: {metadata_path}")

    # ---- STEP 1: QC ----
    flush_print("\n" + "=" * 60)
    flush_print("STEP 1: Quality Control (Pure Python)")
    flush_print("=" * 60)

    all_qc = {}
    for _, row in samples_df.iterrows():
        sid = row['sample']
        r1 = os.path.join(RAW_DATA_ROOT, row['file1'])
        r2 = os.path.join(RAW_DATA_ROOT, row['file2'])

        if not os.path.exists(r1):
            flush_print(f"  SKIP {sid}: R1 not found at {r1}")
            continue

        qc = run_qc_for_sample(sid, r1, r2, qc_dir, max_reads=QC_SAMPLE_SIZE)
        all_qc[sid] = qc

    # Save QC summary
    if all_qc:
        qc_rows = []
        for sid, qc in all_qc.items():
            for read_end in ['R1', 'R2']:
                if read_end in qc:
                    r = qc[read_end]
                    qc_rows.append({
                        'sample': sid,
                        'read': read_end,
                        'total_reads': r['total_reads'],
                        'mean_length': r['mean_read_length'],
                        'gc_pct': r['gc_content_pct'],
                        'mean_quality': r['mean_quality'],
                        'q20_pct': r['q20_pct'],
                        'q30_pct': r['q30_pct']
                    })
        qc_summary = pd.DataFrame(qc_rows)
        qc_summary_path = os.path.join(qc_dir, "qc_summary.csv")
        qc_summary.to_csv(qc_summary_path, index=False)
        flush_print(f"\nQC Summary saved: {qc_summary_path}")

    # ---- STEP 2+3: TRIM → QUANTIFY → CLEANUP (per sample) ----
    flush_print("\n" + "=" * 60)
    flush_print("STEP 2+3: Trim → Quantify → Cleanup (per sample)")
    flush_print("  Trimmed files are deleted after quantification to save disk space.")
    flush_print("  Only cutadapt reports (JSON) and kallisto results are kept.")
    flush_print("=" * 60)

    # Build kallisto index first
    if not build_kallisto_index():
        flush_print("  Cannot proceed without kallisto index.")
        return

    quant_files = {}
    for idx, row in enumerate(samples_df.itertuples(), 1):
        sid = row.sample
        r1 = os.path.join(RAW_DATA_ROOT, row.file1)
        r2 = os.path.join(RAW_DATA_ROOT, row.file2)

        if not os.path.exists(r1):
            flush_print(f"  SKIP {sid}: R1 not found")
            continue

        # Check if kallisto output already exists (skip entirely)
        abundance_file = os.path.join(quant_dir, sid, "abundance.tsv")
        if os.path.exists(abundance_file):
            flush_print(f"  [{idx}/{len(samples_df)}] {sid}: quant exists, skipping")
            quant_files[sid] = abundance_file
            continue

        flush_print(f"\n  [{idx}/{len(samples_df)}] {sid}")

        # 2a. Trim
        r1_t, r2_t = run_trimming(sid, r1, r2, trim_dir)

        # 2b. Quantify
        af = run_kallisto_quant(sid, r1_t, r2_t, quant_dir, strandedness=STRANDEDNESS)
        if af:
            quant_files[sid] = af

        # 2c. Cleanup trimmed files (keep cutadapt JSON report)
        for tf in [r1_t, r2_t]:
            if tf.endswith('.trimmed.fq.gz') and os.path.exists(tf):
                sz_mb = os.path.getsize(tf) / (1024 * 1024)
                try:
                    os.remove(tf)
                    flush_print(f"    Cleaned: {os.path.basename(tf)} ({sz_mb:.0f} MB freed)")
                except OSError:
                    flush_print(f"    WARN: Could not delete {os.path.basename(tf)}")

    flush_print(f"\nTrim+Quantify done: {len(quant_files)}/{len(samples_df)} samples.")

    # ---- STEP 4: MERGE & DIFFERENTIAL EXPRESSION ----
    flush_print("\n" + "=" * 60)
    flush_print("STEP 4: Merge Counts & Differential Expression (PyDESeq2)")
    flush_print("=" * 60)

    sample_ids = list(quant_files.keys())
    counts_df, tpm_df, efflen_df = merge_kallisto_counts(quant_dir, sample_ids)

    if counts_df is None:
        flush_print("  No quantification data. Cannot run DE analysis.")
        return

    # ---- 比对率落盘 + QC 门控 ----
    quant_summary = qc_gate.collect_quant_summary(quant_dir, sample_ids)
    quant_summary_path = os.path.join(quant_dir, "quant_summary.csv")
    quant_summary.to_csv(quant_summary_path, index=False)
    flush_print(f"  Alignment summary: {quant_summary_path}")

    excluded_samples = _run_qc_gate(all_qc, quant_summary, qc_dir)
    if excluded_samples:
        keep = [s for s in counts_df.columns if s not in set(excluded_samples)]
        counts_df = counts_df[keep]
        tpm_df = tpm_df[keep]
        if efflen_df is not None:
            efflen_df = efflen_df[[c for c in efflen_df.columns if c in keep]]
        sample_ids = [s for s in sample_ids if s not in set(excluded_samples)]

    # Save raw matrices
    counts_path = os.path.join(de_dir, "transcript_counts_matrix.csv")
    tpm_path = os.path.join(de_dir, "transcript_tpm_matrix.csv")
    counts_df.round(0).astype(int).to_csv(counts_path)
    tpm_df.to_csv(tpm_path)
    flush_print(f"  Transcript count matrix: {counts_path}")
    flush_print(f"  Transcript TPM matrix:   {tpm_path}")

    # Build tx2gene mapping and aggregate to gene level
    tx2gene = build_tx2gene_from_fasta(TRANSCRIPTOME_FASTA)
    gene_counts, agg_info = aggregate_counts_to_gene(
        counts_df, tx2gene, tpm_df=tpm_df, efflen_df=efflen_df, method=GENE_AGGREGATION,
    )

    gene_counts_path = os.path.join(de_dir, "gene_counts_matrix.csv")
    gene_counts.to_csv(gene_counts_path)
    flush_print(f"  Gene count matrix: {gene_counts_path}")

    # ---- 运行溯源：参考指纹 + 工具版本 + 技术参数 ----
    prov = provenance_mod.capture_run_provenance(
        OUTPUT_BASE,
        experiment_id=str(cfg_info.get("experiment_id")),
        protocol_id=PROTOCOL_ID,
        batch_id=BATCH_ID,
        transcriptome_fasta=TRANSCRIPTOME_FASTA,
        kallisto_exe=KALLISTO_EXE,
        kallisto_index=KALLISTO_INDEX,
        strandedness=STRANDEDNESS,
        gene_aggregation=agg_info.get("gene_aggregation", GENE_AGGREGATION),
        technical_params=TECHNICAL_PARAMS,
        extra={"n_samples": len(sample_ids), "n_genes": int(gene_counts.shape[0])},
    )
    flush_print(f"  Run provenance: {prov.get('_path')}")

    # 为 gene_counts_matrix 添加 Ensembl/NCBI 数据库注释
    try:
        annotate_csv_file(gene_counts_path)
    except Exception as e:
        flush_print(f"  [WARN] Gene annotation skipped: {e}")

    # Run DE analysis per tissue
    flush_print("\n--- Differential Expression Analysis per Tissue ---")

    # Comparison list from config (fallback to pairwise in config loader)
    comparisons = list(COMPARISONS)

    # Write resolved_comparisons.json for traceability
    resolved_comp_path = os.path.join(de_dir, "resolved_comparisons.json")
    comp_records = [
        {"numerator": n, "denominator": d, "comparison_name": c}
        for n, d, c in comparisons
    ]
    with open(resolved_comp_path, 'w', encoding='utf-8') as f:
        json.dump({
            "experiment_id": cfg_info.get("experiment_id"),
            "n_comparisons": len(comparisons),
            "comparisons": comp_records,
        }, f, indent=2, ensure_ascii=False)
    flush_print(f"  Resolved comparisons saved: {resolved_comp_path}")
    for n, d, c in comparisons:
        flush_print(f"    - {c}: {n} vs {d}")

    for tissue in TISSUES:
        flush_print(f"\n{'='*40}")
        flush_print(f"Tissue: {tissue}")
        flush_print(f"{'='*40}")

        # Get samples for this tissue
        tissue_samples = [s for s in gene_counts.columns
                          if s.endswith(f'_{tissue}')]

        if len(tissue_samples) < 2:
            flush_print(f"  Skipping {tissue}: too few samples ({len(tissue_samples)})")
            continue

        tissue_de_dir = os.path.join(de_dir, tissue)
        os.makedirs(tissue_de_dir, exist_ok=True)

        tissue_counts = gene_counts[tissue_samples]
        tissue_meta = metadata_df.loc[tissue_samples]

        flush_print(f"  Samples: {len(tissue_samples)}")
        for grp in tissue_meta['condition'].unique():
            n = (tissue_meta['condition'] == grp).sum()
            flush_print(f"    {grp}: {n}")

        # Run each comparison
        for numerator, denominator, comp_name in comparisons:
            comp_samples = tissue_meta[
                tissue_meta['condition'].isin([numerator, denominator])
            ].index.tolist()

            if len(comp_samples) < 2:
                flush_print(f"  [{comp_name}] Skipping: too few samples")
                continue

            n_num = (tissue_meta.loc[comp_samples, 'condition'] == numerator).sum()
            n_den = (tissue_meta.loc[comp_samples, 'condition'] == denominator).sum()

            if n_num < 1 or n_den < 1:
                flush_print(f"  [{comp_name}] Skipping: need >= 1 per group")
                continue

            flush_print(f"\n  --- {comp_name} ({tissue}) ---")
            flush_print(f"    {numerator}: {n_num}, {denominator}: {n_den}")

            # Warn about low-replicate evidence level
            if n_num < 2 or n_den < 2:
                flush_print(f"    ⚠ Low replicate: results are effect_size_only — no statistical conclusions")

            de_done = False
            has_batch_design = (
                'batch' in tissue_meta.columns
                and tissue_meta.loc[comp_samples, 'batch'].nunique() > 1
            )
            design_formula = 'batch + condition' if has_batch_design else 'condition'

            # Try PyDESeq2 first (needs >= 2 per group)
            if n_num >= 2 and n_den >= 2:
                try:
                    de_results = run_pydeseq2_analysis(
                        tissue_counts[comp_samples],
                        tissue_meta.loc[comp_samples],
                        design_formula=design_formula,
                        contrast=('condition', numerator, denominator),
                    )

                    if de_results is not None:
                        out_file = os.path.join(tissue_de_dir, f"de_{comp_name}.csv")
                        de_results.to_csv(out_file)
                        flush_print(f"  Saved: {out_file}")
                        # 添加数据库注释
                        try:
                            annotate_csv_file(out_file)
                        except Exception as e:
                            flush_print(f"    [WARN] DE annotation skipped: {e}")
                        # 多组学稳定接口导出
                        try:
                            export_multiomics_csv(de_results, tissue_de_dir, comp_name)
                        except Exception as e:
                            flush_print(f"    [WARN] Multi-omics export skipped: {e}")
                        de_done = True
                except Exception as e:
                    flush_print(f"  PyDESeq2 failed for {comp_name}: {e}")
                    if has_batch_design:
                        flush_print(
                            "  [BLOCKED] Explicit multi-batch design cannot be fit; "
                            "unadjusted Wilcoxon fallback is disabled."
                        )
                        continue

            # Fallback: Wilcoxon rank-sum (works with n >= 1 per group)
            if not de_done:
                if has_batch_design:
                    flush_print(
                        "  [BLOCKED] Multi-batch comparison requires an identifiable "
                        "batch-adjusted model; unadjusted fallback is disabled."
                    )
                    continue
                if n_num < 2 or n_den < 2:
                    flush_print(f"  [Fallback] Wilcoxon effect-only (n_case={n_num}, n_control={n_den})")
                else:
                    flush_print(f"  [Fallback] Using Wilcoxon rank-sum for {comp_name} ...")
                try:
                    de_results = run_wilcoxon_analysis(
                        tissue_counts[comp_samples],
                        tissue_meta.loc[comp_samples],
                        numerator=numerator,
                        denominator=denominator,
                    )
                    if de_results is not None:
                        out_file = os.path.join(tissue_de_dir, f"de_{comp_name}_wilcoxon.csv")
                        de_results.to_csv(out_file)
                        flush_print(f"  Saved (Wilcoxon): {out_file}")
                        # 添加数据库注释
                        try:
                            annotate_csv_file(out_file)
                        except Exception as e:
                            flush_print(f"    [WARN] DE annotation skipped: {e}")
                        # 多组学稳定接口导出
                        try:
                            export_multiomics_csv(de_results, tissue_de_dir, comp_name)
                        except Exception as e:
                            flush_print(f"    [WARN] Multi-omics export skipped: {e}")
                except Exception as e:
                    flush_print(f"  Wilcoxon also failed for {comp_name}: {e}")

    # ---- GLOBAL CLEANUP: delete any leftover .fq.gz files ----
    flush_print("\n" + "=" * 60)
    flush_print("Cleanup: scanning for leftover trimmed FASTQ files...")
    flush_print("=" * 60)
    leftover_freed = 0
    leftover_count = 0
    for dirpath, _, filenames in os.walk(trim_dir):
        for fname in filenames:
            if fname.endswith('.fq.gz'):
                fpath = os.path.join(dirpath, fname)
                sz = os.path.getsize(fpath)
                try:
                    os.remove(fpath)
                    leftover_freed += sz
                    leftover_count += 1
                    flush_print(f"  Deleted: {fname} ({sz / (1024*1024):.0f} MB)")
                except OSError as e:
                    flush_print(f"  WARN: Cannot delete {fname} ({e})")
    if leftover_count == 0:
        flush_print("  No leftover .fq.gz files found. Clean!")
    else:
        flush_print(f"  Removed {leftover_count} files, freed {leftover_freed / (1024**3):.2f} GB")

    # ---- SUMMARY ----
    end_time = datetime.now()
    duration = end_time - start_time

    flush_print("\n" + "=" * 60)
    flush_print("PIPELINE COMPLETE")
    flush_print("=" * 60)
    flush_print(f"Duration: {duration}")
    flush_print(f"\nOutput directory: {OUTPUT_BASE}")
    flush_print(f"  01_qc/              - Quality control reports")
    flush_print(f"  02_trimmed/         - cutadapt reports only (trimmed FASTQ auto-deleted)")
    flush_print(f"  03_kallisto_quant/  - kallisto quantification results")
    flush_print(f"  04_de_results/{tissue}/ - DE analysis results per tissue")
    flush_print(f"  sample_metadata.csv - Sample metadata")
    flush_print(f"\nDE comparisons performed per tissue:")
    for num, den, cname in comparisons:
        flush_print(f"  - {cname}: {num} vs {den}")

    if args.integrate:
        flush_print("\n" + "=" * 60)
        flush_print("STEP 5: Multi-batch / Multi-protocol integration")
        flush_print("=" * 60)
        run_integration_analysis(cfg_info["experiment_id"])


if __name__ == "__main__":
    main()
