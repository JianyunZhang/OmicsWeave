"""
Quick FASTQ QC Test - Process ONE sample to verify the pipeline works.
Pure Python, no external tools needed. Compatible with Python 3.7+.
"""
import os
import sys
import gzip
import json
import collections
import time
import argparse
import numpy as np
import pandas as pd

from .experiment_config import load_transcriptomics_config, config_arg_help

# Config-driven (set in main)
RAW_DATA_ROOT = None
OUTPUT_DIR = None
QC_SAMPLE_SIZE = 50000  # Only read first 50k reads for speed


def flush_print(msg):
    """Print with immediate flush."""
    print(msg, flush=True)


def phred_to_quality(qual_string):
    return [ord(c) - 33 for c in qual_string]


def fastq_qc(filepath, max_reads=50000):
    """QC a single FASTQ.GZ file."""
    flush_print(f"    Reading: {os.path.basename(filepath)}")
    flush_print(f"    (Sampling first {max_reads} reads from network drive...)")
    
    total_reads = 0
    total_bases = 0
    gc_count = 0
    n_count = 0
    read_lengths = []
    qual_sum = 0
    qual_count = 0
    q20_count = 0
    q30_count = 0
    
    start = time.time()
    
    with gzip.open(filepath, 'rt', encoding='ascii') as f:
        while total_reads < max_reads:
            header = f.readline()
            if not header:
                break
            seq = f.readline().strip()
            _ = f.readline()  # +
            qual = f.readline().strip()
            
            total_reads += 1
            seq_len = len(seq)
            total_bases += seq_len
            read_lengths.append(seq_len)
            
            gc_count += seq.count('G') + seq.count('C')
            n_count += seq.count('N')
            
            scores = phred_to_quality(qual)
            qual_sum += sum(scores)
            qual_count += len(scores)
            q20_count += sum(1 for s in scores if s >= 20)
            q30_count += sum(1 for s in scores if s >= 30)
            
            if total_reads % 10000 == 0:
                elapsed = time.time() - start
                flush_print(f"      ... {total_reads} reads processed ({elapsed:.1f}s)")
    
    elapsed = time.time() - start
    
    result = {
        'file': os.path.basename(filepath),
        'reads_sampled': total_reads,
        'total_bases': total_bases,
        'mean_read_length': round(np.mean(read_lengths), 1),
        'gc_content_pct': round(gc_count / total_bases * 100, 2) if total_bases > 0 else 0,
        'n_content_pct': round(n_count / total_bases * 100, 4) if total_bases > 0 else 0,
        'mean_quality': round(qual_sum / qual_count, 2) if qual_count > 0 else 0,
        'q20_pct': round(q20_count / qual_count * 100, 2) if qual_count > 0 else 0,
        'q30_pct': round(q30_count / qual_count * 100, 2) if qual_count > 0 else 0,
        'processing_time_sec': round(elapsed, 1)
    }
    
    return result


def main():
    global RAW_DATA_ROOT, OUTPUT_DIR

    parser = argparse.ArgumentParser(
        description="Quick FASTQ QC test for ONE sample (config-driven).",
        epilog=config_arg_help(),
    )
    parser.add_argument("--experiment", default=None, help="Experiment id under experiments/<id>/transcriptomics*/")
    parser.add_argument("--config", default=None, help="Path to transcriptomics run-entry config.yaml or batch.yaml")
    parser.add_argument("--outdir", default=None, help="Override output_dir (transcriptomics_results root)")
    parser.add_argument("--sample", required=True, help="Sample name, e.g. 52_brain (must exist in sample_list.tsv).")
    parser.add_argument("--reads", type=int, default=QC_SAMPLE_SIZE, help="Number of reads to sample for QC.")
    args = parser.parse_args()

    cfg = load_transcriptomics_config(
        experiment_id=args.experiment,
        config_path=args.config,
        outdir_override=args.outdir,
    )

    RAW_DATA_ROOT = cfg.raw_data_root
    os.environ["TRANSCRIPTOMICS_OUTDIR"] = os.path.abspath(str(cfg.output_dir))
    OUTPUT_DIR = os.path.join(os.environ["TRANSCRIPTOMICS_OUTDIR"], "01_qc")
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    df = pd.read_csv(cfg.sample_list_tsv, sep="\t")
    if "sample" not in df.columns:
        raise SystemExit("sample_list.tsv missing required column: sample")
    row = df[df["sample"] == args.sample]
    if row.empty:
        raise SystemExit(f"Sample not found in sample_list.tsv: {args.sample}")
    r = row.iloc[0].to_dict()

    fq1 = str(r.get("fastq_1") or r.get("fq1") or r.get("r1") or "").strip()
    fq2 = str(r.get("fastq_2") or r.get("fq2") or r.get("r2") or "").strip()
    if not fq1:
        raise SystemExit("sample_list.tsv missing fastq_1 (or fq1/r1) for this sample.")

    def _resolve_fastq(p: str) -> str:
        if p.startswith("\\\\") or (len(p) > 2 and p[1] == ":" and (p[2] == "\\" or p[2] == "/")):
            return p
        return os.path.join(RAW_DATA_ROOT, p)

    test_r1 = _resolve_fastq(fq1)
    test_r2 = _resolve_fastq(fq2) if fq2 else None

    flush_print("=" * 60)
    flush_print("  FASTQ Quality Control (Pure Python)")
    flush_print("=" * 60)
    flush_print(f"Experiment: {cfg.experiment_id}")
    flush_print(f"Python:     {sys.version}")
    flush_print(f"Sample:     {args.sample}")
    flush_print(f"FASTQ1:     {test_r1}")
    if test_r2:
        flush_print(f"FASTQ2:     {test_r2}")
    flush_print(f"Output:     {OUTPUT_DIR}")
    flush_print(f"Reads:      {args.reads}")
    flush_print("")

    if not os.path.exists(test_r1):
        raise SystemExit(f"ERROR: Cannot access {test_r1}")
    if test_r2 and (not os.path.exists(test_r2)):
        raise SystemExit(f"ERROR: Cannot access {test_r2}")

    # QC
    flush_print("  [R1] Processing...")
    qc_r1 = fastq_qc(test_r1, max_reads=int(args.reads))
    flush_print("  [R1] Done!")

    qc_r2 = None
    if test_r2:
        flush_print("  [R2] Processing...")
        qc_r2 = fastq_qc(test_r2, max_reads=int(args.reads))
        flush_print("  [R2] Done!")

    # Print results
    flush_print("\n" + "=" * 60)
    flush_print(f"  QC RESULTS: {args.sample}")
    flush_print("=" * 60)

    for label, qc in [("R1", qc_r1), ("R2", qc_r2)]:
        if qc is None:
            continue
        flush_print(f"\n  {label}: {qc['file']}")
        flush_print(f"    Reads sampled:    {qc['reads_sampled']:,}")
        flush_print(f"    Total bases:      {qc['total_bases']:,}")
        flush_print(f"    Mean read length: {qc['mean_read_length']} bp")
        flush_print(f"    GC content:       {qc['gc_content_pct']}%")
        flush_print(f"    N content:        {qc['n_content_pct']}%")
        flush_print(f"    Mean quality:     Q{qc['mean_quality']}")
        flush_print(f"    Q20 bases:        {qc['q20_pct']}%")
        flush_print(f"    Q30 bases:        {qc['q30_pct']}%")
        flush_print(f"    Processing time:  {qc['processing_time_sec']}s")

    # Save
    report = {args.sample: {"R1": qc_r1, "R2": qc_r2}}
    report_path = os.path.join(OUTPUT_DIR, f"{args.sample}_qc_quick.json")
    with open(report_path, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)
    flush_print(f"\n  Report saved: {report_path}")

    # Quality assessment
    flush_print("\n" + "=" * 60)
    flush_print("  QUALITY ASSESSMENT")
    flush_print("=" * 60)

    for label, qc in [("R1", qc_r1), ("R2", qc_r2)]:
        if qc is None:
            continue
        issues = []
        if qc["mean_quality"] < 20:
            issues.append("Low mean quality (< Q20)")
        if qc["q30_pct"] < 70:
            issues.append("Low Q30 percentage (< 70%)")
        if qc["gc_content_pct"] < 30 or qc["gc_content_pct"] > 60:
            issues.append(f"Unusual GC content ({qc['gc_content_pct']}%)")
        if qc["n_content_pct"] > 1:
            issues.append(f"High N content ({qc['n_content_pct']}%)")

        if issues:
            flush_print(f"  {label}: WARNINGS")
            for issue in issues:
                flush_print(f"    - {issue}")
        else:
            flush_print(f"  {label}: PASS - All quality metrics look good!")

    flush_print("\nTest complete.")


if __name__ == "__main__":
    main()
