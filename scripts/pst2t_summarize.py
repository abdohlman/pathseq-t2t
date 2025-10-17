#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
pst2t_summarize.py
Combined summarization + normalization step for PathSeq-T2T.

What this script does (single entry point for the new "summarize" step):
1) Reads PathSeq filtering metrics + T2T flagstats from --filter-stats-dir
2) Reads raw classifier outputs from --classification-stats-dir
   - Kraken2: <sample>.paired.kraken.report.txt and <sample>.unpaired.kraken.report.txt
   - MetaPhlAn4: <sample>.metaphlan.report.txt
3) Writes exactly two types of outputs into --results-dir (user-defined):
   a) <sample>.summary.tsv  (ONE wide row merging filtering + classification key totals)
   b) Normalized tables (RPM using PRIMARY_READS):
        - <sample>.kraken.txt       (if Kraken2 present)
        - <sample>.metaphlan.txt    (if MetaPhlAn present)

Exit codes:
- 0 on success
- nonzero with a helpful message on error
"""

import argparse
import os
import sys
from typing import Optional, Dict

import pandas as pd

# ------------------------ Utilities ------------------------

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def require_dir(path: str, label: str):
    if not path or not os.path.isdir(path):
        eprint(f"ERROR: {label} not found or not a directory: {path}")
        sys.exit(2)


def require_parent(path: str):
    parent = os.path.dirname(path) or '.'
    os.makedirs(parent, exist_ok=True)


# ------------------------ Filtering summary (from pst2t_summarize_filtering.py) ------------------------

def parse_flagstat_primary(path: str) -> Optional[int]:
    """
    Expect 3-column TSV rows, with a row 'primary' in column 3.
    Example line: 5995946    0    primary
    Returns int(primary) or None if not found.
    """
    try:
        with open(path, 'r', encoding='utf-8') as fh:
            for line in fh:
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 3 and parts[2].strip() == 'primary':
                    val = parts[0].strip()
                    if val.isdigit():
                        return int(val)
        return None
    except Exception as ex:
        eprint(f'WARNING: Could not read flagstat file: {path} ({ex})')
        return None


def collect_primary_reads(filter_stats_dir: str, sample: str) -> pd.Series:
    rs: Dict[str, int] = {}

    mapping = {
        'PRIMARY_READS':                    f'{sample}.flagstat.tsv',
        'T2T_ALIGNED_PAIRED_PRIMARY':       f'{sample}.qcfilt_paired.t2t_aln.flagstat.tsv',
        'T2T_UNALIGNED_PAIRED_PRIMARY':     f'{sample}.qcfilt_paired.t2t_unaln.flagstat.tsv',
        'T2T_ALIGNED_UNPAIRED_PRIMARY':     f'{sample}.qcfilt_unpaired.t2t_aln.flagstat.tsv',
        'T2T_UNALIGNED_UNPAIRED_PRIMARY':   f'{sample}.qcfilt_unpaired.t2t_unaln.flagstat.tsv',
    }
    for key, fname in mapping.items():
        path = os.path.join(filter_stats_dir, fname)
        if os.path.isfile(path):
            v = parse_flagstat_primary(path)
            if v is not None:
                rs[key] = v

    rs['T2T_UNALIGNED_TOTAL_READS'] = int(rs.get('T2T_UNALIGNED_PAIRED_PRIMARY', 0)) + int(
        rs.get('T2T_UNALIGNED_UNPAIRED_PRIMARY', 0)
    )

    return pd.Series(rs, dtype='Int64')


PATHSEQ_COLS = [
    'PRIMARY_READS',
    'READS_AFTER_PREALIGNED_HOST_FILTER',
    'READS_AFTER_QUALITY_AND_COMPLEXITY_FILTER',
    'READS_AFTER_HOST_FILTER',
    'READS_AFTER_DEDUPLICATION',
    'FINAL_PAIRED_READS',
    'FINAL_UNPAIRED_READS',
    'FINAL_TOTAL_READS',
]


def read_pathseq_metrics_table(path: str) -> Optional[pd.Series]:
    if not os.path.isfile(path):
        return None

    try:
        df = pd.read_table(path, comment='#', sep='\t')
        if df.empty or 'PRIMARY_READS' not in df.columns:
            df = pd.read_table(path, sep='\t', skiprows=6)
        if 'PRIMARY_READS' not in df.columns:
            with open(path, 'r', encoding='utf-8') as fh:
                lines = fh.readlines()
            hdr_idx = None
            for i, ln in enumerate(lines):
                if ln.startswith('PRIMARY_READS'):
                    hdr_idx = i
                    break
            if hdr_idx is None:
                return None
            df = pd.read_table(path, sep='\t', skiprows=hdr_idx)

        if df.shape[0] >= 1:
            s = df.iloc[0]
            s = pd.to_numeric(s, errors='coerce').fillna(0)
            s = s.astype('Int64')
            return s
        return None
    except Exception as ex:
        eprint(f'WARNING: Could not parse PathSeq metrics: {path} ({ex})')
        return None


def collect_pathseq_metrics(filter_stats_dir: str, sample: str) -> pd.Series:
    rs: Dict[str, int] = {}
    files = {
        'UNALIGNED': os.path.join(filter_stats_dir, f'{sample}.unaligned.filter_metrics.txt'),
        'EXCLUDED':  os.path.join(filter_stats_dir, f'{sample}.excluded.filter_metrics.txt'),
    }

    per_kind: Dict[str, pd.Series] = {}
    for kind, fpath in files.items():
        s = read_pathseq_metrics_table(fpath)
        if s is None:
            s = pd.Series(dtype='Int64')
        per_kind[kind] = s

        for col in PATHSEQ_COLS:
            rs[f'{kind}_{col}'] = int(s.get(col, 0))

    for col in PATHSEQ_COLS:
        rs[f'QCFILTER_{col}'] = int(per_kind['UNALIGNED'].get(col, 0)) + int(per_kind['EXCLUDED'].get(col, 0))

    return pd.Series(rs, dtype='Int64')


FILTER_SUMMARY_ORDER = (
    ['PRIMARY_READS']
    + [f'{kind}_{col}' for kind in ('QCFILTER', 'UNALIGNED', 'EXCLUDED') for col in PATHSEQ_COLS]
    + ['T2T_UNALIGNED_TOTAL_READS', 'T2T_UNALIGNED_PAIRED_PRIMARY', 'T2T_UNALIGNED_UNPAIRED_PRIMARY']
)


def filtering_summary(filter_stats_dir: str, sample: str) -> pd.Series:
    rs_flagstat = collect_primary_reads(filter_stats_dir, sample)
    rs_pathseq = collect_pathseq_metrics(filter_stats_dir, sample)
    merged = pd.concat([rs_flagstat, rs_pathseq])
    row = merged.reindex(FILTER_SUMMARY_ORDER).fillna(0).astype(int)
    row.name = sample
    return row


# ------------------------ Classification summary + normalization ------------------------

KRAKEN2_COLS = [
    'pct_reads',
    'reads_clade',
    'reads_taxon',
    'minimizers_count',
    'minimizers_distinct',
    'rank',
    'tax_id',
    'name',
]

MICROBIAL_TAX_IDS = [2, 4751, 2157, 10239]  # Bacteria, Fungi, Archaea, Viruses

def load_kraken_report(path: str) -> pd.DataFrame:
    if os.path.getsize(path) == 0:
        df = pd.DataFrame(columns=KRAKEN2_COLS)
        for c in ('reads_clade','reads_taxon','pct_reads'):
            df[c] = df[c].astype('float64').fillna(0)
        df['tax_id'] = df.index.astype('Int64')
        return df

    df = pd.read_table(path, header=None, dtype=str)
    if df.shape[1] < 8:
        raise ValueError(f'Malformed Kraken2 report: {path}')
    df = df.iloc[:, :8].copy()
    df.columns = KRAKEN2_COLS
    df['name'] = df['name'].str.strip()
    for c in ['pct_reads','reads_clade','reads_taxon','minimizers_count','minimizers_distinct','tax_id']:
        df[c] = pd.to_numeric(df[c], errors='coerce').fillna(0)
    return df


def merge_kraken_reports(paired_path: str, unpaired_path: str) -> pd.DataFrame:
    present = [os.path.isfile(paired_path), os.path.isfile(unpaired_path)]
    if not any(present):
        raise FileNotFoundError('No Kraken2 reports found.')
    dfs = []
    if present[0]:
        dfp = load_kraken_report(paired_path).copy()
        dfp['reads_clade'] = 2 * dfp['reads_clade']
        dfp['reads_taxon'] = 2 * dfp['reads_taxon']
        dfs.append(dfp)
    if present[1]:
        dfu = load_kraken_report(unpaired_path).copy()
        dfu['reads_clade'] = 2 * dfu['reads_clade']
        dfu['reads_taxon'] = 2 * dfu['reads_taxon']
        dfs.append(dfu)

    df = (
        pd.concat(dfs, ignore_index=True)
        .groupby(['name','tax_id','rank'], as_index=False)[['reads_clade','reads_taxon','minimizers_count','minimizers_distinct']]
        .sum()
    )
    total_reads = df.loc[df['tax_id'].isin([0,1]), 'reads_clade'].sum()
    df['pct_reads'] = (100.0 * df['reads_clade'] / total_reads).round(4) if total_reads > 0 else 0.0
    df = df[['pct_reads','reads_clade','reads_taxon','rank','tax_id','name']]
    return df


MPA_COLUMNS = [
    'clade_name',
    'clade_taxid',
    'relative_abundance',
    'coverage',
    'estimated_number_of_reads_from_the_clade',
]

def load_metaphlan_table(path: str) -> pd.DataFrame:
    df = pd.read_table(path, comment='#', header=None)
    if df.shape[1] < 5:
        raise RuntimeError(f'Unexpected MetaPhlAn format: {path}')
    df.columns = MPA_COLUMNS
    df['clade_name'] = df['clade_name'].astype(str)
    df['clade_taxid'] = pd.to_numeric(df['clade_taxid'], errors='coerce').fillna(-1).astype(int)
    df['relative_abundance'] = pd.to_numeric(df['relative_abundance'], errors='coerce').fillna(0.0)
    df['coverage'] = pd.to_numeric(df['coverage'], errors='coerce').fillna(0.0)
    df['estimated_number_of_reads_from_the_clade'] = pd.to_numeric(
        df['estimated_number_of_reads_from_the_clade'], errors='coerce'
    ).fillna(0.0)
    return df


def summarize_classification(classif_dir: str, sample: str) -> pd.Series:
    paired_path   = os.path.join(classif_dir, f'{sample}.paired.kraken.report.txt')
    unpaired_path = os.path.join(classif_dir, f'{sample}.unpaired.kraken.report.txt')
    mpa_path      = os.path.join(classif_dir, f'{sample}.metaphlan.report.txt')

    summary: Dict[str, int] = {}

    # Kraken2 totals
    try:
        df_k = merge_kraken_reports(paired_path, unpaired_path)
        summary['TOTAL_READS_TESTED']    = int(df_k.loc[df_k['tax_id'].isin([0,1]), 'reads_clade'].sum())
        summary['UNCLASSIFIED_TOTAL_K2'] = int(df_k.loc[df_k['tax_id'] == 0, 'reads_clade'].sum())
        summary['CLASSIFIED_TOTAL_K2']   = int(df_k.loc[df_k['tax_id'] == 1, 'reads_clade'].sum())
        summary['MICROBIAL_TOTAL_K2']    = int(df_k.loc[df_k['tax_id'].isin([2,4751,2157,10239]), 'reads_clade'].sum())
        summary['HUMAN_TOTAL_K2']        = int(df_k.loc[df_k['tax_id'] == 9606, 'reads_clade'].sum())
    except FileNotFoundError:
        pass

    # MetaPhlAn totals (parsed from header + table)
    if os.path.isfile(mpa_path):
        total_reads = 0
        classified_reads = 0
        with open(mpa_path, 'r', encoding='utf-8') as fh:
            for line in fh:
                if not line.startswith('#'):
                    break
                if 'reads processed' in line:
                    parts = line.strip('# \n').split()
                    for i,p in enumerate(parts):
                        if p.isdigit() and parts[i+1] == 'reads' and parts[i+2].startswith('processed'):
                            total_reads = int(p); break
                if line.startswith('#estimated_reads_mapped_to_known_clades:'):
                    try:
                        classified_reads = int(line.split(':')[1])
                    except Exception:
                        pass
        unclassified_reads = max(0, total_reads - classified_reads)

        df_m = load_metaphlan_table(mpa_path)
        bacteria_reads = int(df_m.loc[df_m['clade_name'].str.startswith('k__Bacteria', na=False), 'estimated_number_of_reads_from_the_clade'].sum())
        archaea_reads  = int(df_m.loc[df_m['clade_name'].str.startswith('k__Archaea',  na=False), 'estimated_number_of_reads_from_the_clade'].sum())

        summary.update(dict(
            TOTAL_READS_MPA=total_reads,
            CLASSIFIED_READS_MPA=classified_reads,
            UNCLASSIFIED_READS_MPA=unclassified_reads,
            BACTERIAL_READS_MPA=bacteria_reads,
            ARCHAEAL_READS_MPA=archaea_reads,
        ))

    if not summary:
        eprint('WARNING: No classifier outputs found in classification-stats dir. Summary will only include filtering metrics.')

    return pd.Series(summary, dtype='Int64')


def write_normalized_tables(outdir: str, sample: str, classif_dir: str, primary_reads: int, verbose: bool = False):
    if primary_reads <= 0:
        eprint('WARNING: PRIMARY_READS <= 0; skipping normalization outputs.')
        return

    # Kraken
    k2_p = os.path.join(classif_dir, f'{sample}.paired.kraken.report.txt')
    k2_u = os.path.join(classif_dir, f'{sample}.unpaired.kraken.report.txt')
    if os.path.isfile(k2_p) or os.path.isfile(k2_u):
        try:
            dfk = merge_kraken_reports(k2_p, k2_u)
            df_rpm = dfk.copy()
            df_rpm['reads_clade_per_million'] = 1e6 * df_rpm['reads_clade'] / primary_reads
            df_rpm['reads_taxon_per_million'] = 1e6 * df_rpm['reads_taxon'] / primary_reads
            df_rpm = df_rpm[['name','tax_id','rank','reads_clade','reads_taxon','reads_clade_per_million','reads_taxon_per_million','pct_reads']]
            out_k2 = os.path.join(outdir, f'{sample}.kraken.txt')
            df_rpm.to_csv(out_k2, sep='\t', index=False)
            if verbose: eprint(f'[summarize] wrote {out_k2}')
        except Exception as ex:
            eprint(f'WARNING: Kraken2 normalization failed: {ex}')

    # MetaPhlAn
    mpa_path = os.path.join(classif_dir, f'{sample}.metaphlan.report.txt')
    if os.path.isfile(mpa_path):
        try:
            dfm = load_metaphlan_table(mpa_path)
            dfm_out = dfm.copy()
            dfm_out['estimated_number_of_reads_from_the_clade_per_million'] = 1e6 * dfm_out['estimated_number_of_reads_from_the_clade'] / primary_reads
            dfm_out = dfm_out[['clade_name','clade_taxid','estimated_number_of_reads_from_the_clade','estimated_number_of_reads_from_the_clade_per_million','relative_abundance','coverage']]
            out_mpa = os.path.join(outdir, f'{sample}.metaphlan.txt')
            dfm_out.to_csv(out_mpa, sep='\t', index=False)
            if verbose: eprint(f'[summarize] wrote {out_mpa}')
        except Exception as ex:
            eprint(f'WARNING: MetaPhlAn normalization failed: {ex}')


# ------------------------ Main ------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Combine filtering + classification summaries and write normalized classifier tables (RPM).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument('--filter-stats-dir', required=True, help='Directory with filter_stats outputs')
    p.add_argument('--classification-stats-dir', required=True, help='Directory with raw classifier reports')
    p.add_argument('--results-dir', required=True, help="Directory to write final outputs (new 'results')")
    p.add_argument('--sample-id', required=True, help='Sample ID (basename used by pipeline)')
    p.add_argument('-v','--verbose', action='store_true', help='Verbose logging to stderr')
    return p.parse_args()


def main():
    args = parse_args()
    require_dir(args.filter_stats_dir, 'filter-stats-dir')
    require_dir(args.classification_stats_dir, 'classification-stats-dir')
    os.makedirs(args.results_dir, exist_ok=True)

    sample = args.sample_id

    # Filtering summary (wide 1-row)
    filt_row = filtering_summary(args.filter_stats_dir, sample)
    primary_reads = int(filt_row.get('PRIMARY_READS', 0))

    # Classification summary (key totals)
    class_row = summarize_classification(args.classification_stats_dir, sample)

    # Merge → one wide row
    combined = pd.concat([filt_row, class_row]).fillna(0).astype(int)
    combined.name = sample

    # Write (a) combined filtering/classification summary
    out_summary = os.path.join(args.results_dir, f'{sample}.summary.tsv')
    combined.to_frame().to_csv(out_summary, sep='\t', index=True, index_label='sample_id')
    if args.verbose:
        eprint(f'[summarize] wrote {out_summary}')

    # Write (b) normalized Kraken/MetaPhlAn tables
    write_normalized_tables(args.results_dir, sample, args.classification_stats_dir, primary_reads, verbose=args.verbose)

if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as ex:
        eprint(f'FATAL: {ex}')
        sys.exit(1)
