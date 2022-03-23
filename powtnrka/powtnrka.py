# Copyright (C) 2024   Alvin Hsu
from collections import defaultdict
import csv
from functools import lru_cache
from pathlib import Path
import re
from typing import Iterable, List, Tuple

import numpy as np
import pandas as pd
from pandas.api.types import is_integer_dtype as pd_isint

from .aligner import PowtnrkaAligner, RepeatAligner
from .trimmer import PowtnrkaTrimmer
from .utils import revcomp, FastqGZWriter


def process_heterogeneous(qry_str: str, rc: bool):
    if rc:
        qry_str = revcomp(qry_str)
    return ['het_seq', 'het_len'], {'het_seq': qry_str, 'het_len': len(qry_str)}


IUPAC = {'A': list('A'), 'C': list('C'), 'G': list('G'), 'T': list('T'),
         'M': list('AC'), 'R': list('AG'), 'W': list('AT'),
         'S': list('CG'), 'Y': list('CT'), 'K': list('GT'),
         'V': list('ACG'), 'H': list('ACT'), 'D': list('AGT'), 'B': list('CGT'),
         'N': list('ACGT')}
@lru_cache
def iupac_expand(iupac: str) -> List[str]:
    if len(iupac) == 1:
        return IUPAC[iupac]
    else:
        return [f'{base}{rest}' for base in IUPAC[iupac[0]] for rest in iupac_expand(iupac[1:])]


def process_repeats(repeat_tract: str, aligner: RepeatAligner, rc: bool):
    repeat_seq = aligner.ref_str
    repeat_list = aligner(repeat_tract)
    possible = iupac_expand(repeat_seq)
    counts = {k: 0 for k in possible}
    counts['other'] = defaultdict(lambda: 0)
    if using_aux := (aligner.aux_str is not None):
        counts[f'non_{aligner.aux_str}'] = 0
        aux_options = iupac_expand(aligner.aux_str)
    rep_len = len(repeat_seq)
    indel = False
    for i, rep in enumerate(repeat_list, 1):
        if (len(rep) != rep_len) and (i < len(repeat_list)):
            indel = True
        if rep in counts:
            counts[rep] += 1
            if using_aux:
                # noinspection PyUnboundLocalVariable
                counts[f'non_{aligner.aux_str}'] += int(rep not in aux_options)
        else:
            counts['other'][rep] += 1
    counts['other'] = dict(counts['other'])
    # Add to results
    if not rc:
        retval = {k: counts[k] for k in possible}
        if using_aux:
            aux_label = f'non_{aligner.aux_str}_total'
            retval[aux_label] = counts[f'non_{aligner.aux_str}']
        retval['other_total'] = sum(list(counts['other'].values()))
        retval['other_counts'] = counts['other']
        retval['indel'] = indel
        retval['repeat_list'] = repeat_list
    else:
        retval = {revcomp(k): counts[k] for k in possible}
        if using_aux:
            aux_label = f'non_{revcomp(aligner.aux_str)}_total'
            retval[aux_label] = counts[f'non_{aligner.aux_str}']
        retval['other_total'] = sum(list(counts['other'].values()))
        retval['other_counts'] = {revcomp(k): v for k, v in counts['other'].items()}
        retval['indel'] = indel
        retval['repeat_list'] = [revcomp(x) for x in repeat_list[::-1]]
        possible = [revcomp(x) for x in possible]
    columns = possible
    if using_aux:
        # noinspection PyUnboundLocalVariable
        columns = columns + [aux_label]
    columns = columns + ['other_total', 'other_counts', 'indel', 'repeat_list']
    return columns, retval


def process_constant(ref_str: str, qry_str: str, nick_windows: List[Tuple[int, int]], rc: bool):
    columns = ['alignment', 'any_indel', 'nick_indel', 'read_end']
    any_indel = False
    nick_indel = False
    q_end = len(qry_str.rstrip('-'))
    r_idx = 0
    alignment = []
    for i, (r, q) in enumerate(zip(ref_str, qry_str)):
        r_idx += int(r != '-')  # Increment r_idx if it is not a gap character
        if r == q:
            alignment.append('|')
        elif (r == '-' or q == '-') and i < q_end:
            any_indel = True
            if not nick_indel:
                for start, end in nick_windows:
                    if start <= r_idx < end:
                        nick_indel = True
            alignment.append(' ')
        else:
            alignment.append('.')
    align_str = ''.join(alignment)
    alignment = (f'{ref_str}\n{align_str}\n{qry_str}' if not rc else
                 f'{revcomp(ref_str)}\n{align_str[::-1]}\n{revcomp(qry_str)}')
    retval = {'alignment': alignment,
              'any_indel': any_indel,
              'nick_indel': nick_indel,
              'read_end': (q_end < len(qry_str))}
    return columns, retval


REF_RE = re.compile('([ACGTMRWSYKVHDBN-]+|(?P<rep>\d|\*)(?P=rep)*)')
def process_alignment(alignment: str, powtnrka_aligner: PowtnrkaAligner,
                      repeat_aligners: List[RepeatAligner], window: int, rc: bool):
    ref_aln, qry_aln = alignment.split('\n')
    # Trim NNNN bases for sequencing from beginning
    ref_aln = ref_aln.lstrip('-')
    qry_aln = qry_aln[-len(ref_aln):]
    columns = []
    results = []
    grouped = []
    ref_idx = 0
    segment_idx = 0
    while ref_aln:
        match = REF_RE.match(ref_aln)
        end = match.span()[1]
        ref_str, ref_aln = ref_aln[:end], ref_aln[end:]
        qry_str, qry_aln = qry_aln[:end], qry_aln[end:]
        ref_0 = ref_str[0]
        if ref_0 == '*':
            col, res = process_heterogeneous(qry_str, rc)
            ref_idx += 1
        elif '0' <= ref_0 <= '9':
            rep_idx = int(ref_0)
            col, res = process_repeats(qry_str, repeat_aligners[rep_idx], rc)
            ref_idx += 1
        else:
            segment_len = len(ref_str.replace('-', ''))
            nick_sites = [x - ref_idx for x in powtnrka_aligner.nick_sites]
            nick_windows = [(x - window, x + window) for x in nick_sites
                            if (x + window >= 0) and (x - window < segment_len)]
            col, res = process_constant(ref_str, qry_str, nick_windows, rc)
            ref_idx += segment_len
        if not rc:
            columns = columns + [f'{segment_idx}_{x}' for x in col]
            results = results + [res[x] for x in col]
            grouped.append(res)
        else:
            columns = [f'{segment_idx}_{x}' for x in col] + columns
            results = [res[x] for x in col] + results
            grouped.append(res)
        segment_idx += 1
    return columns, results, grouped


def run_powtnrka(fastq_iter: Iterable[Tuple[str, np.ndarray]], out_dir: Path, out_prefix: str,
                 aligners: Iterable[Tuple[str, PowtnrkaAligner, List[RepeatAligner]]],
                 trimmer: PowtnrkaTrimmer, rc: bool, min_score: float = 60.0,
                 min_qual: float = 30.0, min_len: int = 50, window: int = 10):
    # Prepare output paths
    trimmed_path = out_dir / 'trimmed_reads.fastq.gz'
    low_qual_path = out_dir / 'low_quality_reads.fastq.gz'
    filtered_path = out_dir / 'reads_passing_filter.csv'
    unaligned_path = out_dir / 'unaligned_reads.fastq.gz'
    # Setup for processing
    powtnrka_aligners = {x0: x1 for x0, x1, _ in aligners}
    repeat_aligners = {x0: x2 for x0, _, x2 in aligners}
    full_results = {x0: [] for x0, _, _ in aligners}
    full_headers = {x0: None for x0, _, _ in aligners}
    full_grouped = {x0: [] for x0, _, _ in aligners}
    total_aligned = 0
    with FastqGZWriter(trimmed_path) as tf, \
         FastqGZWriter(low_qual_path) as qf, \
         filtered_path.open('w+') as ff, \
         FastqGZWriter(unaligned_path) as uf:
        ff_header = [x for a, _, _ in aligners for x in [f'{a}_score', f'{a}_alignment']]
        ff_header = ('Read #', 'Aligned?', 'Assigned') + tuple(ff_header)
        ff = csv.writer(ff)
        ff.writerow(ff_header)
        for i, (u_read, u_qual) in enumerate(fastq_iter):
            t_read, t_qual = trimmer(u_read, u_qual)
            if t_qual.mean() < min_qual or len(t_read) < min_len:
                tf.write(t_read, t_qual, title=f'{i} (Q discarded)')
                qf.write(t_read, t_qual, title=f'{i} (Q discarded)')
                continue
            read_aligns = {}
            read_scores = {}
            unaligned = True
            # Initial candidate alignments
            for a_name, aligner, _ in aligners:
                read_aligns[a_name], read_scores[a_name] = aligner(t_read)
                if read_scores[a_name] >= min_score:
                    unaligned = False
            best_a_name = max(read_scores, key=lambda x: read_scores.get(x)) \
                          if not unaligned else 'UNALIGNED'
            ff_row = [i, (not unaligned), best_a_name]
            for a, _, _ in aligners:
                ff_row.append(read_aligns[a])
                ff_row.append(read_scores[a])
            ff.writerow(tuple(ff_row))
            if unaligned:
                tf.write(t_read, t_qual, title=f'{i} (unaligned)')
                uf.write(t_read, t_qual, title=f'{i} (unaligned)')
                continue
            total_aligned += 1
            tf.write(t_read, t_qual, title=f'{i} ({best_a_name})')
            # Polish the best alignment
            columns, results, grouped = process_alignment(read_aligns[best_a_name],
                                                          powtnrka_aligners[best_a_name],
                                                          repeat_aligners[best_a_name],
                                                          window,
                                                          rc)
            if full_headers[best_a_name] is None:
                full_headers[best_a_name] = ['read_id', 'score', 'trimmed', 'raw_align'] + columns
            result_row = [i, read_scores[best_a_name], t_read, read_aligns[best_a_name]] + results
            full_results[best_a_name].append(result_row)
            full_grouped[best_a_name].append(grouped)

    reads_processed = i

    # Write files for each amplicon
    for a_name, columns in full_headers.items():
        df = pd.DataFrame(full_results[a_name], columns=columns).set_index('read_id', drop=True)
        df.to_csv(out_dir / f'{out_prefix}_{a_name}_FULL_alignments.csv')
        grouped = full_grouped[a_name]
        grouped_cols = [x.split('_', maxsplit=1) for x in columns]
        with pd.ExcelWriter(out_dir / f'{out_prefix}_{a_name}_SUMMARY.xlsx') as s_writer, \
                pd.ExcelWriter(out_dir / f'{out_prefix}_{a_name}_DETAILED.xlsx') as d_writer:
            df = pd.DataFrame([('SUMMARY',)], columns=[''])
            df.to_excel(s_writer, sheet_name='SUMMARY', index=False, header=False)
            summary_col = 3
            any_indel = 0   # Gets typecast to pd.Series
            nick_indel = 0  # Gets typecast to pd.Series
            for segment_i in range(len(grouped[0])):
                sub_data = [x[segment_i] for x in grouped]
                sub_columns = [x[1] for x in grouped_cols if (len(x) == 2) and
                                                             (x[0] == str(segment_i))]
                df = pd.DataFrame(sub_data)[sub_columns]
                df.to_excel(d_writer, sheet_name=f'{segment_i}')
                if 'any_indel' in df.columns:  # Constant
                    any_indel = any_indel | df['any_indel']
                    nick_indel = nick_indel | df['nick_indel']
                    s_any, s_nick = df['any_indel'].sum(), df['nick_indel'].sum()
                    p_any, p_nick = 100 * s_any / total_aligned, 100 * s_nick / total_aligned
                    summary_df = pd.DataFrame(data={'CONSTANT': ['# Reads', '% Total'],
                                                    'any_indel': [s_any, p_any],
                                                    'nick_indel': [s_nick, p_nick]})
                    summary_df = summary_df.set_index('CONSTANT', drop=True)
                    summary_df.to_excel(s_writer, sheet_name='SUMMARY', startcol=summary_col)
                    summary_col += 3
                elif 'indel' in df.columns:  # Repeat
                    any_indel = any_indel | df['indel']
                    nick_indel = any_indel | df['indel']
                    s_indel = df['indel'].sum()
                    p_indel = 100 * s_indel / total_aligned
                    summary_df = pd.DataFrame(data={'REPEAT': ['# Reads', '% Total'],
                                                    'indel': [s_indel, p_indel]})
                    summary_df = summary_df.set_index('REPEAT', drop=True)
                    summary_df.to_excel(s_writer, sheet_name='SUMMARY', startcol=summary_col)
                    summary_col += 2
                    repeat_cols = {x: df[x].values for x in df.columns if pd_isint(df[x].dtype)}
                    summary_df = pd.DataFrame(data=repeat_cols)
                    max_repeat_len = summary_df.sum(axis=1).max()
                    bins = [0, 1] + list(range(10, max_repeat_len + 10, 10))
                    bin_labels = (['[0, 0]', '[1, 10)'] +
                                  [f'[{i}, {i + 10})' for i in range(10, max_repeat_len, 10)])
                    repeat_hist = {k: np.histogram(v, bins=bins)[0] for k, v in repeat_cols.items()}
                    summary_df = pd.DataFrame(data=repeat_hist)
                    summary_df['bins'] = [f'# {x}' for x in bin_labels]
                    summary_df = summary_df.set_index('bins', drop=True)
                    summary_df.to_excel(s_writer, sheet_name='SUMMARY', startcol=summary_col)
                    summary_df = 100 * summary_df / total_aligned
                    summary_df['bins'] = [f'% {x}' for x in bin_labels]
                    summary_df = summary_df.set_index('bins', drop=True)
                    summary_df.to_excel(s_writer, sheet_name='SUMMARY', startcol=summary_col,
                                        startrow=len(summary_df) + 1)
                    summary_col += len(summary_df.columns) + 1
                    bins = list(range(0, max_repeat_len + 1))
                    repeat_hist = {k: np.histogram(v, bins=bins)[0] for k, v in repeat_cols.items()}
                    summary_df = pd.DataFrame(data=repeat_hist)
                    summary_df['bins'] = [f'# {x}' for x in bins[:-1]]
                    summary_df = summary_df.set_index('bins', drop=True)
                    summary_df.to_excel(s_writer, sheet_name=f'{segment_i}_histogram')
                else:  # Heterogeneous
                    summary_df = pd.DataFrame(data={'HET': ['mean_len'],
                                                    'heterogeneous': df['het_len'].mean()})
                    summary_df = summary_df.set_index('HET', drop=True)
                    summary_df.to_excel(s_writer, sheet_name='SUMMARY', startcol=summary_col)
                    summary_col += 2
            data = [('reads_processed', reads_processed),
                    ('total_aligned', total_aligned),
                    (f'{a_name}_aligned', len(df)),
                    (f'{a_name}_any_indel', any_indel.sum()),
                    (f'{a_name}_nick_indel', nick_indel.sum()),
                    (f'{a_name}_%_total', 100 * len(df) / total_aligned),
                    (f'{a_name}_%_indel', 100 * any_indel.sum() / total_aligned)]
            summary_df = pd.DataFrame(data, columns=['Label', 'Value'])
            summary_df.to_excel(s_writer, sheet_name='SUMMARY', index=False, header=False)
