# Copyright (C) 2024   Alvin Hsu
import argparse
from functools import lru_cache
import gzip
from pathlib import Path
import re
from typing import Any, Dict, Generator, List, Optional, Tuple

import numpy as np

from .aligner import PowtnrkaAligner, RepeatAligner
from .trimmer import PowtnrkaTrimmer, NullTrimmer


def make_parser() -> argparse.ArgumentParser:
    """
    Makes an ArgumentParser that processes arguments needed to run PowTNRka.
    """
    parser = argparse.ArgumentParser()
    # Files
    parser.add_argument('r1_fastq', type=str,
                        help='The fastq file to analyze.')
    parser.add_argument('out_dir', type=str,
                        help='Sets the output directory for a run of PowTNRka.')
    # Amplicons
    parser.add_argument('--amplicon', '-a', type=str, default=None,
                        help='Reference amplicon string in PowTNRka format.')
    parser.add_argument('--edited', '-e', type=str, default=None,
                        help='Edited amplicon string in PowTNRka format.')
    parser.add_argument('--amplicons_fasta', '-af', type=str, default=None,
                        help='If set, reads amplicons from this file instead of command line.')
    # Guides (for finding nick sites)
    parser.add_argument('--grna', '-g', action='append', default=list(),
                        help='If set, adds a gap incentive where this guide RNA targets. Can be'
                             'specified multiple times.')
    parser.add_argument('--window_size', '-w', type=int, default=10,
                        help='Window size for quantifying indels around nicks.')
    # Quality filtering
    parser.add_argument('--min_trimmed_quality', '-min_q', type=float, default=30,
                        help='Minimum average quality for a read post-trimming')
    parser.add_argument('--min_trimmed_length', '-min_l', type=float, default=50,
                        help='Minimum length for a read post-trimming')
    # Trimming reads
    parser.add_argument('--no_trim', '-no_trim', action='store_true',
                        help='Set this flag to disable read trimming based on Q score.')
    parser.add_argument('--trim_bp_window', '-tw', type=int, default=5,
                        help='Number of bases to compute average Q-score over to trim reads.')
    parser.add_argument('--trim_bp_quality', '-tq', type=float, default=20,
                        help='If the average Q-score falls below this number, the read is trimmed.')
    # Read alignment score
    parser.add_argument('--amplicon_min_alignment_score', '-amas', type=float, default=60,
                        help='Minimum % of amplicon that must be aligned to a read.')
    # Misc options
    parser.add_argument('--gz', '-gz', action='store_true',
                        help='Set this flag to indicate that r1_fastq is a gzipped file.')
    parser.add_argument('--revcomp', '-rc', action='store_true',
                        help='Set this flag to indicate that r1_fastq is reverse complemented '
                             'relative to the reference amplicons.')
    return parser


SPECIAL_RE = re.compile('\([\w|]+\)|\*')


@lru_cache
def make_repeat_aligner(reference: str):
    if '|' in reference:
        ref_str, aux_str = reference.split('|')
        return RepeatAligner(ref_str, aux_str=aux_str)
    return RepeatAligner(reference)


RC_TRANS = str.maketrans('ACGTMRWSYKBDHVN()-',
                         'tgcakywsrmvhdbn)(-')
def revcomp(s: str) -> str:
    return s.translate(RC_TRANS)[::-1].upper()


@lru_cache
def make_aligners(reference: str, grnas: Tuple[str], rc: bool) \
        -> Tuple[PowtnrkaAligner, List[RepeatAligner]]:
    if rc:
        reference = revcomp(reference)
    specials = SPECIAL_RE.findall(reference)
    parts = []
    num_repeats = 0
    repeats = []
    repeat_aligners = []
    while specials:
        special = specials.pop(0)
        idx = reference.find(special)
        parts.append(reference[:idx])
        if special == '*':
            parts.append('*')
        else:
            parts.append(str(num_repeats))
            repeats.append(special[1:-1])
            repeat_aligners.append(make_repeat_aligner(special[1:-1]))
            num_repeats += 1
        reference = reference[idx + len(special):]
    parts.append(reference)
    reference = ''.join(parts)
    nicks = []
    for grna in grnas:
        if (idx := reference.find(grna)) >= 0:
            nicks.append(idx + len(grna) - 3)
        elif (idx := reference.find(revcomp(grna))) >= 0:
            nicks.append(idx + 3)
    return PowtnrkaAligner(reference, num_repeats, repeats, nicks), repeat_aligners


def fasta_reader(fname: str) -> Generator[Tuple[str, str], None, None]:
    assert Path(fname).resolve().is_file(), f'File not found: {fname}'
    with open(fname, 'r') as f:
        name, seq = '', []
        for i, line in enumerate(f):
            if line[0] == '>':
                if i > 0:
                    yield name, ''.join(seq)
                name = line.strip()[1:]
                seq = []
            else:
                seq.append(line.strip())
        yield name, ''.join(seq)


def process_args(args) -> Dict[str, Any]:
    """
    Processes the args parsed by the parser returned by make_parser into a dict of kwargs for use
    by PowTNRka's main function.
    """
    # Files
    assert args.r1_fastq is not None, 'R1_fastq not specified'
    r1_path = Path(args.r1_fastq).resolve()
    assert r1_path.is_file(), f'R1_fastq file not found: {args.r1_fastq}'
    fastq_iter = fastq_reader(args.r1_fastq, gz=args.gz)
    assert args.out_dir is not None, 'out_dir not specified'
    out_path = Path(args.out_dir).resolve()
    out_path.mkdir(parents=True, exist_ok=True)
    out_prefix = out_path.name
    # Amplicons
    assert (args.amplicon is not None) or (args.amplicons_fasta is not None), \
        'Need to specify amplicon with either --amplicon or --amplicons_fasta'
    aligners = []
    grna = tuple(args.grna)
    if len(grna) == 1:
        grna = tuple(grna[0].split(','))
    if args.amplicons_fasta is not None:
        assert (args.amplicon is None) and (args.edited is None), \
            'Amplicons must be specified either in the command or in a file, but not both'
        assert Path(args.amplicons_fasta).resolve().is_file(), \
            f'Amplicons file not found: {args.amplicons_fasta}'
        for name, seq in fasta_reader(args.amplicons_fasta):
            print(name, seq)
            aligners.append((name, *make_aligners(seq, grna, args.revcomp)))
    else:
        aligners.append(('Reference', *make_aligners(args.amplicon, grna, args.revcomp)))
        if args.edited is not None:
            aligners.append(('Edited', *make_aligners(args.edited, grna, args.revcomp)))
    # Trimming reads
    if args.no_trim:
        trimmer = NullTrimmer()
    else:
        trimmer = PowtnrkaTrimmer(args.trim_bp_window, args.trim_bp_quality)
    # Return kwargs for run_powtnrka()
    return dict(fastq_iter=fastq_iter, out_dir=out_path, out_prefix=out_prefix, aligners=aligners,
                trimmer=trimmer, rc=args.revcomp, min_score=args.amplicon_min_alignment_score,
                min_qual=args.min_trimmed_quality, min_len=args.min_trimmed_length)


def fastq_reader(fname: str, gz=False) -> Generator[Tuple[str, np.ndarray], None, None]:
    """
    A Generator that yields tuples of (read, quality) where the quality has been converted from
    ASCII format to an np.ndarray of ints.
    """
    assert Path(fname).resolve().is_file(), f'File not found: {fname}'
    if gz:
        with gzip.open(fname, 'rb') as f:
            for i, line in enumerate(f):
                if i % 4 == 1:
                    read = line.decode().strip()
                elif i % 4 == 3:
                    qual = np.frombuffer(line.strip(), dtype=np.uint8) - 33
                    yield read, qual
    else:
        with open(fname, 'r') as f:
            for i, line in enumerate(f):
                if i % 4 == 1:
                    read = line.strip()
                elif i % 4 == 3:
                    qual = np.frombuffer(line.strip().encode('ascii'), dtype=np.uint8) - 33
                    yield read, qual


class FastqGZWriter(object):
    def __init__(self, path: Path):
        self.path = path
        self.file = None

    def __enter__(self):
        self.file = gzip.open(self.path, 'wb+')
        return self

    def write(self, seq: str, qual: Optional[np.ndarray] = None, title: Optional[str] = None):
        title_line = f'>{title}' if title is not None else '>'
        qual_line = (qual + 33).tobytes().decode('ascii') if qual is not None else len(seq)*'!'
        entry = f'{title_line}\n{seq}\n+\n{qual_line}\n'
        self.file.write(entry.encode('ascii'))

    def __exit__(self, *_, **__):
        self.file.close()
