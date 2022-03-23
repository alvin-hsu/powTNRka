#!/broad/liulabdata/Alvin_Hsu/powTNRka/env/bin/python
# Copyright (C) 2024   Alvin Hsu
import argparse
from datetime import datetime
from pathlib import Path
from subprocess import Popen, PIPE

from jinja2 import Template
import pandas as pd

_SRC_PATH = Path(__file__).resolve().parent
_COL_MAP = {'r1_fastq': 'r1',
            'amplicon': 'a',
            'edited': 'e',
            'amplicons_fasta': 'af',
            'grna': 'g',
            'window_size': 'w',
            'min_trimmed_quality': 'min_q',
            'min_trimmed_length': 'min_l',
            'trim_bp_window': 'tw',
            'trim_bp_quality': 'tq',
            'amplicon_min_alignment_score': 'amas',
            'revcomp': 'rc'}
_SPECIALS = ['r1', 'a', 'e', 'rc', 'gz', 'no_trim']


def _bash_escape(s):
    return "'" + s.replace("'", "'\''") + "'"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('batch_file', type=str,
                        help='File containing list of fastq names')
    parser.add_argument('--email', type=str, default=None,
                        help='Optional: Sends an email to this email address '
                             'when the job has completed.')
    args = parser.parse_args()
    # Validate batch file input
    batch_path = Path(args.batch_file).resolve()
    batch_name = batch_path.stem
    batch_dir = batch_path.parent
    assert batch_path.is_file()
    batch_df = pd.read_csv(batch_path, sep='\t')
    batch_df = batch_df.rename(columns=_COL_MAP)
    df_cols = batch_df.columns
    # Make sure fastq files are specified and present
    assert 'r1' in df_cols, 'Missing `r1` column in batch file'
    for fname in batch_df['r1']:
        assert (batch_dir / fname).is_file(), f'File not found: {fname}'
    # Make sure amplicons are specified
    if ('af' in df_cols) or ('amplicons_fasta' in df_cols):
        assert 'a' not in df_cols, 'amplicon column specified but using file reference'
        assert 'e' not in df_cols, 'edited column specified but using file reference'
    else:
        assert 'a' in df_cols, 'Not using file reference but no reference amplicon specified'
    # Make sure boolean arguments are boolean
    if 'no_trim' in df_cols:
        assert batch_df['no_trim'].dtype == bool, '`no_trim` column must be TRUE or FALSE'
    if 'rc' in df_cols:
        assert batch_df['rc'].dtype == bool, '`rc` column must be TRUE or FALSE'
    # Set up directory structure
    powtnrka_dir = batch_dir / f'powTNRka_on_{batch_name}'
    if powtnrka_dir.is_dir():
        old_time = datetime.fromtimestamp(powtnrka_dir.stat().st_mtime).strftime('%y%m%d_%H%M%S')
        old_dir = batch_dir / f'powTNRka_on_{batch_name}_OLD_{old_time}'
        powtnrka_dir.rename(old_dir)
    powtnrka_dir.mkdir(parents=True)
    job_dir = powtnrka_dir / 'jobs'
    job_dir.mkdir()
    # Add required fields to DataFrame
    out_dirs = []
    for i, (_, row) in enumerate(batch_df.iterrows(), 1):
        r1_path = Path(row['r1'])
        name = r1_path.name.split('.')[0]
        out_dirs.append(str(powtnrka_dir / f'{i}_{name}'))
    batch_df['out_dir'] = out_dirs
    batch_df['gz'] = (batch_df['r1'].str.slice(-2) == 'gz')
    if 'rc' not in df_cols and args.rc:
        batch_df['rc'] = True
    # Set up job files and submit
    with (_SRC_PATH / 'powTNRka.job').open('r') as f:
        job_template = Template(f.read())
    job_ids = []
    for i, row in batch_df.iterrows():
        job_command = [str(_SRC_PATH / 'env' / 'bin' / 'python'),
                       str(_SRC_PATH / 'main.py'),
                       _bash_escape(str(batch_dir / row['r1'])),
                       _bash_escape(row['out_dir'])]
        if 'no_trim' in df_cols and row['no_trim']:
            job_command.append('--no_trim')
        if row['gz']:
            job_command.append('-gz')
        if 'rc' in df_cols and row['rc']:
            job_command.append('-rc')
        for col in df_cols:
            if col in _SPECIALS:
                continue
            job_command.append(f'-{col}')
            job_command.append(_bash_escape(str(row[col])))
        sample_name = Path(row['r1']).stem
        job_name = f'{i + 1}_{sample_name}'
        with (job_dir / f'{job_name}.job').open('w+') as f:
            f.write(job_template.render(work_dir=str(batch_dir),
                                        job_dir=str(job_dir),
                                        name=job_name,
                                        command=' '.join(job_command)))
        retval = -1
        while retval != 0:
            proc = Popen(['qsub', '-terse', str(job_dir / f'{job_name}.job')], stdout=PIPE)
            retval = proc.wait()
        # noinspection PyUnboundLocalVariable
        job_id, _ = proc.communicate()
        job_id = int(job_id.strip())
        job_ids.append(job_id)
        print(f'Submitted powTNRka job {job_name}.job ({job_id})')
    # Set up aggregate file and submit
    job_ids = [str(x) for x in job_ids]
    agg_name = str(job_dir / 'powTNRka_aggregate')
    email_flags = f'#$ -m e\n#$ -M {args.email}' if args.email is not None else ''
    agg_command = [str(_SRC_PATH / 'env' / 'bin' / 'python'),
                   str(_SRC_PATH / 'aggregate.py'),
                   _bash_escape(str(powtnrka_dir))]
    with (_SRC_PATH / 'aggregate.job').open('r') as f:
        agg_template = Template(f.read())
    with open(f'{agg_name}.job', 'w+') as f:
        f.write(agg_template.render(job_ids=job_ids,
                                    work_dir=str(batch_dir),
                                    job_dir=str(job_dir),
                                    email_flags=email_flags,
                                    command=' '.join(agg_command)))
    retval = -1
    while retval != 0:
        proc = Popen(['qsub', '-terse', f'{agg_name}.job'], stdout=PIPE)
        retval = proc.wait()
    job_id, _ = proc.communicate()
    job_id = job_id.decode().strip()
    print(f'Submitted aggregation job {agg_name}.job ({job_id})')


if __name__ == '__main__':
    main()
