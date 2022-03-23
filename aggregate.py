# Copyright (C) 2024   Alvin Hsu
from pathlib import Path

import pandas as pd


def main(args):
    out_path = Path(args.output_dir).resolve()
    out_name = out_path.name
    with pd.ExcelWriter(out_path / f'{out_name}_SUMMARY.xlsx') as xlsx_writer:
        summary_i = 0
        for row_dir in out_path.glob('*'):
            name = row_dir.name
            summary_files = list(row_dir.glob(f'{name}_*_SUMMARY.xlsx'))
            if len(summary_files) == 0:
                continue
            row_df = pd.DataFrame({'name': [name]})
            row_df.to_excel(xlsx_writer, sheet_name='AGGREGATE', header=False, index=False,
                            startrow=summary_i)
            summary_i += 1
            for p in summary_files:
                row_df = pd.read_excel(p, header=None)
                old_cols = list(row_df.columns)
                row_df['directory'] = name
                row_df['amplicon'] = p.stem
                reordered = ['directory', 'amplicon'] + old_cols
                row_df = row_df[reordered]
                row_df.to_excel(xlsx_writer, sheet_name='AGGREGATE', header=False, index=False,
                                startrow=summary_i)
                summary_i += len(row_df)


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('output_dir', type=str,
                        help='Parent directory of individual output directories')
    args = parser.parse_args()
    main(args)
