#!/usr/bin/env python

import sys
import hashlib
import argparse

import pandas as pd

from pathlib import Path



if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Generate a hash from a primer-centric BED file'
    )
    parser.add_argument('primer_file', type=Path)
    args = parser.parse_args()

    bed_fields = ['chrom', 'chromStart', 'chromEnd', 'name', 'poolName', 'orientation', 'seq']
    fields_to_hash = ['chromStart', 'chromEnd', 'orientation', 'seq']

    df = pd.read_csv(
        args.primer_file,
        sep='\t',
        names=bed_fields,
        dtype=dict(
            chrom=str,
            chromStart=int,
            chromEnd=int,
            name=str,
            poolName=str,
            orientation=str,
            seq=str
        )
    )

    # Normalisation
    df_norm = df[[*fields_to_hash]].sort_values('chromStart')  # Sort by start pos
    df_norm['seq'] = df_norm['seq'].str.upper()  # Normalise case
    df_norm = df_norm.applymap(lambda x: x.strip() if isinstance(x, str) else x)  # Strip trailing and leading whitespace
    df_norm_text = df_norm.to_csv(sep='\t', header=False, index=False)  

    hexdigest = hashlib.md5(df_norm_text.encode()).hexdigest()

    print('Primer scheme hash:', file=sys.stderr)
    print(hexdigest[:16])  # Use first 64 bits only
    print('Hash function input:', df_norm_text.strip(), sep='\n', file=sys.stderr)
