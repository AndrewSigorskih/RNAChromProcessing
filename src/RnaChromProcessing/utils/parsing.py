import re
from pathlib import Path

import pandas as pd

from .errors import exit_with_error
from .file_utils import check_file_exists


BEDRC_COLS = ('chr', 'start', 'end', 'name', 'strand', 'biotype', 'source')
GTF_COLS = ('chr', 'type', 'start', 'end', 'strand', 'attrs')

GENE_ID_PAT = re.compile(r'(?<=gene_id \")[^\"]+(?=\";)')
GENE_NAME_PAT = re.compile(r'(?<=gene_name \")[^\"]+(?=\";)')


# Universal functions

def load_config(path: str) -> dict:
    """loads config file in json or yaml format"""
    check_file_exists(path)
    with open(path, 'r') as f:
        if path.endswith('.json'):
            import json
            config = json.load(f)
        elif path.endswith('.yml') or path.endswith('.yaml'):
            import yaml
            config = yaml.safe_load(f)
        else:
            msg = f'Unknown config file extension: {path}. Only json and yml/yaml formats are supported!'
            exit_with_error(msg)
    return config


def read_bedrc(filepath: Path) -> pd.DataFrame:
    """Load gene annotation in bedrc format"""
    data = pd.read_csv(
        filepath, sep='\t', usecols=BEDRC_COLS
    )
    return data


def read_gtf(filepath: Path) -> pd.DataFrame:
    """Load gene annotation in gtf format"""
    skip = 0
    with open(filepath, 'r') as f:
        while next(f).startswith('#'):
            skip += 1
    data = pd.read_csv(
        filepath, sep='\t', header=None, skiprows=skip,
        usecols=[0,2,3,4,6,8], names=GTF_COLS
    )
    data = data[data['type'] == 'gene']
    return data


# Subset selections

def _read_name_or_id(s: str) -> str:
    # will throw exception if both attrs ase absent (invalid annotation)
    capture = re.search(GENE_NAME_PAT, s) or re.search(GENE_ID_PAT, s)
    return capture.group(0)


def fetch_genes_from_gtf(annotation: Path, genes_file: Path) -> pd.DataFrame:
    """Load gene annotation in gtf format, read a list of gene names OR ids in on-by-line format and
    select specified genes from annotation."""
    # load gene names
    with open(genes_file, 'r') as f:
        # will search for gene_name "DDX11L2"; OR gene_id "ENSG00000290825.1"; for all listed genes
        genes_list: str = '|'.join([f'"{line.strip()}";' for line in f])

    # load annotation and extract selected genes
    data = read_gtf(annotation)
    data = data[data['attrs'].str.contains(genes_list)]
    if not data.shape[0]:
        exit_with_error('Could not find any of the listed genes in provided annotation file!')

    # create a column containing gene name/id
    try:
        data['name'] = data['attrs'].apply(_read_name_or_id)
    except AttributeError:
        exit_with_error(
            'GTF annotation should have "gene_name" or "gene_id" attribute '
            'for "gene" type records!'
        )
    return data


def fetch_genes_from_bedrc(annotation: Path, genes_file: Path) -> pd.DataFrame:
    """Load gene annotation in bedrc format, read a list of gene names OR ids in on-by-line format and
    select specified genes from annotation."""
    # load gene names
    with open(genes_file, 'r') as f:
        genes_list = {line.strip() for line in f}

    # load annotation and extract selected genes
    data = read_bedrc(annotation)
    data = data[data['name'].isin(genes_list)]
    if not data.shape[0]:
        exit_with_error('Could not find any of the listed genes in provided annotation file!')
    return data
