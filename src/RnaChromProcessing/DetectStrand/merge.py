import shutil
from logging import getLogger
from os.path import getmtime
from pathlib import Path

import pandas as pd

from ..utils import (
    check_file_exists, exit_with_error
)


logger = getLogger('strand')


def merge_tables(inputs: list[Path], output: Path) -> None:
    # checking inputs
    for file in inputs:
        check_file_exists(file)

    output_dir = output.parent.absolute()
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    # if only one input file provided
    if len(inputs) == 1:
        logger.warning("Only one input file provided, copying to output..")
        shutil.copy(inputs[0], output)
        return
    
    # sort input files by time of mod
    inputs.sort(key=getmtime)

    # read earliest file, then add contents of other files to it
    data = pd.read_csv(inputs[0], sep='\t', index_col=(0, 1))
    colset = set(data.columns)

    files_to_merge = len(inputs) - 1
    for i, file in enumerate(inputs[1:]):
        logger.debug(f"Merging file {i} of {files_to_merge}")
        tab = pd.read_csv(file, sep='\t', index_col=(0, 1))

        if set(tab.columns) != colset:
            exit_with_error(f"Incompatible columns between {file} and {inputs[0]}!")

        for idx in tab.index:
            data.loc[idx, :] = tab.loc[idx, :]

    # save results
    data.to_csv(output, sep='\t')
    logger.info('Done.')
    