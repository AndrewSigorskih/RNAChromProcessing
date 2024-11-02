from dataclasses import dataclass
from typing import List, Literal, TextIO

import numpy as np
import pandas as pd
from pydantic import PositiveInt

from .basicstage import BasicStage, SamplePair
from ...utils import run_command

DNA_COLUMNS = {
    'dna_chr': str, 'dna_start': np.uint32, 'dna_end': np.uint32,
    'id': str, 'dna_score': np.uint16, 'dna_strand': str, 'dna_cigar': str
}

RNA_COLUMNS = {
    'rna_chr': str, 'rna_start': np.uint32, 'rna_end': np.uint32,
    'id': str, 'rna_score': np.uint16, 'rna_strand': str, 'rna_cigar': str
}

TAB_HEADER = (
    'rna_chr', 'rna_start', 'rna_end', 'id', 'rna_strand', 'rna_cigar',
    'dna_chr', 'dna_start', 'dna_end', 'dna_strand', 'dna_cigar'
)


def _process_id(read_id: str) -> str:
    """remove 'SRRXXXXX' part from read ID if exists"""
    if '.' in read_id:
        return read_id.split('.', maxsplit=1)[1]
    return read_id


@dataclass
class BedRow:
    chr: str
    start: str
    end: str
    id: str
    score: str
    strand: str
    cigar: str

    @classmethod
    def from_string(cls, row: str):
        return cls(*row.strip().split())._prepare_id()
    
    @classmethod
    def from_file(cls, file: TextIO):
        return cls.from_string(next(file))

    def _prepare_id(self):
        self.id = _process_id(self.id)
        return self

    def to_rna(self) -> str:
        return '\t'.join(
            (self.chr, self.start, self.end, self.id, self.strand, self.cigar)
        )
    
    def to_dna(self) -> str:
        return '\t'.join(
            (self.chr, self.start, self.end, self.strand, self.cigar)
        )


class Contacts(BasicStage):
    mode: Literal['fast', 'chunks', 'low-mem'] = 'low-mem'
    chunksize: PositiveInt = 10_000_000

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
    
    def run(self, samples: List[SamplePair]) -> List[SamplePair]:
        """make tab-separated contacts file from two bed files"""
        # choose function to run
        if self.mode == 'fast':
            func = self._make_contacts_fast
        elif self.mode == 'chunks':
            func = self._make_contacts_chunks
        elif self.mode == 'low-mem':
            func = self._make_contacts_iter
        # prepare filepaths
        output_samples = self._make_output_samples(samples, new_suff='tab')
        # run function
        self.run_function(func, samples, output_samples)
        # return results
        return output_samples

    def _make_contacts_fast(self,
                            inp_sample: SamplePair,
                            out_sample: SamplePair) -> int:
        """read 2 bed files in memory, produce 1 combined contacts file"""
        dna = pd.read_csv(inp_sample.dna_file, sep='\t', header=None,
                          names=DNA_COLUMNS.keys(), dtype=DNA_COLUMNS)
        rna = pd.read_csv(inp_sample.rna_file, sep='\t', header=None,
                          names=RNA_COLUMNS.keys(), dtype=RNA_COLUMNS)
        dna['id'] = dna['id'].apply(_process_id)
        rna['id'] = rna['id'].apply(_process_id)
        rna = pd.merge(rna, dna, on='id', how='inner')
        rna = rna.drop(['rna_score', 'dna_score'], axis=1)
        rna.to_csv(out_sample.rna_file, sep='\t', index=False, header=True)
        return 0
    
    def _make_contacts_chunks(self,
                              inp_sample: SamplePair,
                              out_sample: SamplePair) -> int:
        # https://stackoverflow.com/questions/58441517/merging-dataframe-chunks-in-pandas
        """read rna bed file in a memory-efficient way, produce 1 combined contacts file"""
         # prepare names for output and temporal files
        dna = pd.read_csv(inp_sample.dna_file, sep='\t', header=None,
                          names=DNA_COLUMNS.keys(), dtype=DNA_COLUMNS)
        dna['id'] = dna['id'].apply(_process_id)
        result_chunks = []
        for chunk in pd.read_csv(inp_sample.rna_file, sep='\t', header=None,
                                 names=RNA_COLUMNS.keys(), dtype=RNA_COLUMNS,
                                 chunksize=self.chunksize):
            chunk['id'] = chunk['id'].apply(_process_id)
            result_chunks.append(pd.merge(dna, chunk, on='id', how='inner'))
        result = pd.concat(result_chunks)
        result = result.drop(['rna_score', 'dna_score'], axis=1)
        result.to_csv(out_sample.rna_file, sep='\t', index=False, header=True)
        return 0
    
    def _make_contacts_iter(self,
                            inp_sample: SamplePair,
                            out_sample: SamplePair) -> int:
        """read 2 bed files line-by-line, produce 1 combined contacts file"""
        # sort bed files by id
        rna_sorted = out_sample.rna_file.with_suffix('.tmp')
        dna_sorted = out_sample.dna_file.with_suffix('.tmp')
        for in_file, sorted_file in zip((inp_sample.dna_file, inp_sample.rna_file),
                                        (dna_sorted, rna_sorted)):
            cmd = f'sort -k 4 {in_file} > {sorted_file}'
            retcode = run_command(cmd, shell=True)
            if retcode: return retcode
        # process sorted files line by line, seek matching ids
        with open(dna_sorted, 'r') as dna_in, \
             open(rna_sorted, 'r') as rna_in, \
             open(out_sample.rna_file, 'w') as rna_out:
            print(*TAB_HEADER, file=rna_out, sep='\t')
            try:
                while True:
                    rna_row, dna_row = BedRow.from_file(rna_in), BedRow.from_file(dna_in)
                    while (rna_row.id != dna_row.id):
                        if rna_row.id < dna_row.id:
                            rna_row = BedRow.from_file(rna_in)
                        elif rna_row.id > dna_row.id:
                            dna_row = BedRow.from_file(dna_in)
                    print(rna_row.to_rna(), dna_row.to_dna(), file=rna_out, sep='\t')
            except StopIteration:  # one or both files ended
                pass
        rna_sorted.unlink()
        dna_sorted.unlink()
        return 0
