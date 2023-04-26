import logging
import os
from tempfile import TemporaryDirectory
from typing import Any, Dict, List

from .default_configs import dedup_default_cfg, rsites_default_cfg, \
                    trim_default_cfg, hisat_default_cfg, contacts_default_cfg
from .stages import AlignedToBed, BamFilter, Contacts,\
                    Dedup, Hisat, Rsites, Trim, StatsCalc
from ..utils import exit_with_error, make_directory, move_exist_ok

logger = logging.getLogger()
SUBDIR_LIST = ('dedup', 'rsites', 'trim', 'hisat', 'bam', 'bed', 'contacts')

class BaseProcessor:
    def __init__(self, cfg: Dict[str, Any]):
        # get basic parameters
        self.cpus: int = cfg.get('cpus', 1)
        self.base_dir: str = os.path.abspath(cfg.get('base_dir', os.getcwd()))
        self.input_dir: str = cfg.get('input_dir', None)
        self.output_dir: str = cfg.get('output_dir', None)
        self.dna_ids: List[str] = cfg.get('dna_ids', [])
        self.rna_ids: List[str] = cfg.get('rna_ids', [])
        self.stats = cfg.get('stats', 'default')
        self.stats_prefix = cfg.get('stats_prefix', 'stats')
        self.keep: List[str] = cfg.get('keep', ['bam', 'contacts'])
        if 'contacts' not in self.keep:
            msg = f'contacts dir not found in "keep" array, resulting contacts files will not be saved!\n{self.keep=}'
            logger.warning(msg)
        # spam errors
        self.validate_inputs()
        # create working directory
        os.chdir(self.base_dir)
        self.work_dir = TemporaryDirectory(dir=self.base_dir)
        self.setup_dirs()
        # get stages-specific configs: get default and update
        self.stages_from_cfgs(cfg)

    def validate_inputs(self):
        """Basic input sanity check"""
        # input dir: check name
        if not self.input_dir:
            exit_with_error('Input directory not specified!')
        self.input_dir = os.path.abspath(self.input_dir)
        # output dir: check name and create if needed
        if not self.output_dir:
            exit_with_error('Output directory not specified!')
        self.output_dir = os.path.abspath(self.output_dir)
        if not os.path.exists(self.output_dir):
            make_directory(self.output_dir)
        # other important inputs
        if (not self.rna_ids) or (not self.dna_ids):
            exit_with_error('Input file ids not specified!')
        if len(self.rna_ids) != len(self.dna_ids):
            exit_with_error('DNA and RNA ids arrays have different length!')
        if not self.keep:
            exit_with_error('Empty keep list!')
    
    def setup_dirs(self):
        """Create dirs for pipeline stages and store names"""
        self.dedup_dir: str = os.path.join(self.work_dir.name, 'dedup')
        self.rsite_dir: str = os.path.join(self.work_dir.name, 'rsites')
        self.trim_dir : str = os.path.join(self.work_dir.name, 'trim')
        self.hisat_dir: str = os.path.join(self.work_dir.name, 'hisat')
        self.bam_dir: str = os.path.join(self.work_dir.name, 'bam')
        self.bed_dir: str = os.path.join(self.work_dir.name, 'bed')
        self.contacts_dir: str = os.path.join(self.work_dir.name, 'contacts')
        for dirname in (
            self.dedup_dir, self.rsite_dir, self.trim_dir,
            self.hisat_dir, self.bam_dir, self.bed_dir, self.contacts_dir
        ):
            make_directory(dirname)
    
    def stages_from_cfgs(self, cfg: Dict[str, Any]):
        """For each stage, create corresponding routine object 
        based on default config and user input"""
        rsites_default_cfg.update(cfg.get('rsites', {}))
        self.rsitefilter: Rsites = Rsites(rsites_default_cfg,
                                          self.input_dir, self.rsite_dir,
                                          self.cpus)
        dedup_default_cfg.update(cfg.get('dedup', {}))
        self.dupremover: Dedup = Dedup(dedup_default_cfg,
                                       self.rsite_dir, self.dedup_dir,
                                       self.cpus)
        trim_default_cfg.update(cfg.get('trim', {}))
        self.trimmer: Trim = Trim(trim_default_cfg,
                                  self.dedup_dir, self.trim_dir,
                                  self.cpus)
        hisat_default_cfg.update(cfg.get('hisat', {}))
        self.aligner: Hisat = Hisat(hisat_default_cfg,
                                    self.trim_dir, self.hisat_dir,
                                    self.cpus)
        self.bamfilter = BamFilter(self.hisat_dir, self.bam_dir,
                                   self.cpus)
        self.bamtobed = AlignedToBed(self.bam_dir, self.bed_dir,
                                     self.cpus)
        contacts_default_cfg.update(cfg.get('contacts', {}))
        self.contactsmerger = Contacts(contacts_default_cfg,
                                       self.bed_dir, self.contacts_dir,
                                       self.cpus)
        
    def save_outputs(self):
        """copy everything needed to out dir"""
        for to_copy in self.keep:
            if to_copy not in SUBDIR_LIST:
                logger.warning(f'Unknown directory to copy: {to_copy}. Skipping..')
                continue
            source_pth = os.path.join(self.work_dir.name, to_copy)
            dest_pth = os.path.join(self.output_dir, to_copy)
            make_directory(dest_pth)
            move_exist_ok(source_pth, dest_pth)

    def run(self):
        """Iteratively run all stages of pipeline and 
        retrieve data"""
        # run pipeline
        logger.info(f'Started processing of {len(self.rna_ids)} pairs of files.')
        os.chdir(self.work_dir.name)
        self.rsitefilter.run(self.dna_ids, self.rna_ids)
        self.dupremover.run(self.dna_ids, self.rna_ids)
        self.trimmer.run(self.dna_ids, self.rna_ids)
        self.aligner.run(self.dna_ids, self.rna_ids)
        self.bamfilter.run(self.dna_ids, self.rna_ids)
        self.bamtobed.run(self.dna_ids, self.rna_ids)
        self.contactsmerger.run(self.dna_ids, self.rna_ids)
        # calculate stats
        StatsCalc(self.output_dir, self.cpus, self.stats,
                  self.stats_prefix, self.dna_ids, self.rna_ids).run()
        # copy outputs
        os.chdir(self.base_dir)
        self.save_outputs()
        logger.info("Done.")
