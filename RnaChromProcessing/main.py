import argparse
import json

from .DataProcessors import BaseProcessor
from .utils import check_file_exists

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('-c', '--config', required=True,
                        help='Configuration file.')
    parser.add_argument('-mode', '--run_mode', required=False,
                        choices=('contacts', 'XRNA'), default='contacts',
                        help='''Action to be performed during run:
* contacts: build contacts data from provided raw FASTQ files
* XRNA: attempt to infer X-RNAs from provided contacts and BAM files.''')
    return parser.parse_args()

class Program:
    def __init__(self) -> None:
        args: argparse.Namespace = parse_args()
        self.mode: str = args.mode
        check_file_exists(args.config)
        with open(args.config, 'r') as f:
            self.config: dict = json.load(f)

    def run(self) -> None:
        if self.mode == 'XRNA':
            raise NotImplementedError("Warning: XRNA collection is not implemented yet!")
            exit() # for when it will be implemented
        # we are collecting contacts
        # actually run base class or any of its pre-defined subclasses here
        # depending on config
        BaseProcessor(self.config).run()


def main() -> None:
    #TODO: actually set up logging
    Program().run()
    

if __name__ == '__main__':
    main()
