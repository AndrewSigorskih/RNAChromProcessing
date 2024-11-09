import argparse
import logging
from pathlib import Path

from .DetectStrand import DetectStrand, merge_tables
from .utils import configure_logger, load_config

logger = logging.getLogger('strand')


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('-c', '--config', required=True,
                        help='Configuration file.')
    parser.add_argument('-v', '--verbose', action='count', default=0,
                        help='''Verbosity level. By default little to none information is printed.
Use -v once to increase information logs about each step, and -vv to 
print every command that is being run.''')
    return parser.parse_args()
    

def main() -> None:
    args = parse_args()
    configure_logger(logger, args.verbose)
    logger.debug(f'Started with arguments: {vars(args)}')

    config = load_config(args.config)
    DetectStrand(**config).run()


def merge() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--inputs', required=True,
                        nargs='+', type=Path, help='''Input file(s).''')
    parser.add_argument('-o', '--output', required=True,
                        type=Path, help='''Output file.''')
    parser.add_argument('-v', '--verbose', action='count', default=0,
                        help='''Verbosity level. By default little to none information is printed.
Use -v to increase information logs about each step.''')
    
    args = parser.parse_args()
    configure_logger(logger, args.verbose)
    logger.debug(f'Started with arguments: {vars(args)}')

    merge_tables(args.inputs, args.output)


if __name__ == '__main__':
    main()
