import argparse
import logging

from pydantic import ValidationError

from .XRNA import XRNAProcessor
from .utils import (
    configure_logger, exit_with_validation_error, load_config
)

logger = logging.getLogger()


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
    try:
        processor = XRNAProcessor(**config)
    except ValidationError as error:
        exit_with_validation_error(error)
    processor.run()
    

if __name__ == '__main__':
    main()
