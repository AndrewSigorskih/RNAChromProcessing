import logging
import sys

from pydantic import ValidationError

logger = logging.getLogger()

class StageFailedError(RuntimeError):
    pass


def exit_with_error(message: str = '') -> None:
    """Print error and exit"""
    if not message:
        message = 'Empty error message!'
    logger.critical(message)
    sys.exit(1)


def exit_with_validation_error(error: ValidationError) -> None:
    msg_lines = ["Config file validation failed:"]
    for sub_error in error.errors():
        msg_lines.append(
            f"\t- {''.join(sub_error['loc'])}: {sub_error['msg']}"
        )
    exit_with_error("\n".join(msg_lines))
