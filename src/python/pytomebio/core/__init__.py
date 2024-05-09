import re
from typing import Callable
from typing import Dict


def key_value_parser(pattern: str, flags: int = 0) -> Callable[[str], Dict[str, str]]:
    """
    Given a regex pattern containing exactly two capture groups, and optional flags, returns a
    function that parses a string and returns a dictionary containing all the captured key/value
    pairs.
    """
    regex = re.compile(pattern, flags)
    return lambda s: dict(regex.findall(s))
