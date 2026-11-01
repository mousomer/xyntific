"""Dependency-light numeric Python floor; prerelease suffixes do not change release numbers."""

import re


def numeric_python_floor(version, requirement):
    floor = re.fullmatch(r">=\s*(\d+(?:\.\d+)*)", requirement)
    release = re.fullmatch(
        r"(\d+(?:\.\d+)*)(?:(?:a|b|rc)\d+|\.dev\d+|\+[-\w.]+)?", version
    )
    if floor is None or release is None:
        raise ValueError("Expected a numeric Python release and an explicit >= floor")
    actual = tuple(int(part) for part in release[1].split("."))
    required = tuple(int(part) for part in floor[1].split("."))
    width = max(len(actual), len(required))
    return actual + (0,) * (width - len(actual)) >= required + (0,) * (
        width - len(required)
    )
