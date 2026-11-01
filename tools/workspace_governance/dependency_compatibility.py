"""Conservative PEP 440 compatibility checks for a shared environment.

This is deliberately not a dependency resolver.  It proves only a small,
useful subset of version-specifier intersections and reports everything else
as unknown.  Callers must refuse an unknown result before mutating a shared
environment rather than treating a lack of proof as compatibility.
"""

from __future__ import annotations

import json
import sys
from collections.abc import Iterable
from dataclasses import dataclass

from packaging.requirements import InvalidRequirement, Requirement
from packaging.version import InvalidVersion, Version


@dataclass(frozen=True)
class Result:
    status: str
    requirement: str
    detail: str


def _bounds(
    requirement: Requirement,
) -> tuple[tuple[Version, bool] | None, tuple[Version, bool] | None] | None:
    """Return a closed/open version interval when it is safe to do so.

    Exclusions, compatible-release clauses, wildcards, direct URLs and
    arbitrary equality deliberately stay outside this proof system.
    """
    if requirement.url or requirement.extras:
        return None
    lower: tuple[Version, bool] | None = None
    upper: tuple[Version, bool] | None = None
    for specifier in requirement.specifier:
        if specifier.operator not in {">", ">=", "<", "<="}:
            return None
        try:
            version = Version(specifier.version)
        except InvalidVersion:
            return None
        if specifier.operator in {">", ">="}:
            candidate = (version, specifier.operator == ">=")
            if (
                lower is None
                or candidate[0] > lower[0]
                or (candidate[0] == lower[0] and not candidate[1])
            ):
                lower = candidate
        else:
            candidate = (version, specifier.operator == "<=")
            if (
                upper is None
                or candidate[0] < upper[0]
                or (candidate[0] == upper[0] and not candidate[1])
            ):
                upper = candidate
    return lower, upper


def _empty(
    bounds: tuple[tuple[Version, bool] | None, tuple[Version, bool] | None],
) -> bool:
    lower, upper = bounds
    return bool(
        lower
        and upper
        and (
            lower[0] > upper[0]
            or (lower[0] == upper[0] and not (lower[1] and upper[1]))
        )
    )


def compare(left: str, right: str) -> Result:
    """Classify two requirements for the same normalized distribution name."""
    try:
        first, second = Requirement(left), Requirement(right)
    except InvalidRequirement as exc:
        return Result("unknown", f"{left!r} / {right!r}", f"cannot parse: {exc}")
    if first.name.lower().replace("_", "-") != second.name.lower().replace("_", "-"):
        return Result(
            "compatible", f"{first.name} / {second.name}", "different distributions"
        )
    if str(first) == str(second):
        return Result("compatible", str(first), "identical declarations")
    first_bounds, second_bounds = _bounds(first), _bounds(second)
    if first_bounds is None or second_bounds is None:
        return Result(
            "unknown",
            f"{first} / {second}",
            "intersection needs resolver semantics outside the conservative interval subset",
        )
    combined = (
        max(
            (bound for bound in (first_bounds[0], second_bounds[0]) if bound),
            default=None,
            key=lambda item: (item[0], not item[1]),
        ),
        min(
            (bound for bound in (first_bounds[1], second_bounds[1]) if bound),
            default=None,
            key=lambda item: (item[0], item[1]),
        ),
    )
    if _empty(combined):
        return Result(
            "contradictory", f"{first} / {second}", "version intervals do not overlap"
        )
    return Result("compatible", f"{first} / {second}", "version intervals overlap")


def classify(existing: Iterable[str], requested: Iterable[str]) -> list[Result]:
    """Compare every same-name declaration and retain unknown results."""
    results: list[Result] = []
    for left in existing:
        for right in requested:
            result = compare(left, right)
            if (
                result.status != "compatible"
                or result.detail != "different distributions"
            ):
                results.append(result)
    return results


def main(argv: list[str] | None = None) -> int:
    raw = json.load(sys.stdin)
    if not isinstance(raw, dict):
        raise TypeError("compatibility input must be an object")
    existing = raw.get("existing", [])
    requested = raw.get("requested", [])
    if not all(
        isinstance(value, list) and all(isinstance(item, str) for item in value)
        for value in (existing, requested)
    ):
        raise ValueError("existing and requested requirements must be string lists")
    results = classify(existing, requested)
    payload = {"results": [result.__dict__ for result in results]}
    json.dump(payload, sys.stdout, sort_keys=True)
    print()
    return 1 if any(result.status != "compatible" for result in results) else 0


if __name__ == "__main__":
    raise SystemExit(main())
