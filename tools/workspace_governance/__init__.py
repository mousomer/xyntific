"""Generic, dependency-light workspace governance pack."""

from .resolver.core import Diagnostic, GovernanceError, GovernanceResolver

__all__ = ["Diagnostic", "GovernanceError", "GovernanceResolver"]
