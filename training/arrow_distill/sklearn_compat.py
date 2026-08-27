"""Narrow compatibility helpers for historical local scikit-learn artifacts."""

from __future__ import annotations

from typing import Any

from sklearn.impute import SimpleImputer


def repair_legacy_simple_imputers(value: Any, seen: set[int] | None = None) -> Any:
    """Populate the 1.9 ``SimpleImputer`` attribute absent in 1.7 pickles.

    The legacy sprint artifacts were fit with scikit-learn 1.7.2. Version 1.9
    reads those models but expects ``_fill_dtype`` during ``transform``. Its
    value is the input dtype recorded by 1.7 as ``_fit_dtype``.
    """
    if seen is None:
        seen = set()
    if id(value) in seen:
        return value
    seen.add(id(value))
    if isinstance(value, SimpleImputer) and not hasattr(value, "_fill_dtype"):
        value._fill_dtype = value._fit_dtype
    if isinstance(value, dict):
        for item in value.values():
            repair_legacy_simple_imputers(item, seen)
    elif isinstance(value, (list, tuple)):
        for item in value:
            repair_legacy_simple_imputers(item, seen)
    elif hasattr(value, "steps"):
        for _, item in value.steps:
            repair_legacy_simple_imputers(item, seen)
    return value
