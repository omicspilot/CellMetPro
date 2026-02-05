"""Sample data for CellMetPro tutorials and testing.

This module provides small synthetic datasets that users can use to
explore CellMetPro's functionality without needing to prepare their
own data first.

Example
-------
>>> from cellmetpro.data import load_sample_expression, load_sample_groups
>>> expression = load_sample_expression()
>>> groups = load_sample_groups()
>>> print(f"Expression: {expression.shape[0]} genes x {expression.shape[1]} cells")
>>> print(f"Groups: {groups['group'].unique()}")
"""

from .samples import (
    create_sample_model,
    get_sample_data_path,
    load_sample_expression,
    load_sample_groups,
    load_sample_reaction_scores,
)

__all__ = [
    "load_sample_expression",
    "load_sample_groups",
    "load_sample_reaction_scores",
    "create_sample_model",
    "get_sample_data_path",
]
