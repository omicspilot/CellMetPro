"""Models module for genome-scale metabolic models.

This module handles loading and managing GEM (Genome-Scale Metabolic)
models for different organisms.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import cobra

# Path to bundled GEM models
GEMS_DIR = Path(__file__).parent / "gems"


def load_human_gem() -> cobra.Model:
    """Load the human genome-scale metabolic model.

    Returns
    -------
    cobra.Model
        Human GEM model.
    """
    raise NotImplementedError


def load_mouse_gem() -> cobra.Model:
    """Load the mouse genome-scale metabolic model.

    Returns
    -------
    cobra.Model
        Mouse GEM model.
    """
    raise NotImplementedError


def load_gem(organism: str) -> cobra.Model:
    """Load a GEM model by organism name.

    Parameters
    ----------
    organism : str
        Organism name: 'human', 'mouse', or path to custom model.

    Returns
    -------
    cobra.Model
        The loaded metabolic model.

    Raises
    ------
    ValueError
        If organism is not supported.
    """
    raise NotImplementedError
