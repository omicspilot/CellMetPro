"""Report generation module for CellMetPro.

This module provides automated report generation including:
- HTML reports with embedded visualizations
- PDF export (optional, requires weasyprint)
- Summary statistics tables
- Customizable templates
"""

from .report import (
    ReportGenerator,
    generate_html_report,
    generate_summary_stats,
)

__all__ = [
    "generate_html_report",
    "generate_summary_stats",
    "ReportGenerator",
]
