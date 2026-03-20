# Contributing to CellMetPro

Thank you for your interest in contributing to CellMetPro! This document provides guidelines for contributing to the project.

## Getting Started

1. **Fork the repository** on GitHub
2. **Clone your fork** locally:
   ```bash
   git clone https://github.com/YOUR-USERNAME/CellMetPro.git
   cd CellMetPro
   ```

3. **Set up development environment**:
   ```bash
   conda create -n cellmetpro python=3.11
   conda activate cellmetpro
   pip install -e ".[dev]"
   ```

## Development Workflow

### 1. Create a Branch

Create a new branch for your feature or bugfix:

```bash
git checkout -b feature/your-feature-name
```

Branch naming conventions:
- `feature/` - New features
- `fix/` - Bug fixes
- `docs/` - Documentation improvements
- `test/` - Test additions or improvements
- `refactor/` - Code refactoring

### 2. Make Your Changes

Follow these guidelines:

#### Code Style

- **Follow PEP 8** guidelines
- **Use type hints** for all function signatures
- **Write docstrings** for all public functions and classes (Google style)
- **Keep functions small and focused**
- **Use descriptive variable names**

Example:

```python
def compare_groups(
    scores: pd.DataFrame,
    group_a: str,
    group_b: str,
    method: str = "wilcoxon",
) -> pd.DataFrame:
    """
    Compare metabolic reaction scores between two cell groups.

    Args:
        scores: Reaction scores matrix (reactions × cells).
        group_a: Label of the first group.
        group_b: Label of the second group.
        method: Statistical test to use ('wilcoxon' or 'ttest').

    Returns:
        DataFrame with reactions, log2FC, p-values, and adjusted p-values.

    Raises:
        ValueError: If group labels are not found in scores.
    """
    ...
```

#### Testing

- **Write tests** for all new code
- **Aim for >80% coverage**
- **Use pytest fixtures** for reusable test data
- **Test edge cases** and error conditions

Example:

```python
import pytest
from cellmetpro.analysis.differential import DifferentialAnalysis

def test_compare_groups_basic(sample_scores, sample_labels):
    """Test basic group comparison returns expected columns."""
    da = DifferentialAnalysis(sample_scores, sample_labels)
    result = da.compare_groups("A", "B")

    assert "reaction" in result.columns
    assert "log2fc" in result.columns
    assert "padj_bh" in result.columns
```

#### Documentation

- **Update docstrings** when changing function behaviour
- **Add examples** for new public API
- **Update README.md** for significant new features
- **Add type hints** to all function signatures

### 3. Run Tests and Checks

Before committing, ensure all checks pass:

```bash
# Format code
black cellmetpro tests

# Lint code
ruff check cellmetpro tests

# Type check
mypy cellmetpro

# Run tests
pytest

# Run tests with coverage
pytest --cov=cellmetpro --cov-report=html
```

### 4. Commit Your Changes

Write clear, descriptive commit messages:

```bash
git add .
git commit -m "feat: add Louvain clustering support

- Add louvain method to ClusteringAnalysis
- Update CLI --method choices
- Add tests for new clustering option"
```

Commit message prefixes:
- `feat:` - New feature
- `fix:` - Bug fix
- `docs:` - Documentation changes
- `test:` - Test additions/changes
- `refactor:` - Code refactoring
- `style:` - Code style changes
- `chore:` - Maintenance tasks

### 5. Push and Create Pull Request

```bash
git push origin feature/your-feature-name
```

Then create a Pull Request on GitHub with:
- **Clear title** describing the change
- **Description** of what changed and why
- **Link to related issues** (if applicable)
- **Breaking changes** noted clearly

## Code Review Process

1. Maintainers will review your PR
2. Address any feedback or requested changes
3. Once approved, your PR will be merged into `dev`, then released to `main`

## Testing Guidelines

### Unit Tests

Test individual functions and classes in isolation:

```python
def test_reaction_score_normalization():
    """Test that scores are normalized to [0, 1]."""
    raw = np.array([0.0, 5.0, 10.0])
    normalized = normalize_scores(raw)
    assert normalized.min() == pytest.approx(0.0)
    assert normalized.max() == pytest.approx(1.0)
```

### Integration Tests

Test multiple components working together:

```python
@pytest.mark.integration
def test_full_differential_pipeline(tmp_path):
    """Test scoring → clustering → differential end-to-end."""
    ...
```

### Test Fixtures

Use fixtures for reusable test data:

```python
@pytest.fixture
def sample_scores():
    """Small reaction scores matrix for testing."""
    rng = np.random.default_rng(42)
    data = rng.random((50, 20))
    return pd.DataFrame(data, index=[f"RXN{i}" for i in range(50)])
```

## Adding New Features

### 1. Discuss First

For major features, open an issue first to discuss:
- Design approach
- API changes
- Breaking changes
- Implementation plan

### 2. Follow Module Structure

New analysis modules should follow the existing pattern:

```python
from cellmetpro.analysis.base import BaseAnalysis

class NewAnalysis(BaseAnalysis):
    """Description of the analysis module."""

    def __init__(self, scores: pd.DataFrame, **kwargs):
        super().__init__(scores)

    def run(self, *args, **kwargs) -> pd.DataFrame:
        """Execute the analysis and return results."""
        ...
```

### 3. Add Tests

Include comprehensive tests:
- Unit tests for new functions
- Integration tests for workflows
- Edge case and error handling tests

### 4. Update Documentation

- Add docstrings
- Update the README feature table if applicable
- Update CHANGELOG.md under an `[Unreleased]` section

## Reporting Issues

### Bug Reports

Include:
- **Clear title and description**
- **Steps to reproduce**
- **Expected vs actual behaviour**
- **Environment details** (OS, Python version, package version)
- **Minimal reproducible example**
- **Error messages and stack traces**

### Feature Requests

Include:
- **Clear use case description**
- **Proposed solution or API**
- **Alternative approaches considered**
- **Potential impact** on existing functionality

## Code of Conduct

- Be respectful and inclusive
- Welcome newcomers
- Accept constructive criticism
- Focus on what is best for the community
- Show empathy towards others

## Questions?

- Open an issue for questions
- Check existing issues and PRs first
- Read the documentation in `README.md`

## License

By contributing, you agree that your contributions will be licensed under the MIT License that covers this project.

---

Thank you for contributing to CellMetPro!
