"""anchovy: single-cell viral consensus and genotype-network pipeline."""

# THE single definition of the version. pyproject.toml reads it from here via
# setuptools' dynamic-version support, so packaging metadata and `anchovy
# --version` cannot disagree. Previously both files carried the literal, with
# pyproject calling itself the single source of truth while a second copy sat
# here -- they happened to agree, but nothing made them.
__version__ = "1.2.0"
