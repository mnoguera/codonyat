#!/usr/bin/env python3
import re

with open("pyproject.toml") as f:
    pyproject = f.read()

# Bump setuptools min
pyproject = pyproject.replace('requires = ["setuptools>=65", "wheel"]',
                               'requires = ["setuptools>=68", "wheel"]')

# Pin biopython upper bound
pyproject = pyproject.replace('"biopython>=1.79"',
                               '"biopython>=1.79,<2.1"')

# Replace classifiers
old = '''classifiers = [
    "Programming Language :: Python :: 3",
    "Operating System :: OS Independent",
]'''
new = '''classifiers = [
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3.11",
    "Programming Language :: Python :: 3.12",
    "Programming Language :: Python :: 3.13",
    "Operating System :: OS Independent",
    "Development Status :: 4 - Beta",
]'''
pyproject = pyproject.replace(old, new)

# Add ruff + pytest config
ruff_config = """

[tool.ruff]
line-length = 88
target-version = "py310"

[tool.ruff.lint]
select = ["E", "F", "W", "I", "N", "UP", "B", "C4"]
ignore = ["E501"]

[tool.pytest.ini_options]
testpaths = ["tests"]
"""
pyproject += ruff_config

with open("pyproject.toml", "w") as f:
    f.write(pyproject)

print("Done")
