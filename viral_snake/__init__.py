"""Viral Snake - A Snakemake pipeline for viral genomics"""

__version__ = "0.1.0"

from .utils import *
from .collections import load_collections
from .cli import main as cli