"""Lanceur : `python main.py <commande>` équivaut à `python -m gravitation <commande>`."""

import sys

from gravitation.cli import main

if __name__ == "__main__":
    sys.exit(main())
