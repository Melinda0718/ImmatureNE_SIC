# -*- coding: utf-8 -*-
"""build_tables.py - build the basic statistical tables (CSV / LaTeX+PDF / Word update).

    python tables/build_tables.py              # full pipeline
    python tables/build_tables.py --no-docx    # CSV + tex/pdf only (debugging)
    python tables/build_tables.py --no-tex     # CSV + docx update only

Data come from the four CSVs in `data/`; table structure comes from the docx files in `tables/`..
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen_tables  # noqa: E402

if __name__ == "__main__":
    gen_tables.build_all(do_tex="--no-tex" not in sys.argv,
                         do_docx="--no-docx" not in sys.argv)
