# Major project report

NITK-style B.Tech. major project report (based on the department thesis template).

Build (needs XeLaTeX; falls back to TeX Gyre Termes if Times New Roman is not installed):

    xelatex report && bibtex report && xelatex report && xelatex report

Edit the names, register numbers, guide and department in the "Project details" block near the top of `report.tex`.
Chapters are in `chapters/`; tables and plot data in `generated/` are produced by `convex/make_tables.py` and `arc/make_tables.py`.
