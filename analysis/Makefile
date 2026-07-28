# analysis/Makefile -- the documented entry point for regenerating the paper's
# figures and tables. Run from the REPOSITORY ROOT:
#
#   make -f analysis/Makefile all
#
# or simply:  Rscript analysis/make.R

R      ?= Rscript
FIGDIR := analysis/figures
TABDIR := analysis/tables

.PHONY: all fig01 fig02 tables clean check

all:
	$(R) analysis/make.R

fig01:
	$(R) analysis/scripts/fig01_roc_grid.R

fig02:
	$(R) analysis/scripts/fig02_deployed_bars.R

tables:
	$(R) analysis/scripts/tab03_provenance.R

## Re-verify every produced .eps against the Taylor & Francis artwork spec.
check:
	$(R) -e 'source("analysis/R/04_check.R"); \
	         for (f in list.files("$(FIGDIR)", "[.]eps$$", full.names=TRUE)) \
	           check_figure(f, min_line_pt = 0.5, min_text_pt = 8)'

## Drop derived artefacts. The cache is a convenience only: removing it and
## re-running must reproduce every figure bit-for-bit.
clean:
	rm -rf analysis/cache $(FIGDIR) $(TABDIR)
