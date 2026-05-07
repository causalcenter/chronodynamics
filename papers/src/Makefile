# Build both arXiv preprints.
#
# Usage:
#   make            # build both PDFs into ../preprints/, then auto-clean
#   make chrono     # build chrono_mass.pdf only
#   make crisis     # build clock_crisis.pdf only
#   make clean      # remove build artefacts (keeps PDFs)
#   make distclean  # also remove generated PDFs from ../preprints/
#   make arxiv      # bundle each paper into a tar.gz suitable for arXiv upload

LATEX     ?= pdflatex
LATEXOPTS ?= -interaction=nonstopmode -halt-on-error -file-line-error

PAPERS    := chrono_mass clock_crisis
OUTDIR    := ../preprints
PDFS      := $(PAPERS:%=$(OUTDIR)/%.pdf)

AUX_EXT   := aux log out toc lof lot bbl blg fls fdb_latexmk synctex.gz \
             nav snm vrb lock run.xml bcf idx ind ilg auxlock figlist makefile

.PHONY: all build chrono crisis clean distclean arxiv

# Top-level target: build all PDFs into ../preprints/, then run clean.
all: build
	@$(MAKE) --no-print-directory clean

build: $(PDFS)

chrono: $(OUTDIR)/chrono_mass.pdf
	@$(MAKE) --no-print-directory clean

crisis: $(OUTDIR)/clock_crisis.pdf
	@$(MAKE) --no-print-directory clean

# Two-pass build picks up cross-references and the embedded thebibliography.
# pdflatex writes the PDF in cwd; we move it to ../preprints/ after the second pass.
$(OUTDIR)/%.pdf: %.tex | $(OUTDIR)
	$(LATEX) $(LATEXOPTS) $<
	$(LATEX) $(LATEXOPTS) $<
	mv $*.pdf $@

$(OUTDIR):
	@mkdir -p $@

clean:
	@for p in $(PAPERS); do \
		for ext in $(AUX_EXT); do \
			rm -f $$p.$$ext; \
		done; \
	done
	rm -f *.aux *.log *.out *.toc *.lof *.lot *.bbl *.blg *.fls \
	      *.fdb_latexmk *.synctex.gz *.nav *.snm *.vrb *.lock \
	      *.run.xml *.bcf *.idx *.ind *.ilg *.auxlock *.figlist \
	      texput.log
	rm -rf _minted-* .texpadtmp/

distclean: clean
	rm -f $(PDFS)
	rm -f arxiv_*.tar.gz

# Bundle each paper into an arXiv-ready tarball: just the .tex (no .bbl needed,
# bibliography is inline). arXiv will run pdflatex twice on its end.
arxiv: $(PAPERS:%=arxiv_%.tar.gz)

arxiv_%.tar.gz: %.tex
	tar -czf $@ $<
