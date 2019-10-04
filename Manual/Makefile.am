MAKEFLAGS = -j1

sphinx_SOURCE = $(srcdir)/source
sphinx_options = -D release=@VERSION@
sphinx_build_command = sphinx-build

if Manual_ENABLED
pkgdata_DATA = \
	sherpa-completion \
	completion.index \
	options.index
man_MANS = Sherpa.1
doc_DATA = sherpamanual.pdf

.NOTPARALLEL:
completion.index options.index: $(srcdir)/source/manpage.rst
	@echo "Creating index for tab completion."
	$(sphinx_build_command) -M completion $(sphinx_SOURCE) . $(sphinx_options)
	@mv completion/* .
	@rm -r completion

.NOTPARALLEL:
Sherpa.1: $(srcdir)/source/manpage.rst $(srcdir)/source/man/**
	$(sphinx_build_command) -M man $(sphinx_SOURCE) . $(sphinx_options)
	@mv man/Sherpa.1 .
	@rm -r man

.NOTPARALLEL:
sherpamanual_html: $(srcdir)/source/index.rst $(srcdir)/source/manual/**.rst
	$(sphinx_build_command) -M html $(sphinx_SOURCE) . $(sphinx_options)
	@mv html sherpamanual_html

.NOTPARALLEL:
sherpamanual.pdf: $(srcdir)/source/index.rst $(srcdir)/source/manual/**.rst
	$(sphinx_build_command) -M latexpdf $(sphinx_SOURCE) . $(sphinx_options)
	@mv latex/sherpamanual.pdf .
	@rm -r latex

.NOTPARALLEL:
sherpamanual.info: $(srcdir)/source/**.rst
	$(sphinx_build_command) -M info $(sphinx_SOURCE) . $(sphinx_options)
	@mv texinfo/sherpamanual.info .
	@rm -r texinfo

.NOTPARALLEL:
install-data-local: sherpamanual.info
	@mkdir -p $(DESTDIR)${infodir}/
	${INSTALL_DATA} sherpamanual.info $(DESTDIR)${infodir}/sherpamanual.info
endif


EXTRA_DIST    = scrape_completion.py \
	source

if Manual_ENABLED
# Distribute Manual so that building it isn't required
EXTRA_DIST    += sherpamanual.info \
	sherpamanual.pdf \
	sherpamanual_html \
	Sherpa.1
endif

# cleanup
CLEANFILES = options.index \
	completion.index \
	sherpamanual.info \
	sherpamanual.pdf \
	Sherpa.1

clean-local:
	-rm -rf sherpamanual_html doctrees
