## For documentation on texinfo, see
## http://www.gnu.org/software/texinfo/manual

info_TEXINFOS = Sherpa.texi

Sherpa_TEXINFOS = \
	Macros.texi \
	Bibliography.texi \
	Introduction.texi \
	Commands.texi \
	Input.texi \
	Parameters.texi \
	Beam.texi \
	Model.texi \
	ISR.texi \
	ME.texi \
	Processes.texi \
	Selector.texi \
	Integration.texi \
	Shower.texi \
	SoftQED.texi \
	MI.texi \
	Tips.texi \
	Customization.texi \
	Run.texi \
	RNG.texi

man1_MANS = Sherpa.1

dist_pkgdata_DATA         = \
	sherpa-completion \
	completion.index

dist_html_DATA = \
	Sherpa.html

update-manpage:
	groff -Tascii -mandoc Sherpa.trf > Sherpa.1

completion.index: $(srcdir)/*.texi $(srcdir)/Sherpa.info
	@echo "Creating index for tab completion."
	@rm -f Sherpa.txt* Sherpa.index
	-@makeinfo --plaintext -o Sherpa.txt Sherpa.texi && \
        awk 'BEGIN{ m=0; }{ if (match($$0,"Appendix B Index")>0) { m=m+1; getline; getline; getline;} if (m==2) { print $$0 } }' Sherpa.txt | cut -d ":" -f 1 > completion.index || touch completion.index
	@rm -f Sherpa.txt*

Sherpa.html: $(srcdir)/*.texi
	if test -x "$$(which texi2html)"; then \
          texi2html Sherpa.texi; \
        else \
          touch Sherpa.html; \
        fi

dvi:

EXTRA_DIST    = @SVNINCLUDE@ Sherpa.1
