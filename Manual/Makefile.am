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

update-manpage:
	groff -Tascii -mandoc Sherpa.trf > Sherpa.1

completion.index: $(srcdir)/*.texi $(srcdir)/Sherpa.info
	@echo "Creating index for tab completion."
	@rm -f Sherpa.txt* Sherpa.index
	-@makeinfo --plaintext -o Sherpa.txt Sherpa.texi && \
	csplit --quiet --prefix Sherpa.txt. Sherpa.txt "%Appendix A Index%1" && \
	csplit --quiet --prefix Sherpa.txt.00. Sherpa.txt.00 "%Appendix A Index%3" && \
	cut Sherpa.txt.00.00 -d ":" -f 1 > completion.index || touch completion.index
	@rm -f Sherpa.txt*

dvi:

EXTRA_DIST    = @SVNINCLUDE@ Sherpa.1
