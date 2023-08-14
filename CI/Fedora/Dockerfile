FROM fedora:35
RUN  yum -y  install  dnf-plugins-core \
                      bc make ninja-build cmake binutils git wget diffutils file sed gawk grep which autoconf automake libtool ccache \
                      gcc-gfortran gcc-c++ clang bzip2   openssl-devel openssl \
                      libzip-devel  zlib zlib-devel pythia8-devel pythia8 python3-lhapdf lhapdf lhapdf-devel  \
                      root-core root root-physics HepMC3-rootIO-devel HepMC3-devel HepMC3-search HepMC3-search-devel  HepMC3-rootIO HepMC3  HepMC HepMC-devel \
                      texinfo latexmk 'tex(latex)' ghostscript 'tex(sectsty.sty)' \
                      openmpi environment-modules  openmpi-devel  \
                      doxygen catch-devel \
                      python3 python3-devel swig python3-matplotlib && yum -y clean all
# Packages coming from HEPrpms ->
RUN  dnf -y copr enable averbyts/HEPrpms &&  yum -y  install  gosam  fastjet fastjet-devel  YODA-devel Rivet-devel Rivet Rivet-devel hztool qd qd-devel \
                      blackhat-devel blackhat recola cernlib cernlib-devel cernlib-static MCFM MG5_aMC openloops \ 
                      f90cache  && yum -y update lhapdf && yum -y clean all
# Packages coming from HEPrpms <-
RUN  yum -y  install  sphinx python3-sphinx python3-sphinxcontrib-bibtex python3-sphinx* \
                      'tex(atbegshi.sty)' 'tex(babel.sty)' 'tex(bigintcalc.sty)' 'tex(bitset.sty)' 'tex(etexcmds.sty)' 'tex(gettitlestring.sty)'   \
                      'tex(iftex.sty)' 'tex(ifvtex.sty)' 'tex(infwarerr.sty)' 'tex(intcalc.sty)' 'tex(kvdefinekeys.sty)' 'tex(kvsetkeys.sty)'  \
                      'tex(ltxcmds.sty)' 'tex(pdfescape.sty)' 'tex(pdftexcmds.sty)' 'tex(uniquecounter.sty)' 'tex(amsfonts.sty)' 'tex(amssymb.sty)'   \
                      'tex(amsbsy.sty)' 'tex(amsgen.sty)' 'tex(amsmath.sty)' 'tex(amsopn.sty)' 'tex(amstext.sty)' 'tex(atveryend.sty)' 'tex(auxhook.sty)' 'tex(alltt.sty)'   \
                      'tex(atbegshi-ltx.sty)' 'tex(atveryend-ltx.sty)' 'tex(fontenc.sty)' 'tex(inputenc.sty)' 'tex(makeidx.sty)' 'tex(textcomp.sty)' 'tex(capt-of.sty)'   \
                      'tex(cmap.sty)' 'tex(epstopdf-base.sty)' 'tex(fancyhdr.sty)' 'tex(fancyvrb.sty)' 'tex(float.sty)' 'tex(fncychap.sty)' 'tex(framed.sty)'   \
                      'tex(geometry.sty)' 'tex(color.sty)' 'tex(graphics.sty)' 'tex(graphicx.sty)' 'tex(keyval.sty)' 'tex(trig.sty)' 'tex(hycolor.sty)' 'tex(hyperref.sty)'   \
                      'tex(nameref.sty)' 'tex(kvoptions.sty)' 'tex(letltxmacro.sty)' 'tex(needspace.sty)' 'tex(hypcap.sty)' 'tex(parskip-2001-04-09.sty)' 'tex(parskip.sty)'  \
                      'tex(refcount.sty)' 'tex(rerunfilecheck.sty)' 'tex(tabulary.sty)' 'tex(tgheros.sty)' 'tex(tgtermes.sty)' 'tex(titlesec.sty)' 'tex(array.sty)'   \
                      'tex(longtable.sty)' 'tex(multicol.sty)' 'tex(upquote.sty)' 'tex(url.sty)' 'tex(varwidth.sty)' 'tex(wrapfig.sty)' 'tex(xcolor.sty)' \
                       rsync --exclude=python3-mistune08 --exclude=python3-sphinxcontrib-openapi && yum -y  clean all  
RUN lhapdf update && lhapdf install PDF4LHC21_40_pdfas
