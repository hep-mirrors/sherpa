FROM archlinux:latest
RUN  pacman -Syu --noconfirm \
                      bc make  cmake binutils git wget  file sed gawk grep which autoconf automake libtool ccache \
                      doxygen texlive-bibtexextra \
                      pkgconf \
                      gcc gcc-fortran \
                      libzip lhapdf \
                      catch2 \
                      doxygen
RUN lhapdf update && lhapdf install PDF4LHC21_40_pdfas
