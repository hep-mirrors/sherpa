#ifndef CXXFLAGS_H
#define CXXFLAGS_H

#define ARCH_LINUX

#define SHERPA_VERSION "1"
#define SHERPA_SUBVERSION "0"

#define SHERPA_BUILD_PATH "/home/ptpp90/stefan/hep/sherpa/rel-1.0.7/SHERPA-MC-1.0"

#ifdef USING__CLHEP
#undef USING__CLHEP
#endif
#ifdef USING__ROOT
#undef USING__ROOT
#endif
#ifdef USING__LHAPDF
#undef USING__LHAPDF
#endif
#define USING__Model
#ifdef USING__Isajet
#undef USING__Isajet
#endif
#ifdef USING__Hdecay
#undef USING__Hdecay
#endif
#define USING__Amisic
#define USING__Hadrons
#define USING__Sherpa

#endif
