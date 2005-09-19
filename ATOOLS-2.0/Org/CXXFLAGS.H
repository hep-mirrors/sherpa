#ifndef CXXFLAGS_H
#define CXXFLAGS_H

#define SHERPA_VERSION "1"
#define SHERPA_SUBVERSION "0"

#define SHERPA_BUILD_PATH "/Users/mnedos/SHERPA"

#ifdef USING__CLHEP
#undef USING__CLHEP
#endif
#define USING__ROOT
#ifdef USING__LHAPDF
#undef USING__LHAPDF
#endif
#define USING__Model
#define USING__Adicic
#define USING__Amisic
#define USING__CSS
#define USING__Ahadic
#define USING__Hadrons
#define USING__Sherpa
#ifdef USING__MCatNLO
#undef USING__MCatNLO
#endif

#endif
