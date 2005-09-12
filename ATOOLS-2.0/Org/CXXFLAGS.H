#ifndef CXXFLAGS_H
#define CXXFLAGS_H

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
