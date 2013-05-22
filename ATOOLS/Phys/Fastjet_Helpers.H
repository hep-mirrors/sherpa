#ifndef ATOOLS_Phys_Fastjet_Helpers_H
#define ATOOLS_Phys_Fastjet_Helpers_H

#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#ifdef USING__FASTJET
#include "fastjet/PseudoJet.hh"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Flavour.H"

namespace ATOOLS {

  bool ToBeClustered(const ATOOLS::Flavour& flav, int bmode);
  bool BTag(const fastjet::PseudoJet& jet, int bmode);
  fastjet::PseudoJet MakePseudoJet(const Flavour& flav, const Vec4D& mom);
}

#endif
#endif
