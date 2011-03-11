#ifndef PDF_Main_Jet_Criterion_H
#define PDF_Main_Jet_Criterion_H

#include "ATOOLS/Phys/NLO_Subevt.H"

namespace PDF {

  double Qij2(const ATOOLS::Vec4D &pi,const ATOOLS::Vec4D &pj,
	      const ATOOLS::Vec4D &pk,
	      const ATOOLS::Flavour &fi=ATOOLS::Flavour(kf_none),
	      const ATOOLS::Flavour &fj=ATOOLS::Flavour(kf_none),
	      const double& dparam=-1.0,const int mode=0);

  double Qij2Min(const ATOOLS::Vec4D_Vector &p,
		 ATOOLS::NLO_subevtlist *const subs);

}// end of namespace PDF

#endif
