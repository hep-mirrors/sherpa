#ifndef PDF_Main_Jet_Criterion_H
#define PDF_Main_Jet_Criterion_H

#include "ATOOLS/Phys/NLO_Subevt.H"
#include "ATOOLS/Org/Getter_Function.H"

namespace PDF {

  struct JetCriterion_Key {
  };// end of struct JetCriterion_Key

  class Jet_Criterion {
  public:

    virtual ~Jet_Criterion();

    virtual double Qij2
    (const ATOOLS::Vec4D &pi,const ATOOLS::Vec4D &pj,
     const ATOOLS::Vec4D &pk,
     const ATOOLS::Flavour &fi=ATOOLS::Flavour(kf_none),
     const ATOOLS::Flavour &fj=ATOOLS::Flavour(kf_none),
     const double& dparam=-1.0,const int mode=0) = 0;

    double Qij2Min(const ATOOLS::Vec4D_Vector &p,
		   ATOOLS::NLO_subevtlist *const subs);

  };// end of class Jet_Criterion

  typedef ATOOLS::Getter_Function
  <Jet_Criterion,JetCriterion_Key> JetCriterion_Getter;

}// end of namespace PDF

#endif
