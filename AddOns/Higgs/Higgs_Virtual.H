#ifndef AddOns__Higgs__Higgs_Virtual_H
#define AddOns__Higgs__Higgs_Virtual_H

#include "PHASIC++/Process/Virtual_ME2_Base.H"

namespace MODEL { class Model_Base; }
namespace AMEGIC { class Basic_Sfuncs; }

namespace HIGGS {

  class Higgs_Virtual: public PHASIC::Virtual_ME2_Base {
  private:

    int    m_int, m_io, m_proc;
    double m_mh, m_gh;

    AMEGIC::Basic_Sfuncs *p_bs;
    std::vector<int>      m_b;

    static MODEL::Model_Base *s_model;

    double  htheta(const double &x);
    Complex lnrat(const double &x,const double &y);

  public:

    Higgs_Virtual(const PHASIC::Process_Info& pi,
		  const ATOOLS::Flavour_Vector& flavs,
		  int mode,int io);

    ~Higgs_Virtual();

    void Calc(const ATOOLS::Vec4D_Vector &p);
    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector &p);

    inline static void SetModel(MODEL::Model_Base *const md) { s_model=md; }

  };

}// end of namespace HIGGS

#endif
