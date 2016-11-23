#ifndef OpenLoops_Virtual_H
#define OpenLoops_Virtual_H

#include "PHASIC++/Process/Virtual_ME2_Base.H"

namespace OpenLoops {

  class OpenLoops_Virtual : public PHASIC::Virtual_ME2_Base {

    int m_ol_id;


  public:

    OpenLoops_Virtual(const PHASIC::Process_Info& pi,
                      const ATOOLS::Flavour_Vector& flavs,
                      int ol_id);

    ~OpenLoops_Virtual() {}
  
    void SwitchMode(const int mode);

    inline bool SetColours(const ATOOLS::Vec4D_Vector& momenta) {
      return true;
    }

    inline double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom) {
      return 4.*M_PI;
    }
  
    void Calc(const ATOOLS::Vec4D_Vector& momenta);

  };

}

#endif
