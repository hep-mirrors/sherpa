#ifndef GoSam_Virtual_H
#define GoSam_Virtual_H

#include "PHASIC++/Process/Virtual_ME2_Base.H"

namespace GoSam {

  class GoSam_Virtual : public PHASIC::Virtual_ME2_Base {

    int m_gs_id;
    bool m_ismapped;

    int m_modebackup;

  public:

    GoSam_Virtual(const PHASIC::Process_Info& pi,
                  const ATOOLS::Flavour_Vector& flavs,
                  int gs_id);

    ~GoSam_Virtual() {}
  
    inline bool SetColours(const ATOOLS::Vec4D_Vector& momenta) {
      return true;
    }

    inline double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom) {
      return 4.*M_PI;
    }
  
    void Calc(const ATOOLS::Vec4D_Vector& momenta);

    bool IsMappableTo(const PHASIC::Process_Info& pi);
  };

}

#endif
