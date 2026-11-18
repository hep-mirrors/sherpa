#ifndef GoSam_Virtual_H
#define GoSam_Virtual_H

#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "PHASIC++/Process/Process_Info.H"
#include "ATOOLS/Phys/Flavour.H"

#include "AddOns/GoSam/GoSam_Interface.H"

namespace GoSam {

  class GoSam_Virtual : public PHASIC::Virtual_ME2_Base {
    size_t m_pn;
    bool m_active;
    int m_proc_id, m_oqcd, m_oqed;
    std::vector<double> m_momenta;
    std::array<double,4> m_result;
    double m_Norm;
    double m_g2qcd, m_g2qed;
  public:
    GoSam_Virtual(const PHASIC::Process_Info& pi,const ATOOLS::Flavour_Vector& flavs,
                  bool active, int proc_id=-1);
    ~GoSam_Virtual() { }

    void Calc(const ATOOLS::Vec4D_Vector& momenta);
    bool SetColours(const ATOOLS::Vec4D_Vector& momenta);
    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom);
#ifdef VALIDATE_GOSAM_INTERFACE
    PHASIC::Virtual_ME2_Base* p_check;
    bool m_errorOnFailedCheck;
    bool m_printAllChecks;
#endif
  };
}


#endif
