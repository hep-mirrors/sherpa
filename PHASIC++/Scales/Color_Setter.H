#ifndef PHASIC_Cluster_Color_Setter_H
#define PHASIC_Cluster_Color_Setter_H

#include "PHASIC++/Main/Color_Integrator.H"
#include "PHASIC++/Process/Single_Process.H"

namespace EXTRAXS { class ME2_Base; }

namespace PHASIC {

  class Color_Setter {
  private:

    typedef std::map<ATOOLS::Flavour_Vector,EXTRAXS::ME2_Base*> Flav_ME_Map;

    PHASIC::Process_Base *p_xs;

    static PHASIC::NLOTypeStringProcessMap_Map m_pmap;
    PHASIC::Process_Vector m_procs;

    Flav_ME_Map m_xsmap;

    int m_cmode;

    bool SetRandomColors(ATOOLS::Cluster_Amplitude *const ampl);
    bool SetSumSqrColors(ATOOLS::Cluster_Amplitude *const ampl);
    bool SetLargeNCColors(ATOOLS::Cluster_Amplitude *const ampl);

  public:

    Color_Setter(const int cmode);

    ~Color_Setter();

    bool Initialize(ATOOLS::Cluster_Amplitude *const ampl);

    Process_Base *GetProcess(ATOOLS::Cluster_Amplitude *const ampl);

    void SetColors(ATOOLS::Cluster_Amplitude *const ampl);

  };//end of class Color_Setter

}// end of namespace PHASIC

#endif
