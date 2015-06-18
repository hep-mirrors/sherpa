#ifndef DIRE_Cluster_Color_Setter_H
#define DIRE_Cluster_Color_Setter_H

#include "PHASIC++/Main/Color_Integrator.H"
#include "PHASIC++/Process/Single_Process.H"

namespace DIRE {

  class Color_Setter {
  private:

    PHASIC::Process_Base *p_xs;

    PHASIC::NLOTypeStringProcessMap_Map m_pmap;
    PHASIC::Process_Vector m_procs;

    int m_cmode;

    bool SetRandomColors(ATOOLS::Cluster_Amplitude *const ampl);

  public:

    Color_Setter(const int mode);

    ~Color_Setter();

    bool SetColors(ATOOLS::Cluster_Amplitude *const ampl);

  };//end of class Color_Setter

}// end of namespace DIRE

#endif
