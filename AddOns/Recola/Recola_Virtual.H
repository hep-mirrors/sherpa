#ifndef Recola_Virtual_H
#define Recola_Virtual_H

#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "PHASIC++/Process/Process_Base.H"

namespace Recola {

  class Recola_Virtual : public PHASIC::Virtual_ME2_Base {

    unsigned int m_recola_id;

  public:

    int m_eventcount;
    std::map<int,PHASIC::Process_Info> m_procmap;
    
    Recola_Virtual(const PHASIC::Process_Info& pi,
                      const ATOOLS::Flavour_Vector& flavs,
                      unsigned int recola_id);

    ~Recola_Virtual() {}

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
