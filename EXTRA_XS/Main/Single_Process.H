#ifndef EXTRA_XS_Main_Single_Process_H
#define EXTRA_XS_Main_Single_Process_H

#include "PHASIC++/Process/Single_Process.H"

namespace EXTRAXS {

  class ME2_Base;
  class Virtual_ME2_Base;
  class ME_Base;

  class Single_Process: public PHASIC::Single_Process {
  private:
    ME2_Base* p_born_me2;
    Virtual_ME2_Base* p_virtual_me2;
    PHASIC::nlo_type::code m_nlotype;

  public:

    // constructors
    Single_Process();

    // destructor
    ~Single_Process();

    // member functions
    bool Initialize();
    bool FillIntegrator(PHASIC::Phase_Space_Handler *const psh);

    double Differential(const ATOOLS::Vec4D_Vector& momenta);
    double Differential2();

    // inline functions
    inline ME2_Base *GetME() const { return p_born_me2; }

  };

}// end of namespace EXTRAXS

#endif



