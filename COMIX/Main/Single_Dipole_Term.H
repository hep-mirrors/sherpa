#ifndef COMIX_Main_Single_Dipole_Term_H
#define COMIX_Main_Single_Dipole_Term_H

#include "COMIX/Main/Single_Process.H"

namespace COMIX {

  class Single_Dipole_Term: public PHASIC::Single_Process {
  protected:

    Single_Process *p_proc;

    Amplitude *p_bg;

    ATOOLS::NLO_subevt *p_sub, *p_msub;

    PHASIC::Process_Integrator *p_rsint;

    double m_w;

  public:

    // constructors
    Single_Dipole_Term(COMIX::Single_Process *const rs,
		       ATOOLS::NLO_subevt *const sub,
		       ATOOLS::NLO_subevt *const msub);

    // destructor
    ~Single_Dipole_Term();

    // member functions
    double Partonic(const ATOOLS::Vec4D_Vector &p,const int mode);
    double Differential(const ATOOLS::Cluster_Amplitude &ampl,
			int mode=0);

    bool Trigger(const ATOOLS::Vec4D_Vector &p);

    bool GeneratePoint();

    bool Combinable(const size_t &idi,const size_t &idj);

    const ATOOLS::Flavour_Vector &CombinedFlavour(const size_t &idij);

    inline void SetScaleSetter(PHASIC::Scale_Setter_Base *const scale)
    { p_scale=scale; }

    inline Amplitude *GetAmplitude() const { return p_bg; }

    inline Single_Process *Process() const { return p_proc; }

    inline NLO_subevt *Sub() const { return p_sub; }

  };// end of class Single_Dipole_Term

}// end of namespace COMIX

#endif



