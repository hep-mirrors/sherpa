#ifndef EXTAMP_RS_Process_H
#define EXTAMP_RS_Process_H

#include "EXTAMP/Process.H"

namespace PHASIC {
  class Tree_ME2_Base;
}

namespace ATOOLS {
  class NLO_subevtlist;
}

namespace EXTAMP {

  class CS_Dipole;
  class Dipole_Wrapper_Process;

  class RS_Process : public Process {

    friend class Dipole_Wrapper_Process;

    typedef std::vector<CS_Dipole*>              Dipole_Vector;
    typedef std::vector<Dipole_Wrapper_Process*> Dipole_Wrappers;

  public:

    RS_Process(const PHASIC::Process_Info& pi);
    
    ~RS_Process();

    double Partonic(const ATOOLS::Vec4D_Vector &p,
		    const int mode);

    /* Inherited virtual from PHASIC::Process_Base. Returns NULL by
       default and needs to be overwritten in such a way as to return
       a list to ATOOLS::NLO_subevts representing the differential
       dipole terms */
    ATOOLS::NLO_subevtlist *GetSubevtList() { return &m_subevents; }
    ATOOLS::NLO_subevtlist *GetRSSubevtList() { return &m_subevents; }

    /* Inherited virtual from PHASIC::Process_Base. Need to overwrite
       this to generically return true because real emission term and
       subtraction terms have to be triggered separately. This is done
       in RS::Process::Partonic. */
    bool Trigger(const ATOOLS::Vec4D_Vector &p) { return true; }

    /* Overloaded from PHASIC::Process_Base such that dipole dummy
       processes are initialized as well */
    void Init(const PHASIC::Process_Info &pi,
	      BEAM::Beam_Spectra_Handler *const beamhandler,
	      PDF::ISR_Handler *const isrhandler,const int mode=0);

    /* Overload the following methods of PHASIC::Process_Base in such
       a way as to forward calls to dipole wrappers as well. */
    void SetScale(const PHASIC::Scale_Setter_Arguments &args);
    void SetKFactor(const PHASIC::KFactor_Setter_Arguments &args);
    void SetSelector(const PHASIC::Selector_Key &key);
    void SetNLOMC(PDF::NLOMC_Base* const nlomc);
    void SetShower(PDF::Shower_Base *const ps);
    void SetGenerator(PHASIC::ME_Generator_Base *const gen);
    void FillProcessMap(PHASIC::NLOTypeStringProcessMap_Map *apmap);
    size_t SetMCMode(const size_t mcmode);

  private:

    PHASIC::Tree_ME2_Base* p_real_me;

    /* Construct a complete vector of CS dipoles. m_dipoles holds
       pointers instead of objects because EXTAMP::Dipole is
       abstract */
    Dipole_Vector ConstructDipoles();
    Dipole_Vector m_dipoles;

    /* Legacy structure: Sherpa needs a representation of dipoles in
       terms of PHASIC::Single_Processes in the MC@NLO
       implementation. */
    Dipole_Wrappers ConstructDipoleWrappers(const Dipole_Vector& dv) const;
    Dipole_Wrappers m_dipole_wrappers;
    void DeleteDipoleWrappers();

    /* Legacy structure: Sherpa needs a representation of dipoles in
       terms of ATOOLS::NLO_subevents. Since this is a vector of
       stupid raw pointers, need some memory management here as
       well. */
    ATOOLS::NLO_subevtlist ConstructSubevents(size_t n_dipoles) const;
    ATOOLS::NLO_subevtlist m_subevents;
    void DeleteSubevents();

    /* Kinematics of dipole terms and update the kinematics of
       NLO_subevents accordingly, thus making the momenta available
       for e.g. the scale setter */
    void CalculateKinematics(const ATOOLS::Vec4D_Vector& p);

    /* Check if any of the dipole misses the minimum alpha
       criterion */
    bool PassesAlphaMin(const Dipole_Vector& dv) const;

    /* If alpha_min cut is not passed, use this method to set all
       subevents to zero */
    void SetSubEventsToZero(ATOOLS::NLO_subevtlist subevents) const;

    /* Populate running coupling multimap of Coupling_Data's (one for
       each NLO_subevent) */
    void ConstructRunningCouplings(MODEL::Coupling_Map& cpls,
				   const ATOOLS::NLO_subevtlist& evts,
				   const Dipole_Vector& dipoles) const;

    static void SmearSubEvents(const Dipole_Vector& dipoles,
			       ATOOLS::NLO_subevtlist& subs,
			       const double& alpha_0,
			       const double& power);
    double m_alpha_0, m_smear_power;

  };

}

#endif
