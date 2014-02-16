#ifndef PHASIC_Main_BBar_Emission_Generator_H
#define PHASIC_Main_BBar_Emission_Generator_H

#include "PHASIC++/Channels/CS_Dipole.H"
#include "PHASIC++/Process/Process_Base.H"

namespace ATOOLS { class Integration_Info; }

namespace PHASIC {

  class BBar_Process;
  class Phase_Space_Handler;

  struct Dipole_Params {
    CS_Dipole *p_dip;
    Process_Vector m_procs;
    ATOOLS::Vec4D_Vector m_p;
    double m_weight;
    inline Dipole_Params(CS_Dipole *const dip=NULL,
			 const Process_Vector &procs=Process_Vector(),
			 const ATOOLS::Vec4D_Vector &p=ATOOLS::Vec4D_Vector(),
			 const double weight=0.0):
      p_dip(dip), m_procs(procs), m_p(p), m_weight(weight) {}
  };// end of struct Dipole_Params

  std::ostream &operator<<(std::ostream &ostr,const Dipole_Params &dp);

  class BBar_Emission_Generator {
  private:

    CSDipole_Vector  m_dipoles;
    CS_Dipole       *p_active;

    ATOOLS::Vec4D_Vector m_p;

    size_t m_nin, m_opt, m_omode;
    double m_weight, m_amin, m_Q2min;

    std::map<Process_Base*,CS_Dipole*> m_dmap;

    std::map<CS_Dipole*,std::map<Process_Base*,Process_Vector> > m_pmap;

    bool AddDipole(Process_Base *const bviproc,CS_Dipole *const dip);

  public:

    BBar_Emission_Generator();

    ~BBar_Emission_Generator();

    // members
    bool InitDipoles(Process_Base *const bviproc,Process_Base *const sproc,
		     Phase_Space_Handler *const psh);

    bool GeneratePoint(const ATOOLS::Vec4D_Vector &p,Cut_Data *const cuts);
    bool GenerateWeight(Cut_Data *const cuts, bool activeonly=false);

    void AddPoint(const double &value);

    void Optimize();
    void EndOptimize();
    void MPISync();

    void WriteOut(std::string pid);
    bool ReadIn(std::string pid);

    void Print();

    Dipole_Params Active(Process_Base *const bviproc) const;

    inline double Weight() const { return m_weight; }

  };// end of class BBar_Emission_Generator

}// end of namespace PHASIC

#endif
