#ifndef PHASIC_Main_Extra_Emission_Generator_H
#define PHASIC_Main_Extra_Emission_Generator_H

#include "PHASIC++/Channels/CS_Dipole.H"

namespace ATOOLS { class Integration_Info; }

namespace PHASIC {

  class Process_Base;
  class POWHEG_Process;
  class Phase_Space_Handler;

  class Extra_Emission_Generator {
  private:

    CSDipole_Vector  m_dipoles;
    CS_Dipole       *p_active;

    size_t m_nin, m_opt, m_omode;
    double m_weight, m_asum, m_amin, m_Q2min;

    ATOOLS::Integration_Info    *p_info;
    std::map<size_t,CS_Dipole*>  m_dmap;

    POWHEG_Process *p_proc;

    bool m_numtrig;

    bool AddDipole(POWHEG_Process *const proc,CS_Dipole *const dip);

  public:

    Extra_Emission_Generator();

    ~Extra_Emission_Generator();

    // members
    bool InitDipoles(POWHEG_Process *const proc,Process_Base *const sproc,
		     Phase_Space_Handler *const psh);

    ATOOLS::Vec4D_Vector GeneratePoint
    (const ATOOLS::Vec4D_Vector &p,Cut_Data *const cuts);

    bool GenerateWeight(const ATOOLS::Vec4D_Vector &p,
			Cut_Data *const cuts, bool activeonly=false);

    double SelectionWeight(const size_t &idx) const;

    void AddPoint(const double &value);

    void Optimize();
    void EndOptimize();
    void MPISync();

    void WriteOut(std::string pid);
    void ReadIn(std::string pid);

    void Print();

    inline double Weight() const { return m_weight; }

    inline size_t Opt() const { return m_opt; }

    inline CS_Dipole *Active() const { return p_active; }

    inline const bool NumericsTrigger() const { return m_numtrig; }

  };// end of class Extra_Emission_Generator

}// end of namespace PHASIC

#endif
