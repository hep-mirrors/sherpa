#ifndef PHASIC_Main_Extra_Emission_Generator_H
#define PHASIC_Main_Extra_Emission_Generator_H

#include "PHASIC++/Channels/CS_Dipole.H"

namespace PHASIC {

  class Process_Integrator;

  class Extra_Emission_Generator {
  private:

    Process_Integrator *p_int;

    CSDipole_Vector m_dipoles;

    size_t m_nin, m_opt, m_omode;

    double m_weight;

    bool m_isqcd;

    bool AddDipole(CS_Dipole *dip);

    bool InitDipoles(Process_Base *const proc);
    bool InitEWDipoles(Process_Base *const proc);

  public:

    // constructor/destructor
    Extra_Emission_Generator(Process_Integrator *const pint, bool isqcd);
    ~Extra_Emission_Generator();

    // members
    bool GeneratePoint(ATOOLS::Cluster_Amplitude *const ampl);
    bool GenerateWeight(const ATOOLS::Cluster_Amplitude *ampl);

    void AddPoint(const double &value);

    void Optimize();
    void EndOptimize();

    void WriteOut(std::string pid);
    void ReadIn(std::string pid);

    void SetAMin(const double &amin);

    inline double Weight() const { return m_weight; }

    inline size_t Opt() const { return m_opt; }

  };// end of class Extra_Emission_Generator

}// end of namespace PHASIC

#endif
