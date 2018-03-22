#ifndef EXTAMP_Process_H
#define EXTAMP_Process_H

#include "PHASIC++/Process/Process_Info.H"
#include "PHASIC++/Process/Single_Process.H"
#include "EXTAMP/CS_Dipole.H"

namespace EXTAMP {

  class Process : public PHASIC::Single_Process {

  public :

    Process(const PHASIC::Process_Info& pi);
    
    ~Process() { };

    /* Inherited from PHASIC::Process_Base. Required for PS
       clustering. Also include a const-correct definition and
       implementation. */
    bool Combinable(const size_t &idi,const size_t &idj);
    
    /* Inherited from PHASIC::Process_Base. Required for PS
       clustering. */
    const ATOOLS::Flavour_Vector &CombinedFlavour(const size_t &idij);

    virtual int  PerformTests();

  protected:

    /* Factor due to spin and color averaging in the initial state,
       i.e. number of unobserved degrees of freedom in initial state
       (polarizations, colors) */
    double NormFac() const { return m_norm; }
    double m_norm;
    EXTAMP::DipoleCase m_dipole_case;
    
    const std::vector<size_t>& PartonIndices() const
    { return m_parton_indices; }

  private :

    /* Given indices of external flavours i_1,i_2,i_3,... and the
       combined index k=(1<<i_1) & (1<<i_2) & (1<<i_3) & ...,
       m_cluster_flav_map[k] contains all flavours that arise from
       diagrammatic clusterings of the external flavours */
    std::map<size_t, ATOOLS::Flavour_Vector>  m_cluster_flav_map;

    /* Fill m_parton_indices */
    void FillPartonIndices();

    /* Holds a list of indices i, such that all m_flavs[i] are
       massless QCD partons */
    std::vector<size_t> m_parton_indices;
    
  };

}

#endif
