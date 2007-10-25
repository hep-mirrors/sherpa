#ifndef Gluon_Decayer_H
#define Gluon_Decayer_H

#include "Dipole_Splitter.H"
#include "Strong_Coupling.H"
#include "Cluster.H"
#include "Blob.H"

namespace AHADIC {
  class Gluon_Decayer {
  private:
    double              m_pt2min;
    Dipole_Splitter   * p_splitter;
    Strong_Coupling   * p_as;

    FlavDecayMap        m_options;
    DipoleList          m_dipoles;

    bool    Shift(Proto_Particle_List * plin=NULL);
    bool    FillDipoleList(Proto_Particle_List * plin=NULL);
    void    PrintDipoleList();

    bool    DecayDipoles();
    DipIter SelectDipole();
  public:
    Gluon_Decayer(Dipole_Splitter *);
    ~Gluon_Decayer();
    bool DecayList(Proto_Particle_List * plin=NULL);
  };
}
#endif
