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

    FlavDecayMap        m_options;
    DipoleList          m_dipoles;

    bool    FillDipoleList(Proto_Particle_List * plin=NULL);
    void    UpdatePPList(Proto_Particle_List * plin=NULL);
    void    PrintDipoleList();

    bool    DecayDipoles();
    DipIter SelectDipole();
    bool    Rescue(DipIter & dip);
    void    AfterSplit(DipIter dip);
    void    MergeDipoles(DipIter & dip1,DipIter & dip2);
    void    SplitIt(DipIter dip,ATOOLS::Vec4D=ATOOLS::Vec4D(0.,0.,0.,0.));
  public:
    Gluon_Decayer(Dipole_Splitter *);
    ~Gluon_Decayer();
    bool DecayList(Proto_Particle_List * plin=NULL);
  };
}
#endif
