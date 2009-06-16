#ifndef AHADIC_Formation_Gluon_Decayer_H
#define AHADIC_Formation_Gluon_Decayer_H

#include "AHADIC++/Tools/Dipole_Splitter.H"
#include "AHADIC++/Tools/Strong_Coupling.H"
#include "AHADIC++/Tools/Cluster.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Math/Histogram.H"
#include <map>

namespace AHADIC {
  class Gluon_Decayer {
  private:
    double              m_pt2min, m_pt2max;
    Dipole_Splitter   * p_splitter;

    FlavDecayMap        m_options;
    DipoleList          m_dipoles;

    bool                                      m_analyse;
    std::map<std::string,ATOOLS::Histogram *> m_histograms;

    bool    FillDipoleList(SP(Proto_Particle_List) plin=NULL);
    void    UpdatePPList(SP(Proto_Particle_List) plin=NULL);
    void    PrintDipoleList();

    bool    DecayDipoles();
    DipIter SelectDipole();
    int     Rescue(DipIter & dip);
    void    AfterSplit(DipIter dip);
    void    MergeDipoles(DipIter & dip1,DipIter & dip2);
    void    SplitIt(DipIter dip,ATOOLS::Vec4D=ATOOLS::Vec4D(0.,0.,0.,0.));
  public:
    Gluon_Decayer(Dipole_Splitter *,bool=false);
    ~Gluon_Decayer();
    bool DecayList(SP(Proto_Particle_List) plin=NULL);
  };
}
#endif
