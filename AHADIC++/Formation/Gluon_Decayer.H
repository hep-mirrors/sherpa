#ifndef AHADIC_Formation_Gluon_Decayer_H
#define AHADIC_Formation_Gluon_Decayer_H

#include "AHADIC++/Tools/Gluon_Splitter.H"
#include "AHADIC++/Tools/Dipole.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Math/Histogram.H"
#include <map>

namespace AHADIC {
  class Gluon_Decayer {
  private:
    Gluon_Splitter  * p_gsplitter;
    double          m_pt2max, m_pt2max_factor, m_pt02;
    FlavDecayMap    m_options;
    DipoleList      m_dipoles;

    bool                                      m_analyse;
    std::map<std::string,ATOOLS::Histogram *> m_histograms;

    void    Init();
    bool    FillDipoleList(Proto_Particle_List * plin=NULL);
    void    UpdatePPList(Proto_Particle_List * plin=NULL);
    void    PrintDipoleList();

    bool    DecayDipoles();
    DipIter SelectDipole();
    int     Rescue(DipIter & dip);
    void    AfterSplit(DipIter dip);
    bool    MergeDipoles(DipIter & dip1,DipIter & dip2);
    void    SplitIt(DipIter dip,ATOOLS::Vec4D=ATOOLS::Vec4D(0.,0.,0.,0.));
    double  PT2Max(Dipole * dip) const;

    long int m_tot, m_s, m_u, m_d, m_di0, m_di1;
  public:
    Gluon_Decayer(bool ana=false);
    ~Gluon_Decayer();
    std::map<std::string,ATOOLS::Histogram *> GetHistos() {return m_histograms;};
    bool DecayList(Proto_Particle_List * plin=NULL);
  };
}
#endif
