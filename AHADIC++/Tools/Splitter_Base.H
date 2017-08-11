#ifndef AHADIC_Tools_Splitter_Base_H
#define AHADIC_Tools_Splitter_Base_H

#include "AHADIC++/Tools/Proto_Particle.H"
#include "MODEL/Main/Strong_Coupling.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Math/Histogram.H"

namespace AHADIC {
  struct  DecaySpecs { double weight,massmin,popweight; };
  typedef std::map<ATOOLS::Flavour,DecaySpecs *> FlavDecayMap;
  typedef FlavDecayMap::iterator                 FDIter;

  struct LightCone {
  public:
    double m_smandel, m_E;
    double m_msplit2, m_msplit, m_mspect2, m_mspect, m_alpha, m_beta;
    ATOOLS::Vec4D m_pA, m_pB;
    LightCone() :
      m_smandel(0.), m_E(0.), 
      m_msplit2(0.), m_msplit(0.), m_mspect2(0.), m_mspect(0.), 
      m_alpha(0.), m_beta(0.),
      m_pA(ATOOLS::Vec4D(0.,0.,0.,0.)), m_pB(ATOOLS::Vec4D(0.,0.,0.,0.)) {}
  };

  struct PoppedPair {
  public:
    ATOOLS::Flavour  m_flav;
    double           m_mpop2;
    double           m_x, m_y, m_z, m_kt2, m_sqq; 
    ATOOLS::Vec4D    m_outmom[2];
    PoppedPair() :
      m_flav(ATOOLS::Flavour(kf_none)), m_mpop2(0.),
      m_x(0.), m_y(0.), m_z(0.), m_kt2(0.), m_sqq(0.) 
    {
      for (size_t i(0);i<2;i++) m_outmom[i] = ATOOLS::Vec4D(0.,0.,0.,0.);
    }
    void ResetKin() {
      m_x = m_y = m_z = m_kt2 = m_sqq = m_mpop2 = 0.;
    }
  };
  std::ostream & operator<<(std::ostream & s, const PoppedPair &);

  class Splitter_Base {
  protected:
    MODEL::Strong_Coupling * p_as;
    double       m_pt2max, m_pt02;
    FlavDecayMap m_options;
    double       m_mmin, m_mmin2, m_mmax, m_mmax2, m_sumwt;

    Proto_Particle * p_split, * p_spect;
    LightCone        m_LC;
    bool             m_swap, m_anti, m_leadspect, m_leadsplit, m_isbeam;
    ATOOLS::Poincare m_boost, m_rotat;
    ATOOLS::Vec4D    m_splitmom,m_spectmom;

    double           m_sumx, m_sumy;
    PoppedPair     * p_pair;  
    std::list<PoppedPair *> m_popped;

    bool m_ana;
    std::string m_anapath;
    std::map<std::string, ATOOLS::Histogram *> m_histograms;
    long int m_pop,m_popu,m_popd,m_pops,m_popud0,m_popsu0,m_popsd0;
    long int m_popuu1,m_popud1,m_popdd1,m_popsu1,m_popsd1,m_popss1;

    virtual void   Init();
    virtual void   Reset();
    virtual void   DefineTags();
    virtual void   ConstructTrafos();
    virtual void   UndoTrafos();
    virtual bool   ConstructLightC();
    virtual bool   AlphaBeta(const double & s,double & alpha,double & beta);
    virtual double SelectY(const double & ymin,const double & ymax,
			   const double & eta,const double & offset=0);
    virtual double SelectZ(const double & delta,const bool & lead);
    virtual bool   SelectFlavour(const double & sqq,const bool & vetodi=false);
    virtual double CalculateSumWT(const double & m_mmax=-1.,
				  const bool & vetodi=false);
    double Lambda(const double & ma2,const double & mb2,const double & mc2);

    virtual bool PoppedMassPossible(const double & m2)=0;

    virtual void   InitAnalysis();
    virtual void   Analysis();
  public:
    Splitter_Base();
    virtual ~Splitter_Base();
  };

  inline double Splitter_Base::
  Lambda(const double & ma2,const double & mb2,const double & mc2) {
    return sqrt(ATOOLS::sqr(ma2-mb2-mc2)-4.*mb2*mc2)/(2.*ma2);
  }
}

#endif
