#ifndef SHRiMPS_Main_Cluster_Algorithm_H
#define SHRiMPS_Main_Cluster_Algorithm_H

#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "PDF/Main/Cluster_Definitions_Base.H"
#include "ATOOLS/Math/Histogram.H"
#include <map>
#include <list>

namespace ATOOLS { 
  class Blob;
  class Particle;
}

namespace PDF { 
  class Cluster_Definitions_Base; 
}

namespace SHRIMPS {
  typedef std::list<ATOOLS::Particle *> ParticleList;

  class Cluster_Algorithm: public ATOOLS::Mass_Selector {
    class JF {
    private:
      double m_ycut;
    public:
      JF(const double & kt2) : 
	m_ycut(kt2/ATOOLS::sqr(ATOOLS::rpa->gen.Ecms())) {}
      ~JF() {};
      inline const double & Ycut() const { return m_ycut; }
      inline double DR()   const { return 1.; }
    };
  private:
    int           m_mode,m_nlad;
    bool          m_resc;
    ATOOLS::Vec4D m_rescvec;
    double        m_minkt2,m_showerfac,m_tmax;

    std::map<std::string, ATOOLS::Histogram * > m_histomap;

    ATOOLS::Cluster_Amplitude     *p_ampl;
    PDF::Cluster_Definitions_Base *p_clus;
    ATOOLS::Mass_Selector         *p_ms;
    JF * p_jf;

    int  ColorConnected(const ATOOLS::ColorID &i,
			const ATOOLS::ColorID &j) const;
    void ProjectOnSinglets(ATOOLS::Blob *const blob,
    			   std::list<ParticleList *> & singlets);
    double PTij2(const ATOOLS::Vec4D & pi,const ATOOLS::Vec4D & pj) const;
    double PTi2(const ATOOLS::Vec4D & pi,const ATOOLS::Vec4D & pbeam) const;
  public:
    Cluster_Algorithm();

    ~Cluster_Algorithm();

    bool Cluster(ATOOLS::Blob *const blob);
    
    // inline functions
    inline ATOOLS::Cluster_Amplitude *Amplitude() { return p_ampl; }
    inline void SetShowerParams(const int & mode,const double & minkt2) {
      m_mode   = mode; 
      m_minkt2 = minkt2;
      if (m_mode<2) m_minkt2 = 0.;
    }
    inline void SetShowerFac(const double & sfac) {
      m_showerfac = sfac; 
    }
    inline void 
    SetClusterDefinitions(PDF::Cluster_Definitions_Base *const cb) 
    { p_clus=cb; }
    inline void SetMinKT2(const double & minkt2) { m_minkt2 = minkt2; }
    inline void SetRescatt(const double & resc)  { m_resc = resc; }
    inline void SetTMax(const double & tmax)  { m_tmax = tmax; }
    inline void SetNLad(const int & nlad)  { m_nlad = nlad; }

    double Mass(const ATOOLS::Flavour &fl) const;
  };// end of class Cluster_Algorithm

}// end of namespace SHRiMPS

#endif
