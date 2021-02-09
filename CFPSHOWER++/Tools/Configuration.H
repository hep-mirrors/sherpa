#ifndef CFPShower_Tools_Configuration_H
#define CFPShower_Tools_Configuration_H

#include "CFPSHOWER++/Shower/Cluster_Definitions.H"
#include "CFPSHOWER++/Tools/Parton.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"

namespace CFPSHOWER {
  class Configuration : public Parton_List {
    ATOOLS::Cluster_Amplitude * p_ampl;
    double        m_t, m_t0, m_tveto;
    ATOOLS::Vec4D m_pos;
  private:
    void Fill(std::map<ATOOLS::Cluster_Leg*,Parton*> & lmap);
    void EstablishRelations();
    void EstablishHistories(std::map<ATOOLS::Cluster_Leg*,Parton*> & lmap,
			    Cluster_Definitions * cluster);
  public:
    Configuration(const double & t, const double & t0);
    Configuration(ATOOLS::Cluster_Amplitude * const ampl,
		  std::map<ATOOLS::Cluster_Leg*,Parton*> & lmap,
		  Cluster_Definitions * cluster);
    ~Configuration();

    
    inline void AddParton(Parton * const parton) { push_back(parton); }
    inline void RemoveParton(Parton * parton)    { remove(parton); }
    inline const ATOOLS::Vec4D & Position()      { return m_pos; }
    inline void SetT(const double & t)           { m_t     = t; }
    inline void SetT0(const double & t0)         { m_t0    = t0; }
    inline void SetT_veto(const double & tveto)  { m_tveto = tveto; }
    inline const double & T() const              { return m_t; }
    inline const double & T0() const             { return m_t0; }
    inline const double & T_veto() const         { return m_tveto; }
    inline void SetPosition(const ATOOLS::Vec4D & pos=
			    ATOOLS::Vec4D(0.,0.,0.,0.)) { m_pos = pos; }
    size_t NPartons();
  };

  std::ostream &operator<<(std::ostream & s,const Configuration & configuration);
  typedef std::vector<Configuration*> Configuration_Vector;
}

#endif
