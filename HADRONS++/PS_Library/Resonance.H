#ifndef HADRONS_PS_Library_Resonance_H
#define HADRONS_PS_Library_Resonance_H

#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Lookup_Tables.H"
#include <list>
#include <string>
#include <vector>

namespace HADRONS {
  enum class resonance_type {
    fixed    = 0,
    running  = 1,
    GS       = 2,
    fixed3   = 10,
    running3 = 11,
    GS3      = 12,
    bespoke  = 20
  };

  struct Res_Params {
    ATOOLS::Flavour              m_inflav;
    std::vector<ATOOLS::Flavour> m_outflavs;
    resonance_type               m_type;
    double m_OSmass, m_OSwidth, m_weight, m_phase;

    Res_Params() {}
    Res_Params(const ATOOLS::Flavour & inflav,
	       const std::vector<ATOOLS::Flavour> & outflavs,
	       const resonance_type & type=resonance_type::running,
	       const double & OSmass=-1.,const double & OSwidth=-1.,
	       const double & weight=1.,const double & phase=0.
	       ) :
    m_inflav(inflav), m_outflavs(outflavs), m_type(type),
    m_OSmass(OSmass), m_OSwidth(OSwidth), m_weight(weight), m_phase(phase) {}
  };
  
  typedef std::list<Res_Params> ResonanceParameters;
  typedef std::map<int,
		   std::map<int,
			    std::map<long unsigned int,
				     std::vector<double> > > > ResDefaults;
  
  
  class Resonance_Base {
  protected:
    ATOOLS::Flavour m_inflav;
    resonance_type  m_type;
    double          m_OSmass, m_OSmass2, m_OSwidth, m_OSwidth2, m_weight, m_phase;
    double          m_threshold, m_threshold2, m_lambda;
    int             m_exponent;
    std::string     m_name;
    std::vector<double> m_decmasses2;

    virtual double CalculateWidth(const double & s) = 0;
  public:
    Resonance_Base(const Res_Params & params);
    ~Resonance_Base() {}
    virtual const ATOOLS::Flavour & Flav()           { return m_inflav; }
    virtual const double Mass(const double & s=0.)   { return m_OSmass; }
    virtual const double Mass2(const double & s=0.)  { return m_OSmass2; }
    virtual const double Width(const double & s=0.)  { return CalculateWidth(s); }
    virtual const double Width2(const double & s=0.) { return ATOOLS::sqr(Width(s)); }
    virtual const double Weight()                    { return m_weight; }
    virtual const double Phase()                     { return m_phase; }
    
    inline double Lambda(const double & M2,const double & m12,const double & m22) {
      return sqrt( ATOOLS::Max(0., ATOOLS::sqr(M2-m12-m22)-4.*m12*m22)/(4.*M2) );
    }
    virtual Complex BreitWigner(const double & s,const double & M2,const double & MW) {
      return M2/Complex( M2-s, -MW );
    }
    virtual Complex AltBreitWigner(const double & s,const double & M2,const double & MW) {
      return s/Complex( M2-s, -MW );
    }
    virtual Complex BreitWigner(const double & s) = 0;
  };

  class Resonance_List : public std::list<std::pair<double, Resonance_Base *> > {
  public:
    Resonance_List(Resonance_Base * start=NULL) {
      if (start) push_back(std::make_pair(1.,start));
    }
    inline ~Resonance_List() { while (!empty()) { delete back().second; pop_back(); } }
    inline void Add(Resonance_Base * R) { push_back(std::make_pair(R->Weight(),R)); }
    inline Complex BreitWigner(const double & s) {
      Complex res(1.,0.);
      for (std::list<std::pair<double, Resonance_Base *> >::iterator rit=begin();
	   rit!=end(); rit++) {
	res *= rit->first * rit->second->BreitWigner(s);
      }
      return res;
    }
  };
  
  class Resonance_Tree {
  private:
    Resonance_Base  * p_start;
    Resonance_List *  p_nodes;
    double            m_norm;
    std::string       m_name;
  public:
    Resonance_Tree(Resonance_Base * start=NULL,const double & norm=1.0,
		   const std::string name=std::string("")) :
      p_start(start), p_nodes(NULL), m_norm(norm), m_name(name) {}
    ~Resonance_Tree() {
      delete p_nodes;
      delete p_start;
    }
    
    Complex operator()(const double & s,const double & s1);
    void    Output();
    inline void Add(Resonance_Base * R) {
      if (R->Weight()<1.e-8) return;
      if (!p_nodes) p_nodes = new Resonance_List();
      p_nodes->Add(R);
    }
    inline void AddNodes(Resonance_List * nodes) { p_nodes = nodes; }
    inline const std::string & Name()            { return m_name; }
    inline Resonance_List * GetNodes()           { return p_nodes; }
  };

  
  class FixedWidth_Resonance : public Resonance_Base {
    double m_MGamma;
    double CalculateWidth(const double & s) { return m_OSwidth; }
  public:
    FixedWidth_Resonance(const Res_Params & params) :
      Resonance_Base(params), m_MGamma(m_OSmass*m_OSwidth) {}
    Complex BreitWigner(const double & s) {
      return Resonance_Base::BreitWigner(s,m_OSmass2,m_MGamma);
    }
  };
  
  class RunningWidth_Resonance : public Resonance_Base {
    double CalculateWidth(const double & s);
  public:
    RunningWidth_Resonance(const Res_Params & params);
    Complex BreitWigner(const double & s); 
  };


  class GS_Resonance : public Resonance_Base {
    double m_m_pi, m_d, m_h0, m_dh0;
    double f(const double & s);
    double h(const double & s);
    double CalculateWidth(const double & s);
  public:
    GS_Resonance(const Res_Params & params);
    Complex BreitWigner(const double & s); 
  };

  
  class RunningWidth3_Resonance : public Resonance_Base {
  private:
    Resonance_Base * p_BW_ij, * p_BW_kj;
    double           m_mi2, m_mj2, m_mk2, m_allm2;
    double           m_delta2_ij, m_sum2_ij, m_delta2_kj, m_sum2_kj;
    double           m_smin_ij, m_smin_kj, m_smin_ik;
    size_t           m_nbins, m_nsbins;
    ATOOLS::OneDim_Table * p_g;
    double           m_gV;
    
    bool   GetParams(std::vector<ATOOLS::Flavour> & flavs,Res_Params & V);
    double dg(const double & Q2,const double & s1,const double & s2);
    void   FillTable(const double & Qmax=-1.);
    double CalculateWidth(const double & s);
  public:
    RunningWidth3_Resonance(const Res_Params & params);
    ~RunningWidth3_Resonance();
    Complex BreitWigner(const double & s);
    void    SetResonance(const size_t & pi,Resonance_Base * res);
    inline const double g(const double & Q2) { return (*p_g)(Q2); }
    inline const double gV()                 { return m_gV; }
  };
}
#endif
