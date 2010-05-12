#ifndef PHASIC_Process_NLO_Helpers_H
#define PHASIC_Process_NLO_Helpers_H

#include "ATOOLS/Phys/Particle_List.H"

namespace PHASIC {

  struct nlo_type {
    enum code {
      lo   = 0,
      born = 1,
      loop = 2,
      vsub = 4,
      real = 8,
      rsub = 16,
      polecheck = 32
    };
  };

  inline const nlo_type::code operator|(const nlo_type::code c1,const nlo_type::code c2)  
  { return (nlo_type::code)((int)c1|(int)c2); }
  inline const nlo_type::code &operator|=(nlo_type::code &c1,const nlo_type::code c2)  
  { return c1=(nlo_type::code)((int)c1|(int)c2); }
  inline const nlo_type::code operator&(const nlo_type::code c1,const nlo_type::code c2)  
  { return (nlo_type::code)((int)c1&(int)c2); }
  inline const nlo_type::code &operator&=(nlo_type::code &c1,const nlo_type::code c2)  
  { return c1=(nlo_type::code)((int)c1&(int)c2); }

  std::ostream &operator<<(std::ostream &str,const nlo_type::code &c);
  std::istream &operator>>(std::istream &str,nlo_type::code &c);

  class NLO_subevt {
  public:
    size_t n;
    const ATOOLS::Flavour* p_fl;
    const ATOOLS::Vec4D*   p_mom;
    bool m_flip;
    double m_result,m_me,m_scale,m_alpha;
    std::string m_ID;

    NLO_subevt() {n=0;p_fl=0;p_mom=0;m_ID=std::string("");m_result=0.;m_scale=0.;m_alpha=0.;}

    NLO_subevt(size_t nn,const ATOOLS::Flavour* fl,const ATOOLS::Vec4D* mom,std::string ID)
    { n=nn; p_fl=fl; p_mom=mom; m_ID=ID; m_result=0.;m_me=0.;m_scale=0.;m_alpha=0.; }

    NLO_subevt(NLO_subevt* nlos)
    { *this = *nlos; m_result=m_me; }

    void Print() {std::cout<<m_ID<<": "<<m_result<<std::endl;}
  
    ATOOLS::Particle_List* CreateParticleList();

    void Mult(const double scal) {
      m_result*=scal;
      m_me*=scal;
    }
    
    NLO_subevt& operator*= (const double scal) {
      m_result*=scal;
      return *this;
    }

    NLO_subevt& operator+= (const double sm) {
      m_result+=sm;
      return *this;
    }

    NLO_subevt& operator=(const NLO_subevt& cp) {
      if (this!=&cp) {
	n        = cp.n;
	p_fl     = cp.p_fl;
	p_mom    = cp.p_mom;
	m_result = cp.m_result;
	m_me     = cp.m_me;
	m_scale  = cp.m_scale;
	m_alpha  = cp.m_alpha;
	m_ID     = cp.m_ID;
      }
      return *this;
    }
    friend std::ostream &operator<< (std::ostream &ostr,const  NLO_subevt *sevt)
    {
      return ostr<<sevt->m_ID<<"\n";
    }
  };

  class NLO_subevtlist : public std::vector<NLO_subevt*> { 
  public:
    double m_muf2;
    NLO_subevtlist& operator*=(const double scal);
    void Print();
    void Mult(const double scal);
  };



}

#endif
