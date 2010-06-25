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

  inline nlo_type::code operator|(const nlo_type::code c1,const nlo_type::code c2)  
  { return (nlo_type::code)((int)c1|(int)c2); }
  inline const nlo_type::code &operator|=(nlo_type::code &c1,const nlo_type::code c2)  
  { return c1=(nlo_type::code)((int)c1|(int)c2); }
  inline nlo_type::code operator&(const nlo_type::code c1,const nlo_type::code c2)  
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
    double m_result,m_me,m_mewgt,m_facscale,m_renscale,m_alpha;
    std::string m_ID;

    NLO_subevt();

    NLO_subevt(size_t nn,const ATOOLS::Flavour* fl,const ATOOLS::Vec4D* mom,std::string ID);

    NLO_subevt(NLO_subevt* nlos);

    void Print();
  
    ATOOLS::Particle_List* CreateParticleList();

    void Mult(const double scal);
    
    void MultMEwgt(const double scal);

    NLO_subevt& operator*= (const double scal);

    NLO_subevt& operator+= (const double sm);

    NLO_subevt& operator=(const NLO_subevt& cp);

    friend std::ostream &operator<< (std::ostream &ostr,const  NLO_subevt *sevt) {
      return ostr<<sevt->m_ID<<"\n";
    }

  };

  class NLO_subevtlist : public std::vector<NLO_subevt*> { 
  public:
    double m_muf2;
    NLO_subevtlist& operator*=(const double scal);
    void Print();
    void Mult(const double scal);
    void MultMEwgt(const double scal);
  };



}

#endif
