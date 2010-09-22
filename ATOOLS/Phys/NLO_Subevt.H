#ifndef ATOOLS_Process_NLO_Subevt_H
#define ATOOLS_Process_NLO_Subevt_H

#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Phys/Decay_Info.H"

namespace ATOOLS {

  class ME_wgtinfo {
  public: 
    int m_nx;
    double m_w0, *p_wx;
    double m_x1, m_x2, m_y1, m_y2, m_mur2;
    ME_wgtinfo();
    ~ME_wgtinfo();
    ME_wgtinfo& operator*=(const double &scal);
    void Flip();
    void AddMEweights(int n);
  };

  struct IDip_ID {
    size_t m_ijt, m_kt;
  public:
    inline IDip_ID(const size_t &ijt=0,const size_t &kt=0):
      m_ijt(ijt), m_kt(kt) {}
    bool operator<(const IDip_ID &di) const;
  };// end of struct IDip_ID

  struct DDip_ID {
    size_t m_i, m_j, m_k;
  public:
    inline DDip_ID(const size_t &i=0,
		   const size_t &j=0,const size_t &k=0):
      m_i(i), m_j(j), m_k(k) {}
    bool operator<(const DDip_ID &di) const;
  };// end of struct DDip_ID

  struct Dip_ID: public IDip_ID, public DDip_ID {
    inline Dip_ID(const size_t &ijt=0,const size_t &kt=0,
		  const size_t &i=0,const size_t &j=0,const size_t &k=0):
      IDip_ID(ijt,kt), DDip_ID(i,j,k) {}
    bool operator<(const Dip_ID &di) const;
  };// end of struct Dip_ID

  std::ostream &operator<<(std::ostream &ostr,const Dip_ID &di);

  class NLO_subevt: public Dip_ID {
  public:

    const Flavour *p_fl;
    const Vec4D   *p_mom;

    const size_t     *p_id;
    DecayInfo_Vector *p_dec;

    NLO_subevt *p_real;
    void       *p_proc;

    size_t m_n, m_nco, m_flip;
    double m_result, m_last[2], m_me, m_mewgt, m_muf2, m_mur2;

    std::string m_pname;

  public:

    inline NLO_subevt(const size_t &n=0,const size_t *id=NULL,
		      const Flavour *fl=NULL,const Vec4D *mom=NULL,
		      const int i=-1,const int j=-1,const int k=-1):
      Dip_ID(0,0,i,j,k), p_fl(fl), p_mom(mom), p_id(id), p_dec(NULL),
      p_real(NULL), p_proc(NULL), m_n(n), m_nco(0), m_flip(0),
      m_result(0.0), m_me(0.0), m_mewgt(0.0), m_muf2(0.0), m_mur2(0.0)
    { m_last[0]=0.0; m_last[1]=0.0; }

    inline NLO_subevt(const NLO_subevt *nlos)
    { *this=*nlos; m_result=0.; m_last[0]=0.; m_last[1]=0.; }

    Particle_List *CreateParticleList() const;

    std::string IDString() const;

    inline void Mult(const double &scal)
    { m_result*=scal; m_last[0]*=scal; m_last[1]*=scal; m_me*=scal; m_mewgt*=scal; }
    inline void MultMEwgt(const double &scal) { m_mewgt*=scal; }

    inline bool IsReal() const { return m_i==m_j; }

    inline NLO_subevt &operator*=(const double &scal)
    { m_result*=scal; m_last[0]*=scal; m_last[1]*=scal; return *this; }
    inline NLO_subevt &operator+=(const double &sm)
    { m_result+=sm; m_last[0]+=sm; m_last[1]+=sm; return *this; }

    template <class Type> inline Type *Proc() const 
    { return static_cast<Type*>(p_proc); }

  };// end of class NLO_subevt

  std::ostream &operator<<(std::ostream &ostr,const NLO_subevt &sevt);

  class NLO_subevtlist: public std::vector<NLO_subevt*> { 
  public:

    void Mult(const double &scal);
    void MultMEwgt(const double &scal);

    NLO_subevtlist& operator*=(const double scal);

  };// end of class NLO_subevtlist

}// end of namespace ATOOLS

#endif
