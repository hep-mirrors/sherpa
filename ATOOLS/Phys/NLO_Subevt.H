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
    inline bool operator==(const IDip_ID &di) const
    { return m_ijt==di.m_ijt && m_kt==di.m_kt; }
  };// end of struct IDip_ID

  std::ostream &operator<<(std::ostream &ostr,const IDip_ID &idi);

  typedef std::set<IDip_ID> IDip_Set;

  typedef std::map<std::string,IDip_Set> StringIDipSet_Map;

  struct DDip_ID {
    size_t m_i, m_j, m_k;
  public:
    inline DDip_ID(const size_t &i=0,
		   const size_t &j=0,const size_t &k=0):
      m_i(i), m_j(j), m_k(k) {}
    bool operator<(const DDip_ID &di) const;
    inline bool operator==(const DDip_ID &di) const
    { return m_i==di.m_i && m_j==di.m_j && m_k==di.m_k; }
  };// end of struct DDip_ID

  std::ostream &operator<<(std::ostream &ostr,const DDip_ID &ddi);

  typedef std::set<DDip_ID> DDip_Set;

  typedef std::map<std::string,DDip_Set> StringDDipSet_Map;

  struct Dip_ID: public IDip_ID, public DDip_ID {
    inline Dip_ID(const size_t &ijt=0,const size_t &kt=0,
		  const size_t &i=0,const size_t &j=0,const size_t &k=0):
      IDip_ID(ijt,kt), DDip_ID(i,j,k) {}
    bool operator<(const Dip_ID &di) const;
  };// end of struct Dip_ID

  std::ostream &operator<<(std::ostream &ostr,const Dip_ID &di);

  struct stp {
    enum id {
      fac     =  0,
      ren     =  1,
      size    =  2
    };
  };// end of struct stp

  class NLO_subevt: public Dip_ID {
  public:

    const Flavour *p_fl;
    const Vec4D   *p_mom;

    const size_t     *p_id;
    DecayInfo_Vector *p_dec;

    NLO_subevt *p_real;
    void       *p_proc;

    size_t m_n, m_flip;
    double m_result, m_last[2], m_me, m_mewgt;

    std::vector<double> m_mu2;

    std::string m_pname;

  public:

    inline NLO_subevt(const size_t &n=0,const size_t *id=NULL,
		      const Flavour *fl=NULL,const Vec4D *mom=NULL,
		      const int i=-1,const int j=-1,const int k=-1):
      Dip_ID(0,0,i,j,k), p_fl(fl), p_mom(mom), p_id(id), p_dec(NULL),
      p_real(NULL), p_proc(NULL), m_n(n), m_flip(0),
      m_result(0.0), m_me(0.0), m_mewgt(0.0), m_mu2(2,0.0)
    { m_last[0]=0.0; m_last[1]=0.0; }

    void CopyXSData(const NLO_subevt *sub);

    Particle_List *CreateParticleList() const;

    std::string IDString(const int mode=0) const;

    inline void Mult(const double &scal)
    { m_result*=scal; m_last[0]*=scal; m_last[1]*=scal;
      m_me*=scal; m_mewgt*=scal; }
    inline void MultME(const double &scal) { m_me*=scal; m_mewgt*=scal; }
    inline void MultMEwgt(const double &scal) { m_mewgt*=scal; }

    inline bool IsReal() const { return m_i==m_j; }

    inline NLO_subevt &operator*=(const double &scal)
    { m_result*=scal; m_last[0]*=scal; m_last[1]*=scal; return *this; }
    inline NLO_subevt &operator+=(const double &sm)
    { m_result+=sm; m_last[0]+=sm; m_last[1]+=sm; return *this; }

    template <class Type> inline Type *Proc() const 
    { return static_cast<Type*>(p_proc); }

    inline void Reset()
    { m_result=m_last[0]=m_last[1]=m_me=m_mewgt=0.0; }

  };// end of class NLO_subevt

  std::ostream &operator<<(std::ostream &ostr,const NLO_subevt &sevt);

  class NLO_subevtlist: public std::vector<NLO_subevt*> { 
  public:

    void Mult(const double &scal);
    void MultME(const double &scal);
    void MultMEwgt(const double &scal);

    NLO_subevtlist& operator*=(const double scal);

  };// end of class NLO_subevtlist

}// end of namespace ATOOLS

#endif
