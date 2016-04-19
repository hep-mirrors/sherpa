#ifndef ATOOLS_Phys_Cluster_Leg_H
#define ATOOLS_Phys_Cluster_Leg_H

#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/CXXFLAGS.H"

namespace ATOOLS {

  struct ColorID {
    int m_i, m_j;
    inline ColorID(const int &i=-1,const int &j=-1): m_i(i), m_j(j) {}
    inline ColorID Conj() const { return ColorID(m_j,m_i); }
    inline bool operator==(const ColorID &c) const 
    { return c.m_i==m_i && c.m_j==m_j; }
    inline bool Singlet() { return m_i==m_j; }
  };// end of class ColorID

  typedef std::vector<ColorID> ColorID_Vector;

  std::ostream &operator<<(std::ostream &ostr,const ColorID &col);

  class Cluster_Amplitude;
  class Cluster_Leg;

  class ClusterLeg_PVector: public std::vector<Cluster_Leg*> {
  public:

    ClusterLeg_PVector();

    ~ClusterLeg_PVector();

  };// end of class ClusterLeg_PVector

  class Cluster_Leg {
  private:

    Cluster_Amplitude *p_ampl;

    size_t  m_id, m_st, m_n, m_d, m_k;
    Vec4D   m_p;
    Flavour m_fl;
    ColorID m_c;

    double  m_kt2[2];

    static ClusterLeg_PVector s_legs;

  public:

    inline Cluster_Leg(Cluster_Amplitude *const ampl,
		       const Cluster_Leg &ref)
    { *this=ref; p_ampl=ampl; }
    inline Cluster_Leg(Cluster_Amplitude *const ampl,
		       const Vec4D &p,const Flavour &fl,
		       const ColorID &c=ColorID()):
      p_ampl(ampl), m_id(0), m_st(0), m_n(0), m_d(0), m_k(0),
      m_p(p), m_fl(fl), m_c(c)
    { m_kt2[0]=m_kt2[1]=-1.0; }

    inline ~Cluster_Leg() {}

    static Cluster_Leg *New(Cluster_Amplitude *const ampl,
			    const Cluster_Leg &ref);
    static Cluster_Leg *New(Cluster_Amplitude *const ampl,
			    const Vec4D &p,const Flavour &fl,
			    const ColorID &c=ColorID());
    void Delete();

    // inline functions
    inline void SetId(const size_t &id)   { m_id=id; }
    inline void SetStat(const size_t &st) { m_st=st; }
    inline void SetNMax(const size_t &n)  { m_n=n;   }
    inline void SetK(const size_t &k)     { m_k=k;   }
      
    inline void SetCol(const ColorID &c)   { m_c=c; }
    inline void SetMom(const Vec4D &p)     { m_p=p; }
    inline void SetFlav(const Flavour &fl) { m_fl=fl; }

    inline void SetKT2(const int i,const double &kt2) { m_kt2[i]=kt2; }

    inline Cluster_Amplitude *Amplitude() const { return p_ampl; }

    inline Vec4D   Mom() const  { return m_p;  }
    inline Flavour Flav() const { return m_fl; }
    inline ColorID Col() const  { return m_c;  }

    inline size_t Id() const   { return m_id; }
    inline size_t Stat() const { return m_st; }
    inline size_t NMax() const { return m_n;  }
    inline size_t K() const    { return m_k;  }

    inline double KT2(const int i) const { return m_kt2[i]; }

  };// end of class Cluster_Leg

  typedef std::vector<Cluster_Leg*> ClusterLeg_Vector;

  std::ostream &operator<<(std::ostream &ostr,const Cluster_Leg &leg);

}// end of namespace ATOOLS

#endif
