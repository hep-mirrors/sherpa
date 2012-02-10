#ifndef ATOOLS_Phys_Cluster_Leg_H
#define ATOOLS_Phys_Cluster_Leg_H

#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/CXXFLAGS.H"

#ifdef USING__Threading
#include <pthread.h>
#endif

namespace ATOOLS {

  struct ColorID {
    int m_i, m_j;
    inline ColorID(const int &i=-1,const int &j=-1): m_i(i), m_j(j) {}
    inline ColorID Conj() { return ColorID(m_j,m_i); }
    inline bool operator==(const ColorID &c) const 
    { return c.m_i==m_i && c.m_j==m_j; }
    inline bool Singlet() { return m_i==m_j; }
  };// end of class ColorID

  std::ostream &operator<<(std::ostream &ostr,const ColorID &col);

  class Cluster_Amplitude;
  class Cluster_Leg;

  class ClusterLeg_PVector: public std::vector<Cluster_Leg*> {
  private:

#ifdef USING__Threading
    pthread_mutex_t m_mtx;
#endif

  public:

    ClusterLeg_PVector();

    ~ClusterLeg_PVector();

    inline void MtxLock()
    {
#ifdef USING__Threading
      pthread_mutex_lock(&m_mtx);
#endif
    }

    inline void MtxUnLock()
    {
#ifdef USING__Threading
      pthread_mutex_unlock(&m_mtx);
#endif
    }

  };// end of class ClusterLeg_PVector

  class Cluster_Leg {
  private:

    Cluster_Amplitude *p_ampl;

    size_t  m_id, m_st, m_n, m_d, m_k;
    Vec4D   m_p;
    Flavour m_fl;
    ColorID m_c;

    double  m_ktstart, m_ktmax, m_ktveto;

    bool    m_connected;
    std::set<Cluster_Leg *> m_spectators;

    static ClusterLeg_PVector s_legs;

  public:

    inline Cluster_Leg(Cluster_Amplitude *const ampl,
		       const Cluster_Leg &ref)
    { *this=ref; p_ampl=ampl; }
    inline Cluster_Leg(Cluster_Amplitude *const ampl,
		       const Vec4D &p,const Flavour &fl,
		       const ColorID &c=ColorID()):
      p_ampl(ampl), m_id(0), m_st(0), m_n(0), m_d(0), m_k(0),
      m_p(p), m_fl(fl), m_c(c),
      m_ktstart(0.), m_ktmax(0.), m_ktveto(0.), m_connected(true)
    {}

    inline ~Cluster_Leg() {}

    static Cluster_Leg *New(Cluster_Amplitude *const ampl,
			    const Cluster_Leg &ref);
    static Cluster_Leg *New(Cluster_Amplitude *const ampl,
			    const Vec4D &p,const Flavour &fl,
			    const ColorID &c=ColorID());
    inline std::set<Cluster_Leg *> & GetSpectators() { return m_spectators; }
    inline void AddToSpectators(Cluster_Leg * cl) { m_spectators.insert(cl); }
    inline int NumberOfSpectators() const { return m_spectators.size(); }

    void Delete();

    // inline functions
    inline void SetId(const size_t &id)   { m_id=id; }
    inline void SetStat(const size_t &st) { m_st=st; }
    inline void SetNMax(const size_t &n)  { m_n=n;   }
    inline void SetDMax(const size_t &d)  { m_d=d;   }
    inline void SetK(const size_t &k)     { m_k=k;   }
      
    inline void SetCol(const ColorID &c)   { m_c=c; }
    inline void SetMom(const Vec4D &p)     { m_p=p; }
    inline void SetFlav(const Flavour &fl) { m_fl=fl; }

    inline void SetKTStart(const double & kt) { m_ktstart = kt; }
    inline void SetKTMax(const double & kt)   { m_ktmax   = kt; }
    inline void SetKTVeto(const double & kt)  { m_ktveto  = kt; }

    inline void SetConnected(const bool & c)  { m_connected = c; }

    inline Cluster_Amplitude *Amplitude() const { return p_ampl; }

    inline Vec4D   Mom() const  { return m_p;  }
    inline Flavour Flav() const { return m_fl; }
    inline ColorID Col() const  { return m_c;  }

    inline size_t Id() const   { return m_id; }
    inline size_t Stat() const { return m_st; }
    inline size_t NMax() const { return m_n;  }
    inline size_t DMax() const { return m_d;  }
    inline size_t K() const    { return m_k;  }

    inline const double KTStart() const { return m_ktstart; }
    inline const double KTMax() const   { return m_ktmax; }
    inline const double KTVeto() const  { return m_ktveto; }

    inline const bool Connected() const { return m_connected; }
  };// end of class Cluster_Leg

  typedef std::vector<Cluster_Leg*> ClusterLeg_Vector;

  std::ostream &operator<<(std::ostream &ostr,const Cluster_Leg &leg);

}// end of namespace ATOOLS

#endif
