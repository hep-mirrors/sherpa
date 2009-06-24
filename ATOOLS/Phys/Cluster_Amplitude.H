#ifndef ATOOLS_Phys_Cluster_Amplitude_H
#define ATOOLS_Phys_Cluster_Amplitude_H

#include "ATOOLS/Phys/Cluster_Leg.H"

#include <vector>
#include <map>

namespace ATOOLS {

  typedef std::map<int,int> CI_Map;

  class ClusterAmplitude_PVector: 
    public std::vector<Cluster_Amplitude*> {
  public:

    ~ClusterAmplitude_PVector();

  };// end of class ClusterAmplitude_PVector

  class Cluster_Amplitude {
  private:

    Cluster_Amplitude *p_prev, *p_next;

    ClusterLeg_Vector m_legs;

    size_t m_oew, m_oqcd, m_swap, m_nin, m_new;
    double m_mur2, m_muf2, m_kt2qcd, m_x1, m_x2, m_rbmax;
    CI_Map m_cmap;

    void *p_jf;

    Mass_Selector *p_ms;

    static ClusterAmplitude_PVector s_ampls;

    friend class ClusterAmplitude_PVector;

  protected:

    Cluster_Amplitude(Cluster_Amplitude *const prev=NULL);

    ~Cluster_Amplitude();

  public:

    static Cluster_Amplitude* New(Cluster_Amplitude *const prev=NULL);

    void Delete();

    void CopyFrom(const Cluster_Amplitude *const master);

    Cluster_Amplitude *Copy() const;
    Cluster_Amplitude *CopyNext() const;
    Cluster_Amplitude *CopyAll() const;

    void CreateLeg(const Vec4D &p,const Flavour &fl,
		   const ColorID &col=ColorID(),
		   const size_t &id=std::string::npos);
    void CombineLegs(Cluster_Leg *const i,Cluster_Leg *const j,
		     const Flavour &fl,
		     const ColorID &col=ColorID());

    Cluster_Amplitude *InitNext();
    void SetNext(Cluster_Amplitude *const next);

    void DeletePrev();
    void DeleteNext();
    void Print() const;

    size_t NQCD() const;
    size_t NEW() const;

    void SwapInOrder();

    Cluster_Leg *IdLeg(const size_t &id) const;
    size_t       IdIndex(const size_t &id) const;

    // inline functions
    inline const ClusterLeg_Vector &Legs() const { return m_legs; }
    inline ClusterLeg_Vector       &Legs()       { return m_legs; }

    inline Cluster_Leg *Leg(const size_t &i) const { return m_legs[i]; }

    inline Cluster_Amplitude *Prev() const { return p_prev; }
    inline Cluster_Amplitude *Next() const { return p_next; }

    inline void SetNIn(const size_t &nin)  { m_nin=nin; }
    inline void SetIdNew(const size_t &id) { m_new=id;  }

    inline void SetOrderEW(const size_t &oew)   { m_oew=oew;   }
    inline void SetOrderQCD(const size_t &oqcd) { m_oqcd=oqcd; }

    inline void SetMuF2(const double &muf2) { m_muf2=muf2; }
    inline void SetMuR2(const double &mur2) { m_mur2=mur2; }

    inline void SetX1(const double &x1) { m_x1=x1; }
    inline void SetX2(const double &x2) { m_x2=x2; }

    inline void SetKT2QCD(const double &kt2) { m_kt2qcd=kt2; }
    inline void SetRBMax(const double &rb)   { m_rbmax=rb;   }

    inline void SetJF(void *const jf) { p_jf=jf; }

    inline void SetMS(Mass_Selector *const ms) { p_ms=ms; }

    inline size_t NIn() const   { return m_nin; }
    inline size_t IdNew() const { return m_new; }

    inline size_t OrderEW() const  { return m_oew;  }
    inline size_t OrderQCD() const { return m_oqcd; }

    inline double MuF2() const { return m_muf2; }
    inline double MuR2() const { return m_mur2; }

    inline double X1() const { return m_x1; }
    inline double X2() const { return m_x2; }

    inline double KT2QCD() const { return m_kt2qcd; }
    inline double RBMax() const  { return m_rbmax;  }

    inline bool InSwaped() const { return m_swap; }

    inline CI_Map       &ColorMap()       { return m_cmap; }
    inline const CI_Map &ColorMap() const { return m_cmap; }

    inline Mass_Selector *MS() const { return p_ms; }

    template <class Type> inline Type *JF() const 
    { return static_cast<Type*>(p_jf); }

    inline void SwapLegs(const size_t &i,const size_t &j)
    { std::swap<Cluster_Leg*>(m_legs[i],m_legs[j]); }

  };// end of class Cluster_Amplitude

  std::ostream &operator<<
    (std::ostream &ostr,const Cluster_Amplitude &ampl);

}// end of namespace ATOOLS

#endif
