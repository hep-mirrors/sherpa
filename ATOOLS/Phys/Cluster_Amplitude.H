#ifndef ATOOLS_Phys_Cluster_Amplitude_H
#define ATOOLS_Phys_Cluster_Amplitude_H

#include "ATOOLS/Phys/Cluster_Leg.H"
#include "ATOOLS/Phys/Decay_Info.H"

#include <vector>
#include <map>

namespace ATOOLS {

  typedef std::map<int,int> CI_Map;

  class ClusterAmplitude_PVector: 
    public std::vector<Cluster_Amplitude*> {
  private:

#ifdef USING__Threading
    pthread_mutex_t m_mtx;
#endif

  public:

    ClusterAmplitude_PVector();

    ~ClusterAmplitude_PVector();

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

  };// end of class ClusterAmplitude_PVector

  class Cluster_Amplitude {
  private:

    Cluster_Amplitude *p_prev, *p_next;

    ClusterLeg_Vector m_legs;
    DecayInfo_Vector  m_decs;

    size_t m_oew, m_oqcd, m_nin, m_new, m_ncl, m_kin, m_nlo;
    double m_mur2, m_muf2, m_Q2, m_mu2, m_kt2, m_z, m_phi, m_lkf;
    CI_Map m_cmap;

    void *p_jf, *p_proc, *p_procs, *p_iinfo, *p_dinfo;

    Mass_Selector *p_ms;

    static ClusterAmplitude_PVector s_ampls;

    friend class ClusterAmplitude_PVector;

  protected:

    Cluster_Amplitude(Cluster_Amplitude *const prev=NULL);

    ~Cluster_Amplitude();

  public:

    static Cluster_Amplitude* New(Cluster_Amplitude *const prev=NULL);

    void Delete();

    void CopyFrom(const Cluster_Amplitude *const master,
		  const int mode=0);

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
    Cluster_Amplitude *InitPrev();
    void SetNext(Cluster_Amplitude *const next);
    void UnsetNext();

    void DeletePrev();
    void DeleteNext();

    void IdSort();

    size_t NQCD() const;
    size_t NEW() const;

    Cluster_Leg *IdLeg(const size_t &id) const;
    size_t       IdIndex(const size_t &id) const;

    void OrderLegs();

    Cluster_Leg *Splitter() const;

    // inline functions
    inline const ClusterLeg_Vector &Legs() const { return m_legs; }
    inline ClusterLeg_Vector       &Legs()       { return m_legs; }

    inline Cluster_Leg *Leg(const size_t &i) const { return m_legs[i]; }

    inline Cluster_Amplitude *Prev() const { return p_prev; }
    inline Cluster_Amplitude *Next() const { return p_next; }

    inline const DecayInfo_Vector &Decays() const { return m_decs; }
    inline DecayInfo_Vector       &Decays()       { return m_decs; }

    inline void SetNIn(const size_t &nin)  { m_nin=nin; }
    inline void SetIdNew(const size_t &id) { m_new=id;  }
    inline void SetNewCol(const size_t &c) { m_ncl=c;   }

    inline void SetKin(const size_t &kin)   { m_kin=kin; }
    inline void SetNLO(const size_t &nlo)   { m_nlo=nlo; }

    inline void SetOrderEW(const size_t &oew)   { m_oew=oew;   }
    inline void SetOrderQCD(const size_t &oqcd) { m_oqcd=oqcd; }

    inline void SetMuF2(const double &muf2) { m_muf2=muf2; }
    inline void SetMuR2(const double &mur2) { m_mur2=mur2; }
    inline void SetMu2(const double &mu2)   { m_mu2=mu2;   }

    inline void SetQ2(const double &Q2) { m_Q2=Q2; }

    inline void SetKT2(const double &kt2) { m_kt2=kt2; }
    inline void SetZ(const double &z)     { m_z=z;     }
    inline void SetPhi(const double &phi) { m_phi=phi; }

    inline void SetJF(void *const jf)   { p_jf=jf;   }
    inline void SetProc(void *const p)  { p_proc=p;  }
    inline void SetProcs(void *const p) { p_procs=p; }
    inline void SetIInfo(void *const i) { p_iinfo=i; }
    inline void SetDInfo(void *const d) { p_dinfo=d; }

    inline void SetMS(Mass_Selector *const ms) { p_ms=ms; }

    inline size_t NIn() const    { return m_nin; }
    inline size_t IdNew() const  { return m_new; }
    inline size_t NewCol() const { return m_ncl; }

    inline size_t Kin() const  { return m_kin; }
    inline size_t NLO() const  { return m_nlo; }

    inline size_t OrderEW() const  { return m_oew;  }
    inline size_t OrderQCD() const { return m_oqcd; }

    inline double MuF2() const { return m_muf2; }
    inline double MuR2() const { return m_mur2; }
    inline double Mu2() const  { return m_mu2;  }

    inline double Q2() const { return m_Q2; }

    inline double KT2() const { return m_kt2; }
    inline double Z() const   { return m_z;   }
    inline double Phi() const { return m_phi; }

    inline void SetLKF(const double &lkf) { m_lkf=lkf; }

    inline double LKF() const { return m_lkf; }

    inline CI_Map       &ColorMap()       { return m_cmap; }
    inline const CI_Map &ColorMap() const { return m_cmap; }

    inline Mass_Selector *MS() const { return p_ms; }

    template <class Type> inline Type *JF() const 
    { return static_cast<Type*>(p_jf); }
    
    static void SetColours(Cluster_Leg *const lij,
                           Cluster_Leg *const li,
                           Cluster_Leg *const lj);
    template <class Type> inline Type *Proc() const 
    { return static_cast<Type*>(p_proc); }
    template <class Type> inline Type *Procs() const 
    { return static_cast<Type*>(p_procs); }
    template <class Type> inline Type *IInfo() const 
    { return static_cast<Type*>(p_iinfo); }
    template <class Type> inline Type *DInfo() const 
    { return static_cast<Type*>(p_dinfo); }

  };// end of class Cluster_Amplitude

  std::ostream &operator<<
    (std::ostream &ostr,const Cluster_Amplitude &ampl);

  typedef std::vector<Cluster_Amplitude*> ClusterAmplitude_Vector;

}// end of namespace ATOOLS

#endif
