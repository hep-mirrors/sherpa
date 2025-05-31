#ifndef MALARIC__Tools__Amplitude_H
#define MALARIC__Tools__Amplitude_H

#include "MALARIC/Tools/Parton.H"
#include "MALARIC/Tools/Splitting.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"

namespace MALARIC {

  class Amplitude: public Parton_Vector {
  private:

    Splitting m_s;
    int       m_mec;
    double    m_t, m_t0;

    ATOOLS::Cluster_Amplitude *p_ampl;

  public:

    Amplitude(ATOOLS::Cluster_Amplitude *const a);

    ~Amplitude();

    void Add(Parton *const p);
    void Remove(Parton *const p);

    ATOOLS::Cluster_Amplitude *GetAmplitude() const;

    inline ATOOLS::Cluster_Amplitude *
    ClusterAmplitude() const { return p_ampl; }

    inline void SetSplit(const Splitting &s) { m_s=s; }

    inline const Splitting &Split() const { return m_s; }

    inline void SetT(const double &t)  { m_t=t;  }
    inline void SetT0(const double &t) { m_t0=t; }

    inline double T() const  { return m_t;  }
    inline double T0() const { return m_t0; }

    inline void SetMEC(const int mec) { m_mec=mec;    }
    inline int  MEC() const           { return m_mec; }

    inline void SetJF(void *const jf) { p_ampl->SetJF(jf); }
    template <class Type> inline Type *JF() const 
    { return p_ampl->JF<Type>(); }

  };// end of class Amplitude

  std::ostream &operator<<(std::ostream &s,const Amplitude &a);

  typedef std::vector<Amplitude*> Amplitude_Vector;

}// end of namespace MALARIC

#endif
