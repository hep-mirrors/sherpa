#ifndef COMIX_Phasespace_PS_Vertex_H
#define COMIX_Phasespace_PS_Vertex_H

#include "METOOLS/Explicit/Vertex.H"
#include "COMIX/Phasespace/PS_Info.H"

using namespace METOOLS;

namespace ATOOLS { class NLO_subevt; }

namespace COMIX {

  typedef std::vector<PS_Info*> PS_Info_Vector;

  class PS_Vertex: public Vertex {
  protected:

    double m_alpha, m_oldalpha, m_weight;
    double m_np, m_sum, m_sum2;
    double m_mnp, m_msum, m_msum2;

    int m_type;

    ATOOLS::NLO_subevt *p_dip;

  public:

    // constructor
    PS_Vertex(const Vertex_Key &key);
    
    // member functions
    void Evaluate();

    void AddPoint(const double &weight);
    void Reset();

    void MPISync();

    // inline functions
    inline double Alpha() const    { return m_alpha;    }
    inline double OldAlpha() const { return m_oldalpha; }

    inline double Weight() const { return m_weight; }

    inline void SetAlpha(const double &alpha)    { m_alpha=alpha;    }
    inline void SetOldAlpha(const double &alpha) { m_oldalpha=alpha; }

    inline void SetWeight(const double &weight) { m_weight=weight; }

    inline double N() const    { return m_np;   }
    inline double Sum() const  { return m_sum;  }
    inline double Sum2() const { return m_sum2; }

    inline void SetN(const double &n)       { m_np=n;      }
    inline void SetSum(const double &sum)   { m_sum=sum;   }
    inline void SetSum2(const double &sum2) { m_sum2=sum2; }

    inline double Mean() const { return m_sum/m_np; }
    inline double Variance() const    
    { return (m_sum2-m_sum*m_sum/m_np)/(m_np-1.0); }
    inline double Sigma() const  
    { return sqrt(Variance()/m_np); }

    inline void SetType(const int &type) { m_type=type; }

    inline int Type() const { return m_type; }

    inline void GetMPIVars(double *v)
    { v[0]=m_mnp; v[1]=m_msum; v[2]=m_msum2; }
    inline void AddMPIVars(const double *v)
    { m_mnp+=v[0]; m_msum+=v[1]; m_msum2+=v[2]; }
    inline void SetMPIVars(const double *v)
    { m_mnp=v[0]; m_msum=v[1]; m_msum2=v[2]; }

    inline void SetDip(ATOOLS::NLO_subevt *const d) { p_dip=d; }

    inline ATOOLS::NLO_subevt *Dip() const { return p_dip; }

  };// end of class PS_Vertex

}// end of namespace COMIX

#endif
