#ifndef PHASIC_Main_CS_Dipole_H
#define PHASIC_Main_CS_Dipole_H

#include "ATOOLS/Phys/Cluster_Amplitude.H"

namespace PHASIC {

  class Vegas;
  class Process_Base;

  class CS_Dipole {
  protected:

    ATOOLS::Cluster_Leg m_lij, m_lk, m_li, m_lj;

    Vegas *p_vegas;
    double m_rn[3];

    double m_alpha, m_oldalpha, m_weight, m_np, m_sum, m_sum2, m_max;
    double m_amin, m_amax;

    std::string m_id;

    double Lambda(const double &s,const double &sb,
		  const double &sc) const;

  public:

    CS_Dipole(const ATOOLS::Cluster_Leg &lij,
	      const ATOOLS::Cluster_Leg &lk,
	      const ATOOLS::Cluster_Leg &li,
	      const ATOOLS::Cluster_Leg &lj);
    virtual ~CS_Dipole();

    virtual bool GeneratePoint
    (ATOOLS::Cluster_Amplitude *const ampl,const double *rns) = 0;
    virtual double GenerateWeight
    (const ATOOLS::Cluster_Amplitude *ampl) = 0;

    virtual CS_Dipole *Copy() const = 0;

    void Map(ATOOLS::Cluster_Amplitude *const ampl);
    void ReMap(ATOOLS::Cluster_Amplitude *const ampl);

    void InitVegas(Process_Base *const proc);

    double Phi(ATOOLS::Vec4D pijt,ATOOLS::Vec4D pkt,
               ATOOLS::Vec4D pi,const bool ii) const;

    void AddPoint(const double &weight);
    void Reset();

    // inline functions
    inline const ATOOLS::Cluster_Leg &LI()  const { return m_li;  }
    inline const ATOOLS::Cluster_Leg &LJ()  const { return m_lj;  }
    inline const ATOOLS::Cluster_Leg &LIJ() const { return m_lij; }
    inline const ATOOLS::Cluster_Leg &LK()  const { return m_lk;  }

    inline Vegas *GetVegas() const { return p_vegas; }

    inline void SetAMin(const double &amin) { m_amin=amin; }
    inline void SetAMax(const double &amax) { m_amax=amax; }

    inline double Alpha() const    { return m_alpha;    }
    inline double OldAlpha() const { return m_oldalpha; }

    inline double Weight() const { return m_weight; }

    inline void SetAlpha(const double &alpha)    { m_alpha=alpha;    }
    inline void SetOldAlpha(const double &alpha) { m_oldalpha=alpha; }

    inline void SetWeight(const double &weight) { m_weight=weight; }

    inline double N() const    { return m_np;   }
    inline double Sum() const  { return m_sum;  }
    inline double Sum2() const { return m_sum2; }
    inline double Max() const  { return m_max;  }

    inline void SetN(const double &n)       { m_np=n;      }
    inline void SetSum(const double &sum)   { m_sum=sum;   }
    inline void SetSum2(const double &sum2) { m_sum2=sum2; }
    inline void SetMax(const double &max)   { m_max=max;   }

    inline double Mean() const { return m_sum/m_np; }
    inline double Variance() const    
    { return (m_sum2-m_sum*m_sum/m_np)/(m_np-1.0); }
    inline double Sigma() const  
    { return sqrt(Variance()/m_np); }

    inline const std::string &Id() const { return m_id; }

  };// end of struct Dipole

  std::ostream &operator<<(std::ostream &ostr,const CS_Dipole &dip);

  typedef std::vector<CS_Dipole*> CSDipole_Vector;
  typedef std::vector<CSDipole_Vector> CSDipole_Matrix;

}// end of namespace PHASIC

#endif
