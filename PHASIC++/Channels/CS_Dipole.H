#ifndef PHASIC_Main_CS_Dipole_H
#define PHASIC_Main_CS_Dipole_H

#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Phys/NLO_Subevt.H"
#include "ATOOLS/Org/Info_Key.H"

namespace PDF { class ISR_Handler; }

namespace PHASIC {

  class Vegas;
  class Cut_Data;
  class Multi_Channel;
  class Phase_Space_Handler;

  class CS_Dipole {
  protected:

    ATOOLS::NLO_subevt m_sub;

    Vegas *p_vegas;
    double m_rn[3];

    Multi_Channel    *p_fsmc, *p_ismc;
    ATOOLS::Info_Key  m_isrspkey, m_isrykey;
    PDF::ISR_Handler *p_isr;

    double m_alpha, m_oldalpha, m_weight, m_rbweight;
    double m_np, m_sum, m_sum2;
    double m_mnp, m_msum, m_msum2;
    double m_amin, m_q2min;

    std::string m_id;

    std::map<size_t,size_t> m_brmap, m_rbmap;

    size_t m_type, m_ijt, m_kt, m_idx;
    bool   m_on, m_bmcw;

    ATOOLS::Flavour m_fli, m_flj, m_flk, m_flij;

    double Lambda(const double &s,const double &sb,
		  const double &sc) const;

    double GetS(const double &Q2,const double &y,
		const double &mi2,const double &mj2,
		const double &mk2) const;
    double GetZ(const double &Q2,const double &sij,
		const double &y,const double &zt,
		const double &mi2,const double &mk2) const;
    double GetKT2(const double &Q2,const double &y,const double &z,
		  const double &mi2,const double &mj2,
		  const double &mk2) const;
    double ConstructLN(const double &Q2,const double &sij,
		       const double &mij2,const double &mk2,
		       const ATOOLS::Vec4D &Q,ATOOLS::Vec4D &pk,
		       ATOOLS::Vec4D &l,ATOOLS::Vec4D &n) const;

  public:

    CS_Dipole(ATOOLS::NLO_subevt *const sub,
	      Phase_Space_Handler *const psh,const bool bmcw=0);

    virtual ~CS_Dipole();

    virtual ATOOLS::Vec4D_Vector GeneratePoint
    (const ATOOLS::Vec4D_Vector &p,
     Cut_Data *const cuts,const double *rns) = 0;
    virtual double GenerateWeight
    (const ATOOLS::Vec4D_Vector &p,Cut_Data *const cuts) = 0;

    virtual bool ValidPoint(const ATOOLS::Vec4D_Vector& p) = 0;

    void InitVegas(const std::string &pid);

    double Phi(ATOOLS::Vec4D pijt,ATOOLS::Vec4D pkt,
               ATOOLS::Vec4D pi,const bool ii) const;

    double Alpha(const int mode) const;

    void AddPoint(const double &value,const double &ewgt,const int mode);
    void Optimize();
    void EndOptimize();
    void MPICollect(std::vector<double> &sv,size_t &i);
    void MPIReturn(std::vector<double> &sv,size_t &i);
    void MPISync();
    void Reset();

    bool IsMapped(CS_Dipole *const dip) const;

    void WriteOut(const std::string &pid,
		  std::vector<std::string> &info) const;
    void ReadIn(const std::string &pid,
		const std::vector<std::string> &info);

    // inline functions
    inline const ATOOLS::NLO_subevt *GetSubEvt() const { return &m_sub; }

    inline void SetAMin(const double &amin)   { m_amin=amin;   }
    inline void SetQ2Min(const double &q2min) { m_q2min=q2min; }

    inline double Alpha() const    { return m_alpha;    }
    inline double OldAlpha() const { return m_oldalpha; }

    inline double Weight() const   { return m_weight;   }
    inline double RBWeight() const { return m_rbweight; }

    inline void SetAlpha(const double &alpha)    { m_alpha=alpha;    }
    inline void SetOldAlpha(const double &alpha) { m_oldalpha=alpha; }

    inline void SetIdx(const size_t &idx) { m_idx=idx; }

    inline double N() const { return m_np; }

    inline double Mean() const { return m_sum/m_np; }
    inline double Variance() const    
    { return (m_sum2-m_sum*m_sum/m_np)/(m_np-1.0); }
    inline double Sigma() const  
    { return sqrt(Variance()/m_np); }

    inline const std::string &Id() const { return m_id; }

    inline size_t Type() const { return m_type; }
    inline size_t Idx() const  { return m_idx;  }

    inline void SetOn(const bool on) { m_on=on; }
    inline bool On() const { return m_on; }

  };// end of struct Dipole

  std::ostream &operator<<(std::ostream &ostr,const CS_Dipole &dip);

  typedef std::vector<CS_Dipole*> CSDipole_Vector;
  typedef std::vector<CSDipole_Vector> CSDipole_Matrix;

}// end of namespace PHASIC

#endif
