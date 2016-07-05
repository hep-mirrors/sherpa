#ifndef PHASIC__Process__KP_Terms_H
#define PHASIC__Process__KP_Terms_H

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/Flavour_Kernels.H"
#include "PHASIC++/Process/Massive_Kernels.H"
#include "MODEL/Main/Coupling_Data.H"

namespace PDF { class PDF_Base; }

namespace PHASIC {

  struct xp_coeff {
    double xp;
    double kpc;
  };

  class KP_Terms {
  private:
    int m_imode;

    Process_Base    *p_proc;
    Flavour_Kernels *p_flkern;
    Massive_Kernels *p_masskern;

    MODEL::Coupling_Data *p_cpl;

    ATOOLS::Flavour_Vector m_flavs;

    bool   m_massive, m_cemode, m_negativepdf;
    double m_a, m_loga, m_kpca[8], m_kpcb[8], m_cpldef;

    std::vector<std::vector<double> > m_dsij;

    std::vector<size_t> m_plist;

    std::vector<xp_coeff> m_xpa, m_xpb;

  public:

    KP_Terms(Process_Base *const proc,const int mode);

    ~KP_Terms();

    void SetCoupling(const MODEL::Coupling_Map *cpls);
    
    void Calculate(const ATOOLS::Vec4D_Vector &mom,
		   const double &x0,const double &x1,
		   const double &eta0,const double &eta1,
		   const double &weight);

    double Get(const double &x0,const double &x1,
               const double &eta0,const double &eta1,
               const ATOOLS::Flavour_Vector &flav,
               const int mode,
               const double &scalefac2=1.0);

    void SetNC(const double &nc);
    void SetAlpha(const double &aff,const double &afi=-1.0,
		  const double &aif=-1.0,const double &aii=-1.0);

    void FillMEwgts(ATOOLS::ME_Weight_Info &wgt);

    void SetDSij(const std::vector<std::vector<double> > &ds);

    inline double Coupling() const { return m_cpldef*p_cpl->Factor(); }

    inline Flavour_Kernels *FlavKern() const { return p_flkern;   }
    inline Massive_Kernels *MassKern() const { return p_masskern; }

    inline void SetIMode(int imode) { m_imode=imode; }
  };// end of class KP_Terms

};// end of namespace PHASIC

#endif
