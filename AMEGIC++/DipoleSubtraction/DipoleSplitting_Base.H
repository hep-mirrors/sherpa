#ifndef AMEGIC_DipoleSubtraction_DipoleSplitting_Base_H
#define AMEGIC_DipoleSubtraction_DipoleSplitting_Base_H

#include "ATOOLS/Math/Vector.H"
#include "AMEGIC++/Main/Pol_Info.H"
#include "MODEL/Main/Coupling_Data.H"
#include <vector>
#include <limits>

namespace MODEL {
  class Model_Base;
}

namespace AMEGIC {

#define nan numeric_limits<double>::quiet_NaN()

  class dpt {
  public:
    enum dipoltype{  
      none = 0,
      f_f  = 1,
      f_i  = 2,
      i_f  = 3,
      i_i  = 4,
      f_fm = 5,
      f_im = 6,
      i_fm = 7
    };
  };

  class spt {
  public:
    enum splittingtype{  
      none = 0,
      qg  = 1,
      gq  = 2,
      qq  = 3,
      gg  = 4,
      sqg = 5,
      gsq = 6,
      gtg = 7,
      ggt = 8
    };
  };

  class DipoleSplitting_Base {
  protected:
    double m_alpha,m_amin,m_kappa;
    int m_nf;
    ATOOLS::Vec4D m_pt1,m_pt2;

    dpt::dipoltype m_type;
    ATOOLS::Vec4D m_pi,m_pj,m_pk,m_ptij,m_ptk;
    int m_i,m_j,m_k,m_tij,m_tk,m_m;
    std::vector<ATOOLS::Vec4D> m_dpollist;
    std::vector<ATOOLS::Vec4D> m_mom;
    std::vector<double> m_pfactors;
//     spt::splittingtype m_ft;
    int m_ft;
    double m_spfdef;
    double m_g1,m_g2;
    MODEL::Coupling_Data *p_cpl;
    void CalcVectors(ATOOLS::Vec4D&,ATOOLS::Vec4D&,double);
  public:
    DipoleSplitting_Base();
    virtual ~DipoleSplitting_Base() {}
    
    virtual void SetMomenta(const ATOOLS::Vec4D*);
    void SetCoupling(const MODEL::Coupling_Map *cpls);
    std::vector<ATOOLS::Vec4D>* GetMomenta() {return &m_mom;}
    ATOOLS::Vec4D Getptij() {return m_ptij;}
    ATOOLS::Vec4D Getptk()  {return m_ptk;}
    virtual double GetF()   {return 0.;}
    virtual double GetDPSF()   {return 0.;}
    virtual void CalcDiPolarizations() {}
    virtual void SetAlpha(double alp) {m_alpha=alp;}
    double GetR(const ATOOLS::Vec4D*,const ATOOLS::Vec4D*);
    std::vector<ATOOLS::Vec4D>* GetDiPolarizations() {return &m_dpollist;}
    std::vector<double>* GetFactors() {return &m_pfactors;}
    dpt::dipoltype GetType() { return m_type;}
    double Vijf(int type,int mode);
    double Vif(int type,int mode);
    double Vie1(int type);
    double Vie2(int type);
    double SPFac() { return m_spfdef*p_cpl->Factor(); }
    double G2() { return m_g2; }
  };
}
#endif
