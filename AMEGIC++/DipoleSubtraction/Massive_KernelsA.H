#ifndef AMEGIC_DipoleSubtraction_Massive_KernelsA_H
#define AMEGIC_DipoleSubtraction_Massive_KernelsA_H

#include "MODEL/Main/Coupling_Data.H"
#include <vector>

namespace MODEL {
  class Model_Base;
}

namespace AMEGIC {

  class Massive_KernelsA {
  protected:
    size_t m_nf,m_nmf;
    double m_cpldef,m_g1,m_g2,m_g3,m_K1,m_K2,m_K3, m_aterm;
    double m_alpha,m_loga,m_kappa;
    double p_VS[3];
    double p_VNS[3];
    double p_Gamma[3];
    std::vector<double> m_massflav;
    MODEL::Coupling_Data *p_cpl;

    double Lambda(double x,double y, double z);
    void CalcVS(double,double,double);
    void CalcVNS(int,double,double,double,bool);
    void CalcVNSq(double s,double mj,double mk);
    void CalcVNSg(double s,double mk,bool);
    void CalcGamma(int,double,double);
    void CalcAterms(int t,double mu2,double s,double mj,double mk);
  public:
    Massive_KernelsA();
    
    void SetCoupling(const MODEL::Coupling_Map *cpls);
    void SetAlpha(double a); 
    void SetKappa(double k) { m_kappa=k; } 
    double Coupling() { return m_cpldef*p_cpl->Factor(); }
    size_t Nf()       { return m_nf; }
    size_t Nmf()      { return m_nmf; }
    double FMass(int i) { return m_massflav[i]; }
    void Calculate(int,double,double,double,double,bool,bool,bool);
    double I_Fin();
    double I_E1();
    double I_E2();
    double Kt1(int type,double x);
    double Kt2(int type);
    double Kt3(int type,double x);
    double Kt4(int type,double x);

    double t1(int type,int spin,double muq,double x);
    double t1p(int type,int spin,double muq,double x);
    double t2(int type,int spin,double muq);
    double t3(int type,int spin,double muq,double x);
    double t4(int type,int spin,double muq,double x);    
    double t5(int type,double x,double xp);
    double t6(int type,double xp);
    double t7(int type,double x,double xp);    
    double at1(int type,int spin,double muq,double x);
    double at2(int type,int spin,double muq);
    double at3(int type,int spin,double muq,double x);
    double at4(int type,int spin,double muq,double x);
  };
}
#endif
