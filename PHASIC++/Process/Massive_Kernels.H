#ifndef PHASIC__Process__Massive_Kernels_H
#define PHASIC__Process__Massive_Kernels_H

#include <vector>
#include <cstddef>

namespace PHASIC {

  class Massive_Kernels {
  protected:
    size_t m_nf,m_nmf;
    double m_NC,m_CA,m_CF,m_TR,m_g1,m_g2,m_g3,m_K1,m_K2,m_K3, m_aterm;
    double m_alpha_ff,m_alpha_ii,m_alpha_if,m_alpha_fi,m_loga,m_kappa;
    double p_VS[3];
    double p_VNS[3];
    double p_Gamma[3];
    std::vector<double> m_massflav;

    double Lambda(double x,double y, double z);
    void CalcVS(double,double,double);
    void CalcVNS(int,double,double,double,bool);
    void CalcVNSq(double s,double mj,double mk);
    void CalcVNSg(double s,double mk,bool);
    void CalcGamma(int,double,double);
    void CalcAterms(int t,double mu2,double s,double mj,double mk);
  public:
    Massive_Kernels();
    
    void SetAlpha(double aff, double afi, double aif, double aii); 
    void SetNC(const double &n);

    void SetKappa(double k) { m_kappa=k; } 
    size_t Nf()       { return m_nf; }
    size_t Nmf()      { return m_nmf; }
    double FMass(int i) { return m_massflav[i]; }
    void Calculate(int,double,double,double,double,bool,bool,bool,bool);
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
