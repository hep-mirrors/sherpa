#ifndef PHASIC__Process__Massive_Kernels_H
#define PHASIC__Process__Massive_Kernels_H

#include "ATOOLS/Phys/NLO_Types.H"
#include <vector>
#include <cstddef>

namespace PHASIC {

  class Flavour_Kernels {
  protected:
  public:
    Flavour_Kernels();
    ~Flavour_Kernels();
  };

  class I_Kernels : public Flavour_Kernels {
  protected:
  public:
    I_Kernels();
    ~I_Kernels();
  };

  class KP_Kernels : public Flavour_Kernels {
  protected:
  public:
    KP_Kernels();
    ~KP_Kernels();
  };

  class Massive_Kernels {
  protected:
    const ATOOLS::sbt::subtype m_stype;
    size_t m_nf,m_nmf,m_subtype;
    double m_NC,m_CA,m_CF,m_TR,m_TRbyCA,m_CFbyCA,m_TRbyCF;
    double m_g1t,m_g2t,m_g3t,m_K1t,m_K2t,m_K3t;
    double m_beta0qcd, m_beta0qed;
    double m_alpha_ff,m_alpha_fi,m_alpha_if,m_alpha_ii,m_kappa;
    double m_logaff,m_logafi,m_logaif,m_logaii;
    double m_VNS,m_gKterm,m_aterm;
    double p_VS[3],p_Gammat[2];
    std::vector<double> m_massflav, m_cpls;
    int    m_Vsubmode,m_collVFF;

    inline double Lambda(double x, double y, double z)
    { return x*x+y*y+z*z-2.*(x*y+x*z+y*z); }

    void CalcVS(ATOOLS::ist::itype type, double s, double mj, double mk);
    void CalcVNS(ATOOLS::ist::itype type, double s, double mj, double mk,
                 bool ini);
    void CalcVNSq(double s, double mj, double mk);
    void CalcVNSg(double s, double mk, bool ini);
    void CalcVNSs(double s, double mj, double mk);
    void CalcGamma(ATOOLS::ist::itype type, double mu2, double s, double m);
    void CalcgKterm(ATOOLS::ist::itype type, double mu2, double s, double mj,
                    bool mode);
    void CalcAterms(ATOOLS::ist::itype type, double mu2, double s,
                    double mj,double mk,bool inij,bool inik);
    void CalcAq(double mu2, double s,double mj,double mk);
    void CalcAg(double mu2, double s,double mk);
    void CalcAs(double mu2, double s,double mj,double mk);
  public:
    Massive_Kernels(ATOOLS::sbt::subtype st,
                    const size_t &nf, const size_t nmf);
    ~Massive_Kernels() {}

    void SetAlpha(double aff, double afi, double aif, double aii); 
    void SetNC(const double &n);
    void SetSubType(const int type);

    inline void   SetKappa(double k)                { m_kappa=k; }

    inline void   SetNf(size_t nf)                  { m_nf=nf; }
    inline void   SetNmf(size_t nmf)                { m_nmf=nmf; }
    inline size_t Nf()                              { return m_nf; }
    inline size_t Nmf()                             { return m_nmf; }
    inline double FMass(size_t i)                   { return m_massflav[i]; }
    inline void   SetVSubtractionMode(int vsm)      { m_Vsubmode=vsm; }
    inline void   SetCollinearVFFSplitting(int csm) { m_collVFF=csm; }

    void Calculate(ATOOLS::ist::itype t, double mu2, double s,
                   double mj, double mk,
                   bool inij, bool inik, bool mode);

    double I_Fin();
    double I_E1();
    double I_E2();

    // Kbar
    double Kb1(int type,double x);
    double Kb2(int type);
    double Kb3(int type,double x);
    double Kb4(int type,double x);

    // K_FS
    double KFS1(int type,double x);
    double KFS2(int type);
    double KFS3(int type,double x);
    double KFS4(int type,double x);

    // Kcal
    double Kc1(int type,int typej,double muja2,double x);
    double Kc2(int type,int typej,double muja2);
    double Kc3(int type,int typej,double muja2,double x);
    double Kc4(int type,int typej,double muja2,double x);

    // Kbar_M
    double KbM1(int type,double muak2,double x);
    double KbM2(int type,double muak2);
    double KbM3(int type,double muak2,double x);
    double KbM4(int type,double muak2,double x);

    // Ktilde
    double Kt1(int type,double x);
    double Kt2(int type);
    double Kt3(int type,double x);
    double Kt4(int type,double x);

    // P
    double P1(int type,double x);
    double P2(int type);
    double P3(int type,double x);
    double P4(int type,double x);

    double t1(int type,int spin,double muq,double x);
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

    inline void   SetCpls(std::vector<double> cpls) { m_cpls=cpls; }
    inline double Cpl(size_t i)                     { return m_cpls[i]; }

    inline double NC() const                        { return m_NC; }
    inline double CF() const                        { return m_CF; }
    inline double CA() const                        { return m_CA; }
    inline double TR() const                        { return m_TR; }
    inline double Beta0QCD() const                  { return m_beta0qcd; }
    inline double Beta0QED() const                  { return m_beta0qed; }

    /*
      function key: (old)
      Kb<x> = {\bar K}
      t<x>  = (1/(1-x)_+ + delta(1-x))
      Kt<x> = {\tilde K}
      P<x>  = P(x)/T^2

      <x>:  (distribution type)
       1    g(x)
       2    h
       3    k(x)
       4    G(eta)

      int type: (flavour combination a a')
       1     q q
       2     g q
       3     q g
       4     g g
    */
  };
}
#endif
