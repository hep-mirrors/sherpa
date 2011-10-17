#ifndef PHASIC__Process__Flavour_Kernels_H
#define PHASIC__Process__Flavour_Kernels_H

namespace PHASIC {

  class Flavour_Kernels {
  protected:

    int m_nf;

    double m_NC, m_CA, m_CF, m_TR, m_g1, m_g2;
    double m_alpha, m_loga;

  public:

    Flavour_Kernels();

    void SetAlpha(const double &a); 
    void SetNC(const double &n);

    double Kb1(int type,double x);
    double Kb2(int type);
    double Kb3(int type,double x);
    double Kb4(int type,double x);
    double t1(double x);
    double t2();
    double t4(double x);
    double Kt1(int type,double x);
    double Kt2(int type);
    double Kt3(int type,double x);
    double Kt4(int type,double x);
    double P1(int type,double x);
    double P2(int type);
    double P3(int type,double x);
    double P4(int type,double x);
    double ft(int type);
    double Nf()       { return m_nf; }
    double hfnc1(double x);

    /*
      function key:
      Kb<x> = {\bar K}
      t<x>  = (1/(1-x)_+ + delta(1-x))
      Kt<x> = {\tilde K}
      P<x>  = P(x)/T^2

      <x>:  (distribution type) 
       1    g(x)
       2    h
       3    k(x)
       4    G(eta)

      int type: (flavour combination)
       1     q q
       2     q g
       3     g qb
       4     g g   
    */
  };

};

#endif
