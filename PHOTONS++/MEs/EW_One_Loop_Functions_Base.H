#ifndef PHOTONS_MEs_EW_One_Loop_Functions_Base_H
#define PHOTONS_MEs_EW_One_Loop_Functions_Base_H

#include "ATOOLS/Math/MyComplex.H"
#include "METOOLS/Loops/Divergence_Array.H"
#include "METOOLS/Loops/Master_Integrals.H"

namespace PHOTONS{
  class EW_One_Loop_Functions_Base {
  protected:
    std::string        m_name;
    int                m_ferm, m_nonferm, m_ew, m_OL_comp;
    // EW parameters
    double             m_alpha, m_e;
    double             m_g;
    Complex            m_sw, m_cw, m_sw2, m_cw2, m_sw4, m_cw4;
    double             MW2;
    double             MZ2;
    double             MH2;
    double             mt2;
    Complex            muW2;
    Complex            muZ2;
    Complex            muH2;
    Complex            mut2;
    double             RW;
    double             RZ;
    double             wh;
    double             zh;

    double             m_s;
    double             m_mu2;

    // Generic constants
    METOOLS::DivArrC   One, Zero;

    // Helper functions
    METOOLS::DivArrC B_f(const double& p2, const Complex& m12, const Complex& m22);
    METOOLS::DivArrC B_fp(const double& p2, const Complex& m12, const Complex& m22);
    METOOLS::DivArrC B_ff(const double& p2, const Complex& m12, const Complex& m22);
    METOOLS::DivArrC wWU();
    METOOLS::DivArrC wWF();
    double Kallen(const double& x, const double& y, const double& z);

    // Self energies
    METOOLS::DivArrC Sigma_ZA(const double& p2);
    METOOLS::DivArrC Sigma_ZZ(const double& p2);
    METOOLS::DivArrC Sigma_W(const double& p2);
    METOOLS::DivArrC Sigma_H(const double& p2);
    // Derivatives of self energies
    METOOLS::DivArrC dSigma_ZZ(const double& p2);
    METOOLS::DivArrC dSigma_W(const double& p2);
    METOOLS::DivArrC dSigma_H(const double& p2);
    METOOLS::DivArrC dSigma_GamGam(const double& p2);
    METOOLS::DivArrC dSigma_GamGam0();

    // Fermion self energies
    METOOLS::DivArrC Sigma_ferm_R(const double& p2, const double& m2, const double& m2p,
			       const double& Qf, const double& If);
    METOOLS::DivArrC Sigma_ferm_L(const double& p2, const double& m2, const double& m2p,
			       const double& Qf, const double& If);
    METOOLS::DivArrC Sigma_ferm_S(const double& p2, const double& m2, const double& m2p,
			       const double& Qf, const double& If);
    // Derivatives of fermion self energies
    METOOLS::DivArrC dSigma_ferm_R(const double& p2, const double& m2, const double& m2p,
				const double& Qf, const double& If);
    METOOLS::DivArrC dSigma_ferm_L(const double& p2, const double& m2, const double& m2p,
				const double& Qf, const double& If);
    METOOLS::DivArrC dSigma_ferm_S(const double& p2, const double& m2, const double& m2p,
				const double& Qf, const double& If);

  public:
    EW_One_Loop_Functions_Base(const double& s, const double& mu2, const int& ew = 1);
    virtual ~EW_One_Loop_Functions_Base();

    void Print_Ren_Constants_Finite();
    void Print_Ren_Constants_UV();
    void Print_Ren_Constants_IR();
    // Counterterms
    METOOLS::DivArrC dZW();
    METOOLS::DivArrC dZH();
    METOOLS::DivArrC dZZZ();
    METOOLS::DivArrC dZAA();
    METOOLS::DivArrC dZZA();
    METOOLS::DivArrC dZAZ();
    METOOLS::DivArrC dMW2();
    METOOLS::DivArrC dMZ2();
    METOOLS::DivArrC dMH2();
    METOOLS::DivArrC dZfermL(const double& m2, const double& m2p,
			     const double& Qf, const double& If);
    METOOLS::DivArrC dZfermR(const double& m2, const double& m2p,
			     const double& Qf, const double& If);
    METOOLS::DivArrC dm(const double& m2, const double& m2p,
			const double& Qf, const double& If);


    METOOLS::DivArrC dZe();
    METOOLS::DivArrC dcw();

    inline double Get_sw() { return real(m_sw); }
    inline double Get_cw() { return real(m_cw); }
    inline double Get_e() { return m_e; }
    inline double Get_MW2() { return MW2; }
  };
}

#endif
