#ifndef PHOTONS_MEs_RVTools_Z_Decay_RV_Diagrams_H
#define PHOTONS_MEs_RVTools_Z_Decay_RV_Diagrams_H

#include "ATOOLS/Math/MyComplex.H"
#include "METOOLS/Loops/Divergence_Array.H"
#include "METOOLS/Loops/Master_Integrals.H"
#include "PHOTONS++/MEs/EW_One_Loop_Functions_Base.H"
#include "ATOOLS/Phys/Particle.H"

namespace PHOTONS{
  // Declare all RV diagram classes
  class Z_Decay_RV_Diagrams {
  protected:
    double m_mu2, m_p12, m_p22, m_pP2; 
    double m_m,m_m2,m_m3,m_m4;
    double m_s,m_s_2,m_s_3;
    double m_s12;
    double m_s1k;
    double m_s2k;
    double m_x,m_x_2,m_x_3,m_x_4,m_x_5,m_x_6,m_x_7,m_x_8,m_x_9,m_x_10,m_x_11;
    double m_z12,m_z12_2,m_z12_3,m_z12_4,m_z12_5,m_z12_6,m_z12_7,m_z12_8,m_z12_9;
    double m_z1k,m_z1k_2,m_z1k_3,m_z1k_4,m_z1k_5,m_z1k_6,m_z1k_7,m_z1k_8,m_z1k_9;
    double m_z2k,m_z2k_2,m_z2k_3,m_z2k_4,m_z2k_5,m_z2k_6,m_z2k_7,m_z2k_8,m_z2k_9;
    Complex m_cL, m_cR;
    ATOOLS::Vec4D m_p1, m_p2, m_pP;
    METOOLS::DivArrC Zero,One;//, D;

    void Init_Powers();
    virtual void Init_Coefficients() = 0;

  public:
    Z_Decay_RV_Diagrams(const double& m, const double& s,
			const Vec4D& p1, const Vec4D& p2, const Vec4D& pP,
			const Complex& cL, const Complex& cR, const double& mu2);
    virtual ~Z_Decay_RV_Diagrams();    
  };

  class Z_Decay_RV_Bubbles : public Z_Decay_RV_Diagrams {
  private:
    Complex * Bubble1Abb;
    Complex * Bubble2Abb;

    void Init_Bubble_1_Coefficients();
    void Init_Bubble_2_Coefficients();
    void Init_Coefficients();

  public:
    Z_Decay_RV_Bubbles(const double& m, const double& s,
		       const Vec4D& p1, const Vec4D& p2, const Vec4D& pP,
		       const Complex& cL, const Complex& cR, const double& mu2);
    virtual ~Z_Decay_RV_Bubbles();
    METOOLS::DivArrC RV_Bubble_Insertion_1(const int& ME, const int& LR,
					   const ATOOLS::Vec4C& epsV, const ATOOLS::Vec4C& epsP);
    METOOLS::DivArrC RV_Bubble_Insertion_2(const int& ME, const int& LR,
					   const ATOOLS::Vec4C& epsV, const ATOOLS::Vec4C& epsP);
  };

  class Z_Decay_RV_Vertices : public Z_Decay_RV_Diagrams {
  private:
    Complex * ZVertex1Abb;
    Complex * ZVertex2Abb;
    Complex * PVertex1Abb;
    Complex * PVertex2Abb;

    void Init_Z_Vertex_1_Coefficients();
    void Init_Z_Vertex_2_Coefficients();
    void Init_P_Vertex_1_Coefficients();
    void Init_P_Vertex_2_Coefficients();
    void Init_Coefficients();

  public:
    Z_Decay_RV_Vertices(const double& m, const double& s,
			const Vec4D& p1, const Vec4D& p2, const Vec4D& pP,
			const Complex& cL, const Complex& cR, const double& mu2);
    virtual ~Z_Decay_RV_Vertices();
    METOOLS::DivArrC RV_Z_Vertex_1(const int& ME, const int& LR,
				   const ATOOLS::Vec4C& epsV, const ATOOLS::Vec4C& epsP);
    METOOLS::DivArrC RV_Z_Vertex_2(const int& ME, const int& LR,
				   const ATOOLS::Vec4C& epsV, const ATOOLS::Vec4C& epsP);
    METOOLS::DivArrC RV_P_Vertex_1(const int& ME, const int& LR,
				   const ATOOLS::Vec4C& epsV, const ATOOLS::Vec4C& epsP);
    METOOLS::DivArrC RV_P_Vertex_2(const int& ME, const int& LR,
				   const ATOOLS::Vec4C& epsV, const ATOOLS::Vec4C& epsP);
  };

  class Z_Decay_RV_Box_1 : public Z_Decay_RV_Diagrams {
  private:
    Complex * Box1Abb;

    void Init_Box_1_Coefficients();
    void Init_Coefficients();

  public:
    Z_Decay_RV_Box_1(const double& m, const double& s,
		     const Vec4D& p1, const Vec4D& p2, const Vec4D& pP,
		     const Complex& cL, const Complex& cR, const double& mu2);
    virtual ~Z_Decay_RV_Box_1();
    METOOLS::DivArrC RV_Box_1(const int& ME, const int& LR,
			      const ATOOLS::Vec4C& epsV, const ATOOLS::Vec4C& epsP);
  };

  class Z_Decay_RV_Box_2 : public Z_Decay_RV_Diagrams {
  private:
    Complex * Box2Abb;

    void Init_Box_2_Coefficients();
    void Init_Coefficients();

  public:
    Z_Decay_RV_Box_2(const double& m, const double& s,
		     const Vec4D& p1, const Vec4D& p2, const Vec4D& pP,
		     const Complex& cL, const Complex& cR, const double& mu2);
    virtual ~Z_Decay_RV_Box_2();
    METOOLS::DivArrC RV_Box_2(const int& ME, const int& LR,
			      const ATOOLS::Vec4C& epsV, const ATOOLS::Vec4C& epsP);
  };
  class Z_Decay_RV_CT : public Z_Decay_RV_Diagrams {
  private:
    
    Complex COne;
    METOOLS::DivArrC B;
    EW_One_Loop_Functions_Base * p_EW;
    METOOLS::DivArrC Z_Vertex_CT_L();
    METOOLS::DivArrC Z_Vertex_CT_R();
    void Init_Coefficients();

  public:
    Z_Decay_RV_CT(const double& m, const double& s,
		     const Vec4D& p1, const Vec4D& p2, const Vec4D& pP,
		     const Complex& cL, const Complex& cR, const double& mu2);
    virtual ~Z_Decay_RV_CT();
    METOOLS::DivArrC RV_Fermion_CT_1(const int& ME, const int& LR, const ATOOLS::Vec4C& epsP);
    METOOLS::DivArrC RV_Vertex_CT_1(const int& ME, const int& LR, const ATOOLS::Vec4C& epsP);
    METOOLS::DivArrC RV_Fermion_CT_2(const int& ME, const int& LR, const ATOOLS::Vec4C& epsP);
    METOOLS::DivArrC RV_Vertex_CT_2(const int& ME, const int& LR, const ATOOLS::Vec4C& epsP);
    METOOLS::DivArrC RV_B_1(const int& ME, const int& LR, const ATOOLS::Vec4C& epsP);
    METOOLS::DivArrC RV_B_2(const int& ME, const int& LR, const ATOOLS::Vec4C& epsP);
  };
}

#endif
