#ifndef PHOTONS_MEs_RVTools_Z_Decay_Two_Loop_Functions_H
#define PHOTONS_MEs_RVTools_Z_Decay_Two_Loop_Functions_H

#include "ATOOLS/Math/MyComplex.H"
#include "METOOLS/Loops/Divergence_Array.H"
#include "METOOLS/Loops/Master_Integrals.H"
#include "PHOTONS++/MEs/EW_One_Loop_Functions_Base.H"
#include "PHOTONS++/MEs/RVTools/Z_Decay_RV_Diagrams.H"
#include "ATOOLS/Phys/Particle.H"

namespace PHOTONS{
  class Z_Decay_Two_Loop_Functions {
  private:
    double m_mu2, m_e, m_If, m_Qf, m_sw, m_cw, m_sw2, m_cw2, m_p12, m_p22, m_pP2; 
    double m_m,m_m2,m_m3,m_m4,m_m6,m_m8,m_m10,m_m12,m_m14,m_m16,m_m18,m_m20,m_m22;
    double m_s,m_s_2,m_s_3,m_s_4,m_s_5,m_s_6,m_s_7,m_s_8,m_s_9;
    double m_s12,m_s12_2,m_s12_3,m_s12_4,m_s12_5,m_s12_6,m_s12_7,m_s12_8,m_s12_9;
    double m_s1k,m_s1k_2,m_s1k_3,m_s1k_4,m_s1k_5,m_s1k_6,m_s1k_7,m_s1k_8,m_s1k_9;
    double m_s2k,m_s2k_2,m_s2k_3,m_s2k_4,m_s2k_5,m_s2k_6,m_s2k_7,m_s2k_8,m_s2k_9;
    Complex m_cL, m_cR;
    ATOOLS::Vec4D m_p1, m_p2, m_pP;
    int spins[2]; // spin of fermion and antifermion
    EW_One_Loop_Functions_Base* p_EW;
    METOOLS::DivArrC Zero,One;//, D;

    Z_Decay_RV_Bubbles * p_bubbles;
    Z_Decay_RV_Vertices * p_vertices;
    Z_Decay_RV_Box_1 * p_box_1;
    Z_Decay_RV_Box_2 * p_box_2;
    Z_Decay_RV_CT * p_ct;

    METOOLS::DivArrC p_bubblecoeff_1[8][2];
    METOOLS::DivArrC p_bubblecoeff_2[8][2];
    METOOLS::DivArrC p_Z_vertexcoeff_1[8][2];
    METOOLS::DivArrC p_Z_vertexcoeff_2[8][2];
    METOOLS::DivArrC p_P_vertexcoeff_1[8][2];
    METOOLS::DivArrC p_P_vertexcoeff_2[8][2];
    METOOLS::DivArrC p_boxcoeff_1[8][2];
    METOOLS::DivArrC p_boxcoeff_2[8][2];
    METOOLS::DivArrC p_fermionct_coeff_1[8][2];
    METOOLS::DivArrC p_fermionct_coeff_2[8][2];
    METOOLS::DivArrC p_vertexct_coeff_1[8][2];
    METOOLS::DivArrC p_vertexct_coeff_2[8][2];
    METOOLS::DivArrC p_B_coeff_1[8][2];
    METOOLS::DivArrC p_B_coeff_2[8][2];

  public:
    Z_Decay_Two_Loop_Functions(const double& m, const double& s,
      const Vec4D& p1, const Vec4D& p2, const Vec4D& pP,
      const Complex& cL, const Complex& cR, const double& mu2);
    virtual ~Z_Decay_Two_Loop_Functions();    

    METOOLS::DivArrC Propagator_CT(const int& anti, const Vec4D& prop);
    METOOLS::DivArrC CT_L();
    METOOLS::DivArrC CT_R();
    METOOLS::DivArrC CT_QED_L();
    METOOLS::DivArrC CT_QED_R();

    void Calculate_RV_Coeffs(const ATOOLS::Vec4C& epsV, const ATOOLS::Vec4C& epsP);
    // Retrieve tabulated coefficients
    METOOLS::DivArrC Get_RV_Bubble_Insertion_1(const int& ME, const int& LR) {return p_bubblecoeff_1[ME][LR];}
    METOOLS::DivArrC Get_RV_Bubble_Insertion_2(const int& ME, const int& LR) {return p_bubblecoeff_2[ME][LR];}
    METOOLS::DivArrC Get_RV_Z_Vertex_1(const int& ME, const int& LR) {return p_Z_vertexcoeff_1[ME][LR];}
    METOOLS::DivArrC Get_RV_Z_Vertex_2(const int& ME, const int& LR) {return p_Z_vertexcoeff_2[ME][LR];}
    METOOLS::DivArrC Get_RV_P_Vertex_1(const int& ME, const int& LR) {return p_P_vertexcoeff_1[ME][LR];}
    METOOLS::DivArrC Get_RV_P_Vertex_2(const int& ME, const int& LR) {return p_P_vertexcoeff_2[ME][LR];}
    METOOLS::DivArrC Get_RV_Box_1(const int& ME, const int& LR) {return p_boxcoeff_1[ME][LR];}
    METOOLS::DivArrC Get_RV_Box_2(const int& ME, const int& LR) {return p_boxcoeff_2[ME][LR];}
    METOOLS::DivArrC Get_RV_Fermion_CT_1(const int& ME, const int& LR) {return p_fermionct_coeff_1[ME][LR];}
    METOOLS::DivArrC Get_RV_Fermion_CT_2(const int& ME, const int& LR) {return p_fermionct_coeff_2[ME][LR];}
    METOOLS::DivArrC Get_RV_Vertex_CT_1(const int& ME, const int& LR) {return p_vertexct_coeff_1[ME][LR];}
    METOOLS::DivArrC Get_RV_Vertex_CT_2(const int& ME, const int& LR) {return p_vertexct_coeff_2[ME][LR];}
    METOOLS::DivArrC Get_RV_B_1(const int& ME, const int& LR) {return p_B_coeff_1[ME][LR];}
    METOOLS::DivArrC Get_RV_B_2(const int& ME, const int& LR) {return p_B_coeff_2[ME][LR];}      
  };
}

#endif
