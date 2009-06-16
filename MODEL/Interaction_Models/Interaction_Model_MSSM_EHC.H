#ifndef MODEL_Interaction_Models_Interaction_Model_MSSM_EHC_H
#define MODEL_Interaction_Models_Interaction_Model_MSSM_EHC_H

#include "MODEL/Interaction_Models/Interaction_Model_THDM.H"
#include "MODEL/Interaction_Models/Interaction_Model_Inos.H"
#include "MODEL/Interaction_Models/Interaction_Model_sLepton_EW.H"
#include "MODEL/Interaction_Models/Interaction_Model_sQCD.H"
#include "MODEL/Interaction_Models/Interaction_Model_sQuark_EW.H"
#include "MODEL/Interaction_Models/Interaction_Model_sLepton_sQuark.H"
#include "MODEL/Interaction_Models/Interaction_Model_Base.H"
#include "ATOOLS/Math/Kabbala.H"

namespace MODEL {  
  class Interaction_Model_MSSM_EHC : public Interaction_Model_Base {
    Interaction_Model_THDM           * p_mothdm;
    Interaction_Model_Inos           * p_moinos;
    Interaction_Model_sLepton_EW     * p_moslepton;
    Interaction_Model_sQCD           * p_mosqcd;
    Interaction_Model_sQuark_EW      * p_mosquark;
    Interaction_Model_sLepton_sQuark * p_moslesqu;
    ATOOLS::Kabbala ghgg,g1,g2,g3,sintW,costW,PL,PR,M_I,root2,vev;
  public:
    Interaction_Model_MSSM_EHC(MODEL::Model_Base *,std::string,std::string);
    void c_FFV(std::vector<Single_Vertex>&,int &);
    void c_FFS(std::vector<Single_Vertex>&,int &);
    void c_VVV(std::vector<Single_Vertex>&,int &);
    void c_VVS(std::vector<Single_Vertex>&,int &);
    void c_SSV(std::vector<Single_Vertex>&,int &);
    void c_SSS(std::vector<Single_Vertex>&,int &);
    void c_VVVV(std::vector<Single_Vertex>&,int &);
    void c_SSVV(std::vector<Single_Vertex>&,int &);
    void c_SSSS(std::vector<Single_Vertex>&,int &);
    void c_FFT(std::vector<Single_Vertex>&,int&)  {}
    void c_VVT(std::vector<Single_Vertex>&,int&)  {} 
    void c_SST(std::vector<Single_Vertex>&,int&)  {} 
    void c_VVVT(std::vector<Single_Vertex>&,int&)  {} 
    void c_FFVT(std::vector<Single_Vertex>&,int&)  {} 
    void c_SSST(std::vector<Single_Vertex>&,int&)  {}
    ~Interaction_Model_MSSM_EHC();
  };
}

#endif


