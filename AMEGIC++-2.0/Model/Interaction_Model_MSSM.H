#ifndef Interaction_Model_MSSM_h
#define Interaction_Model_MSSM_h

#include "Interaction_Model_THDM.H"
#include "Interaction_Model_Inos.H"
#include "Interaction_Model_sLepton_EW.H"
#include "Interaction_Model_sQCD.H"
#include "Interaction_Model_sQuark_EW.H"
#include "Interaction_Model_Base.H"
#include "Kabbala.H"

namespace AMEGIC {  
  class Interaction_Model_MSSM : public Interaction_Model_Base {
    Interaction_Model_THDM       * p_mothdm;
    Interaction_Model_Inos       * p_moinos;
    Interaction_Model_sLepton_EW * p_moslepton;
    Interaction_Model_sQCD       * p_mosqcd;
    Interaction_Model_sQuark_EW  * p_mosquark;

  public:
    Interaction_Model_MSSM(MODEL::Model_Base *,std::string,std::string);
    void c_FFV(Single_Vertex *,int &);
    void c_FFS(Single_Vertex *,int &);
    void c_VVV(Single_Vertex *,int &);
    void c_VVS(Single_Vertex *,int &);
    void c_SSV(Single_Vertex *,int &);
    void c_SSS(Single_Vertex *,int &);
    void c_VVVV(Single_Vertex *,int &);
    void c_SSVV(Single_Vertex *,int &);
    void c_SSSS(Single_Vertex *,int &);
    void c_FFT(Single_Vertex  *,int&)  {}
    void c_VVT(Single_Vertex  *,int&)  {} 
    void c_SST(Single_Vertex  *,int&)  {} 
    void c_VVVT(Single_Vertex *,int&)  {} 
    void c_FFVT(Single_Vertex *,int&)  {} 
    void c_SSST(Single_Vertex *,int&)  {}
  };
}

#endif


