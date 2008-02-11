#include "Interaction_Model_MSSM_EHC.H"
#include "MathTools.H"
#include "Message.H"
#include <stdio.h>


using namespace MODEL;
using namespace ATOOLS;
using namespace std;

Interaction_Model_MSSM_EHC::Interaction_Model_MSSM_EHC(MODEL::Model_Base * _model,
					       std::string _cplscheme,std::string _yukscheme) :
  Interaction_Model_Base(_model,_cplscheme,_yukscheme)
{ 
  p_moehc     = new Interaction_Model_EHC_S(p_model,_cplscheme,_yukscheme); 
  p_mothdm    = new Interaction_Model_THDM(p_model,_cplscheme,_yukscheme); 
  p_moinos    = new Interaction_Model_Inos(p_model,_cplscheme,_yukscheme); 
  p_moslepton = new Interaction_Model_sLepton_EW(p_model,_cplscheme,_yukscheme); 
  p_mosqcd    = new Interaction_Model_sQCD(p_model,_cplscheme,_yukscheme); 
  p_mosquark  = new Interaction_Model_sQuark_EW(p_model,_cplscheme,_yukscheme); 
  p_moslesqu  = new Interaction_Model_sLepton_sQuark(p_model,_cplscheme,_yukscheme); 
}

void Interaction_Model_MSSM_EHC::c_FFV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  p_mothdm->c_FFV(vertex,vanz);
  p_moinos->c_FFV(vertex,vanz);
  p_moslepton->c_FFV(vertex,vanz);
  p_mosqcd->c_FFV(vertex,vanz);
  p_mosquark->c_FFV(vertex,vanz);
}

void Interaction_Model_MSSM_EHC::c_VVV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  p_mothdm->c_VVV(vertex,vanz);
  p_moinos->c_VVV(vertex,vanz);
  p_moslepton->c_VVV(vertex,vanz);
  p_mosqcd->c_VVV(vertex,vanz);
  p_mosquark->c_VVV(vertex,vanz);
}

void Interaction_Model_MSSM_EHC::c_VVVV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  p_moehc->c_VVVV(vertex,vanz);
  p_mothdm->c_VVVV(vertex,vanz);
  p_moinos->c_VVVV(vertex,vanz);
  p_moslepton->c_VVVV(vertex,vanz);
  p_mosqcd->c_VVVV(vertex,vanz);
  p_mosquark->c_VVVV(vertex,vanz);
}

void Interaction_Model_MSSM_EHC::c_FFS(std::vector<Single_Vertex>& vertex,int& vanz)  
{ 
  p_mothdm->c_FFS(vertex,vanz);
  p_moinos->c_FFS(vertex,vanz);
  p_moslepton->c_FFS(vertex,vanz);
  p_mosqcd->c_FFS(vertex,vanz);
  p_mosquark->c_FFS(vertex,vanz);
}

void Interaction_Model_MSSM_EHC::c_VVS(std::vector<Single_Vertex>& vertex,int& vanz)  
{ 
  p_moehc->c_VVS(vertex,vanz);
  p_mothdm->c_VVS(vertex,vanz);
  p_moinos->c_VVS(vertex,vanz);
  p_moslepton->c_VVS(vertex,vanz);
  p_mosqcd->c_VVS(vertex,vanz);
  p_mosquark->c_VVS(vertex,vanz);
}

void Interaction_Model_MSSM_EHC::c_SSS(std::vector<Single_Vertex>& vertex,int& vanz)  
{ 
  p_mothdm->c_SSS(vertex,vanz);
  p_moinos->c_SSS(vertex,vanz);
  p_moslepton->c_SSS(vertex,vanz);
  p_mosqcd->c_SSS(vertex,vanz);
  p_mosquark->c_SSS(vertex,vanz);
}

void Interaction_Model_MSSM_EHC::c_SSV(std::vector<Single_Vertex>& vertex,int& vanz) 
{ 
  p_mothdm->c_SSV(vertex,vanz);
  p_moinos->c_SSV(vertex,vanz);
  p_moslepton->c_SSV(vertex,vanz);
  p_mosqcd->c_SSV(vertex,vanz);
  p_mosquark->c_SSV(vertex,vanz);
}

void Interaction_Model_MSSM_EHC::c_SSVV(std::vector<Single_Vertex>& vertex,int& vanz) 
{ 
  p_mothdm->c_SSVV(vertex,vanz);
  p_moinos->c_SSVV(vertex,vanz);
  p_moslepton->c_SSVV(vertex,vanz);
  p_mosqcd->c_SSVV(vertex,vanz);
  p_mosquark->c_SSVV(vertex,vanz);
}

void Interaction_Model_MSSM_EHC::c_SSSS(std::vector<Single_Vertex>& vertex,int& vanz) 
{ 
  p_mothdm->c_SSSS(vertex,vanz);
  p_moinos->c_SSSS(vertex,vanz);
  p_moslepton->c_SSSS(vertex,vanz);
  p_mosqcd->c_SSSS(vertex,vanz);
  p_mosquark->c_SSSS(vertex,vanz);
  p_moslesqu->c_SSSS(vertex,vanz);
}

Interaction_Model_MSSM_EHC::~Interaction_Model_MSSM_EHC()
{
  delete p_moehc;
  delete p_mothdm;
  delete p_moinos;
  delete p_moslepton;
  delete p_mosqcd;
  delete p_mosquark;
  delete p_moslesqu;
}
