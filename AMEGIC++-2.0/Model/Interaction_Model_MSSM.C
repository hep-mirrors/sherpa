#include "Interaction_Model_MSSM.H"
#include "MathTools.H"
#include "Message.H"
#include <stdio.h>


using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;

Interaction_Model_MSSM::Interaction_Model_MSSM(MODEL::Model_Base * _model,
					       std::string _cplscheme,std::string _yukscheme) :
  Interaction_Model_Base(_model,_cplscheme,_yukscheme)
{ 
  p_mothdm    = new Interaction_Model_THDM(p_model,_cplscheme,_yukscheme); 
  p_moinos    = new Interaction_Model_Inos(p_model,_cplscheme,_yukscheme); 
  p_moslepton = new Interaction_Model_sLepton_EW(p_model,_cplscheme,_yukscheme); 
  p_mosqcd    = new Interaction_Model_sQCD(p_model,_cplscheme,_yukscheme); 
  p_mosquark  = new Interaction_Model_sQuark_EW(p_model,_cplscheme,_yukscheme); 
  p_moslesqu  = new Interaction_Model_sLepton_sQuark(p_model,_cplscheme,_yukscheme); 
}

void Interaction_Model_MSSM::c_FFV(Single_Vertex* vertex,int& vanz)
{
  p_mothdm->c_FFV(vertex,vanz);
  p_moinos->c_FFV(vertex,vanz);
  p_moslepton->c_FFV(vertex,vanz);
  p_mosqcd->c_FFV(vertex,vanz);
  p_mosquark->c_FFV(vertex,vanz);
}

void Interaction_Model_MSSM::c_VVV(Single_Vertex* vertex,int& vanz)
{
  p_mothdm->c_VVV(vertex,vanz);
  p_moinos->c_VVV(vertex,vanz);
  p_moslepton->c_VVV(vertex,vanz);
  p_mosqcd->c_VVV(vertex,vanz);
  p_mosquark->c_VVV(vertex,vanz);
}

void Interaction_Model_MSSM::c_VVVV(Single_Vertex* vertex,int& vanz)
{
  p_mothdm->c_VVVV(vertex,vanz);
  p_moinos->c_VVVV(vertex,vanz);
  p_moslepton->c_VVVV(vertex,vanz);
  p_mosqcd->c_VVVV(vertex,vanz);
  p_mosquark->c_VVVV(vertex,vanz);
}

void Interaction_Model_MSSM::c_FFS(Single_Vertex* vertex,int& vanz)  
{ 
  p_mothdm->c_FFS(vertex,vanz);
  p_moinos->c_FFS(vertex,vanz);
  p_moslepton->c_FFS(vertex,vanz);
  p_mosqcd->c_FFS(vertex,vanz);
  p_mosquark->c_FFS(vertex,vanz);
}

void Interaction_Model_MSSM::c_VVS(Single_Vertex* vertex,int& vanz)  
{ 
  p_mothdm->c_VVS(vertex,vanz);
  p_moinos->c_VVS(vertex,vanz);
  p_moslepton->c_VVS(vertex,vanz);
  p_mosqcd->c_VVS(vertex,vanz);
  p_mosquark->c_VVS(vertex,vanz);
}

void Interaction_Model_MSSM::c_SSS(Single_Vertex* vertex,int& vanz)  
{ 
  p_mothdm->c_SSS(vertex,vanz);
  p_moinos->c_SSS(vertex,vanz);
  p_moslepton->c_SSS(vertex,vanz);
  p_mosqcd->c_SSS(vertex,vanz);
  p_mosquark->c_SSS(vertex,vanz);
}

void Interaction_Model_MSSM::c_SSV(Single_Vertex* vertex,int& vanz) 
{ 
  p_mothdm->c_SSV(vertex,vanz);
  p_moinos->c_SSV(vertex,vanz);
  p_moslepton->c_SSV(vertex,vanz);
  p_mosqcd->c_SSV(vertex,vanz);
  p_mosquark->c_SSV(vertex,vanz);
}

void Interaction_Model_MSSM::c_SSVV(Single_Vertex* vertex,int& vanz) 
{ 
  p_mothdm->c_SSVV(vertex,vanz);
  p_moinos->c_SSVV(vertex,vanz);
  p_moslepton->c_SSVV(vertex,vanz);
  p_mosqcd->c_SSVV(vertex,vanz);
  p_mosquark->c_SSVV(vertex,vanz);
}

void Interaction_Model_MSSM::c_SSSS(Single_Vertex* vertex,int& vanz) 
{ 
  p_mothdm->c_SSSS(vertex,vanz);
  p_moinos->c_SSSS(vertex,vanz);
  p_moslepton->c_SSSS(vertex,vanz);
  p_mosqcd->c_SSSS(vertex,vanz);
  p_mosquark->c_SSSS(vertex,vanz);
  p_moslesqu->c_SSSS(vertex,vanz);
}

Interaction_Model_MSSM::~Interaction_Model_MSSM()
{
  delete p_mothdm;
  delete p_moinos;
  delete p_moslepton;
  delete p_mosqcd;
  delete p_mosquark;
}
