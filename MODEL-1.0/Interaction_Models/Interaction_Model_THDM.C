#include "Interaction_Model_THDM.H"
#include "MathTools.H"
#include "Message.H"
#include <stdio.h>


using namespace MODEL;
using namespace ATOOLS;
using namespace std;

Interaction_Model_THDM::Interaction_Model_THDM(MODEL::Model_Base * _model,
					       std::string _cplscheme,std::string _yukscheme) :
  Interaction_Model_Base(_model,_cplscheme,_yukscheme)
{ 
  p_moew    = new Interaction_Model_EW(p_model,_cplscheme,_yukscheme); 
  p_moqcd   = new Interaction_Model_QCD(p_model,_cplscheme,_yukscheme); 
  p_mohiggs = new Interaction_Model_Higgs_THDM(p_model,_cplscheme,_yukscheme); 
}

void Interaction_Model_THDM::c_FFV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  p_moew->c_FFV(vertex,vanz);
  p_moqcd->c_FFV(vertex,vanz);
  p_mohiggs->c_FFV(vertex,vanz);
}

void Interaction_Model_THDM::c_VVV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  p_moew->c_VVV(vertex,vanz);
  p_moqcd->c_VVV(vertex,vanz);
  p_mohiggs->c_VVV(vertex,vanz);
}
void Interaction_Model_THDM::c_VVVV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  p_moew->c_VVVV(vertex,vanz);
  p_moqcd->c_VVVV(vertex,vanz);
  p_mohiggs->c_VVVV(vertex,vanz);
}

void Interaction_Model_THDM::c_FFS(std::vector<Single_Vertex>& vertex,int& vanz)  
{ 
  p_mohiggs->c_FFS(vertex,vanz); 
}

void Interaction_Model_THDM::c_VVS(std::vector<Single_Vertex>& vertex,int& vanz)  
{ 
  p_mohiggs->c_VVS(vertex,vanz); 
}

void Interaction_Model_THDM::c_SSS(std::vector<Single_Vertex>& vertex,int& vanz)  
{ 
  p_mohiggs->c_SSS(vertex,vanz); 
}

void Interaction_Model_THDM::c_SSV(std::vector<Single_Vertex>& vertex,int& vanz) 
{ 
  p_mohiggs->c_SSV(vertex,vanz); 
}

void Interaction_Model_THDM::c_SSVV(std::vector<Single_Vertex>& vertex,int& vanz) 
{ 
  p_mohiggs->c_SSVV(vertex,vanz); 
}

void Interaction_Model_THDM::c_SSSS(std::vector<Single_Vertex>& vertex,int& vanz) 
{ 
  p_mohiggs->c_SSSS(vertex,vanz);
}

Interaction_Model_THDM::~Interaction_Model_THDM()
{
  delete p_moew;
  delete p_moqcd;
  delete p_mohiggs;
}
