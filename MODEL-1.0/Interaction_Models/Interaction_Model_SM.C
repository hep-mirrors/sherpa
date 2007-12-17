#include "Interaction_Model_SM.H"
#include "MathTools.H"
#include "Message.H"
#include <stdio.h>


using namespace MODEL;
using namespace ATOOLS;
using namespace std;

Interaction_Model_SM::Interaction_Model_SM(MODEL::Model_Base * _model,
					   std::string _cplscheme,std::string _yukscheme) :
  Interaction_Model_Base(_model,_cplscheme,_yukscheme)
{ 
  p_moqcd  = new Interaction_Model_QCD(p_model,_cplscheme,_yukscheme); 
  p_moew   = new Interaction_Model_EW(p_model,_cplscheme,_yukscheme); 
}

void Interaction_Model_SM::c_FFV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  p_moqcd->c_FFV(vertex,vanz);
  p_moew->c_FFV(vertex,vanz);
}

void Interaction_Model_SM::c_VVV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  p_moqcd->c_VVV(vertex,vanz);
  p_moew->c_VVV(vertex,vanz);
}
void Interaction_Model_SM::c_VVVV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  p_moqcd->c_VVVV(vertex,vanz);
  p_moew->c_VVVV(vertex,vanz);
}

void Interaction_Model_SM::c_FFS(std::vector<Single_Vertex>& vertex,int& vanz)  { p_moew->c_FFS(vertex,vanz); }
void Interaction_Model_SM::c_VVS(std::vector<Single_Vertex>& vertex,int& vanz)  { p_moew->c_VVS(vertex,vanz); }
void Interaction_Model_SM::c_SSS(std::vector<Single_Vertex>& vertex,int& vanz)  { p_moew->c_SSS(vertex,vanz); }
void Interaction_Model_SM::c_SSVV(std::vector<Single_Vertex>& vertex,int& vanz) { p_moew->c_SSVV(vertex,vanz); }
void Interaction_Model_SM::c_SSSS(std::vector<Single_Vertex>& vertex,int& vanz) { p_moew->c_SSSS(vertex,vanz); }

Interaction_Model_SM::~Interaction_Model_SM()
{
  delete p_moqcd;
  delete p_moew;
}
