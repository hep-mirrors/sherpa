#include "Interaction_Model_SM_AGC.H"
#include "MathTools.H"
#include "Message.H"
#include <stdio.h>


using namespace MODEL;
using namespace ATOOLS;
using namespace std;

Interaction_Model_SM_AGC::Interaction_Model_SM_AGC(MODEL::Model_Base * _model,
					   std::string _cplscheme,std::string _yukscheme) :
  Interaction_Model_Base(_model,_cplscheme,_yukscheme)
{ 
  p_moqcd  = new Interaction_Model_QCD(p_model,_cplscheme,_yukscheme); 
  p_moaew  = new Interaction_Model_AEW(p_model,_cplscheme,_yukscheme); 
}

void Interaction_Model_SM_AGC::c_FFV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  p_moqcd->c_FFV(vertex,vanz);
  p_moaew->c_FFV(vertex,vanz);
}

void Interaction_Model_SM_AGC::c_VVV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  p_moqcd->c_VVV(vertex,vanz);
  p_moaew->c_VVV(vertex,vanz);
}
void Interaction_Model_SM_AGC::c_VVVV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  p_moqcd->c_VVVV(vertex,vanz);
  p_moaew->c_VVVV(vertex,vanz);
}

void Interaction_Model_SM_AGC::c_FFS(std::vector<Single_Vertex>& vertex,int& vanz)  { p_moaew->c_FFS(vertex,vanz); }
void Interaction_Model_SM_AGC::c_VVS(std::vector<Single_Vertex>& vertex,int& vanz)  { p_moaew->c_VVS(vertex,vanz); }
void Interaction_Model_SM_AGC::c_SSS(std::vector<Single_Vertex>& vertex,int& vanz)  { p_moaew->c_SSS(vertex,vanz); }
void Interaction_Model_SM_AGC::c_SSVV(std::vector<Single_Vertex>& vertex,int& vanz) { p_moaew->c_SSVV(vertex,vanz); }
void Interaction_Model_SM_AGC::c_SSSS(std::vector<Single_Vertex>& vertex,int& vanz) { p_moaew->c_SSSS(vertex,vanz); }

Interaction_Model_SM_AGC::~Interaction_Model_SM_AGC()
{
  delete p_moqcd;
  delete p_moaew;
}
