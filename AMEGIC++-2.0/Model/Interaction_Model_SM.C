#include "Interaction_Model_SM.H"
#include "MathTools.H"
#include "Message.H"
#include <stdio.h>


using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;

Interaction_Model_SM::Interaction_Model_SM(MODEL::Model_Base * _model,
					   std::string _cplscheme,std::string _yukscheme) :
  Interaction_Model_Base(_model,_cplscheme,_yukscheme)
{ 
  p_moqcd  = new Interaction_Model_QCD(p_model,_cplscheme,_yukscheme); 
  p_moew   = new Interaction_Model_EW(p_model,_cplscheme,_yukscheme); 
}

void Interaction_Model_SM::c_FFV(Single_Vertex* vertex,int& vanz)
{
  p_moqcd->c_FFV(vertex,vanz);
  p_moew->c_FFV(vertex,vanz);
}

void Interaction_Model_SM::c_VVV(Single_Vertex* vertex,int& vanz)
{
  p_moqcd->c_VVV(vertex,vanz);
  p_moew->c_VVV(vertex,vanz);
}
void Interaction_Model_SM::c_VVVV(Single_Vertex* vertex,int& vanz)
{
  p_moqcd->c_VVVV(vertex,vanz);
  p_moew->c_VVVV(vertex,vanz);
}

void Interaction_Model_SM::c_FFS(Single_Vertex* vertex,int& vanz)  { p_moew->c_FFS(vertex,vanz); }
void Interaction_Model_SM::c_VVS(Single_Vertex* vertex,int& vanz)  { p_moew->c_VVS(vertex,vanz); }
void Interaction_Model_SM::c_SSS(Single_Vertex* vertex,int& vanz)  { p_moew->c_SSS(vertex,vanz); }
void Interaction_Model_SM::c_SSVV(Single_Vertex* vertex,int& vanz) { p_moew->c_SSVV(vertex,vanz); }
void Interaction_Model_SM::c_SSSS(Single_Vertex* vertex,int& vanz) { p_moew->c_SSSS(vertex,vanz); }
