#include "Interaction_Model_SMHL.H"
#include "MathTools.H"
#include "Message.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;

Interaction_Model_SMHL::Interaction_Model_SMHL(MODEL::Model_Base * _model,
					       std::string _cplscheme,std::string _yukscheme) :
  Interaction_Model_Base(_model,_cplscheme,_yukscheme)
{ 
  p_moqcd     = new Interaction_Model_QCD(p_model,_cplscheme,_yukscheme); 
  p_moew      = new Interaction_Model_EW(p_model,_cplscheme,_yukscheme); 
  p_mohl      = new Interaction_Model_HL(p_model,_cplscheme,_yukscheme); 
}

void Interaction_Model_SMHL::c_FFV(Single_Vertex* vertex,int& vanz)
{
  p_moqcd->c_FFV(vertex,vanz);
  p_moew->c_FFV(vertex,vanz);
}

void Interaction_Model_SMHL::c_VVV(Single_Vertex* vertex,int& vanz)
{
  p_moqcd->c_VVV(vertex,vanz);
  p_moew->c_VVV(vertex,vanz);
}
void Interaction_Model_SMHL::c_VVVV(Single_Vertex* vertex,int& vanz)
{
  p_moqcd->c_VVVV(vertex,vanz);
  p_moew->c_VVVV(vertex,vanz);
}

void Interaction_Model_SMHL::c_FFS(Single_Vertex* vertex,int& vanz) {p_moew->c_FFS(vertex,vanz);}
void Interaction_Model_SMHL::c_VVS(Single_Vertex* vertex,int& vanz) 
{
  p_moew->c_VVS(vertex,vanz);
  p_mohl->c_VVS(vertex,vanz);
}
void Interaction_Model_SMHL::c_SSS(Single_Vertex* vertex,int& vanz) {p_moew->c_SSS(vertex,vanz);}

void Interaction_Model_SMHL::c_SSSS(Single_Vertex* vertex,int& vanz) { p_moew->c_SSSS(vertex,vanz); }
void Interaction_Model_SMHL::c_SSVV(Single_Vertex* vertex,int& vanz) { p_moew->c_SSVV(vertex,vanz); }


Interaction_Model_SMHL::~Interaction_Model_SMHL()
{
  delete  p_moqcd;
  delete  p_moew;
  delete  p_mohl;
}

















