#include "Interaction_Model_THDM.H"
#include "MathTools.H"
#include "Message.H"
#include <stdio.h>


using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;

Interaction_Model_THDM::Interaction_Model_THDM(MODEL::Model_Base * _model,
					       std::string _cplscheme,std::string _yukscheme) :
  Interaction_Model_Base(_model,_cplscheme,_yukscheme)
{ 
  p_mosm    = new Interaction_Model_SM(p_model,_cplscheme,_yukscheme); 
  p_mohiggs = new Interaction_Model_Higgs(p_model,_cplscheme,_yukscheme); 
}

void Interaction_Model_THDM::c_FFV(Single_Vertex* vertex,int& vanz)
{
  p_mosm->c_FFV(vertex,vanz);
  p_mohiggs->c_FFV(vertex,vanz);
}

void Interaction_Model_THDM::c_VVV(Single_Vertex* vertex,int& vanz)
{
  p_mosm->c_VVV(vertex,vanz);
  p_mohiggs->c_VVV(vertex,vanz);
}
void Interaction_Model_THDM::c_VVVV(Single_Vertex* vertex,int& vanz)
{
  p_mosm->c_VVVV(vertex,vanz);
  p_mohiggs->c_VVVV(vertex,vanz);
}

void Interaction_Model_THDM::c_FFS(Single_Vertex* vertex,int& vanz)  
{ 
  p_mohiggs->c_FFS(vertex,vanz); 
}

void Interaction_Model_THDM::c_VVS(Single_Vertex* vertex,int& vanz)  
{ 
  p_mohiggs->c_VVS(vertex,vanz); 
}

void Interaction_Model_THDM::c_SSS(Single_Vertex* vertex,int& vanz)  
{ 
  p_mohiggs->c_SSS(vertex,vanz); 
}

void Interaction_Model_THDM::c_SSV(Single_Vertex* vertex,int& vanz) 
{ 
  p_mohiggs->c_SSV(vertex,vanz); 
}

void Interaction_Model_THDM::c_SSVV(Single_Vertex* vertex,int& vanz) 
{ 
  p_mohiggs->c_SSVV(vertex,vanz); 
}

void Interaction_Model_THDM::c_SSSS(Single_Vertex* vertex,int& vanz) 
{ 
  p_mohiggs->c_SSSS(vertex,vanz);
}
