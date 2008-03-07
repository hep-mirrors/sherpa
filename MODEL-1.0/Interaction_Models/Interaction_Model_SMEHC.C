#include "Interaction_Model_SMEHC.H"
#include "MathTools.H"
#include "Message.H"

using namespace MODEL;
using namespace ATOOLS;
using namespace std;

DECLARE_GETTER(Interaction_Model_SMEHC_Getter,"SM+EHC",
	       Interaction_Model_Base,Interaction_Model_Arguments);

Interaction_Model_Base *Interaction_Model_SMEHC_Getter::
operator()(const Interaction_Model_Arguments &args) const
{
  return new Interaction_Model_SMEHC
    (args.p_model,args.m_cplscheme,args.m_yukscheme);
}

void Interaction_Model_SMEHC_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{ 
  str<<"Standard Model + Effective Higgs Coupling"; 
}

Interaction_Model_SMEHC::Interaction_Model_SMEHC(MODEL::Model_Base * _model,
					       std::string _cplscheme,std::string _yukscheme) :
  Interaction_Model_Base("SM+EHC",_model,_cplscheme,_yukscheme)
{ 
  m_loops=true;
  p_mosm      = new Interaction_Model_SM(p_model,_cplscheme,_yukscheme); 
  p_mohl      = new Interaction_Model_EHC(p_model,_cplscheme,_yukscheme); 
}

void Interaction_Model_SMEHC::c_FFV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  p_mosm->c_FFV(vertex,vanz);
}

void Interaction_Model_SMEHC::c_VVV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  p_mosm->c_VVV(vertex,vanz);
}
void Interaction_Model_SMEHC::c_VVVV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  p_mosm->c_VVVV(vertex,vanz);
  p_mohl->c_VVVV(vertex,vanz);
}
void Interaction_Model_SMEHC::c_FFS(std::vector<Single_Vertex>& vertex,int& vanz) 
{
  p_mosm->c_FFS(vertex,vanz);
}
void Interaction_Model_SMEHC::c_VVS(std::vector<Single_Vertex>& vertex,int& vanz) 
{
  p_mosm->c_VVS(vertex,vanz);
  p_mohl->c_VVS(vertex,vanz);
}
void Interaction_Model_SMEHC::c_SSS(std::vector<Single_Vertex>& vertex,int& vanz) 
{
  p_mosm->c_SSS(vertex,vanz);
}
void Interaction_Model_SMEHC::c_SSSS(std::vector<Single_Vertex>& vertex,int& vanz) 
{ 
  p_mosm->c_SSSS(vertex,vanz); 
}
void Interaction_Model_SMEHC::c_SSVV(std::vector<Single_Vertex>& vertex,int& vanz) 
{ 
  p_mosm->c_SSVV(vertex,vanz); 
}


Interaction_Model_SMEHC::~Interaction_Model_SMEHC()
{
  delete  p_mosm;
  delete  p_mohl;
}

















