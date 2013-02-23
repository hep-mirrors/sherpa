#include "MODEL/Interaction_Models/Interaction_Model_ADD.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"

using namespace MODEL;
using namespace ATOOLS;
using namespace std;

DECLARE_GETTER(Interaction_Model_ADD,"ADD",
	       Interaction_Model_Base,Interaction_Model_Arguments);

Interaction_Model_Base *ATOOLS::Getter
<Interaction_Model_Base,Interaction_Model_Arguments,Interaction_Model_ADD>::
operator()(const Interaction_Model_Arguments &args) const
{
  return new Interaction_Model_ADD
    (args.p_model,args.m_cplscheme,args.m_yukscheme);
}

void ATOOLS::Getter<Interaction_Model_Base,Interaction_Model_Arguments,
		    Interaction_Model_ADD>::PrintInfo
(std::ostream &str,const size_t width) const
{ 
  str<<"The ADD model of large extra dimensions"; 
}

Interaction_Model_ADD::Interaction_Model_ADD(MODEL::Model_Base * _model,
					   std::string _cplscheme,std::string _yukscheme) :
  Interaction_Model_Base("ADD",_model,_cplscheme,_yukscheme)
{ 
  m_tensors=true;
  p_moqcd     = new Interaction_Model_QCD(p_model,_cplscheme,_yukscheme); 
  p_moew      = new Interaction_Model_EW(p_model,_cplscheme,_yukscheme); 
  p_moqcdgrav = new Interaction_Model_QCD_Grav(p_model,_cplscheme,_yukscheme); 
  p_moewgrav  = new Interaction_Model_EW_Grav(p_model,_cplscheme,_yukscheme); 
}

void Interaction_Model_ADD::c_FFV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  p_moqcd->c_FFV(vertex,vanz);
  p_moew->c_FFV(vertex,vanz);
}

void Interaction_Model_ADD::c_VVV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  p_moqcd->c_VVV(vertex,vanz);
  p_moew->c_VVV(vertex,vanz);
}
void Interaction_Model_ADD::c_VVVV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  p_moqcd->c_VVVV(vertex,vanz);
  p_moew->c_VVVV(vertex,vanz);
}

void Interaction_Model_ADD::c_FFS(std::vector<Single_Vertex>& vertex,int& vanz) {p_moew->c_FFS(vertex,vanz);}
void Interaction_Model_ADD::c_VVS(std::vector<Single_Vertex>& vertex,int& vanz) {p_moew->c_VVS(vertex,vanz);}
void Interaction_Model_ADD::c_SSS(std::vector<Single_Vertex>& vertex,int& vanz) {p_moew->c_SSS(vertex,vanz);}

void Interaction_Model_ADD::c_SSSS(std::vector<Single_Vertex>& vertex,int& vanz) { p_moew->c_SSSS(vertex,vanz); }
void Interaction_Model_ADD::c_SSVV(std::vector<Single_Vertex>& vertex,int& vanz) { p_moew->c_SSVV(vertex,vanz); }

void Interaction_Model_ADD::c_FFT(std::vector<Single_Vertex>& vertex,int& vanz)
{
  p_moewgrav->c_FFT(vertex,vanz);
}
void Interaction_Model_ADD::c_VVT(std::vector<Single_Vertex>& vertex,int& vanz)
{
  p_moqcdgrav->c_VVT(vertex,vanz);
  p_moewgrav->c_VVT(vertex,vanz);
}
void Interaction_Model_ADD::c_SST(std::vector<Single_Vertex>& vertex,int& vanz)
{
  p_moewgrav->c_SST(vertex,vanz);
}
void Interaction_Model_ADD::c_VVVT(std::vector<Single_Vertex>& vertex,int& vanz)
{
  p_moqcdgrav->c_VVVT(vertex,vanz);
  p_moewgrav->c_VVVT(vertex,vanz);
}
void Interaction_Model_ADD::c_SSST(std::vector<Single_Vertex>& vertex,int& vanz)
{
  p_moewgrav->c_SSST(vertex,vanz);
}
void Interaction_Model_ADD::c_FFVT(std::vector<Single_Vertex>& vertex,int& vanz)
{
  p_moqcdgrav->c_FFVT(vertex,vanz);
  p_moewgrav->c_FFVT(vertex,vanz);
}

Interaction_Model_ADD::~Interaction_Model_ADD()
{
  delete  p_moqcd;
  delete  p_moew;
  delete  p_moqcdgrav;
  delete  p_moewgrav;
}

















