#include "MODEL/Interaction_Models/Interaction_Model_SM_U1_B.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include <stdio.h>


using namespace MODEL;
using namespace ATOOLS;
using namespace std;

DECLARE_GETTER(Interaction_Model_SM_U1_B_Getter,"SM+U1_B",
	       Interaction_Model_Base,Interaction_Model_Arguments);

Interaction_Model_Base *Interaction_Model_SM_U1_B_Getter::
operator()(const Interaction_Model_Arguments &args) const
{
  return new Interaction_Model_SM_U1_B
    (args.p_model,args.m_cplscheme,args.m_yukscheme);
}

void Interaction_Model_SM_U1_B_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{ 
  str<<"The Standard Model + leptophobic Z'"; 
}

Interaction_Model_SM_U1_B::
Interaction_Model_SM_U1_B(MODEL::Model_Base * _model,
				string _cplscheme,string _yukscheme) :
  Interaction_Model_Base("SM+U1_B",_model,_cplscheme,_yukscheme)
{ 
  p_moew  = new Interaction_Model_EW(p_model,_cplscheme,_yukscheme); 
  p_moqcd = new Interaction_Model_QCD(p_model,_cplscheme,_yukscheme); 

  PL    = Kabbala(string("P_L"),1.);
  PR    = Kabbala(string("P_R"),1.);
  M_I   = Kabbala(string("i"),Complex(0.,1.));
  g1_B  = Kabbala(string("g'_1"),ScalarConstant(std::string("g'_1")));
}

void Interaction_Model_SM_U1_B::c_FFV(vector<Single_Vertex>& vertex,int& vanz)
{
  p_moew->c_FFV(vertex,vanz);
  p_moqcd->c_FFV(vertex,vanz);

  msg_Out()<<METHOD<<":\n";

  Flavour Zprime(kf_Z0_2), quark1, quark2;
  Kabbala kcpl0,kcpl1;
  if (!Zprime.IsOn()) return;
  for (short int i=1;i<=6;i+=2) {
    quark1 = Flavour((kf_code)(i));
    if (!quark1.IsOn() || !quark1.Strong()) continue;
    kcpl0 = kcpl1 = -g1_B*M_I;
    if (ATOOLS::IsZero(kcpl0.Value())) continue;
    vertex[vanz].in[0]   = quark1;
    vertex[vanz].in[1]   = Zprime;
    vertex[vanz].in[2]   = quark1;
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = vertex[vanz].cpl[0];
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    
    
    vertex[vanz].Color.push_back(Color_Function(cf::D));     
    vertex[vanz].Color.back().SetParticleArg(0,2);     
    vertex[vanz].Color.back().SetStringArg('0','2');     
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
    vertex[vanz].Lorentz.back()->SetParticleArg(1);     

    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;

    msg_Out()<<"   add "<<quark1<<"-"<<Zprime<<"-"<<quark1<<"-vertex.\n";
  }
  for (short int i=2;i<=6;i+=2) {
    quark1 = Flavour((kf_code)(i));
    if (!quark1.IsOn() || !quark1.Strong()) continue;
    for (short int j=i;j<=6;j+=2) {
      quark2 = Flavour((kf_code)(j));
      if (!quark2.IsOn() || !quark1.Strong()) continue;
      kcpl0 = kcpl1 = -g1_B*M_I*K_UpMix(j/2-1,i/2-1);
      if (ATOOLS::IsZero(kcpl0.Value())) continue;
      vertex[vanz].in[0]   = quark1;
      vertex[vanz].in[1]   = Zprime;
      vertex[vanz].in[2]   = quark2;
      vertex[vanz].cpl[0]  = kcpl0;
      vertex[vanz].cpl[1]  = vertex[vanz].cpl[0];
      vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    
    
      vertex[vanz].Color.push_back(Color_Function(cf::D));     
      vertex[vanz].Color.back().SetParticleArg(0,2);     
      vertex[vanz].Color.back().SetStringArg('0','2');     
      vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
      vertex[vanz].Lorentz.back()->SetParticleArg(1);     
      
      vertex[vanz].on      = 1;
      vertex.push_back(Single_Vertex());vanz++;
      
      msg_Out()<<"   add "<<quark1<<"-"<<Zprime<<"-"<<quark2<<"-vertex.\n";
    }
  }
}

void Interaction_Model_SM_U1_B::c_VVV(vector<Single_Vertex>& vertex,int& vanz)
{
  p_moew->c_VVV(vertex,vanz);
  p_moqcd->c_VVV(vertex,vanz);
}
void Interaction_Model_SM_U1_B::c_VVVV(vector<Single_Vertex>& vertex,int& vanz)
{
  p_moew->c_VVVV(vertex,vanz);
  p_moqcd->c_VVVV(vertex,vanz);
}

void Interaction_Model_SM_U1_B::c_FFS(vector<Single_Vertex>& vertex,int& vanz)
{ 
  p_moew->c_FFS(vertex,vanz);
  p_moqcd->c_FFS(vertex,vanz);
}

void Interaction_Model_SM_U1_B::c_VVS(vector<Single_Vertex>& vertex,int& vanz)
{ 
  p_moew->c_VVS(vertex,vanz);
  p_moqcd->c_VVS(vertex,vanz);
}

void Interaction_Model_SM_U1_B::c_SSS(vector<Single_Vertex>& vertex,int& vanz)
{ 
  p_moew->c_SSS(vertex,vanz);
  p_moqcd->c_SSS(vertex,vanz);
}

void Interaction_Model_SM_U1_B::c_SSVV(vector<Single_Vertex>& vertex,int& vanz)
{ 
  p_moew->c_SSVV(vertex,vanz);
  p_moqcd->c_SSVV(vertex,vanz);
}

void Interaction_Model_SM_U1_B::c_SSSS(vector<Single_Vertex>& vertex,int& vanz)
{ 
  p_moew->c_SSSS(vertex,vanz);
  p_moqcd->c_SSSS(vertex,vanz);
}

Kabbala Interaction_Model_SM_U1_B::K_UpMix(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("\\tilde U_{")+string(hi)+string(hj)+string("}"),
		 ComplexMatrixElement(std::string("UpMix"),i,j));
} 
  
Interaction_Model_SM_U1_B::~Interaction_Model_SM_U1_B()
{
  delete p_moew;
  delete p_moqcd;
}
