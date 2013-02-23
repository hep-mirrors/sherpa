#include "MODEL/Interaction_Models/Interaction_Model_Base.H"
#include "MODEL/Interaction_Models/Interaction_Model_SM.H"
#include "ATOOLS/Math/Kabbala.H"

namespace MODEL {
  class Interaction_Model_SM_TauPi : public Interaction_Model_Base {
    Interaction_Model_SM * p_sm;
  public:
    Interaction_Model_SM_TauPi(MODEL::Model_Base *,std::string,std::string);
    void c_FFV(std::vector<Single_Vertex>& s,int & i) { p_sm->c_FFV(s, i); }
    void c_FFS(std::vector<Single_Vertex>& s,int & i);
    void c_VVV(std::vector<Single_Vertex>& s,int & i) { p_sm->c_VVV(s, i); }
    void c_VVS(std::vector<Single_Vertex>& s,int & i) { p_sm->c_VVS(s, i); }
    void c_SSV(std::vector<Single_Vertex>& s,int & i) { p_sm->c_SSV(s, i); }
    void c_SSS(std::vector<Single_Vertex>& s,int & i) { p_sm->c_SSS(s, i); }
    void c_VVVV(std::vector<Single_Vertex>& s,int & i) { p_sm->c_VVVV(s, i); }
    void c_SSVV(std::vector<Single_Vertex>& s,int & i) { p_sm->c_SSVV(s, i); }
    void c_SSSS(std::vector<Single_Vertex>& s,int & i) { p_sm->c_SSSS(s, i); }
    void c_FFT(std::vector<Single_Vertex>& s,int& i) { p_sm->c_FFT(s, i); }
    void c_VVT(std::vector<Single_Vertex>& s,int& i) { p_sm->c_VVT(s, i); }
    void c_SST(std::vector<Single_Vertex>& s,int& i) { p_sm->c_SST(s, i); }
    void c_VVVT(std::vector<Single_Vertex>& s,int& i) { p_sm->c_VVVT(s, i); }
    void c_FFVT(std::vector<Single_Vertex>& s,int& i) { p_sm->c_FFVT(s, i); }
    void c_SSST(std::vector<Single_Vertex>& s,int& i) { p_sm->c_SSST(s, i); }
    ~Interaction_Model_SM_TauPi();
};
}


#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"
#include <stdio.h>


using namespace MODEL;
using namespace ATOOLS;
using namespace std;

DECLARE_GETTER(Interaction_Model_SM_TauPi,"SM+TauPi",
	       Interaction_Model_Base,Interaction_Model_Arguments);

Interaction_Model_Base *ATOOLS::Getter
<Interaction_Model_Base,Interaction_Model_Arguments,Interaction_Model_SM_TauPi>::
operator()(const Interaction_Model_Arguments &args) const
{
  return new Interaction_Model_SM_TauPi
    (args.p_model,args.m_cplscheme,args.m_yukscheme);
}

void ATOOLS::Getter<Interaction_Model_Base,Interaction_Model_Arguments,
		    Interaction_Model_SM_TauPi>::PrintInfo
(std::ostream &str,const size_t width) const
{ 
  str<<"The Standard Model plus TauPi";
}

Interaction_Model_SM_TauPi::Interaction_Model_SM_TauPi(MODEL::Model_Base * _model,
					   std::string _cplscheme,std::string _yukscheme) :
  Interaction_Model_Base("SM+TauPi",_model,_cplscheme,_yukscheme)
{ 
  p_sm = new Interaction_Model_SM(p_model,_cplscheme,_yukscheme);
}

void Interaction_Model_SM_TauPi::c_FFS(std::vector<Single_Vertex>& vertex,int& vanz)
{
  p_sm->c_FFS(vertex,vanz);

  Flavour tau(kf_tau);
  Flavour nutau(kf_nutau);
  Flavour pi_plus(kf_pi_plus);
  Kabbala GF = Kabbala(std::string("G_F"), ScalarConstant(std::string("GF")));
  Kabbala f_pi = Kabbala(std::string("F_PI"), ScalarConstant(std::string("F_PI")));
  Kabbala PL    = Kabbala(string("P_L"),1.);
  Kabbala PR    = Kabbala(string("P_R"),1.);
  Kabbala M_I   = Kabbala(string("i"),Complex(0.,1.));
  Kabbala root2 = Kabbala(string("\\sqrt{2}"),sqrt(2.));

  Kabbala kcpl0 = -M_I*f_pi*GF/root2; //*V_ud
  Kabbala kcpl1 = Kabbala(std::string("0"),0.);

  vertex[vanz].in[0] = tau;
  vertex[vanz].in[1] = pi_plus.Bar();
  vertex[vanz].in[2] = nutau;

  vertex[vanz].cpl[0]  = kcpl0;
  vertex[vanz].cpl[1]  = kcpl1;
  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();

  vertex[vanz].Color.push_back(Color_Function(cf::None));

  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("TAUPI",LF_Key()));

  vertex[vanz].on      = 1;
  vertex.push_back(Single_Vertex());vanz++;
}

Interaction_Model_SM_TauPi::~Interaction_Model_SM_TauPi()
{
  delete p_sm;
}
