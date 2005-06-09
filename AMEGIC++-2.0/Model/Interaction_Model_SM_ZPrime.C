#include "Interaction_Model_SM_ZPrime.H"
#include "MathTools.H"
#include "Message.H"
#include "Run_Parameter.H"
#include <stdio.h>


using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;


Interaction_Model_SM_ZPrime::Interaction_Model_SM_ZPrime(
MODEL::Model_Base * _model, std::string _cplscheme,std::string _yukscheme)
: Interaction_Model_Base(_model,_cplscheme,_yukscheme)
{ // The Standard Model is part of this extended model
  p_moSM  = new Interaction_Model_SM(p_model,_cplscheme,_yukscheme);

  // set up constants for the model
  double Ecms2 = sqr(rpa.gen.Ecms());

  sintW = Kabbala(std::string("\\sin\\theta_W"),
		  sqrt(ScalarConstant(std::string("sin2_thetaW"))));
  costW = Kabbala(std::string("\\cos\\theta_W"),
		  sqrt(1.-ScalarConstant(std::string("sin2_thetaW"))));

  // coupling constants
  g1    = Kabbala(string("g_1"),
		  sqrt(4.*M_PI*ScalarFunction(std::string("alpha_QED"),Ecms2)));
  gP = Kabbala(string("g_1/\\cos\\theta_W"), g1.Value()/costW.Value());

  PL    = Kabbala(string("P_L"),1.);
  PR    = Kabbala(string("P_R"),1.);
  M_I   = Kabbala(string("i"),Complex(0.,1.));

  // the parameter specifying the LR model
  // - sqrt(2.) will describe a totally LR-symm model
  // - sqrt(2./3.) describes an E6-inspired model
  alphaLR = Kabbala(std::string("\\alpha_{LR}"), sqrt(2./3.));
};

void Interaction_Model_SM_ZPrime::c_FFV(std::vector<Single_Vertex>& vertex,int& vanz)
{
// create the vertices for the standard model
  p_moSM->c_FFV(vertex,vanz);


// create FFV vertices with Z' if it's on
  Flavour flZPrime(kf::ZPrime);
  if (flZPrime.IsOn()) {

    // parse through all fermions than couple to Z' and create vertices
    int PossibleFermions[12] = {1,2,3,4,5,6,11,12,13,14,15,16};
    for (int i=0; i<12; i++) {

      // initialize the currently parsed fermion
      int FermionNumber = PossibleFermions[i];
      Flavour flFermion = Flavour(kf::code(FermionNumber));
      Kabbala B = Kabbala(string("B_{")+flFermion.TexName()+string("}"),
                          flFermion.BaryonNumber());
      Kabbala L = Kabbala(string("L_{")+ flFermion.TexName()+string("}"),
                          flFermion.LeptonNumber());
      Kabbala Y3R = Kabbala(string("YR_{")+flFermion.TexName()+string("}"),
                            flFermion.IsoWeak());

      if (flFermion.IsOn()) {
        // create the vertex for that particular fermion and a Z'.
	// Right-handed neutrinos will not take part in any interaction.
        Kabbala kcpl0;
        if ((FermionNumber==12)||(FermionNumber==14)||(FermionNumber==16))
          {kcpl0 = Kabbala("0.0", 0.);}
        else {kcpl0 = -M_I * gP * (Y3R * alphaLR + (L-B)/(alphaLR*2));};
	Kabbala kcpl1 = -M_I * gP * (L-B) / (alphaLR*2);
	
	// set couplings and particle info for current vertex
	vertex[vanz].in[0] = flFermion;
	vertex[vanz].in[1] = flZPrime;
        vertex[vanz].in[2] = Flavour(kf::code(FermionNumber));
	vertex[vanz].cpl[0] = kcpl0.Value();
	vertex[vanz].cpl[1] = kcpl1.Value();
	vertex[vanz].cpl[2] = 0.;
	vertex[vanz].cpl[3] = 0.;
        vertex[vanz].Str = (kcpl0*PR+kcpl1*PL).String(); 

	// Color Function for vertex
	vertex[vanz].ncf       = 1;
	if (flFermion.Strong()) {
	  vertex[vanz].Color     = new Color_Function(cf::D);
	  vertex[vanz].Color->SetParticleArg(0,2);
	  vertex[vanz].Color->SetStringArg('0','2');
	} 
        else 
          vertex[vanz].Color = new Color_Function(cf::None);

	// Lorenz function for vertex
	vertex[vanz].nlf     = 1;
        vertex[vanz].Lorentz = new Lorentz_Function(lf::Gamma);
	vertex[vanz].Lorentz->SetParticleArg(1);

	vertex[vanz].on     = 1;
	vertex.push_back(Single_Vertex());vanz++; 
	}; 
    };
  };
}


// no other couplings of the ZPrime are built in, yet. All the following
// methods simply call the Standard Model to create its vertices.
void Interaction_Model_SM_ZPrime::c_VVV(std::vector<Single_Vertex>& vertex,int& vanz)
{ // ZPrime does not couple on any gauge boson.
  // Reason: None of the quantum numbers Z' does couple to (B, L, Y)
  //         are carried by any gauge boson.
  p_moSM->c_VVV(vertex,vanz); }
void Interaction_Model_SM_ZPrime::c_VVVV(std::vector<Single_Vertex>& vertex,int& vanz)
{ // No Z' interactions here - same reason as in c_VVV
  p_moSM->c_VVVV(vertex,vanz); }

// no interaction with the Higgs particles implemented => no interaction here
void Interaction_Model_SM_ZPrime::c_FFS(std::vector<Single_Vertex>& vertex,int& vanz)  { p_moSM->c_FFS(vertex,vanz); }
void Interaction_Model_SM_ZPrime::c_VVS(std::vector<Single_Vertex>& vertex,int& vanz)  { p_moSM->c_VVS(vertex,vanz); }
void Interaction_Model_SM_ZPrime::c_SSS(std::vector<Single_Vertex>& vertex,int& vanz)  { p_moSM->c_SSS(vertex,vanz); }
void Interaction_Model_SM_ZPrime::c_SSVV(std::vector<Single_Vertex>& vertex,int& vanz) { p_moSM->c_SSVV(vertex,vanz); }
void Interaction_Model_SM_ZPrime::c_SSSS(std::vector<Single_Vertex>& vertex,int& vanz) { p_moSM->c_SSSS(vertex,vanz); }

Interaction_Model_SM_ZPrime::~Interaction_Model_SM_ZPrime()
{
  delete p_moSM;
}
