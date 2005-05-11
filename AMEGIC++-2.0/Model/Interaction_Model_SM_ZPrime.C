#include "Interaction_Model_SM_ZPrime.H"
#include "MathTools.H"
#include "Message.H"
#include "Run_Parameter.H"
#include <stdio.h>


using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;


Interaction_Model_SM_ZPrime::Interaction_Model_SM_ZPrime(MODEL::Model_Base * _model,
					   std::string _cplscheme,std::string _yukscheme) :
  Interaction_Model_Base(_model,_cplscheme,_yukscheme)
{ // The Standard Model is part of this extended model
  p_moSM  = new Interaction_Model_SM(p_model,_cplscheme,_yukscheme);

  // set up constants for the model
  double Ecms2 = sqr(rpa.gen.Ecms());

  // sin and cos of the Weinberg angle
  sintW = Kabbala(std::string("\\sin\\theta_W"),
		  sqrt(ScalarConstant(std::string("sin2_thetaW"))));
  costW = Kabbala(std::string("\\cos\\theta_W"),
		  sqrt(1.-ScalarConstant(std::string("sin2_thetaW"))));


  // coupling constant of the gamma
  g1    = Kabbala(string("g_1"),
		  sqrt(4.*M_PI*ScalarFunction(std::string("alpha_QED"),Ecms2)));
  // coupling constant of Z'
  gP = Kabbala(string("g_1/\\cos\\theta_W"), g1.Value()/costW.Value());

  PL    = Kabbala(string("P_L"),1.);
  PR    = Kabbala(string("P_R"),1.);
// unneeded constants from original EW-Model
//  M_I   = Kabbala(string("i"),Complex(0.,1.));
//  root2 = Kabbala(string("\\sqrt{2}"),sqrt(2.));
//  vev   = Kabbala(string("v_{EW}"),ScalarConstant(std::string("vev")));
};

void Interaction_Model_SM_ZPrime::c_FFV(Single_Vertex* vertex,int& vanz)
{
// create the vertices for the standard model
  p_moSM->c_FFV(vertex,vanz);


// create the ZPrime-Flavour
  Flavour flZPrime(kf::ZPrime);

  // only generate vertices if ZPrime is on
  if (flZPrime.IsOn()) {
    /* Particle information
    PRINT_INFO("ZPrime is on");
    cout << kf::ZPrime; PRINT_INFO(" << Particle Number: ");
    cout << flZPrime.Charge(); PRINT_INFO(" << Charge");
    cout << flZPrime.Mass(); PRINT_INFO(" << Mass");
    cout << flZPrime.Spin(); PRINT_INFO(" << Spin"); */

    // list of all fermions that can undergo p -> Z' + p reactions
    int PossibleFermions[12] = {1,2,3,4,5,6,11,12,13,14,15,16};
    // scan the list
    for (int i=0; i<12; i++) {
      int FermionNumber = PossibleFermions[i];
      Flavour flFermion = Flavour(kf::code(FermionNumber));
      PRINT_INFO("Creating Vertices for Particle");
      cout << FermionNumber << " : " << flFermion <<"\n";

      /*
      // if the fermion is on it can undergo the p -> p + Z' reaction
      if (flFermion.IsOn()) {
        // create the vertex for that particular fermion and a Z'
        Kabbala kcpl0: // *** add value here
	Kabbala kcpl1; // *** add value here

	vertex[vanz].in[0] = vertex[vanz].in[2] = flFermion;
	vertex[vanz].in[1] = flZPrime;
	vertex[vanz].cpl[0] = kcpl0.Value();
	vertex[vanz].cpl[1] = kcpl1.Value();
	vertex[vanz].cpl[2] = 0.;
	vertex[vanz].cpl[3] = 0.;
	// PR = "Parton right handed" ???
	vertex[vanz].Str = (kcpl0*PR+kcpl1*PL).String(); 
	}; */
    };
  };
}


void Interaction_Model_SM_ZPrime::c_VVV(Single_Vertex* vertex,int& vanz)
{
  p_moSM->c_VVV(vertex,vanz);
}
void Interaction_Model_SM_ZPrime::c_VVVV(Single_Vertex* vertex,int& vanz)
{
  p_moSM->c_VVVV(vertex,vanz);
}

void Interaction_Model_SM_ZPrime::c_FFS(Single_Vertex* vertex,int& vanz)  { p_moSM->c_FFS(vertex,vanz); }
void Interaction_Model_SM_ZPrime::c_VVS(Single_Vertex* vertex,int& vanz)  { p_moSM->c_VVS(vertex,vanz); }
void Interaction_Model_SM_ZPrime::c_SSS(Single_Vertex* vertex,int& vanz)  { p_moSM->c_SSS(vertex,vanz); }
void Interaction_Model_SM_ZPrime::c_SSVV(Single_Vertex* vertex,int& vanz) { p_moSM->c_SSVV(vertex,vanz); }
void Interaction_Model_SM_ZPrime::c_SSSS(Single_Vertex* vertex,int& vanz) { p_moSM->c_SSSS(vertex,vanz); }

Interaction_Model_SM_ZPrime::~Interaction_Model_SM_ZPrime()
{
  delete p_moSM;
}
