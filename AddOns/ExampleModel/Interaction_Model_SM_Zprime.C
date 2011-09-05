#include "MODEL/Interaction_Models/Interaction_Model_SM.H"

#include "ATOOLS/Math/Kabbala.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include <stdio.h>

using namespace std;
using namespace ATOOLS;

namespace MODEL {
  
  class Interaction_Model_SM_Zprime : public Interaction_Model_SM {
    Kabbala g1, g2, sintW, costW, PL, PR, M_I, alphaLR;
  public:
    Interaction_Model_SM_Zprime(Model_Base *,string,string);

    // overwrite SM FFV vertex filler
    void c_FFV(vector<Single_Vertex>&, int &);
  };
  
  Interaction_Model_SM_Zprime::Interaction_Model_SM_Zprime(Model_Base * model, 
                                                           string cplscheme,
                                                           string yukscheme) :
  Interaction_Model_SM(model, cplscheme, yukscheme)
  { 
    m_code = "SM+Zprime";
    
    // set up constants for the model
    sintW = Kabbala(string("\\sin\\theta_W"),
                    sqrt(ScalarConstant(string("sin2_thetaW"))));
    costW = Kabbala(string("\\cos\\theta_W"),
                    sqrt(1.-ScalarConstant(string("sin2_thetaW"))));
    
    // coupling constants
    double Ecms2 = sqr(rpa->gen.Ecms());
    g1  = Kabbala(string("g_1"),
                  sqrt(4.*M_PI*ScalarFunction(string("alpha_QED"),Ecms2)));
    g2  = Kabbala(string("g_1/\\cos\\theta_W"), g1.Value()/costW.Value());
    
    PL  = Kabbala(string("P_L"),1.);
    PR  = Kabbala(string("P_R"),1.);
    M_I = Kabbala(string("i"),Complex(0.,1.));
    
    // the parameter specifying the LR model
    // - sqrt(2.) will describe a totally LR-symm model
    // - sqrt(2./3.) describes an E6-inspired model
    alphaLR = Kabbala(string("\\alpha_{LR}"), sqrt(2./3.));
  }
  
  void Interaction_Model_SM_Zprime::c_FFV(vector<Single_Vertex>& vertex,int& vanz)
  {
    // create the vertices for the standard model
    Interaction_Model_SM::c_FFV(vertex,vanz);
    
    
    // create FFV vertices with Z' if it's on
    kf_code kfZp=32;
    Flavour flZprime(kfZp);
    if (flZprime.IsOn()) {
      
      // parse through all fermions than couple to Z' and create vertices
      int PossibleFermions[12] = {1,2,3,4,5,6,11,12,13,14,15,16};
      for (int i=0; i<12; i++) {
        
        // initialize the currently parsed fermion
        int FermionNumber = PossibleFermions[i];
        Flavour flFermion = Flavour((kf_code)(FermionNumber));
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
          else {kcpl0 = -M_I * g2 * (Y3R * alphaLR + (L-B)/(alphaLR*2));};
          Kabbala kcpl1 = -M_I * g2 * (L-B) / (alphaLR*2);
          
          // set couplings and particle info for current vertex
          vertex[vanz].in[0] = flFermion;
          vertex[vanz].in[1] = flZprime;
          vertex[vanz].in[2] = Flavour((kf_code)(FermionNumber));
          vertex[vanz].cpl[0] = kcpl0;
          vertex[vanz].cpl[1] = kcpl1;
          vertex[vanz].Str = (kcpl0*PR+kcpl1*PL).String(); 
          
          // Color Function for vertex
          
          if (flFermion.Strong()) {
            vertex[vanz].Color.push_back(Color_Function(cf::D));;
            vertex[vanz].Color.back().SetParticleArg(0,2);
            vertex[vanz].Color.back().SetStringArg('0','2');
          } 
          else 
            vertex[vanz].Color.push_back(Color_Function(cf::None));;
          
          // Lorenz function for vertex
          
          vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
          vertex[vanz].Lorentz.back()->SetParticleArg(1);
          
          vertex[vanz].on     = 1;
          vertex.push_back(Single_Vertex());vanz++; 
        }; 
      };
    };
  }
}


// Now follows some magic to make the model known to Sherpa and print a summary

using namespace MODEL;

DECLARE_GETTER(Interaction_Model_SM_Zprime_Getter,"SM+Zprime",
               Interaction_Model_Base,Interaction_Model_Arguments);

Interaction_Model_Base *Interaction_Model_SM_Zprime_Getter::
    operator()(const Interaction_Model_Arguments &args) const
{
  return new Interaction_Model_SM_Zprime
      (args.p_model,args.m_cplscheme,args.m_yukscheme);
}

void Interaction_Model_SM_Zprime_Getter::PrintInfo
    (ostream &str,const size_t width) const
{ 
  str<<"Standard Model + Z'"; 
}
