#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/Run_Parameter.H"


#include "YFS/Main/Griffin_Interface.H"

using namespace PHASIC;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;


namespace Griffin {



Griffin_Interface::Griffin_Interface(){
  RegisterDefaults();
}

Griffin_Interface::~Griffin_Interface() {}



void Griffin_Interface::RegisterDefaults()
{
  Settings& s = Settings::GetMainSettings();
  s["Griffin_Order"].SetDefault(1);
  s["Griffin_Delta_Alpha"].SetDefault(-1);
}



bool Griffin_Interface::Initialize(const Process_Info& pi)
{
  Settings& s = Settings::GetMainSettings();
  double GF = s["GF"].Get<double>();
  m_mode = s["Griffin_Order"].Get<int>();
  double delap = s["Griffin_Delta_Alpha"].Get<double>();
  double mz2 = sqr(Flavour(kf_Z).Mass());
  if(IsEqual(delap,-1)){
    delap = 1-(*aqed)(0)/(*aqed)(mz2);
  }
  m_inital = pi.m_ii.GetExternal()[0];
  m_final = pi.m_fi.GetExternal()[0];


  griffin_input.set(MZ, Flavour(kf_Z).Mass());
  griffin_input.set(MW, Flavour(kf_Wplus).Mass());
  griffin_input.set(al, s_model->ScalarConstant("alpha_QED"));
  griffin_input.set(als, s_model->ScalarConstant("alpha_S"));
  griffin_input.set(gamZ, Flavour(kf_Z).Width());
  griffin_input.set(GamW, Flavour(kf_Wplus).Width());
  griffin_input.set(MH, Flavour(kf_h0).Mass());
  griffin_input.set(MT, Flavour(kf_t).Mass());
  griffin_input.set(MB, Flavour(kf_b).Mass()); // MSbar mass at scale mu=MZ for mb(mb)=4.20
  griffin_input.set(MC, Flavour(kf_c).Mass()); // MSbar mass at scale mu=MZ
  griffin_input.set(G_MS, Flavour(kf_s).Mass());
  griffin_input.set(MU, Flavour(kf_u).Mass());
  griffin_input.set(MD, Flavour(kf_d).Mass());
  griffin_input.set(ML, Flavour(kf_tau).Mass());
  griffin_input.set(G_MM, Flavour(kf_mu).Mass());
  griffin_input.set(ME, Flavour(kf_e).Mass());
  griffin_input.set(Delal, delap);
  griffin_input.set(Gmu, GF);

  PrintLogo(msg->Info());
  std::cout << "Griffin initialization complete...\n";
  cout << endl << "Complex-pole masses: MW=" << griffin_input.get(MWc) << ", MZ=" 
    << griffin_input.get(MZc) << endl
    << "Delta Alpha = "<< delap << endl;



  return true;
}



double Griffin_Interface::EvaluateLoop(const Vec4D_Vector& momenta)
{
  if(momenta.size()!=4){
    THROW(fatal_error, "Griffin is for 2->2 only");
  }
    FA_SMNLO FAi(ELE, griffin_input), FAf(MUO, griffin_input);
    SW_SMNLO SWi(ELE, griffin_input), SWf(MUO, griffin_input);
  if (m_mode == 0){
    FA_SMLO FAi(ELE, griffin_input), FAf(MUO, griffin_input);
    SW_SMLO SWi(ELE, griffin_input), SWf(MUO, griffin_input);
  }
  else if (m_mode == 1){
    FA_SMNLO FAi(ELE, griffin_input), FAf(MUO, griffin_input);
    SW_SMNLO SWi(ELE, griffin_input), SWf(MUO, griffin_input);
  }
  else if (m_mode==2){
    FA_SMNNLO FAi(ELE, griffin_input), FAf(MUO, griffin_input);
    SW_SMNNLO SWi(ELE, griffin_input), SWf(MUO, griffin_input); 
  }
  double s = (momenta[2]+momenta[3]).Abs2();
  double cost = (momenta[2]+momenta[3]).CosTheta();
  // double cost = 0;
  // if(sqrt(s)<40) return 0;
  // cost=1;
  matel M(m_inital, m_final, VEC, VEC, FAi, FAf, SWi, SWf, s, cost, griffin_input);

  if(m_mode==0) matel M(m_inital, m_final, VEC, VEC, FAi, FAf, SWi, SWf, s, cost, griffin_input);
  else if(m_mode==1) mat_SMNLO M(m_inital, m_final, VEC, VEC, FAi, FAf, SWi, SWf, s, cost, griffin_input);
  else if (m_mode==2) mat_SMNNLO M(m_inital, m_final, VEC, VEC, FAi, FAf, SWi, SWf, s, cost, griffin_input);
  M.setkinvar(s, cost);
  Cplx resvv, resva, resav, resaa;

  M.setform(VEC, VEC);
  resvv = M.result();
  M.setform(AXV, VEC);
  resav = M.result();
  M.setform(VEC, AXV);
  resva = M.result();
  M.setform(AXV, AXV);
  resaa = M.result();
  double res;
  res = real((1+cost*cost)*(resvv*conj(resvv) + resav*conj(resav)
               + resva*conj(resva) + resaa*conj(resaa)) +
    + 4*cost*(resvv*conj(resaa) + resva*conj(resav)).real()
    -2*(1+cost*cost)*(resvv*conj(resav) + resva*conj(resaa)).real()
    -4*cost*(resvv*conj(resva) + resav*conj(resaa)));
  // return res*3*s*m_rescale_alpha/32/M_PI;
  // return m_rescale_alpha*3*res*s/(32*M_PI);
  return res*s*s*m_alpha;
}

void Griffin_Interface::EvaluateBorn(const Vec4D_Vector& momenta)
{
  const int NN = momenta.size();
  double fpp[NN][4];

}





void Griffin_Interface::PrintLogo(std::ostream &s){
  s<<"======================================================"<<endl;
  s<<"                                                      "<<endl;
  s<<"======================================================"<<endl;

  s<<"     ______ ____   ____ ______ ______ ____ _   __"<<endl;
  s<<"    / ____// __ \\ /  _// ____// ____//  _// | / /"<<endl;
  s<<"   / / __ / /_/ / / / / /_   / /_    / / /  |/ / "<<endl;
  s<<"  / /_/ // _, _/_/ / / __/  / __/  _/ / / /|  /  "<<endl;
  s<<"  \\____//_/ |_|/___//_/    /_/    /___//_/ |_/   "<<endl;

  s<<"======================================================"<<endl;
  s<<"                 version 1.0                          "<<endl;
  s<<"           Lisong Chen and Ayres Freitas             "<<endl;
  s<<"          https://arxiv.org/abs/2211.16272            "<<endl;
  s<<"======================================================"<<endl;
  rpa->gen.AddCitation
    (1,"Griffin is published under \\cite{Chen:2022dow}.");
}


}
