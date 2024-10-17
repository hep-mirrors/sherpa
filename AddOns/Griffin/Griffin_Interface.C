#include "AddOns/Griffin/Griffin_Interface.H"


Griffin::Griffin_Interface::Griffin_Interface() :
      ME_Generator_Base("Griffin") {RegisterDefaults();}

Griffin::Griffin_Interface::~Griffin_Interface() {}


// define statics
int Griffin::Griffin_Interface::m_inital=11;
int Griffin::Griffin_Interface::m_final=13;
Griffin::griffinorder::code Griffin::Griffin_Interface::m_order=griffinorder::nnlo;
griffin::inval Griffin::Griffin_Interface::m_griffin;

void Griffin::Griffin_Interface::RegisterDefaults()
{
   Scoped_Settings s{ Settings::GetMainSettings()["GRIFFIN"] };
   s["Order"].SetDefault(griffinorder::nnlo);
   s["Delta_Alpha"].SetDefault(-1);
}


bool Griffin::Griffin_Interface::Initialize(MODEL::Model_Base *const model,
          BEAM::Beam_Spectra_Handler *const beam,
          PDF::ISR_Handler *const isr,
        YFS::YFS_Handler *const yfs)
{
  PrintLogo(msg->Info());
  rpa->gen.AddCitation(
      1, "The Griffin library is described in \\cite{Chen:2022dow}.");
  Settings& s = Settings::GetMainSettings();
  Scoped_Settings ss{ Settings::GetMainSettings()["GRIFFIN"] };
  double GF = s["GF"].Get<double>();
  m_order = ss["Order"].Get<griffinorder::code>();
  double delap = ss["Delta_Alpha"].Get<double>();
  double mz2 = sqr(Flavour(kf_Z).Mass());
  if(IsEqual(delap,-1)){
    delap = 1-(*aqed)(0)/(*aqed)(mz2);
  }
  PRINT_VAR(s_model->ScalarConstant("alpha_QED"));
  m_griffin.set(MZ, Flavour(kf_Z).Mass());
  m_griffin.set(MW, Flavour(kf_Wplus).Mass());
  m_griffin.set(al, s_model->ScalarConstant("alpha_QED"));
  m_griffin.set(als, s_model->ScalarConstant("alpha_S"));
  m_griffin.set(GamZ, Flavour(kf_Z).Width());
  m_griffin.set(GamW, Flavour(kf_Wplus).Width());
  m_griffin.set(MH, Flavour(kf_h0).Mass());
  m_griffin.set(MT, Flavour(kf_t).Mass());
  m_griffin.set(MB, Flavour(kf_b).Mass()); // MSbar mass at scale mu=MZ for mb(mb)=4.20
  m_griffin.set(MC, Flavour(kf_c).Mass()); // MSbar mass at scale mu=MZ
  m_griffin.set(MS, Flavour(kf_s).Mass());
  m_griffin.set(MU, Flavour(kf_u).Mass());
  m_griffin.set(MD, Flavour(kf_d).Mass());
  m_griffin.set(ML, Flavour(kf_tau).Mass());
  m_griffin.set(MM, Flavour(kf_mu).Mass());
  m_griffin.set(ME, Flavour(kf_e).Mass());
  m_griffin.set(Delal, delap);
  m_griffin.set(Gmu, GF);
  return true;
}

void Griffin::Griffin_Interface::RegisterProcess(const PHASIC::Process_Info &pi){
  // Since griffin is only 2->2 m_initial[0]=bar(m_initial[1])
  m_inital = pi.m_ii.GetExternal()[0];
  m_final = pi.m_fi.GetExternal()[0];
}

void Griffin::Griffin_Interface::RegisterProcess(const External_ME_Args& args){
  // Since griffin is only 2->2 m_initial[0]=bar(m_initial[1])
  m_inital = args.m_inflavs[0];
  m_final =  args.m_outflavs[0];
}

void Griffin::Griffin_Interface::EvaluateLoop(const Vec4D_Vector& momenta, METOOLS::DivArrD& virt){
  if(momenta.size()!=4){
    THROW(fatal_error, "Griffin library is for 2->2 scattering only");
  }
  // EvaluateLO(momenta, virt);
  // PRINT_VAR(virt.Finite());
  EvaluateNNLO(momenta, virt);
  // PRINT_VAR(virt.Finite());
  // if(m_order==griffinorder::lo) EvaluateLO(momenta, virt);
  // else if(m_order==griffinorder::nlo) EvaluateNLO(momenta, virt);
  // else EvaluateNNLO(momenta,virt);
  // PRINT_VAR(virt.Finite());
}

void Griffin::Griffin_Interface::EvaluateLO(const Vec4D_Vector& momenta, DivArrD &res){
  double s = (momenta[0]+momenta[1]).Abs2();
  double t = (momenta[0]-momenta[2]).Abs2();
  double cost = 1.+2.*t/s;
  if(cost > 1. || cost < -1){
    msg_Error()<<"CosTheta out of range in "<<METHOD<<endl;
  }
  // double cost = (momenta[2]).CosTheta(momenta[0]);
  // PRINT_VAR(costt/cost);
  // if(sqrt(s)<40) {
  //   res.Finite() = 0;
  // }
  FA_SMLO FAi(m_inital, m_griffin), FAf(m_final, m_griffin);
  // SW_SMLO SWi(m_inital, m_griffin), SWf(m_final, m_griffin);
  // cout << "sineff^i (LO+) = " << SWi.result() << endl;
  // cout << "sineff^f (LO+) = " << SWf.result() << endl;
  Cplx resvv, resva, resav, resaa;
  double sw=s_model->ComplexConstant("csin2_thetaW").real();
  // double  FAi = 0.0345, FAf = 0.0345;
  matel M(m_inital, m_final, VEC, VEC, FAi.result().real(), FAf.result().real(), sw, sw, s, cost, m_griffin);
  M.setkinvar(s, cost);

  M.setform(VEC, VEC);
  resvv = M.result();
  M.setform(AXV, VEC);
  resav = M.result();
  M.setform(VEC, AXV);
  resva = M.result();
  M.setform(AXV, AXV);
  resaa = M.result();
  res.Finite() = real((1+cost*cost)*(resvv*conj(resvv) + resav*conj(resav)
               + resva*conj(resva) + resaa*conj(resaa)) +
    + 4*cost*(resvv*conj(resaa) + resva*conj(resav)).real());
}

void Griffin::Griffin_Interface::EvaluateNLO(const Vec4D_Vector& momenta, DivArrD &res){
  double s = (momenta[0]+momenta[1]).Abs2();
  double t = (momenta[0]-momenta[2]).Abs2();
  double cost = 1.+2.*t/s;
  // double cost = 0;
  // if(sqrt(s)<40) {
  //   res.Finite() = 0;
  // }
  FA_SMNLO FAi(m_inital, m_griffin), FAf(m_final, m_griffin);
  // SW_SMNLO SWi(m_inital, m_griffin), SWf(m_final, m_griffin);
  // cout << "FAi = " << FAi.result() << endl;
  // cout << "FAi = " << FAf.result() << endl;
  double sw=s_model->ComplexConstant("csin2_thetaW").real();

  mat_SMNLO M(m_inital, m_final, VEC, VEC, FAi.result().real(), FAf.result().real(), sw, sw*1.001, s, cost, m_griffin);

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
  res.Finite() = real((1+cost*cost)*(resvv*conj(resvv) + resav*conj(resav)
               + resva*conj(resva) + resaa*conj(resaa)) +
    + 4*cost*(resvv*conj(resaa) + resva*conj(resav)).real());
    // -2*(1+cost*cost)*(resvv*conj(resav) + resva*conj(resaa)).real()
    // -4*cost*(resvv*conj(resva) + resav*conj(resaa)));
  // return res*3*s*s*m_rescale_alpha/32/M_PI;
  // res.Finite() *= 3*s/32/M_PI;
}

void Griffin::Griffin_Interface::EvaluateNNLO(const Vec4D_Vector& momenta, DivArrD &res){
  double s = (momenta[0]+momenta[1]).Abs2();
  double t = (momenta[0]-momenta[2]).Abs2();
  double cost = 1.+2.*t/s;
  // double cost = 0;
  if(sqrt(s)<40) {
    res.Finite() = 0;
  }
  FA_SMNNLO FAi(m_inital, m_griffin), FAf(m_final, m_griffin);
  // SW_SMNNLO SWi(m_inital, m_griffin), SWf(m_final, m_griffin);
  double sw=s_model->ComplexConstant("csin2_thetaW").real();
  mat_SMNNLO M(m_inital, m_final, VEC, VEC, FAi.result().real(), FAf.result().real(), sw, sw, s, cost, m_griffin);

  M.setkinvar(s, cost);
  Cplx resvv, resva, resav, resaa;
  Cplx res1, res2;
  res1 = M.result();
  res2 = M.resoffZ();

  M.setform(VEC, VEC);
  resvv = M.result();
  M.setform(AXV, VEC);
  resav = M.result();
  M.setform(VEC, AXV);
  resva = M.result();
  M.setform(AXV, AXV);
  resaa = M.result();
  res.Finite() = real((1+cost*cost)*(resvv*conj(resvv) + resav*conj(resav)
               + resva*conj(resva) + resaa*conj(resaa)) +
    + 4*cost*(resvv*conj(resaa) + resva*conj(resav)).real());
    // -2*(1+cost*cost)*(resvv*conj(resav) + resva*conj(resaa)).real()
    // -4*cost*(resvv*conj(resva) + resav*conj(resaa)));
  // return res*3*s*s*m_rescale_alpha/32/M_PI;
  // res.Finite() *= 3*s/32/M_PI;
  // PRINT_VAR(res.Finite());
}

void Griffin::Griffin_Interface::EvaluateBorn(const Vec4D_Vector& momenta, double& bornres)
{
  FA_SMLO FAi(m_inital, m_griffin), FAf(m_final, m_griffin);
  SW_SMLO SWi(m_inital, m_griffin), SWf(m_final, m_griffin);
  double s = (momenta[0]+momenta[1]).Abs2();
  double cost = (momenta[0]+momenta[1]).CosTheta();
  matel M(m_inital, m_final, VEC, VEC, FAi, FAf, SWi, SWf, s, cost, m_griffin);
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
  bornres = res;
}

int Griffin::Griffin_Interface::PerformTests()
{
  return 1;
}


std::istream &Griffin::operator>>(std::istream &str, griffinorder::code &mode)
{
  std::string tag;
  str>>tag;
  mode=griffinorder::nnlo;
  if      (tag.find("LO")!=std::string::npos)    mode=griffinorder::lo;
  else if (tag.find("NLO")!=std::string::npos)   mode=griffinorder::nlo;
  else if (tag.find("NNLO")!=std::string::npos)    mode=griffinorder::nnlo;
  else THROW(fatal_error, "Unknown GRIFFIN: Order");
  return str;
}

std::ostream &Griffin::operator<<(std::ostream &str,const griffinorder::code &ym)
{
  if      (ym==griffinorder::lo)     return str<<"LO";
  else if (ym==griffinorder::nlo)    return str<<"NLO";
  else if (ym==griffinorder::nnlo)   return str<<"NNLO";
  return str<<"unknown";
}

void Griffin::Griffin_Interface::PrintLogo(std::ostream &s){
  s<<"======================================================"<<endl;
  s<<"                                                      "<<endl;
  s<<"======================================================"<<endl;

  s<<"     ______ ____   ____ ______ ______ ____ _   __"<<endl;
  s<<"    / ____// __ \\ /  _// ____// ____//  _// | / /"<<endl;
  s<<"   / / __ / /_/ / / / / /_   / /_    / / /  |/ / "<<endl;
  s<<"  / /_/ // _, _/_/ / / __/  / __/  _/ / / /|  /  "<<endl;
  s<<"  \\____//_/ |_|/___//_/    /_/    /___//_/ |_/   "<<endl;

  s<<"======================================================"<<endl;
  s<<"                 version 1.1                          "<<endl;
  s<<"           Lisong Chen and Ayres Freitas             "<<endl;
  s<<"          https://arxiv.org/abs/2211.16272            "<<endl;
  s<<"======================================================"<<endl;
}




DECLARE_GETTER(Griffin::Griffin_Interface,"Griffin",PHASIC::ME_Generator_Base,PHASIC::ME_Generator_Key);

ME_Generator_Base *ATOOLS::Getter<PHASIC::ME_Generator_Base,PHASIC::ME_Generator_Key,
                                  Griffin::Griffin_Interface>::
operator()(const PHASIC::ME_Generator_Key &key) const
{
  return new Griffin::Griffin_Interface();
}

void ATOOLS::Getter<PHASIC::ME_Generator_Base,PHASIC::ME_Generator_Key,Griffin::Griffin_Interface>::
PrintInfo(std::ostream &str,const std::size_t width) const
{
  str<<"Interface to the Griffin loop ME generator";
}
