#include "AddOns/Griffin/Griffin_Interface.H"

using namespace MODEL;

Griffin::Griffin_Interface::Griffin_Interface() :
      ME_Generator_Base("Griffin") {RegisterDefaults();}

Griffin::Griffin_Interface::~Griffin_Interface() {}


// define statics
int Griffin::Griffin_Interface::m_inital=11;
int Griffin::Griffin_Interface::m_final=13;
double Griffin::Griffin_Interface::m_norm=1;
double Griffin::Griffin_Interface::m_B=0;
ew_scheme::code Griffin::Griffin_Interface::m_ewscheme=ew_scheme::alphamZ;
Griffin::griffinorder::code Griffin::Griffin_Interface::m_order=griffinorder::nnlo;
griffin::inval Griffin::Griffin_Interface::m_griffin;
// griffin::SMvalGMwMz Griffin::Griffin_Interface::m_griffinGMwMz;

void Griffin::Griffin_Interface::RegisterDefaults()
{
   Scoped_Settings s{ Settings::GetMainSettings()["GRIFFIN"] };
   s["Order"].SetDefault(griffinorder::nnlo);
   s["Delta_Alpha"].SetDefault(0.06);
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
  ew_scheme::code ewscheme = s["EW_SCHEME"].Get<ew_scheme::code>();
  // switch(ewscheme) {
  //   case ew_scheme::alpha0:
  //     delap=0.06;
  //     break;
  //   case ew_scheme::alphamZ:
  //     delap= 0.06;
  //     break;
  //   case 3:
  //     delap=0.06;
  //     break;
  //   default:
  //     THROW(not_implemented, "Wrong EW Scheme for Griffin. Use alpha0,alphamZ,or Gmu.")
  // }
  m_ewscheme = s["EW_SCHEME"].Get<ew_scheme::code>();
  // if(m_ewscheme==ew_scheme::alphamZ){
    m_griffin.set(0, Flavour(kf_Wplus).Mass());
    m_griffin.set(1, Flavour(kf_Z).Mass());
    m_griffin.set(2, Flavour(kf_h0).Mass());
    m_griffin.set(3, Flavour(kf_e).Mass());
    m_griffin.set(4, Flavour(kf_mu).Mass());
    m_griffin.set(5, Flavour(kf_tau).Mass());
    m_griffin.set(6, Flavour(kf_d).Mass());
    m_griffin.set(7, Flavour(kf_s).Mass());
    m_griffin.set(8, Flavour(kf_b).Mass()); // MSbar mass at scale mu=MZ for mb(mb)=4.20
    m_griffin.set(9, Flavour(kf_u).Mass());
    m_griffin.set(10, Flavour(kf_c).Mass()); // MSbar mass at scale mu=MZ
    m_griffin.set(11, Flavour(kf_t).Mass());
    if(m_ewscheme!=3) m_griffin.set(12, s_model->ScalarConstant("alpha_QED"));
    m_griffin.set(13, s_model->ScalarConstant("alpha_S"));
    m_griffin.set(Delal, delap);
    // m_griffin.set(15, delap);
    m_griffin.set(16, GF);
    m_griffin.set(17, Flavour(kf_Wplus).Width());
    m_griffin.set(18, Flavour(kf_Z).Width());
   //  m_griffin.set(MZ, Flavour(kf_Z).Mass());
   //  m_griffin.set(MW, Flavour(kf_Wplus).Mass());
   // m_griffin.set(al, 1/137.03599976);
   //  m_griffin.set(als, s_model->ScalarConstant("alpha_S"));
   //  m_griffin.set(GamZ, Flavour(kf_Z).Width());
   //  m_griffin.set(GamW, Flavour(kf_Wplus).Width());
   //  m_griffin.set(MH, Flavour(kf_h0).Mass());
   //  m_griffin.set(MT, Flavour(kf_t).Mass());
   //  m_griffin.set(MB, Flavour(kf_b).Mass());  // MSbar mass at scale mu=MZ
   //  m_griffin.set(Delal, delap);
   //  m_griffin.set(Gmu, GF);
  // }
  if(m_ewscheme==3){
    // m_griffinGMwMz = SMvalGMwMz(m_griffin);
    SMvalGMwMz myinput2(m_griffin);

    cout << endl;
    cout << "complex-pole mass: mw = " << myinput2.get(MWc) << endl;
    cout << "PDG mass:          mw = " << myinput2.get(MW) << endl;
    cout << "complex-pole mass: mZ = " << myinput2.get(MZc) << endl;
    cout << "PDG mass:          mZ = " << myinput2.get(MZ) << endl;
    cout << "alpha(0):          al = " << myinput2.get(al) << endl;
    cout << "                   1/al = " << 1/myinput2.get(al) << endl;
    cout << endl;
  }
  return true;
}

void Griffin::Griffin_Interface::RegisterProcess(const PHASIC::Process_Info &pi){
  // Since griffin is only 2->2 m_initial[0]=bar(m_initial[1])
  m_inital = pi.m_ii.GetExternal()[0];
  m_final = pi.m_fi.GetExternal()[0];
  if( m_inital == m_final ) m_norm = 2;  // For BhaBha
}

void Griffin::Griffin_Interface::RegisterProcess(const External_ME_Args& args){
  // Since griffin is only 2->2 m_initial[0]=bar(m_initial[1])
  m_inital = args.m_inflavs[0];
  m_final =  args.m_outflavs[0];
  if( m_inital == m_final ) m_norm = 2;  // For BhaBha
}

void Griffin::Griffin_Interface::EvaluateLoop(const Vec4D_Vector& momenta, METOOLS::DivArrD& virt){
  if(momenta.size()!=4){
    THROW(fatal_error, "Griffin library is for 2->2 scattering only");
  }
  double s = (momenta[0]+momenta[1]).Abs2();
  m_griffin.set(12, (*aqed)((s)));
  EvaluateLO(momenta, virt);
  if(m_order==griffinorder::nlo) EvaluateNLO(momenta, virt);
  else if(m_order==griffinorder::nloe) EvaluateNLOE(momenta, virt);
  else EvaluateNNLO(momenta,virt);
  virt.Finite() -= m_B;
}

void Griffin::Griffin_Interface::EvaluateLO(const Vec4D_Vector& momenta, DivArrD &res){
  double s = (momenta[0]+momenta[1]).Abs2();
  double t = (momenta[0]-momenta[2]).Abs2();
  double cost = 1.+2.*t/s;
  Cplx resvv, resva, resav, resaa;
  double sw=s_model->ComplexConstant("csin2_thetaW").real();
  // cost = cos((momenta[2]).Theta());
  if(cost > 1. || cost < -1){
    msg_Error()<<"CosTheta out of range in "<<METHOD<<endl;
  }
  FA_SMLO FAi(m_inital, m_griffin), FAf(m_final, m_griffin);
  SW_SMLO SWi(m_inital, m_griffin), SWf(m_final, m_griffin);

  matel M(m_inital, m_final, VEC, VEC, FAi.result().real(), FAf.result().real(), SWi.result().real(),  SWf.result().real(), s, cost, m_griffin);
        M.setkinvar(s, cost);

    M.setform(VEC, VEC);
    resvv = M.result();
    M.setform(AXV, VEC);
    resav = M.result();
    M.setform(VEC, AXV);
    resva = M.result();
    M.setform(AXV, AXV);
    resaa = M.result();
    m_B = real((1+cost*cost)*(resvv*conj(resvv) + resav*conj(resav)
                 + resva*conj(resva) + resaa*conj(resaa)) +
      + 4*cost*(resvv*conj(resaa) + resva*conj(resav)).real());
}

void Griffin::Griffin_Interface::EvaluateNLO(const Vec4D_Vector& momenta, DivArrD &res){
  double s = (momenta[0]+momenta[1]).Abs2();
  double t = (momenta[0]-momenta[2]).Abs2();
  double cost = 1.+2.*t/s;

  // if(m_ewscheme==3){
      FA_SMNLO FAi(m_inital, m_griffin), FAf(m_final, m_griffin);
      SW_SMNLO SWi(m_inital, m_griffin), SWf(m_final, m_griffin);
      mat_SMNLO M(m_inital, m_final, VEC, VEC, FAi.result().real(), FAf.result().real(),  SWi.result().real(),  SWf.result().real(), s, cost, m_griffin);

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
    // }
    // else{
    //   FA_SMNLO FAi(m_inital, m_griffinGMwMz), FAf(m_final, m_griffinGMwMz);
    //   SW_SMNLO SWi(m_inital, m_griffinGMwMz), SWf(m_final, m_griffinGMwMz);
    //   mat_SMNLO M(m_inital, m_final, VEC, VEC, FAi.result().real(), FAf.result().real(),  SWi.result().real(),  SWf.result().real(), s, cost, m_griffinGMwMz);

    //   M.setkinvar(s, cost);
    //   Cplx resvv, resva, resav, resaa;

    //   M.setform(VEC, VEC);
    //   resvv = M.result();
    //   M.setform(AXV, VEC);
    //   resav = M.result();
    //   M.setform(VEC, AXV);
    //   resva = M.result();
    //   M.setform(AXV, AXV);
    //   resaa = M.result();
    //   res.Finite() = real((1+cost*cost)*(resvv*conj(resvv) + resav*conj(resav)
    //                + resva*conj(resva) + resaa*conj(resaa)) +
    //     + 4*cost*(resvv*conj(resaa) + resva*conj(resav)).real());

    // }
}

void Griffin::Griffin_Interface::EvaluateNLOE(const Vec4D_Vector& momenta, DivArrD &res){
  double s = (momenta[0]+momenta[1]).Abs2();
  double t = (momenta[0]-momenta[2]).Abs2();
  double m1 = momenta[0].Mass();
  double m2 = momenta[2].Mass();
  double cost = 1+2.*t/s;
  FA_SMLO FAi0(m_inital, m_griffin), FAf0(m_final, m_griffin);
  SW_SMLO SWi0(m_inital, m_griffin), SWf0(m_final, m_griffin);
  FA_SMNLO FAi(m_inital, m_griffin), FAf(m_final, m_griffin);
  SW_SMNLO SWi(m_inital, m_griffin), SWf(m_final, m_griffin);
  double sw=s_model->ComplexConstant("csin2_thetaW").real();
  matel M0(m_inital, m_final, VEC, VEC, FAi0.result().real(), FAf0.result().real(), sw, sw, s, cost, m_griffin);
  mat_SMeNLO M(m_inital, m_final, VEC, VEC, FAi.result().real(), FAf.result().real(), sw, sw, s, cost, m_griffin);

  M0.setkinvar(s, cost);
  M.setkinvar(s, cost);
  Cplx resvv, resva, resav, resaa;
  Cplx res0vv, res0va, res0av, res0aa;
  M0.setform(VEC, VEC);
  res0vv = M0.result();
  M0.setform(AXV, VEC);
  res0av = M0.result();
  M0.setform(VEC, AXV);
  res0va = M0.result();
  M0.setform(AXV, AXV);
  res0aa = M0.result();
  M.setform(VEC, VEC);
  resvv = M.result();
  M.setform(AXV, VEC);
  resav = M.result();
  M.setform(VEC, AXV);
  resva = M.result();
  M.setform(AXV, AXV);
  resaa = M.result();
  res.Finite() = real((1+cost*cost)*
                ((res0vv+2*(resvv-res0vv))*conj(res0vv) 
                 + (res0av+2*(resav-res0av))*conj(res0av)
                 + (res0va+2*(resva-res0va))*conj(res0va) 
                 + (res0aa+2*(resaa-res0aa))*conj(res0aa)) +
                + 4*cost*(res0vv*conj(res0aa) + (resvv-res0vv)*conj(res0aa)
                          + res0vv*conj(resaa-res0aa) 
                        + res0va*conj(res0av) + (resva-res0va)*conj(res0av)
                          + res0va*conj(resav-res0av)));
    // -2*(1+cost*cost)*(resvv*conj(resav) + resva*conj(resaa)).real()
    // -4*cost*(resvv*conj(resva) + resav*conj(resaa)));
  // return res*3*s*s*m_rescale_alpha/32/M_PI;
  // res.Finite() -= m_B;
  // res.Finite() *= 3*s*sqrt(s);
  // res.Finite() *= 2;
}

void Griffin::Griffin_Interface::EvaluateNNLO(const Vec4D_Vector& momenta, DivArrD &res){
  double s = (momenta[0]+momenta[1]).Abs2();
  double t = (momenta[0]-momenta[2]).Abs2();
  double cost = 1.+2.*t/s;
  // double cost = 0;
  // if(sqrt(s)<40) {
  //   res.Finite() = 0;
  // }
  // if(m_ewscheme==3){
     // FA_SMNNLO FAi(m_inital, m_griffinGmu), FAf(m_final, m_griffinGmu);
     // SW_SMNNLO SWi(m_inital, m_griffinGmu), SWf(m_final, m_griffinGmu);
     // mat_SMNNLO M(m_inital, m_final, VEC, VEC,  FAi.result().real(), FAf.result().real(),  SWi.result().real(),  SWf.result().real(), s, cost, m_griffinGmu);

     //  M.setkinvar(s, cost);
     //  Cplx resvv, resva, resav, resaa;

     //  M.setform(VEC, VEC);
     //  resvv = M.result();
     //  M.setform(AXV, VEC);
     //  resav = M.result();
     //  M.setform(VEC, AXV);
     //  resva = M.result();
     //  M.setform(AXV, AXV);
     //  resaa = M.result();
     //  res.Finite() = real((1+cost*cost)*(resvv*conj(resvv) + resav*conj(resav)
     //               + resva*conj(resva) + resaa*conj(resaa)) +
     //    + 4*cost*(resvv*conj(resaa) + resva*conj(resav)).real());
  // }
  // else{
    FA_SMNNLO FAi(m_inital, m_griffin), FAf(m_final, m_griffin);
    SW_SMNNLO SWi(m_inital, m_griffin), SWf(m_final, m_griffin);
    // PRINT_VAR( SWi.result());
    // PRINT_VAR( SWf.result());
     double sw=s_model->ComplexConstant("csin2_thetaW").real();
     // mat_SMNNLO M(m_inital, m_final, VEC, VEC,  FAi.result().real(), FAf.result().real(),  SWi.result().real(),  SWf.result().real(), s, cost, m_griffin);
     mat_SMNNLO M(m_inital, m_final, VEC, VEC,  FAi.result().real(), FAf.result().real(),  sw, sw, s, cost, m_griffin);

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
  // }
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
  if      (tag == "NNLO")    mode = griffinorder::nnlo;
  else if (tag == "NLOE")    mode = griffinorder::nloe;
  else if (tag == "NLO")     mode = griffinorder::nlo;
  else if (tag == "LO")      mode = griffinorder::lo;
  else THROW(fatal_error, "Unknown GRIFFIN: Order = " + tag);
  return str;
}

std::ostream &Griffin::operator<<(std::ostream &str,const griffinorder::code &ym)
{
  if      (ym==griffinorder::lo)     return str<<"LO";
  else if (ym==griffinorder::nlo)    return str<<"NLO";
  else if (ym==griffinorder::nnlo)   return str<<"NNLO";
  else if (ym==griffinorder::nloe)   return str<<"Exact NLO";
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
