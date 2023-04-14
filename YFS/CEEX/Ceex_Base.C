#include "YFS/CEEX/Ceex_Base.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Random.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "EXTAMP/External_ME_Interface.H"
#include "PHASIC++/Process/External_ME_Args.H"

using namespace YFS;


Amplitude::Amplitude() {
  for (int h0 = 0; h0 <= 1; ++h0) {
    for (int h1 = 0; h1 <= 1; ++h1) {
      for (int h2 = 0; h2 <= 1; ++h2) {
        for (int h3 = 0; h3 <= 1; ++h3) {
          m_A[h0][h1][h2][h3] = Complex(0, 0);
        }
      }
    }
  }
}


Ceex_Base::Ceex_Base(const Flavour_Vector &flavs)
{
  RegisterDefaults();
  Scoped_Settings s{ Settings::GetMainSettings()["CEEX"] };
  Settings& ss = Settings::GetMainSettings();
  m_onlyz = s["ONLYZ"].Get<int>();
  m_onlyg = s["ONLYG"].Get<int>();
  m_checkxs = s["CHECK_XS"].Get<int>();
  string widthscheme = ss["WIDTH_SCHEME"].Get<string>();
  m_fixedwidth = (widthscheme == "Fixed");
  m_flavs = flavs;
  if (flavs.size() != 4) {
    THROW(fatal_error, "CEEX is only for 2->2");
  }
  if (flavs[2].IsNeutrino() && flavs[3].IsNeutrino()) {
    m_onlyz = true;
  }

  m_Q1Q2I = flavs[0].Charge() * flavs[1].Charge();
  m_QIQF  = flavs[0].Charge() * flavs[2].Charge();
  m_Q1Q2F = flavs[2].Charge() * flavs[3].Charge();
  m_MZ = Flavour(kf_Z).Mass();
  m_gZ = Flavour(kf_Z).Width();
  double mw = Flavour(kf_Wplus).Mass();
  double MH = Flavour(kf_h0).Mass();
  double  GH  = Flavour(kf_h0).Width();
  double  GW  = Flavour(kf_Wplus).Width();
  double  GZ  = Flavour(kf_Z).Width();
  m_I   = Complex(0., 1.);

  p_xyz = new XYZFunc(flavs, {});
  double F_L = 0.;
  double F_R = 0.;

  // m_alpha  = MODEL::s_model->ScalarConstant("alpha_QED");
  // m_alpha = (*aqed)(0);
  // m_cW = sqrt(std::abs(cw2));
  m_sin2tw = MODEL::s_model->ComplexConstant("csin2_thetaW");
  // for (int i = 0; i < MODEL::s_model->Vertices().size(); ++i)
  // {
  //   cout<<MODEL::s_model->Vertices()[i]<<endl;
  // }
  m_sW = m_sin2tw;
  m_e = sqrt(4.*M_PI * m_alpha);
  m_cW = 1. - m_sW;
  m_norm = sqrt(16. * m_sW * (1. - m_sW));
  m_qe       = m_flavs[0].Charge();
  m_qf       = m_flavs[2].Charge();
  m_Q1Q2I = m_flavs[0].Charge() * m_flavs[1].Charge();
  m_ae       = 2.*m_flavs[0].IsoWeak();
  m_af       = 2.*m_flavs[2].IsoWeak();
  m_ve       = (m_ae - 4.*m_qe * m_sin2tw) / m_norm;
  m_vf       = (m_af - 4.*m_qf * m_sin2tw) / m_norm;
  m_ae /= m_norm;
  m_af /= m_norm;
  m_mass_I = flavs[0].Mass();
  m_mass_F = flavs[2].Mass();
  // full EW couplings
  m_I_L = -m_I * sqrt(4 * M_PI * m_alpha) / (2.*m_sW * m_cW) * (2.*flavs[0].IsoWeak()
          - 2.*flavs[0].Charge() * m_sW * m_sW);

  m_I_R = -m_I * sqrt(4 * M_PI * m_alpha) / (2.*m_sW * m_cW) * (-2.*flavs[0].Charge() * m_sW * m_sW);

  m_F_L = -m_I * sqrt(4 * M_PI * m_alpha) / (2.*m_sW * m_cW) * (2.*flavs[2].IsoWeak()
          - 2.*flavs[2].Charge() * m_sW * m_sW);

  m_F_R = -m_I * sqrt(4 * M_PI * m_alpha) / (2.*m_sW * m_cW) * (-2.*flavs[2].Charge() * m_sW * m_sW);
  m_cL = -m_I * sqrt(4 * M_PI * m_alpha) / (2.*m_sW * m_cW) * (2.*m_flavs[1].IsoWeak()
         - 2.*m_flavs[1].Charge() * m_sW * m_sW) / m_norm;
  m_cR = -m_I * sqrt(4 * M_PI * m_alpha) / (2.*m_sW * m_cW) * (-2.*m_flavs[1].Charge() * m_sW * m_sW) / m_norm;
  m_zeta = {1, 1, 0, 0};
  m_eta = {0, 0, 1, 0};
  m_b  = {0.0,  0.8723e0, -0.7683e0, 0.3348e0};
}

void Ceex_Base::RegisterDefaults()
{
  Scoped_Settings s{ Settings::GetMainSettings()["CEEX"] };
  s["ONLYZ"].SetDefault(0);
  s["ONLYG"].SetDefault(0);
  s["CHECK_XS"].SetDefault(0);
}


void Ceex_Base::Init(const Vec4D_Vector &p)
{
  m_momenta = p;
  // m_momenta[0] = {45.6,0,0,45.5999999971368};
  // m_momenta[1] = {45.6,0,0,-45.5999999971368};
  // m_momenta[2] = {44.8916194574982,-39.5820913415418,-1.83384213624145,-21.0983743897405};
  // m_momenta[3] = {43.8356556626458,39.5793416016701,1.84299500194237,18.751870080043};
  // // m_momenta[0] = {45.5,0,0,-45.499999997130551};
  // m_momenta[1] = {45.5,0,0,45.499999997130551};
  // m_momenta[2] = {45.5, 37.33205, 13.2138, 22.404277};
  // m_momenta[3] = {45.5, -37.33205, -13.2138, -22.404277};

  // m_momenta[2] ={45.499134036866245,-36.445336068203119,-11.005958010401864,24.915182705605929};
  // m_momenta[3] = {45.499725780646209,36.445299766383002,11.005944820779755,-24.916322233717317};
  Poincare cms = Poincare(m_momenta[0] + m_momenta[1]);
  Poincare Rot = Poincare(Vec4D(0., 0., 0., 1.));
  for (size_t i(0); i < p.size(); ++i) {
    cms.Boost(m_momenta[i]);
    // cms.Boost(m_bornmomenta[i]);
  }
  m_crude = 2.0 / (4.0 * M_PI);
  for (size_t i(0); i < m_isrphotons.size(); ++i) {
    m_crude /= pow(2 * M_PI, 3);
  }
  // m_momenta[0] = {125, 0., 0.0, 124.99999999895552};
  // m_momenta[1] = {125, 0., 0.0, -124.99999999895552};
  // m_momenta[2] = {72.065607801457745,-23.712397199988942,38.623097764918803,56.030519628725635};
  // m_momenta[3] = {69.046456665460752,23.649546194782079,-38.407278737009577,52.277928869207031};
  // // // m_momenta[2] = {101.24991567373615, -4.4008060185383222,4.1892855613621593,101.06738831498141};
  // m_momenta[3] = { 19.985103367247920 ,4.399561991,-4.1974654430772027,-19.037286012953118 };
  // m_momenta[3] = {45.499740238089927,36.445324799172084,11.005952439021524,-24.916308653701236};
  // m_momenta[2] ={45.499179692335453,-36.445359187080896,-11.005964933257927,24.915229204002642};
  // m_norm  =  2.0/(4.0*M_PI);
  // m_s  = (m_momenta[0] + m_momenta[1]).Abs2();
  m_sp = (m_momenta[2] + m_momenta[3]).Abs2();
  m_sQ = (m_momenta[2] + m_momenta[3]).Abs2();
  // m_sp = sqr(90.4535);
  // m_sp = 8181.8441939847635;
  m_Sc = Complex(m_s, 0);
  m_Tc = Complex((p[0] - p[3]).Abs2(), 0);
  m_Uc = Complex((p[0] - p[2]).Abs2(), 0);
  m_T = 0;
  // m_Sini[0].clear();
  // m_Sini[1].clear();
  // MakeProp();
}

void Ceex_Base::LoadME() {
  Coupling_Map cpls;
  MODEL::s_model->GetCouplings(cpls);

  m_pi.m_mincpl[1] = m_pi.m_mincpl[1] + 1;
  m_pi.m_maxcpl[1] = m_pi.m_maxcpl[1] + 1;
  m_pi.m_fi.m_nlotype = ATOOLS::nlo_type::loop;
  // p_loop_me = PHASIC::Virtual_ME2_Base::GetME2(m_pi);
  m_pi.m_fi.m_nlotype = ATOOLS::nlo_type::loop;
  const string tag = "Recola_Born";
  PHASIC::External_ME_Args args(m_pi.m_ii.GetExternal(),
                                m_pi.m_fi.GetExternal(),
                                m_pi.m_maxcpl);

  p_lo_me = PHASIC::Tree_ME2_Base::GetME2(tag, args);
  if (!p_lo_me)
    THROW(not_implemented, "Couldn't find virtual ME for this process.");
  p_lo_me->SetCouplings(cpls);
  // p_lo_me->SetSubType(ATOOLS::sbt::qed);


}


void Ceex_Base::MakeProp()
{
  if (m_fixedwidth) {
    m_propZ =  1. / Complex(m_sp - sqr(m_MZ), m_gZ * m_MZ);
  }
  else {
    m_propZ =   1. / Complex(m_sp - sqr(m_MZ), m_gZ * m_sp / m_MZ);
  }
  m_propG = 1. / Complex(m_sp,0);
  if (m_onlyz)  m_propG = 0;
  if (m_onlyg)  m_propZ = 0;
  m_prop =  m_propZ + m_propG;
}

Complex Ceex_Base::CouplingZ(double  j, int mode) {
  if (m_onlyg) return 0.;
  // m_zcpl = 0.;
  Complex zcpl;
  // if(mode==1) j=-j;
  if (mode == 1) {
    //mode 1 T
    zcpl = m_ve * m_vf - dcmplx(j) * m_ae * m_vf + dcmplx(j) * m_ve * m_af - m_ae * m_af;
    // zcpl = (m_ve-j*m_ae)*(m_vf+j*m_af);
  }
  else if (mode == 0) {
    zcpl = m_ve * m_vf - dcmplx(j) * m_ae * m_vf - dcmplx(j) * m_ve * m_af + m_af * m_ae;
    // zcpl = (m_ve-j*m_ae)*(m_vf-j*m_af);
  }
  else msg_Error() << METHOD << "\n wrong mode\n";

  if (zcpl == 0.) {
    msg_Error() << "Z coupling is Zero!\n";
  }
    // zcpl = m_ve * m_vf - dcmplx(j) * m_ae * m_vf + dcmplx(j) * m_ve * m_af - m_af * m_ae;
    // zcpl = m_ve * m_vf + dcmplx(j) * m_ae * m_vf +dcmplx(j) * m_ve * m_af + m_ae * m_af;

  // return m_zcpl/sqrt(2);
  // return m_zcpl*m_norm/sqrt(2);
  return zcpl;
}


Complex Ceex_Base::CouplingG() {
  // z-> (Ve -Hel1*Ae)*(Vf +Hel1*Af)
  // if (m_onlyz) return 0;
  m_gcpl = Complex(m_Q1Q2I, 0);
  return m_gcpl;
}


Complex Ceex_Base::T(const Vec4D &p1, const Vec4D &p2, int h1, int h2) {
  // return S(p1,p2,p1.Mass(),p2.Mass(),h1,h2);
  // h1 = -h1;
  // h2 = -h2;
  // // if(p1.Mass()==0) h1=-h1;
  // // if(p2.Mass()==0) h2=-h2;
  Complex s(0, 0);
  double sq1 = Xi(p1, p2); //sqrt(m_zeta*p1/(m_zeta*p2));
  double sq2 = Xi(p2, p1); //sqrt(m_zeta*p2/(m_zeta*p1));
  if (h1 == -h2) {
    if (h1 == 1) s = Splus(p1, p2);
    else if (h1 == -1) s = Sminus(p1, p2);
    return s;
  }
  else if (h1 == h2) {
    // s =  p1.Mass()*sq2+p2.Mass()*sq1;
    s =  p1.Mass() * sq2;
  }
  else {
    msg_Error() << METHOD << "Wrong helicities\n";
  }
  return s;
}


Complex Ceex_Base::Tp(const Vec4D &p1, const Vec4D &p2, int h1, int h2) {
  // return S(p1,p2,-p1.Mass(),-p2.Mass(),h1,h2);
  h1 = -h1;
  h2 = -h2;
  // // if(p1.Mass()==0) h1=-h1;
  // // if(p2.Mass()==0) h2=-h2;
  Complex s(0, 0);
  double sq1 = Xi(p1, p2); //sqrt(m_zeta*p1/(m_zeta*p2));
  double sq2 = Xi(p2, p1); //sqrt(m_zeta*p2/(m_zeta*p1));
  if (h1 == -h2) {
    if (h1 == 1) s = Splus(p1, p2);
    else s = Sminus(p1, p2);
    return s;
  }
  else if (h1 == h2) {
    // if(h1==-1) s = -p2.Mass()*sq1-p1.Mass()*sq2;
    s = -p2.Mass() * sq1;
  }
  else {
    msg_Error() << METHOD << "Wrong helicities\n";
  }
  return s;
}

Complex Ceex_Base::U(const Vec4D &p1, const Vec4D &p2, int h1, int h2) {
  // return S(p1,p2,-p1.Mass(),-p2.Mass(),h1,h2);
  // return S(p1,p2,h1,h2);
  // h1 = -h1;
  // h2 = -h2;
  // if(p1.Mass()==0) h1=-h1;
  // if(p2.Mass()==0) h2=-h2;
  Complex s(0, 0);
  double sq1 = sqrt(m_zeta * p1 / (m_zeta * p2));
  double sq2 = sqrt(m_zeta * p2 / (m_zeta * p1));
  if (h1 == -h2) {
    if (h1 == -1) s = Splus(p1, p2);
    else s = Sminus(p1, p2);
  }
  else if (h1 == h2) {
    // s =  p2.Mass()*sq1+p1.Mass()*sq2;
    // s =  p2.Mass()*sq1;
    s =  -p2.Mass() * sq1;
  }
  else {
    msg_Error() << METHOD << "Wrong helicities\n";
  }
  return s;
}

Complex Ceex_Base::Up(const Vec4D &p1, const Vec4D &p2, int h1, int h2) {
  // return S(p1,p2,-p1.Mass(),-p1.Mass(), h1,h2);
  // h2=-h2;
  // return S(p1,p2,h1,h2);
  // h1 = -h1;
  // if(p1.Mass()==0) h1=-h1;
  // if(p2.Mass()==0) h2=-h2;
  Complex s(0, 0);
  double sq1 = sqrt(m_zeta * p1 / (m_zeta * p2));
  double sq2 = sqrt(m_zeta * p2 / (m_zeta * p1));
  if (h1 == -h2) {
    if (h1 == 1) s = Splus(p1, p2);
    else s = Sminus(p1, p2);
    return s;
  }
  else if (h1 == h2) {
    // s =  -p1.Mass()*sq2-p2.Mass()*sq1;
    s =  p1.Mass() * sq2;
    // s =  -p2.Mass()*sq1;
  }
  else {
    msg_Error() << METHOD << "Wrong helicities\n";
  }
  return s;
}

Complex Ceex_Base::S(const Vec4D &p1, const Vec4D &p2, int h1, int h2) {
  Complex s(0, 0);
  // if(p1.Mass()==0) h1=-h1;
  // if(p2.Mass()==0) h2=-h2;
  double sq1 = sqrt(m_zeta * p2 / (m_zeta * p1));
  double sq2 = sqrt(m_zeta * p1 / (m_zeta * p2));
  if (h1 == -h2) {
    if (h1 > 0) s = Splus(p1, p2);
    else s = Sminus(p1, p2);
  }
  else if (h1 == h2) {
    s = p1.Mass() * sq1 + p2.Mass() * sq2;
  }
  else {
    msg_Error() << METHOD << "Wrong helicities\n";
  }
  return s;
}

Complex Ceex_Base::S(const Vec4D &p1, const Vec4D &p2, double m1, double m2, int h1, int h2) {
  Complex s(0, 0);
  Vec4D phat1, phat2;
  // phat1 = p1-m_zeta*m1*m1/(2*m_zeta*p1);
  // phat2 = p2-m_zeta*m2*m2/(2*m_zeta*p2);
  if (IsEqual(m1, 0)) h1 = -h1;
  if (IsEqual(m2, 0)) h2 = -h2;
  double sq1 = sqrt(m_zeta * p2 / (m_zeta * p1));
  double sq2 = sqrt(m_zeta * p1 / (m_zeta * p2));
  if (h1 == -h2) {
    if (h1 > 0) s = Splus(p1, p2);
    else s = Sminus(p1, p2);
  }
  else if (h1 == h2) {
    s = m1 * sq1 + m2 * sq2;
  }
  else {
    msg_Error() << METHOD << "Wrong helicities\n";
  }
  return s;
}


Complex Ceex_Base::Splus(const Vec4D &p, const Vec4D &q) {
  Complex sp = -Complex(q[2], q[3]) * sqrt((p[0] - p[1]) / (q[0] - q[1]));
  sp += Complex(p[2], p[3]) * sqrt((q[0] - q[1]) / (p[0] - p[1]));
  return sp;
}

Complex Ceex_Base::Sminus(const Vec4D &p, const Vec4D &q) {
  return -conj(Splus(p, q));
  // Complex sp = Complex(q[2], -q[3]) * sqrt((p[0] - p[1]) / (q[0] - q[1]));
  // sp -= Complex(p[2], -p[3]) * sqrt((q[0] - q[1]) / (p[0] - p[1]));
  // return sp;
}


Complex Ceex_Base::BoxGG() {
  Complex MG = 1e-30;
  Complex t1 = log(m_Tc / m_Uc) * (log(MG * MG / m_Sc) + m_I * M_PI);
  Complex t2 = 0.5 * m_I * m_Sc * (m_Uc - m_Tc) / sqr(m_Uc) * (0.5 * m_I * sqr(log(-m_Tc / m_Sc))
               + m_I * M_PI * log(-m_Tc / m_Sc));
  Complex t3 = -0.5 * m_I * m_Sc / m_Uc * (log(-m_Tc / m_Sc) + m_I * M_PI);
  return t1 + t2 + t3;
}

Complex Ceex_Base::BoxGZ() {
  // Box Gamma-Z, From W. Brown, R. Decker, E. Pashos, Phys. Rev. Lett., 52 (1984)
  Complex MG = 1e-30;
  Complex mb2 = Complex(sqr(m_MZ), -m_MZ * m_gZ);
  Complex t1 = log(m_Tc / m_Uc) * (log(MG * MG / sqrt(m_Tc * m_Uc)));
  Complex t2 = -2.*log(m_Tc / m_Uc) * log((mb2 - m_Sc) / mb2)
               + DiLog((mb2 + m_Uc) / mb2) - DiLog((mb2 + m_Tc) / mb2);

  Complex t3 = (mb2 - m_Sc) * (m_Uc - m_Tc - mb2) / (sqr(m_Uc)) * (
                 log(-m_Tc / m_Sc) * log((mb2 - m_Sc) / mb2)
                 + DiLog((mb2 + m_Tc) / mb2) - DiLog((mb2 - m_Sc) / mb2));
  Complex t4 = sqr(mb2 - m_Sc) / (m_Sc * m_Uc) * log((mb2 - m_Sc) / mb2)
               + (mb2 - m_Sc) / m_Uc * log(-m_Tc / mb2);

  return t1 + t2 + t3 + t4;

}


Complex Ceex_Base::BoxSubtract() {
  Complex MG = 1e-30;
  Complex sub = log(m_Tc / m_Uc) * log(MG * MG / sqrt(m_Uc * m_Tc))
                + 0.5 * log(m_Tc / m_Uc);
  return sub;
}

Complex Ceex_Base::Soft(Vec4D k, Vec4D p1, Vec4D p2, int hel) {
  Complex coeff1 = sqrt(2) * sqrt(m_zeta * p1 / (m_zeta * k));
  Complex coeff2 = sqrt(2) * sqrt(m_zeta * p2 / (m_zeta * k));
  Complex soft;
  Vec4D phat1 = p1 - m_zeta * sqr(p1.Mass()) / (2 * m_zeta * p1);
  Vec4D phat2 = p2 - m_zeta * sqr(p2.Mass()) / (2 * m_zeta * p2);
  if (hel == 1)  soft = coeff1 * Splus(k, phat1) / (2 * k * phat1) - coeff2 * Splus(k, phat2) / (2 * k * phat1);
  else if (hel == -1)  soft = -coeff1 * Sminus(k, phat1) / (2 * k * phat1) + coeff2 * Sminus(k, phat2) / (2 * k * phat2);
  else msg_Error() << METHOD << "\n Wrong photon hel\n";
  return soft;
}

Complex Ceex_Base::bsigma(Vec4D p1, Vec4D p2, int hel) {
  // eq 230 hep-ph/0006359
  Complex coeff1 = sqrt(2) * Xi(p2, p1); // sqrt(m_zeta * p2 / (m_zeta * p1));
  Vec4D phat1 = p1 - m_zeta * sqr(p1.Mass()) / (2 * m_zeta * p1);
  Vec4D phat2 = p2 - m_zeta * sqr(p2.Mass()) / (2 * m_zeta * p2);
  if (hel < 0) return coeff1 * Sminus(p1, p2);
  else return coeff1 * Splus(p1, p2);
}

Complex Ceex_Base::Sfactor(const Vec4D &p1, const Vec4D &p2, const Vec4D &k, int hel) {
  Complex s, b1, b2;
  if (Sminus(k, p1) != -conj(Splus(k, p1))) {
    msg_Error() << "Wrong soft factors in " << METHOD << std::endl;
  }
  if (hel == -1) {
    b1 = sqrt(2) * Xi(p1, k) * Sminus(k, p1);
    b2 = sqrt(2) * Xi(p2, k) * Sminus(k, p2);
  }
  else {
    b1 = sqrt(2) * Xi(p1, k) * Splus(k, p1);
    b2 = sqrt(2) * Xi(p2, k) * Splus(k, p2);
  }
  m_b1 = Splus(k, p1);
  m_b2 = Splus(k, p2);
  m_Soft = -0.5 * (b1 / (p1 * k) - b2 / (p2 * k));
  s = -0.5 * (b1 / (p1 * k) - b2 / (p2 * k));
  m_pp1 = p1;
  m_pp2 = p2;
  m_kk = k;
  // Complex s = -0.5*m_e*(bsigma(k,p2,hel)/(k*p2)-bsigma(k,p1,hel)/(k*p1));
  // S = -0.5*sqrt(2)*m_e*Xi(k,p1);
  // dcmplx bf1 =  Sqr2* iProd1(sigma,ph,p1)*XiProd(p1,ph); //!!! =GPS_bfact(sigma,ph,p1);
  // dcmplx bf2 =  Sqr2* iProd1(sigma,ph,p2)*XiProd(p2,ph); //!!! =GPS_bfact(sigma,ph,p2);
  // Soft = -bf1/(2*pk1) +bf2/(2*pk2);



  // the amplitude level eikoanl when squared agrees
  // with the full eikonal up to a constant factor (2\pi)^3
  // double eik = Eikonal(k,p1,p2);
  // double SS = (s*conj(s)).real();
  // double ratio = eik/SS;
  // if(!IsEqual(ratio,1/pow(2*M_PI,3),1e-3)){
  //   msg_Error()<<METHOD<<"Photon Eikonal is wrong in CEEX\n"
  //             <<"Ratio should be "<<1/pow(2*M_PI,3)<<endl
  //             <<"Ratio is "<<ratio<<endl
  //             <<"p1 = "<<p1<<endl
  //             <<"p2 = "<<p2<<endl
  //             <<"k = "<<k<<endl
  //             <<"Photon helicity = "<<hel<<endl;
  // }
  // if(hel==-1) return conj(m_e*s);
  return -(m_e * s);
  // return -sqrt(4 * M_PI * m_alpha) * (-bsigma(k, p1, hel)/(2.*k*p1)+bsigma(k,p2,hel)/(2.*k*p2));
}

void Ceex_Base::CalculateSfactors() {
  m_Sprod = Complex(1, 0);
  for (int i = 0; i < m_isrphotons.size(); ++i)
  {
    // Vec4D k = m_isrphotons[i];
    // m_Sini[0].push_back(Complex(1,0)*Sfactor(m_bornmomenta[0], m_bornmomenta[1], m_isrphotons[i], m_PhoHel[i]));
    // m_Sini[1].push_back(-conj(m_Sini[0][i]));
    // m_Sprod *=  Sfactor(m_bornmomenta[0], m_bornmomenta[1], m_isrphotons[i], m_PhoHel[i]);
    m_Sprod *=  Sfactor(m_bornmomenta[0], m_bornmomenta[1], m_isrphotons[i], m_PhoHel[i]);
    // if(m_PhoHel[i] == 1 ) m_Sprod *=  Sfactor(m_bornmomenta[0], m_bornmomenta[1], m_isrphotons[i], 1);
    // else m_Sprod *= -conj( Sfactor(m_bornmomenta[0], m_bornmomenta[1], m_isrphotons[i], 1));
  }
  // m_Sprod*=-m_e;
}


void Ceex_Base::MakePhotonHel() {
  m_PhoHel.clear();
  for (int i = 0; i < m_isrphotons.size(); ++i)
  {
    // m_norm /= pow(2.*M_PI, 3);
    if (ran->Get() < 0.5) {
      m_PhoHel.push_back(1);
    }
    else m_PhoHel.push_back(-1);
  }
}


void Ceex_Base::UGamma(const Vec4D &p1, const Vec4D &p2, const Vec4D &k, int sigma, Amplitude &AmpU) {

  AmpU.m_U[0][0] = UGamma(p1, p2, k, sigma, 1, 1);
  AmpU.m_U[0][1] = UGamma(p1, p2, k, sigma, 1, -1);
  AmpU.m_U[1][0] = UGamma(p1, p2, k, sigma, -1, 1);
  AmpU.m_U[1][1] = UGamma(p1, p2, k, sigma, -1, -1);

}

void Ceex_Base::VGamma(const Vec4D &p1, const Vec4D &p2, const Vec4D &k, int sigma, Amplitude &AmpV) {

  AmpV.m_U[0][0] = VGamma(p1, p2, k, sigma, 1, 1);
  AmpV.m_U[0][1] = VGamma(p1, p2, k, sigma, 1, -1);
  AmpV.m_U[1][0] = VGamma(p2, p1, k, sigma, -1, 1);
  AmpV.m_U[1][1] = VGamma(p1, p2, k, sigma, -1, -1);

}


Complex Ceex_Base::UGamma(const Vec4D &p1, const Vec4D &p2, const Vec4D &k, int h1, int i, int j, bool negMass) {
  // eq 222 hep-ph/0006359
  // h1 photon helicity
  // m_Ugamma = CMatrix(2); // ubar eps u
  // h1=1 U^+, h1=-1 U^-
  // matrix indexing  | (+1+1) (+1-1) |
  // i,j= (.,.)       | (-1+1) (-1-1) |
  if (h1 == -1) return -conj(UGamma(p2, p1, k, 1, j, i));
  double sqr2 = sqrt(2);
  double m1 = p1.Mass();
  double m2 = p2.Mass();
  // if(abs(m1) < 0.1*Flavour(11).Mass()) m1 = 1e-60;
  // if(abs(m2) < 0.1*Flavour(11).Mass()) m2 = 1e-60;
  if (negMass) {
    m1 = -m1;
    m2 = -m2;
  }
  Vec4D p1hat = p1 - m_zeta * m1 * m1 / (2 * m_zeta * p1);
  Vec4D p2hat = p2 - m_zeta * m2 * m2 / (2 * m_zeta * p2);
  // m2 = 1e-60;
  // if (m1 < 1e-5 && m1 != Flavour(11).Mass()) m1 = 0.;
  // if (m2 < 1e-5 && m2 != Flavour(11).Mass()) m2 = 0.;
  Complex coeff1 = sqrt(2) * Xi(p1, k); // sqrt(m_zeta * p1 / (m_zeta * k));
  Complex coeff2 = sqrt(2) * Xi(p2, k); // sqrt(m_zeta * p2 / (m_zeta * k));
  if ( i == 1 && j == 1) return sqr2 * Xi(p2, k) * Splus(k, p1hat);
  else if (i == 1 && j == -1) return 0;
  else if (i == -1 && j == 1) return sqr2 * (m2 * Xi(p1, p2) - m1 * Xi(p2, p1));
  else if (i == -1 && j == -1) return sqr2 * Xi(p1, k) * Splus(k, p2hat);
  else {
    msg_Out() << "h1 = " << i << "\n h2 = " << j;
    THROW(fatal_error, "Wrong helicities in CEEX");
  }
}

Complex Ceex_Base::VGamma(const Vec4D &p1, const Vec4D &p2, const Vec4D &k, int h1, int i, int j) {
  // eq 222 hep-ph/0006359
  // h1 photon helicity
  // m_Vgamma = CMatrix(2);
  // ugamm with m-> -m and hel -> -hel
  // h1 = 1;
  // h1=1 V^+, h1=-1 V^-
  // matrix indexing  | (+1+1) (+1-1) |
  // i,j= (.,.)       | (-1+1) (-1-1) |
  return UGamma(p1, p2, k, h1, -i, -j, false);
}

Complex Ceex_Base::iProd(const int L, const Vec4D &p, const Vec4D &q) {
  Complex  Prod;
  if (     L ==  1 ) {
    Prod =
      -sqrt( (p[0] - p[1]) / (q[0] - q[1]) ) * Complex(q[2], q[3])
      + sqrt( (q[0] - q[1]) / (p[0] - p[1]) ) * Complex(p[2], p[3]);
  } else  if ( L == -1 ) {
    Prod =
      -sqrt( (q[0] - q[1]) / (p[0] - p[1]) ) * Complex(p[2], p[3])
      + sqrt( (p[0] - p[1]) / (q[0] - q[1]) ) * Complex(q[2], q[3]);
    Prod = conj(Prod);
  }   else {
    cout << "##### KKceex::iProd1: Wrong L= " << L << endl;
  }
  return Prod;
}


void Ceex_Base::Calculate() {
  m_beta10  = 0.0;
  m_beta01  = 0.0;
  m_beta00  = 0.0;
  MakePhotonHel();
  MakeProp();
  // m_isrphotons.clear();
  // // Vec4D k = {108.59830848948251,6.2820803938209269e-2,-0.21601945041252191,-108.59807547061062};
  // Vec4D k = {2.40913477274814,-0.00116820730070409,0.0011352582029645,2.40913422202764};
  // m_isrphotons.push_back(k);
  // m_isrphotons.push_back({0.063590107,0.0039179472,-0.010288124,-0.062629912});
  // k = {23.367440069979207,-5.2352790624640613e-4,8.3505392052715729e-3,23.367438572049004 };
  // m_isrphotons.push_back(k);
  CalculateSfactors();
  m_cfac = m_Sprod;//   *(m_s/m_sp);
  InfraredSubtractedME_0_0();
  // InfraredSubtractedME_0_1();
  for (int i = 0; i < m_isrphotons.size(); ++i) {
    // InfraredSubtractedME_1_0(m_isrphotons[i], m_PhoHel[i]);
  }
  // InfraredSubtractedME_2_0();
  // m_beta10 /= pow(8*M_PI*M_PI*M_PI,0.5*(3*m_isrphotons.size()+1.));// pow(2.*M_PI,0.5*(3*m_isrphotons.size()+1));
  double BornCru = 4. / 3.*m_born;
  double flux = (m_momenta[2] + m_momenta[3]).Abs2() / (m_momenta[0] + m_momenta[1]).Abs2();
  Complex amp = m_beta00 + m_beta01 + m_beta10 + m_beta20;
  double amp2 = (amp * conj(amp)).real();
  // double amp2 = ((m_beta01)*conj((m_beta01))).real();
  // m_result = 0.5*(amp2).real();//*m_sp/m_s;
  // m_result *= (m_mome)
  // m_result = m_born*(1.+(amp2).real());//*m_s/m_sp;
  // m_result *= 0.5;/*m_s/m_sp;//*m_crude;
  m_result = amp2/2;
  // PRINT_VAR(m_sp);
  // // PRINT_VAR(m_s);
  // PRINT_VAR(m_born);
  // PRINT_VAR(m_born/m_result);
}


Complex Ceex_Base::BornAmplitude(const Vec4D_Vector &k) {
  Complex amp;
  double hel1, hel2, hel3, hel4;
  int mode;
  for (int h0 = 1; h0 <= 2; ++h0) {
    for (int h1 = 1; h1 <= 2; ++h1) {
      for (int h2 = 1; h2 <= 2; ++h2) {
        for (int h3 = 1; h3 <= 2; ++h3) {
          hel1 = 3 - 2 * h0;
          hel2 = 3 - 2 * h1;
          hel3 = 3 - 2 * h2;
          hel4 = 3 - 2 * h3;
          if (hel1 == -hel2 ) {
            m_T = T(k[2], k[0], hel3, hel1) * Tp(k[1], k[3], hel2, hel4);
            m_U = Up(k[2], k[1], hel3, hel2) * U(k[0], k[3], hel1, hel4);
            m_ampborn[h0][h1][h2][h3] = (CouplingZ(hel1, 1) + CouplingG()) * m_T + (CouplingZ(hel1, 1) + CouplingG()) * m_U;
            m_bornAmp.m_A[h0][h1][h2][h3] = (CouplingZ(hel1, 1) + CouplingG()) * m_T + (CouplingZ(hel1, 1) + CouplingG()) * m_U;
            amp += (CouplingZ(hel1, 0) * m_propZ + CouplingG() * m_propG) * m_U + (CouplingZ(hel1, 1) * m_propZ + CouplingG() * m_propG) * m_T;
          }
        }
      }
    }
  }
  return amp;
}


int Ceex_Base::MapHel(int &h) {
  if (h == 0) return -1;
  else if (h == 1) return 1;
  else msg_Error() << METHOD << "Wrong helicity!\n";
  return -10;
}

void Ceex_Base::BornAmplitude(const Vec4D_Vector &k, Amplitude &M) {
  Complex amp;
  int hel1, hel2, hel3, hel4;
  for (int h0 = 0; h0 <= 1; ++h0) {
    for (int h1 = 0; h1 <= 1; ++h1) {
      for (int h2 = 0; h2 <= 1; ++h2) {
        for (int h3 = 0; h3 <= 1; ++h3) {
          hel1 = 1 - 2 * h0;
          hel2 = 1 - 2 * h1;
          hel3 = 1 - 2 * h2;
          hel4 = 1 - 2 * h3;
          if (hel1 == -hel2 ) {
            m_T  = T(k[2], k[0], hel3, hel1) * Tp(k[1], k[3], hel2, hel4);
            m_U  = Up(k[2], k[1], hel3, hel2) * U(k[0], k[3], hel1,  hel4);
            m_Tamp[h0][h1][h2][h3] = m_T;
            m_Uamp[h0][h1][h2][h3] = m_U;
          }
          else {
            m_Tamp[h0][h1][h2][h3] = 0.;
            m_Uamp[h0][h1][h2][h3] = 0.;

          }
        }
      }
    }
  }
  for (int j = 0; j <= 1; j++) {
    double h = 1. - 2.*j;
    m_UC[j] = CouplingZ(h, 0) * m_propZ + CouplingG() * m_propG;
    m_TC[j] = CouplingZ(h, 1) * m_propZ + CouplingG() * m_propG;
  }
  for (int h0 = 0; h0 <= 1; h0++) {
    for (int h1 = 0; h1 <= 1; h1++) {
      for (int h2 = 0; h2 <= 1; h2++) {
        for (int h3 = 0; h3 <= 1; h3++) {
          M.m_A[h0][h1][h2][h3] = m_TC[h0] * m_Tamp[h0][h1][h2][h3] + m_UC[h0] * m_Uamp[h0][h1][h2][h3];
          m_bornAmp.m_A[h0][h1][h2][h3] = M.m_A[h0][h1][h2][h3];
        }
      }
    }
  }
}






Complex Ceex_Base::BornAmplitude(const Vec4D_Vector &k, int h0, int h1, int h2, int h3) {
  Complex amp;
  if (h0 == -h1 ) {
    m_T = T(k[2], k[0], h2, h0) * Tp(k[1], k[3], h1, h3);
    m_U = Up(k[2], k[1], h2, h1) * U(k[0], k[3], h0, h3);
    amp = (CouplingZ(h0, 0) * m_propZ + CouplingG() * m_propG) * m_U + (CouplingZ(h0, 1) * m_propZ + CouplingG() * m_propG) * m_T;
  }
  return amp;
}


Complex Ceex_Base::BornAmplitude(Vec4D p1, Vec4D p2, Vec4D p3, Vec4D p4, int h0, int h1, int h2, int h3)
{
  Vec4D_Vector tmp;
  tmp.push_back(p1);
  tmp.push_back(p2);
  tmp.push_back(p3);
  tmp.push_back(p4);
  return BornAmplitude(tmp, h0, h1, h2, h3);
}

void Ceex_Base::InfraredSubtractedME_0_0() {
  m_T *= 0;
  m_U *= 0;
  m_bornsum *= 0.;
  double hel1, hel2, hel3, hel4;
  m_beta00 = 0;
  // m_Sprod= Complex(1,0);
  Vec4D Px = m_momenta[0] + m_momenta[1];
  double crude = m_born / (4 * M_PI);
  double norm = 1;
  // m_crude = (m_Sprod * conj(m_Sprod)).real();
  // PRINT_VAR(crude);
  int mode;
  double svarx = Px.Abs2();
  Amplitude AmpBorn;
  BornAmplitude(m_momenta, AmpBorn);
  // BornAmplitude(m_momenta, m_bornAmp);
  // Complex amp = BornAmplitude(m_momenta);
  // double amp2 = (amp * conj(amp)).real();
  SumAmplitude(m_beta00, m_bornAmp, m_e*m_e);

  // m_beta00 =  (4 * M_PI * m_alpha) * amp;
}

void Ceex_Base::InfraredSubtractedME_0_1() {
  // boxes first, only if FSR resummed
  Complex SubBox = BoxSubtract();
  Complex coef = m_alpha * m_qe * m_qf / M_PI;
  m_BoxGGtu = coef * BoxGG() - SubBox;
  m_BoxGZtu = coef * BoxGZ() - SubBox;
  // tmp swap u and t for crossing
  Complex u1, t1;
  t1 = m_Tc;
  u1 = m_Uc;
  m_Uc = t1;
  m_Tc = u1;
  m_BoxGGut = coef * (-BoxGG()) - SubBox;
  m_BoxGZut = coef * (-BoxGZ()) - SubBox;
  // swap back
  m_Uc = u1;
  m_Tc = t1;
  //
  // vertex corrections
  // using factorised form
  // PRINT_VAR(m_I * M_PI - 1.);
  Complex LE = log(m_Sc / sqr(m_mass_I)) - Complex(1., M_PI);
  Complex LF = log(m_Sc / sqr(m_mass_F)) - Complex(1., M_PI);
  Complex deltI = 0.5 * sqr(m_qe) * m_alpha / Complex(M_PI, 0) * LE;
  Complex deltF = 0.5 * sqr(m_qf) * m_alpha / Complex(M_PI, 0) * LF;
  deltF = 1.;
  m_T *= 0;
  m_U *= 0;
  m_vertexI = deltI;
  Complex ampbox, boxgg, boxgz, ampvertex;
  int hel1, hel2, hel3, hel4;
  for (int h0 = -1; h0 <= 1; h0 += 2) {
    for (int h1 = -1; h1 <= 1; h1 += 2) {
      for (int h2 = -1; h2 <= 1; h2 += 2) {
        for (int h3 = -1; h3 <= 1; h3 += 2) {
          if (h0 == -h1) {
            if (h0 == h2) {
              boxgg = m_BoxGGtu;
              boxgz = m_BoxGZtu;
            }
            else if (h0 == -h2) {
              boxgg = m_BoxGGut;
              boxgz = m_BoxGZut;
            }
            ampbox += T(m_momenta[2], m_momenta[0], h2, h0) * Tp(m_momenta[1], m_momenta[3], h1, h3) * CouplingG() * m_propG * boxgg
                      + T(m_momenta[2], m_momenta[0], h2, h0) * Tp(m_momenta[1], m_momenta[3], h1, h3) * CouplingZ(h0, 1) * m_propZ * boxgz
                      + Up(m_momenta[2], m_momenta[1], h2, h1) * U(m_momenta[0], m_momenta[3], h0, h3) * CouplingG() * m_propG * boxgg
                      + Up(m_momenta[2], m_momenta[1], h2, h1) * U(m_momenta[0], m_momenta[3], h0, h3) * CouplingZ(h0, 0) * m_propZ * boxgz;

            // not factored as future EW-corr wont
            m_T = T(m_momenta[2], m_momenta[0], h2, h0) * Tp(m_momenta[1], m_momenta[3], h1, h3);
            m_U = Up(m_momenta[2], m_momenta[1], h2, h1) * U(m_momenta[0], m_momenta[3], h0, h3);
          }
        }
      }
    }
  }
  for (int h0 = 0; h0 <= 1; ++h0) {
    for (int h1 = 0; h1 <= 1; ++h1) {
      for (int h2 = 0; h2 <= 1; ++h2) {
        for (int h3 = 0; h3 <= 1; ++h3) {
          ampvertex += m_TC[h0] * m_Tamp[h0][h1][h2][h3] + m_UC[h0] * m_Uamp[h0][h1][h2][h3];
        }
      }
    }
  }
  Amplitude ampew;
  Complex amp;

  BornAmplitude(m_momenta, ampew);
  SumAmplitude(amp, ampew, (deltI) * (deltF));
  double amp2 = (amp * conj(amp)).real();
  double flux = m_sp / (m_momenta[2] + m_momenta[3]).Abs2();
  m_beta01 = (4 * M_PI * m_alpha) * amp;
}

void Ceex_Base::InfraredSubtractedME_1_0(Vec4D &k, int helk) {
  Vec4D_Vector arm1, arm2, p;
  Vec4D p1, p2, p3, p4;
  p = m_momenta;
  // helk = -1;
  arm1 = {k, m_bornmomenta[1], p[2], p[3]};
  arm2 = {m_bornmomenta[0], k, p[2], p[3]};
  Complex amp1ISR(0, 0), sum1(0, 0), sum2(0, 0);
  double p1k = k * m_bornmomenta[0];
  double p2k = k * m_bornmomenta[1];
  Amplitude AmpBornU, AmpBornV, AmpV, AmpU;
  BornAmplitude(arm1, AmpBornU);
  BornAmplitude(arm2, AmpBornV);
  UGamma(m_bornmomenta[0], k, k, helk, AmpU);
  VGamma(k, m_bornmomenta[1], k, helk, AmpV);

  AddU(m_beta10, AmpBornU, AmpU, m_e / p1k / 2.);
  AddV(m_beta10, AmpBornV, AmpV, -m_e / p2k / 2.);

  m_Amp1U = AmpU;
  m_Amp1V = AmpV;

  m_beta10 /= pow(2 * M_PI, 3);

  m_beta10 /= Sfactor(m_bornmomenta[0], m_bornmomenta[1], k, helk); //m_Sprod;
  // m_beta10*=sqrt((m_bornmomenta[2]+m_bornmomenta[3]+k).Abs2());
  // m_beta10/=sqrt((m_bornmomenta[2]+m_bornmomenta[3]).Abs2());
}


void Ceex_Base::InfraredSubtractedME_2_0()
{
  int hel1, hel2, hel3, hel4;
  Complex b2double(0.0);
  Complex b2single(0.0);
  Complex b2rest(0.0);
  Complex b20(0.0);
  for (int i = 0; i < m_isrphotons.size(); ++i) {
    for (int j = 0; j < i; ++j) {
      Vec4D k1 = m_isrphotons[i];
      Vec4D k2 = m_isrphotons[j];
      b2double = BetaDouble_2_0(k1, k2, m_PhoHel[i], m_PhoHel[j]);
      b2rest   = BetaRest_2_0(k1, k2, m_PhoHel[i], m_PhoHel[j]);
      b2single = BetaSingle_2_0(k1, k2, m_PhoHel[i], m_PhoHel[j]);
      b2single += BetaSingle_2_0(k2, k1, m_PhoHel[j], m_PhoHel[i]);
      // b2rest = BetaRest_2_0(k1,k2,m_PhoHel[i], m_PhoHel[j]);
      b20 += b2double + b2single + b2rest;
      // b20 /= pow(2*M_PI,3)*sqrt(2);
      b20 /= pow(2 * M_PI, 6);
      b20 /= Sfactor(m_bornmomenta[0], m_bornmomenta[1], k1, m_PhoHel[i]);
      b20 /= Sfactor(m_bornmomenta[0], m_bornmomenta[1], k2, m_PhoHel[j]);

    }
  }
  m_beta20 = (b20);
}


Complex Ceex_Base::BetaDouble_2_0(Vec4D &k1, Vec4D &k2, int h1, int h2)
{
  Vec4D_Vector p;
  Vec4D p1, p2, p3, p4;
  p = m_bornmomenta;
  Complex amp1ISR(0, 0), sum1(0, 0), sum2(0, 0);
  double m1 = m_momenta[0].Mass();
  double m2 = m_momenta[1].Mass();
  double m3 = m_momenta[2].Mass();
  double m4 = m_momenta[3].Mass();

  double denA = 2.*p[0] * (k1 + k2) - 2 * k1 * k2;
  double denB = 2.*p[1] * (k1 + k2) - 2 * k1 * k2;

  double deltA = 2.*k1 * k2 / denA;
  double deltB = 2.*k1 * k2 / denB;
  k1 = m_isrphotons[0];
  k2 = m_isrphotons[1];
  Complex s1a = -m_qe * bsigma(k1, p[0], h1) / (2 * k1 * p[0]);
  Complex s2a = -m_qe * bsigma(k2, p[0], h2) / (2 * k2 * p[0]);
  Complex s1b = m_qe * bsigma(k1, p[1], h1) / (2 * k1 * p[1]);
  Complex s2b = m_qe * bsigma(k2, p[1], h2) / (2 * k2 * p[1]);
  return 4 * M_PI * m_alpha * (s1a * s2a * deltA + s1b * s2b * deltB) * BornAmplitude(p);
}


Complex Ceex_Base::BetaSingle_2_0(Vec4D &k1, Vec4D &k2, int hk1, int hk2) {
  int hel1, hel2, hel3, hel4;
  // r_if = 2k_i*p_f
  // r_ij = 2k_i*k_j
  Complex Single;
  Vec4D p1, p2, p3, p4;
  Vec4D_Vector p;
  Amplitude AmpU1, AmpU2, AmpV1, AmpV2;
  Amplitude AmpBornU1, AmpBornU2, AmpBornV1, AmpBornV2;
  p = m_bornmomenta;
  Complex amp1ISR(0, 0), sum1(0, 0), sum2(0, 0);
  double m1 = m_momenta[0].Mass();
  double m2 = m_momenta[1].Mass();
  double m3 = m_momenta[2].Mass();
  double m4 = m_momenta[3].Mass();
  Vec4D Q;
  Q = m_bornmomenta[0] + m_bornmomenta[1] - k1 - k2;
  double r11 = 2 * k1 * k1;
  double r12 = 2 * k1 * k2;

  double r1a = 2 * k1 * p[0];
  double r1b = 2 * k1 * p[1];

  double r2a = 2 * k2 * p[0];
  double r2b = 2 * k2 * p[1];
  Complex s1a = -m_qe * bsigma(k1, p[0], hk1) / (2 * k1 * p[0]);
  Complex s1b = m_qe * bsigma(k1, p[1], hk1) / (2 * k1 * p[1]);

  Complex s2a = -m_qe * bsigma(k2, p[0], hk2) / (2 * k2 * p[0]);
  Complex s2b = m_qe * bsigma(k2, p[1], hk2) / (2 * k2 * p[1]);
  // BornAmplitude(arm1,AmpBornU);
  // BornAmplitude(arm2, AmpBornV);
  // UGamma(m_bornmomenta[0],k, k, helk,AmpU);
  // VGamma(k, m_bornmomenta[1], k, helk,AmpV);

  // AddU(m_beta10, AmpBornU, AmpU, m_e/p1k/2.);
  // AddV(m_beta10, AmpBornV, AmpV, -m_e/p2k/2.);
  Vec4D_Vector armU1 = {k1, m_bornmomenta[1], p[2], p[3]};
  Vec4D_Vector armU2 = {k2, m_bornmomenta[1], p[2], p[3]};
  Vec4D_Vector armV1 = {m_bornmomenta[0], k1, p[2], p[3]};
  Vec4D_Vector armV2 = {m_bornmomenta[0], k2, p[2], p[3]};

  BornAmplitude(armU2, AmpBornU2);
  UGamma(m_bornmomenta[0], k2, k1, hk2, AmpU2);
  Complex f1 = -s2a / (r12 - r1a - r2a);
  AddU(Single, AmpBornU2, AmpU2, f1);

  BornAmplitude(armV1, AmpBornV1);
  VGamma(k1, k2, m_bornmomenta[1], hk1, AmpV1);
  f1 = s1b / (r12 - r1a - r2a);
  AddV(Single, AmpBornV1, AmpV1, f1);

  BornAmplitude(armU1, AmpBornU1);
  UGamma(m_bornmomenta[0], k1, k1, hk1, AmpU1);
  f1 = -s2a * (1 / (r12 - r1a - r2a) - 1. / (-r2b));
  AddU(Single, AmpBornU1, AmpU1);



  BornAmplitude(armV2, AmpBornV2);
  VGamma(k2, k2, m_bornmomenta[1], hk1, AmpU1);
  f1 = s1b * (1 / (r12 - r1a - r2a) - 1. / (-r2b));
  AddU(Single, AmpBornU1, AmpU1);



  // BornAmplitude(arm2, AmpBornV);
  // hk1 = 3-2*hk1;
  // hk2 = 3-2*hk2;

  return Single;
}

Complex Ceex_Base::BetaRest_2_0(Vec4D &k1, Vec4D &k2, int hk1, int hk2) {
  int hel1, hel2, hel3, hel4;
  Vec4D p1, p2, p3, p4;
  Vec4D_Vector p;
  p = m_bornmomenta;
  Complex amp1ISR(0, 0), sum1(0, 0), sum2(0, 0);

  hk1 = 0.5 * (3 - hk1);
  hk2 = 0.5 * (3 - hk2);
  double r11 = 2 * k1 * k1;
  double r12 = 2 * k1 * k2;

  double r1a = 2 * k1 * p[0];
  double r1b = 2 * k1 * p[1];

  double r2a = 2 * k2 * p[0];
  double r2b = 2 * k2 * p[1];
  Complex s1a = -m_qe * bsigma(k1, p[0], hk1) / (2 * k1 * p[0]);
  Complex s2a = -m_qe * bsigma(k2, p[0], hk2) / (2 * k2 * p[0]);

  Complex s1b = m_qe * bsigma(k1, p[1], hk1) / (2 * k1 * p[1]);
  Complex s2b = m_qe * bsigma(k2, p[1], hk2) / (2 * k2 * p[1]);

  Complex betarest(0.0);
  Complex  t1, t2, t3, t4;

  Complex Amp, Single;
  Amplitude AmpU2, AmpBornU2;
  Amp = BornAmplitude(m_bornmomenta);
  UGamma(m_bornmomenta[0], k2, k1, hk2, AmpU2);
  Complex f1 = -s2a / (r12 - r1a - r2a);
  AddU(Single, AmpBornU2, AmpU2, f1);


  return betarest;///(Sfactor(m_bornmomenta[0], m_bornmomenta[1], k1, hk1)*Sfactor(m_bornmomenta[0], m_bornmomenta[1], k2, hk2));
}


Complex Ceex_Base::T_mass(const Vec4D &p1, const Vec4D &p2,  double m1, double m2, int h1, int h2) {
  // h1 = -h1;
  // h2 = -h2;
  // if(p1.Mass()==0) h1=-h1;
  // if(p2.Mass()==0) h2=-h2;
  Complex s(0, 0);
  double sq1 = sqrt(m_zeta * p2 / (m_zeta * p1));
  double sq2 = sqrt(m_zeta * p1 / (m_zeta * p2));
  if (h1 == -h2) {
    if (h1 == 1) s = Splus(p1, p2);
    else if (h1 == -1) s = Sminus(p1, p2);
    return s;
  }
  else if (h1 == h2) {
    s =  m1 * sq1; //+p2.Mass()*sq2;;
  }
  else {
    msg_Error() << METHOD << "Wrong helicities\n";
  }
  return s;
}


Complex Ceex_Base::Tp_mass(const Vec4D &p1, const Vec4D &p2,  double m1, double m2, int h1, int h2) {
  // return S(p1,p2,h1,h2);
  h1 = -h1;
  h2 = -h2;
  // if(p1.Mass()==0) h1=-h1;
  // if(p2.Mass()==0) h2=-h2;
  Complex s(0, 0);
  double sq1 = sqrt(m_zeta * p2 / (m_zeta * p1));
  double sq2 = sqrt(m_zeta * p1 / (m_zeta * p2));
  if (h1 == -h2) {
    if (h1 > 0) s = Splus(p1, p2);
    else s = conj(Splus(p2, p1));
    return s;
  }
  else if (h1 == h2) {
    s = -m2 * sq2;
  }
  else {
    msg_Error() << "Wrong helicities\n";
  }
  return s;
}

Complex Ceex_Base::U_mass(const Vec4D &p1, const Vec4D &p2, double m1, double m2, int h1, int h2) {
  // return S(p1,p2,h1,h2);
  h1 = -h1;
  h2 = -h2;
  // if(p1.Mass()==0) h1=-h1;
  // if(p2.Mass()==0) h2=-h2;
  Complex s(0, 0);
  double sq1 = sqrt(m_zeta * p2 / (m_zeta * p1));
  double sq2 = sqrt(m_zeta * p1 / (m_zeta * p2));
  if (h1 == -h2) {
    if (h1 > 0) s = Splus(p1, p2);
    else s = Sminus(p1, p2);
  }
  else if (h1 == h2) {
    s =  -m2 * sq2;; //+p2.Mass()*sq2;
  }
  else {
    msg_Error() << "Wrong helicities\n";
  }
  return s;
}

Complex Ceex_Base::Up_mass(const Vec4D &p1, const Vec4D &p2, double m1, double m2, int h1, int h2) {
  // return S(p1,p2,h1,h2);
  // h1 = -h1;
  // if(p1.Mass()==0) h1=-h1;
  // if(p2.Mass()==0) h2=-h2;
  Complex s(0, 0);
  double sq1 = sqrt(m_zeta * p2 / (m_zeta * p1));
  double sq2 = sqrt(m_zeta * p1 / (m_zeta * p2));
  if (h1 == -h2) {
    if (h1 > 0) s = Splus(p1, p2);
    else s = Sminus(p1, p2);
    return s;
  }
  else if (h1 == h2) {
    s =  m1 * sq1; //+p2.Mass()*sq2;
  }
  else {
    msg_Error() << METHOD << "Wrong helicities\n";
  }
  return s;
}


Complex Ceex_Base::BornAmplitude_mass(const Vec4D_Vector &k, double m1, double m2, double m3, double m4) {
  Complex amp;
  double hel1, hel2, hel3, hel4;
  int mode;
  for (int h0 = 1; h0 <= 2; ++h0) {
    for (int h1 = 1; h1 <= 2; ++h1) {
      for (int h2 = 1; h2 <= 2; ++h2) {
        for (int h3 = 1; h3 <= 2; ++h3) {
          hel1 = 3 - 2 * h0;
          hel2 = 3 - 2 * h1;
          hel3 = 3 - 2 * h2;
          hel4 = 3 - 2 * h3;
          if (hel1 == -hel2 ) {
            // msg_Out()<<"h0,h1,h2,h3 "<<h0<<h1<<h2<<h3<<"\n";
            m_T = T_mass(k[2], k[0], m3, m1, hel3, hel1) * Tp_mass(k[1], k[3], m2, m4, hel2, hel4);
            m_U = Up_mass(k[2], k[1], m3, m2, hel3, hel2) * U_mass(k[0], k[3], m1, m4, hel1, hel4);
            // m_Tamp[h0][h1][h2][h3] = m_T;
            // m_Uamp[h0][h1][h2][h3] = m_U;
            // amp += (CouplingZ(hel4, 1) * m_propZ + CouplingG() * m_propG) * m_U + (CouplingZ(hel4, 1) * m_propZ + CouplingG() * m_propG) * m_T;
            amp += (CouplingZ(hel1, 0) * m_propZ + CouplingG() * m_propG) * m_U + (CouplingZ(hel1, 1) * m_propZ + CouplingG() * m_propG) * m_T;
            // PRINT_VAR(m_T);
            // PRINT_VAR(m_U);
            // PRINT_VAR(CouplingZ(hel1, 1) * m_propZ + CouplingG() * m_propG);
            // PRINT_VAR(CouplingZ(hel1, 0) * m_propZ + CouplingG() * m_propG);
            // PRINT_VAR(CouplingZ(hel1, 0));
            // PRINT_VAR(CouplingZ(hel1, 1));
            // PRINT_VAR(1./m_propG);
            // PRINT_VAR(1./m_propZ);
          }
        }
      }
    }
  }
  // PRINT_VAR(m_ve);
  // PRINT_VAR(m_vf);
  // PRINT_VAR(m_ae);
  // PRINT_VAR(m_af);
  // PRINT_VAR(m_momenta);
  // exit(1);
  return amp;
}

Complex Ceex_Base::BornAmplitude_mass(const Vec4D_Vector &k, double m1, double m2, double m3, double m4, int h0, int h1, int h2, int h3) {
  Complex amp;
  if (h0 == -h1 ) {
    // msg_Out()<<"h0,h1,h2,h3 "<<h0<<h1<<h2<<h3<<"\n";

    m_T = T_mass(k[2], k[0], m3, m1, h2, h0) * Tp_mass(k[1], k[3], m2, m4, h1, h3);
    m_U = -Up_mass(k[2], k[1], m3, m2, h2, h1) * U_mass(k[0], k[3], m1, m4, h0, h3);
    amp = (CouplingZ(h0, 1) * m_propZ + CouplingG() * m_propG) * m_U + (CouplingZ(h0, 0) * m_propZ + CouplingG() * m_propG) * m_T;
    // PRINT_VAR(m_T);
    // PRINT_VAR(m_U);
    // PRINT_VAR(CouplingZ(h0, 1) * m_propZ + CouplingG() * m_propG);
    // PRINT_VAR(CouplingZ(h0, 0) * m_propZ + CouplingG() * m_propG);
    // PRINT_VAR(CouplingZ(h0, 0));
    // PRINT_VAR(CouplingZ(h0, 1));
    // PRINT_VAR(1./m_propG);
    // PRINT_VAR(1./m_propZ);
  }
  // p_lo_me->Calc(m_momenta);
  // p_lo_me->CalcAmp(m_momenta,hel);
  return amp;
}



void Ceex_Base::SumAmplitude(Complex &sum, const Amplitude &Amp, const Complex fac) {
  for (int h0 = 0; h0 <= 1; ++h0) {
    for (int h1 = 0; h1 <= 1; ++h1) {
      for (int h2 = 0; h2 <= 1; ++h2) {
        for (int h3 = 0; h3 <= 1; ++h3) {
          sum += fac * Amp.m_A[h0][h1][h2][h3];
        }
      }
    }
  }
}


void Ceex_Base::SumAmplitude(Complex &sum, const Amplitude &Amp1, const Amplitude &Amp2, const Complex fac) {
  for (int h0 = 0; h0 <= 1; ++h0) {
    for (int h1 = 0; h1 <= 1; ++h1) {
      for (int h2 = 0; h2 <= 1; ++h2) {
        for (int h3 = 0; h3 <= 1; ++h3) {
          for (int  j = 0; j <= 1; j++) {
            sum += fac * Amp2.m_A[j][h1][h2][h3] * Amp1.m_U[j][h1];
          }
        }
      }
    }
  }
}

double Ceex_Base::Eikonal(Vec4D k, Vec4D p1, Vec4D p2) {
  return -m_alpha / (4 * M_PI * M_PI) * (p1 / (p1 * k) - p2 / (p2 * k)).Abs2();
}


void Ceex_Base::AddU(Complex &sum, const Amplitude &Born, const Amplitude &U, const Complex fac) {
  for (int h0 = 0; h0 <= 1; ++h0) {
    for (int h1 = 0; h1 <= 1; ++h1) {
      for (int h2 = 0; h2 <= 1; ++h2) {
        for (int h3 = 0; h3 <= 1; ++h3) {
          Complex CSum(0, 0);
          for (int  j = 0; j <= 1; j++) {
            CSum += fac * Born.m_A[j][h1][h2][h3] * U.m_U[j][h1];
          }
          sum += CSum;
        }
      }
    }
  }
}


void Ceex_Base::AddV(Complex &sum, const Amplitude &Born, const Amplitude &V, const Complex fac) {
  for (int h0 = 0; h0 <= 1; ++h0) {
    for (int h1 = 0; h1 <= 1; ++h1) {
      for (int h2 = 0; h2 <= 1; ++h2) {
        for (int h3 = 0; h3 <= 1; ++h3) {
          Complex CSum(0, 0);
          for (int  j = 0; j <= 1; j++) {
            CSum += fac * (V.m_V[h2][j]) * Born.m_A[h0][j][h2][h3];
          }
          sum += CSum;
        }
      }
    }
  }
}

double Ceex_Base::Xi(const Vec4D p, const Vec4D q) {
  return sqrt((m_zeta * p) / (q * m_zeta));
}

void Ceex_Base::Reset() {
  m_result = 0;
}
