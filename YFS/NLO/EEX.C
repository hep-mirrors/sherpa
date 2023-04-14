#include "YFS/NLO/EEX.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "EXTRA_XS/Main/ME2_Base.H"
#include "PHASIC++/Process/Tree_ME2_Base.H"
#include "PHASIC++/Process/External_ME_Args.H"

using namespace YFS;
using namespace ATOOLS;

double smax = 1000;
int xbins = 50;
int ybins = 50;
int ymin = 0;
double ymax = 20;
Real_ff::Real_ff(int order):
  m_order(order)
{
  Scoped_Settings s{ Settings::GetMainSettings()["YFS"] };
  m_fsrmode  = s["FSR"].Get<int>();
  m_realFSR  = m_realISR = 0.0;
  s["VMIN1"].SetDefault(1e-9);
  s["VMIN2"].SetDefault(1e-9);
  s["VMIN3"].SetDefault(1e-9);
  s["VIRT_NLL"].SetDefault(1);
  s["USE_FAC"].SetDefault(1);
  s["USE_NNLO_QED"].SetDefault(0);
  m_vlim1 = s["VMIN1"].Get<double>();
  m_vlim2 = s["VMIN2"].Get<double>();
  m_vlim3 = s["VMIN2"].Get<double>();
  m_usefulleik = s["BETA_FULL"].SetDefault(1).Get<int>();
  m_formfactor = s["FULL_FORM"].Get<int>();
  m_usenllv = s["VIRT_NLL"].Get<int>();
  m_use_fac = s["USE_FAC"].Get<bool>();
  m_use_nnlo = s["USE_NNLO_QED"].Get<bool>();
#ifdef YFS_DEBUG_REAL
  m_histograms[string("beta00")]  = new Histogram(0, ymin, ymax, ybins);
  m_histograms[string("beta10")]  = new Histogram(0, ymin, ymax, ybins);
  m_histograms[string("beta02")]  = new Histogram(0, ymin, ymax, ybins);
  m_histograms[string("beta01")]  = new Histogram(0, ymin, 100, ybins);
  m_histograms[string("beta11")]  = new Histogram(0, ymin, 100, ybins);

#endif

}

Real_ff::~Real_ff()
{
#ifdef YFS_DEBUG_REAL
  Histogram * histo;
  string name;
  PRINT_VAR(m_histograms.size());
  for (map<string, Histogram *>::iterator hit = m_histograms.begin();
       hit != m_histograms.end(); hit++) {
    histo = hit->second;
    name  = string("./YFS-hist/") + hit->first + string(".dat");
    PRINT_VAR(name);
    // histo->MPISync();
    histo->Finalize() ;
    histo->Output(name);
    delete histo;
  }
#endif

}

void Real_ff::CalculateVirt() {
  // if (m_use_model_alpha) m_alpha = s_model->ScalarConstant("alpha_QED");
  // else m_alpha  = (*aqed)(0); // cant be in constructor
  m_alpi = m_alpha / M_PI;
  m_q1 = m_beam1;
  m_q2 = m_beam2;
//  if(m_fsrmode!=0){
//    m_q1 = D.m_newmomenta[0];
//    m_q2 = D.m_newmomenta[1];

//  }
//  else{
//    m_q1 = p[2];
//    m_q2 = p[3];
// }

  // double betaI1 = CalculateBeta(p[0]);
  // double betaI2 = CalculateBeta(p[1]);
  // double betaF1 = CalculateBeta(m_q1);
  // double betaF2 = CalculateBeta(m_q2);

  // double tI = (1.+betaI1*betaI2)/(betaI1+betaI2);
  // double logargI =  (1.+betaI1)*(1.+betaI2)/((1.-betaI1)*(1.-betaI2));
  // logargI = (p[0]+p[1]).Abs2()/sqr(p[0].Mass());
  // double tF = (1.+betaF1*betaF2)/(betaF1+betaF2);
  // double logargF = (m_q1+m_q2).Abs2()/sqr(m_q2.Mass());
  // m_gammaI  = 2*m_alpi*(log(logargI)-1.); // See Marek's phd thesis A.2.1
  // m_gammaF  = 2*m_alpi*(log(logargF)-1.);
  // m_delI2 = sqr(m_alpi*(log(logargI)))/2.;
  // m_delF2 = sqr(m_alpi*(log(logargF)))/2.;
  // m_delI2 = sqr(m_alpi*log(logargI))/2;
  // if(m_fsrmode==0) m_gammaF = m_delF2 = 0;
  // if(m_fsrmode==2) m_gammaI = m_delI2 = 0;
  // PRINT_VAR(m_gammaI);
  m_beta00 = m_born;//everything is divided by born
  m_beta01f = m_born * (1. + m_gammaF / 2.);
  m_beta01i = m_born * (1. + m_gammaI / 2.);
  m_beta02f = m_born * (1. + m_gammaF / 2. + m_delF2);
  m_beta02i =  m_born * (1. + m_gammaI / 2. + m_delI2);
  m_beta01 = m_born * (1. + m_gammaI / 2.) * (1. + m_gammaF / 2.);
  if (m_gammaF < 0) {
    PRINT_VAR(m_gammaF / 2.);
    PRINT_VAR(sqrt((m_q1 + m_q2).Abs2()));
    PRINT_VAR(m_q1.Mass());
    PRINT_VAR(m_q1);
    abort();
  }
  m_beta02 = m_born * (1. + m_gammaI / 2. + m_delI2) * (1. + m_gammaF / 2. + m_delF2);
  // PRINT_VAR(m_beta02/m_born);
  m_beta03 = m_born * (1. + m_gammaI / 2. + m_delI2 + pow(m_gammaI, 3.) / 48) * (1. + m_gammaF / 2. + m_delF2 + pow(m_gammaF, 3) / 48.);
  m_setvirt = true;
  if (m_use_nnlo) GetBeta02();
}


void Real_ff::GetBeta02() {
  // in limit mZ >> ml (Berends, Neerven, Burgers 1988 "Radiative Corrections at LEP energies", eq. 2.20/2.21/2.22)
  double Z3 = 1.202056903159594; // Zeta(3)
  double oldbeta02 = m_beta02;
  // m_L = m_L-1;
  // m_beta02 = m_beta00 + sqr(m_alpha/M_PI)*(1./2.*sqr(m_L)+(sqr(M_PI)/4.)*m_L);;//+13./4.+sqr(M_PI)*(17./16.-log(2.))-2*pow(M_PI,4.)/15.-9./2.*Z3)*m_beta00;
  m_beta02 = m_beta01 + sqr(m_alpha / M_PI) * (1. / 2.*sqr(m_L) + (-13. / 16. - sqr(M_PI) / 4. + 3.*Z3) * m_L + 13. / 4. + sqr(M_PI) * (17. / 16. - log(2.)) - 2 * pow(M_PI, 4.) / 15. - 9. / 2.*Z3) * m_beta00;
  // m_beta02 = m_beta00 + m_gamma/2*m_beta00 + sqr(m_alpha/M_PI)*(1./2.*sqr(m_L));//+(-13./16.-sqr(M_PI)/4.+3.*Z3)*m_L+13./4.+sqr(M_PI)*(17./16.-log(2.))-2*pow(M_PI,4.)/15.-9./2.*Z3)*m_beta00;
  m_beta30 = m_beta02 + pow(m_gammaI, 3.) / 48;
}


double Real_ff::AddVirtual() {
  if(m_looptool){
    m_beta03 = m_born * (1. + m_delI2 + pow(m_gammaI, 3.) / 48) * (1. + m_delF2 + pow(m_gammaF, 3) / 48.);
    return m_beta03;
  }
  switch (m_order) {
  case 0:
    return m_beta00;
    break;
  case 1:
    return m_beta01;
    break;
  case 2:
    return m_beta02;
    break;
  case 3:
    return m_beta03;
    break;
  default:
    return 0;
  }
}


double Real_ff::AddVirtual(int order) {
  switch (order) {
  case 0:
    return m_beta00;
    break;
  case 1:
    return m_beta01;
    break;
  case 2:
    return  m_beta02;//m_born*(m_delI2)*(m_delF2);;
    break;
  case 3:
    return m_beta03;//m_born*(m_delI2 + pow(m_gammaI,3.)/48)*(1.+ m_delF2+ pow(m_gammaF,3)/48.);;
    break;
  default:
    return 0;
  }
}


void Real_ff::SetIncoming(YFS::Dipole_Vector::iterator dipole, Vec4D_Vector &born, Vec4D_Vector &k, std::vector<double> y, std::vector<double> z) {
  YFS::Dipole D=dipole[0];
  m_photons = D.GetPhotons();
  m_y = y;
  m_z = z;
  m_sfsr.clear();
  m_beta10f.clear();
  m_beta11f.clear();
  m_beta12f.clear();
  m_beta21f.clear();
  m_alpi = m_alpha / M_PI;
  m_q1 = D.m_oldmomenta[0];
  m_q2 = D.m_oldmomenta[1];
  m_beam1 = D.m_oldmomenta[0]; // called beam not to break for now. Is actually born final state momentum
  m_beam2 = D.m_oldmomenta[1];
  m_p1p2 = m_beam1 * m_beam2;
  Vec4D sumk;
  for (auto kk : k) sumk += kk;
  m_spp  = (D.m_newmomenta[0] + D.m_newmomenta[1]).Abs2();
  m_u  = 1 - (m_spp) / sqr(rpa->gen.Ecms());

  double beta1 = CalculateBeta(m_q1);
  double beta2 = CalculateBeta(m_q2);
  double t1 = (1. + beta1 * beta2) / (beta1 + beta2);
  double logarg =  2.*(1. + beta1) * (1. + beta2) / ((1. - beta1) * (1. - beta2));
  double QF2 = D.m_QiQj;
  logarg = (D.m_newmomenta[0] + D.m_newmomenta[1]).Abs2() / sqr(m_beam1.Mass());
  m_gamma  = 2 * m_alpi * (log(logarg) - 1.); // See Mareks phd thesis A.2.1
  m_gammap = m_alpi * t1 * (log(logarg / 2.));
  m_gammaF = m_gamma;
  m_mass = (m_q1.Mass() + m_q2.Mass()) / 2.;
  m_mass2 = sqr(m_mass);
  for (auto &kk : k) m_sfsr.push_back(Eikonal(kk));
  ATOOLS::Poincare poin(m_beam1 + m_beam2);
  for (auto &kk : m_photons) poin.Boost(kk);
  m_gammaF  = 2 * m_alpi * (log(logarg) - 1.);
  m_delF2 = sqr(m_alpi * (log(logarg))) / 2.;
  if (m_fsrmode == 0) m_gammaF = m_delF2 = 0;
  if (m_fsrmode == 2) m_gammaI = m_delI2 = 0;
}

void Real_ff::SetIncoming(YFS::Dipole *p_dipole) {
  m_sp   = p_dipole->Sprime();
  m_v = 1. - m_sp / sqr(rpa->gen.Ecms());
  // if (m_use_model_alpha) m_alpha = s_model->ScalarConstant("alpha_QED");
  // else m_alpha  = (*aqed)(0); // cant be in constructor
  m_alpi = m_alpha / M_PI;
  m_beam1 = p_dipole->GetBornMomenta(0);
  m_beam2 = p_dipole->GetBornMomenta(1);
  m_born = m_beta00;
  m_sisr.clear();
  m_beta10i.clear();
  m_beta11i.clear();
  m_beta12i.clear();
  m_beta21i.clear();
  m_p1p2 = m_beam1 * m_beam2;

  m_mass   = p_dipole->Mass();
  m_mass2  = sqr(m_mass);
  m_L      = log(m_s / m_mass2);
  m_gamma  = 2.*m_alpi * (m_L - 1.);
  m_gammap = 2.*m_alpi * (m_L);
  m_gammaI = m_gamma;
  m_photons = p_dipole->GetPhotons();
  double logargI =  (m_beam1 + m_beam2).Abs2() / sqr(m_beam1.Mass());
  m_gammaI  = 2 * m_alpi * (log(logargI) - 1.); // See Marek's phd thesis A.2.1
  m_delI2 = sqr(m_alpi * (log(logargI))) / 2.;
}

void Real_ff::ResetReal() {
  m_real = m_born * (1. + m_realFSR) * (1. + m_realISR);
}

double Real_ff::Calc(const ATOOLS::Vec4D_Vector &p) {
  msg_Error() << "Virtual Method called" << std::endl;
  return 0;
}


double Real_ff::operator()(const ATOOLS::Vec4D_Vector& momenta) {
  return 0;
  msg_Error() << "Virtual Method called" << std::endl;
}

void Real_ff::Calculate() {
  // CalculateVirt();
  m_beta10.clear();
  // m_born=1;
  m_beta11.clear();
  m_beta12.clear();
  m_Sfac.clear();
  m_Hfac.clear();
  m_d10vec.clear();
  m_beta1 = 0;
  m_real = 0;
  m_beta20 = m_beta21 = 0;
  if (m_order == 1 && !m_realtool) {
    for (auto k: m_photons) {
      if (m_fsrmode == 0 ) {
        // PRINT_V
        m_real += Beta10(k);///m_Sfac[i];///Eikonal(m_photons[i]) ;///m_Sfac[i];
        // m_beta10i.push_back(Beta10(m_photons[i], i));
      }

      else if (m_fsrmode != 0)  {
        m_real += Beta10(k); ///m_Sfac[i];//Eikonal(m_photons[i]);//+m_beta11;
        // m_beta10f.push_back(Beta10(m_photons[i], i));
      }
    }
  }
  if (m_order == 2) {
    m_beta20 = m_beta21 = 0;
    m_real = 0;
    for (int j = 1; j < m_photons.size(); j++) {
      for (int i = 0; i < j; i++) {
        // if(i==j) continue;
        if (m_fsrmode == 0) {
          m_real += Beta20(m_photons[i], m_photons[j], i, j); ///(Eikonal(m_photons[i])*Eikonal(m_photons[j]));
          m_beta21i.push_back(Beta20(m_photons[i], m_photons[j], i, j));
        }
        else {
          m_real += Beta20(m_photons[i], m_photons[j], i, j); ///(Eikonal(m_photons[i])*Eikonal(m_photons[j]));
          m_beta21f.push_back(Beta20(m_photons[i], m_photons[j], i, j));
        }
      }
    }
    for (int i = 0; i < m_photons.size(); i++) {
      if (m_fsrmode == 0) {
        m_real += Beta11(m_photons[i], i); ///Eikonal(m_photons[i]);
        m_beta11i.push_back(Beta11(m_photons[i], i));
      }
      else {
        m_real += Beta11(m_photons[i], i); ///Eikonal(m_photons[i]);
        m_beta11f.push_back(Beta11(m_photons[i], i));
      }
    }
  }
  if (m_order == 3) {
    m_real = 0;
    m_beta30 = 0;
    Beta3();
    m_real += m_beta30;
    for (int i = 0; i < m_photons.size(); i++) {
      if (m_fsrmode == 0) {
        m_real += Beta12(m_photons[i], i); ///Eikonal(m_photons[i]);
        m_beta12i.push_back(Beta12(m_photons[i], i));

      }
      else {
        m_real += Beta12(m_photons[i], i); ///Eikonal(m_photons[i]);
        m_beta12f.push_back(Beta12(m_photons[i], i));
      }
    }
    for (int j = 1; j < m_photons.size(); j++) {
      for (int i = 0; i < j; ++i) {
        if (m_fsrmode == 0) {
          m_real += Beta21(m_photons[i], m_photons[j], i, j); ///(Eikonal(m_photons[i])*Eikonal(m_photons[j]));
          m_beta21i.push_back(Beta21(m_photons[i], m_photons[j], i, j));
        }
        else {
          // m_real+= Beta21(m_photons[i],m_photons[j], i, j);///(Eikonal(m_photons[i])*Eikonal(m_photons[j]));
          m_beta21f.push_back(Beta21(m_photons[i], m_photons[j], i, j));
        }
      }
    }
  }
  if (m_order == 3) {
    m_real = 0;
    m_beta30 = 0;
    Beta3();
    m_real += m_beta30;
    for (int i = 0; i < m_photons.size(); i++) {
      if (m_fsrmode == 0) {
        m_real += Beta12(m_photons[i], i); ///Eikonal(m_photons[i]);
        m_beta12i.push_back(Beta12(m_photons[i], i));

      }
      else {
        m_real += Beta12(m_photons[i], i); ///Eikonal(m_photons[i]);
        m_beta12f.push_back(Beta12(m_photons[i], i));
      }
    }
    for (int j = 1; j < m_photons.size(); j++) {
      for (int i = 0; i < j; ++i) {
        if (m_fsrmode == 0) {
          m_real += Beta20(m_photons[i], m_photons[j], i, j); ///(Eikonal(m_photons[i])*Eikonal(m_photons[j]));
          m_beta21i.push_back(Beta21(m_photons[i], m_photons[j], i, j));
        }
        else {
          m_real += Beta20(m_photons[i], m_photons[j], i, j); ///(Eikonal(m_photons[i])*Eikonal(m_photons[j]));
          m_beta21f.push_back(Beta21(m_photons[i], m_photons[j], i, j));
        }
      }
    }
  }

  DEBUG_FUNC("Real + Virtual = " << m_real << std::endl);
}


double Real_ff::W(double a, double b) {
  double den = pow(1. - a, 2.) + pow(1. - b, 2.);
  double t1  = m_mass2 / (2 * m_beam1 * m_beam2) * (1. - a) * (1 - b) / den;
  double t2  = a / b + b / a;

  return 1. - t1 * t2;
}

double Real_ff::wm0(double a, double b) {
  double del;
  if (m_fsrmode == 0) del = m_mass2 / m_s;
  else del = m_mass2 / m_spp;
  return 1. - 2.*del - del * (a / b + b / a);
}

double Real_ff::wmd(double a, double b) {
  double del;
  if (m_fsrmode == 0) del = m_mass2 / m_s;
  else del = m_mass2 / m_spp;
  return  1. + del * (a / b + b / a) * (a * a + b * b) / (pow(1. - a, 2) + pow(1. - b, 2));
}

double Real_ff::wmd(double del, double a, double b) {
  return  1. + del * (a / b + b / a) * (a * a + b * b) / (pow(1. - a, 2) + pow(1. - b, 2));
}


double Real_ff::Delta(int order, double z) {
  double delta;
  switch (order) {
  case 0: delta = 0.5 * m_gamma;
  case 1: delta = 0.5 * m_gamma - 0.25 * m_gamma * log(z);
  case 2: delta = 0.5 * m_gamma - 0.25 * m_gamma * log(z) + m_gamma * m_gamma * (0.125 - 0.125 * log(z) + sqr(log(z)) / 24.);
  }
  if (z < 0) delta = 0;
  if (IsBad(delta)) {
    msg_Error() << "YFS m_virtual = " << delta << std::endl
                << " for z = " << z << std::endl;
  }
  m_virtual = delta;
  return delta;
}

double Real_ff::Alpha(Vec4D k) {
  return k * m_beam2 / (m_beam1 * m_beam2);
}


double Real_ff::Beta(Vec4D k) {
  return k * m_beam1 / (m_beam1 * m_beam2);
}



double Real_ff::Eikonal(Vec4D k) {
  return -m_alpi / (4.*M_PI) * (m_beam1 / (m_beam1 * k) - m_beam2 / (m_beam2 * k)).Abs2();

}

double Real_ff::chi(double u, double a, double b) {
  return 0.25 * sqr(1. - u) * (sqr(1. - a) + sqr(1. - b));
}

double Real_ff::chi3(double u, double a1, double b1, double a2, double b2) {
  return 0.125 * sqr(1 - u) * (sqr(1. - a1) + sqr(1. - a1)) * (sqr(1 - b1) + sqr(1. - b2));
}

double Real_ff::EikonalInterferance(Vec4D k) {
  return m_alpi / (4.*M_PI) * 2.*m_beam1 * m_beam2 / ((k * m_beam1) * (k * m_beam2));
}

void Real_ff::D1(Vec4D k, double a, double b, double wm, int order) {
  // double del1 = 0;
  double kk = -1;
  if (m_fsrmode >= 1) kk = 1;
  double E = EikonalInterferance(k);
  double t1 = 0.5 * (pow(1. - a, 2.) + pow(1. - b, 2.)) * m_beta00;
  double z = (1. - a) * (1. - b);
  double del2 = 0;
  double del1 = m_gamma / 2. + kk * m_gamma / 4.*log(z);
  // +m_alpi*(log(a)*log(1-b)+log(b)*log(1-a));
  //               +DiLog(a)+DiLog(b)
  //               -1./2*sqr(log(1-a))  -1./2*sqr(log(1-b))
  //               +3./2*log(1-a)     +3./2*log(1-b)
  //               +1./2*a*(1-a)/(1+sqr(1-a))
  //               +1./2*b*(1-b)/(1+sqr(1-b)));
  // if(m_fsrmode == 0){
  del2 = sqr(m_gamma) / 8.*(1. - log(z)) + sqr(m_gamma) / 24.*sqr(log(z));
  //   if(m_usenllv) del1 = m_gammaI/2.-m_gammaI/4.*log(z);
  // }
  // else{
  //   del1 = m_gammaF/2.+m_gammaF/4.;
  // }
  if (m_fsrmode >= 1) del2  = 0 ;
  m_D10 = t1;//*(1+0.5*m_gamma);//*(1.+0.5*log(1-a-b)));
  m_D11 = m_D10 * (1. + del1);
  m_D12 = m_D10 * (1. + del1 + del2);
  if (IsBad(del1) || IsBad(del2)) {
    PRINT_VAR(del1);
    PRINT_VAR(del2);
    PRINT_VAR(z);
    PRINT_VAR(m_gamma);
  }
}



double Real_ff::Beta10(Vec4D k) {
  double d, t1, wm, S, w0;
  double a = k * m_beam1 / (m_beam1 * m_beam2);
  double b = k * m_beam2 / (m_beam1 * m_beam2);
  double at = a;
  double bt = b;
  S = Eikonal(k);
  double virtfac = 1;
  // S = 2./(a*b)*wmd(a,b);
  if (m_fsrmode >= 1) {
    // a = m_y[i];
    // b = m_z[i];
    at = a / (1. + a + b);
    bt = b / (1. + a + b);
    if (m_use_fac) virtfac = 1 + m_gammaI / 2;
    // S = 2./(a*b)*wm0(a,b);
    // double hfac = S*wmd(m_mass2/m_s,at,bt);
    // m_Sfac.push_back(S);
    D1(k, at, bt, wm, m_order);
    m_D10 *= S;//*wmd(a,b);//*(1+m_gammaI/2);;
    m_D11 *= S;
    m_D12 *= S;
    m_d10vec.push_back(m_D10);
    double beta10 = (m_D10 - S * m_beta00) * virtfac; //*(1+m_gamma/2);
    // return 0.5;
    // return m_D10/S;
    return beta10 / S;
  }
  m_D10 = 0;
  if (m_use_fac) virtfac = 1 + m_gammaF / 2;
  D1(k, a, b, wm, m_order);
  m_Sfac.push_back(S);
  m_D10 *= S;//*wmd(a,b);
  m_D11 *= S;
  m_D12 *= S;

  m_d10vec.push_back(m_D10);
  double beta10 = (m_D10 - S * m_beta00) * virtfac; //*(1+m_gamma/2);
  return beta10 / S; ///S;
}

double Real_ff::Beta11(Vec4D k, int i) {
  double d, t1, wm, S, w0;
  double a = k * m_beam1 / (m_beam1 * m_beam2);
  double b = k * m_beam2 / (m_beam1 * m_beam2);
  // double a = m_y[i];
  // double b = m_z[i];
  double at = a;
  double bt = b;
  wm = W(a, b);
  S = Eikonal(k);
  // S = 2./(a*b)*wmd(a,b);

  double virtfac = 1;

  if (m_fsrmode >= 1) {
    at = a / (1. + a + b);
    bt = b / (1. + a + b);
    D1(k, at, bt, wm, m_order);
    m_D10 *= S * wmd(a, b); //*(1+m_gammaI/2);
    m_D11 *= S * (1 + m_gammaI / 2);
    m_D12 *= S;
    m_d10vec.push_back(m_D10);
    if (m_use_fac) virtfac = 1 + m_gammaI / 2;
    double beta11 = (m_D11 - S * m_beta01);
    // if(m_formfactor == 2 ) beta11 = (m_D11 - S*m_beta00);
    return beta11 / S;
  }
  D1(k, a, b, wm, m_order);
  m_D10 *= S;//*wmd(a,b);//*(1+m_gammaF/2);
  m_D11 *= S;//*wmd(a,b);//*(1+m_gammaF/2);
  m_D12 *= S;
  // m_Sfac.push_back(S);
  m_d10vec.push_back(m_D10);
  if (m_use_fac) virtfac = 1 + m_gammaF / 2;
  m_D11 *= 1 + m_gammaF / 2;
  double beta11 = (m_D11 - S * m_beta01);
  // if(m_formfactor == 2 ) beta11 = (m_D11 - S*m_beta00);

  return beta11 / S;
}


double Real_ff::Beta12(Vec4D k, int i) {
  double d, t1, wm, S, w0;
  double a = k * m_beam1 / (m_beam1 * m_beam2);
  double b = k * m_beam2 / (m_beam1 * m_beam2);
  // a = m_y[i];
  // b = m_z[i];
  double at = a;
  double bt = b;
  wm = W(a, b);
  S = Eikonal(k);
  // S = 2/(a*b)*wm0(a,b);
  double virtfac = 1;

  if (m_fsrmode >= 1) {
    at = a / (1. + a + b);
    bt = b / (1. + a + b);
    // S = 2/(a*b)*wm0(a,b);
    // m_Sfac.push_back(S);
    D1(k, at, bt, wm, m_order);
    m_D10 *= S;
    m_D11 *= S;
    m_D12 *= S;
    // m_Sfac.push_back(S);
    // m_d10vec.push_back(m_D10);
    if (m_use_fac) virtfac = 1 + m_gamma / 2;
    // Factorised corrections from Stasek
    // Shown to be closer to analytical
    // bety12 = (1+deli1+deli2)*(dist11 -m_Beta00*(1+delf1)*sfacj)
    double beta12 = (1 + m_gammaI / 2 + m_delI2) * (m_D11 - S * m_beta00 * (1 + m_gammaF / 2.)); //*virtfac;
    // double beta12 = (1+m_gammaF/2+sqr(m_gammaF/2))*(m_D12 - S*m_beta02);//*virtfac;
    return beta12 / S;
  }
  D1(k, a, b, wm, m_order);
  m_D10 *= S;//*wmd(at,bt);//*(1+m_gammaF/2);
  m_D11 *= S;//*wmd(at,bt);//*(1+m_gammaF/2);
  m_D12 *= S;//*wmd(at,bt);
  // m_Sfac.push_back(S);
  if (m_use_fac) virtfac = 1 + m_gammaF / 2 + m_delF2;
  m_d10vec.push_back(m_D10);
  m_D12 *= virtfac;
  double beta12 = (m_D12 - S * m_beta02); //*virtfac;
  return beta12 / S;
}


void Real_ff::D2(int i, int j)
{
  double A = -1;
  if (m_fsrmode >= 1) A = 1;

  Vec4D k1 = m_photons[i];
  Vec4D k2 = m_photons[j];
  double a1 = k1 * m_beam1 / m_p1p2;
  double b1 = k1 * m_beam2 / m_p1p2;
  double a2 = k2 * m_beam1 / m_p1p2;
  double b2 = k2 * m_beam2 / m_p1p2;
  double a1p = a1 / (1. + A * a2);
  double b1p = b1 / (1. + A * b2);
  double a2p = a2 / (1. + A * a1);
  double b2p = b2 / (1. + A * b1);
  double a1pp, a2pp, b1pp, b2pp;
  double p0, p1;

  if (m_fsrmode != 0) {
    double aa1 = m_y[i];
    double bb1 = m_z[i];
    double aa2 = m_y[j];
    double bb2 = m_z[j];
    a1 = aa1 / (1 + aa1 + bb1);
    b1 = bb1 / (1 + aa1 + bb1);
    a2 = aa2 / (1 + aa2 + bb2);
    b2 = bb2 / (1 + aa2 + bb2);

    a1pp = a1 / (1. + a2);
    b1pp = b1 / (1. + b2);
    a2pp = a2 / (1. + a1);
    b2pp = b2 / (1. + b1);

    a1p = a1pp / (1 + a1pp + b1pp);
    a2p = a2pp / (1 + a2pp + b2pp);

    b1p = b1pp / (1 + a1pp + b1pp);
    b2p = b2pp / (1 + a2pp + b2pp);
    double v1 = a1 + b1;
    double v2 = a2 + b2;
    if (v2 >= v1) {
      p0 = chi(a1, a2pp, b2pp);
      p1 = chi(b1, a2pp, b2pp);
    }
    else {
      p0 = chi(a2, a1pp, b1pp);
      p1 = chi(b2, a1pp, b1pp);
    }
    m_g1 = (p0 + p1);
    // PRINT_VAR(m_g1);
    // PRINT_VAR(p0);
    // PRINT_VAR(p1);
    // PRINT_VAR(m_born);
    // PRINT_VAR(m_beam1);
    // PRINT_VAR(m_beam2);
    // PRINT_VAR(m_p1p2);
    // PRINT_VAR(a1);
    // PRINT_VAR(b1);
    // PRINT_VAR(a2);
    // PRINT_VAR(a2);
    double z1 = (1. - a1) * (1. - b1);
    double z1z2 = (1. - a1 - a2) * (1. - b1 - b2);
    double delv = m_gamma / 2. + A * m_gamma / 6.*(log(z1) + log(z1z2));
    if (z1 <= 0 || z1z2 <= 0 ) delv = 0;
    // if(m_fsrmode>=1) delv=0;
    m_D20 += m_born * m_g1;
    m_D21 += m_born * m_g1; //*(1+delv);
    return;
  }
  else {
    // a1 = m_y[i];///(1+m_y[i]+m_z[i]);
    // b1 = m_z[i];///(1+m_y[i]+m_z[i]);
    // a2 = m_y[j];///(1+m_y[j]+m_z[j]);
    // b2 = m_z[j];///(1+m_y[j]+m_z[j]);
    a1p = a1 / (1. + A * a2);
    b1p = b1 / (1. + A * b2);
    a2p = a2 / (1. + A * a1);
    b2p = b2 / (1. + A * b1);
  }

  double v1 = a1 + b1;
  double v2 = a2 + b2;
  if (v2  > v1) {
    p0 = chi(a2, a1p, b1p);
    p1 = chi(b2, a1p, b1p);
  }
  else {
    p0 = chi(a1, a2p, b2p);
    p1 = chi(b1, a2p, b2p);
  }
  m_g1 = p0 + p1;
  // m_gg2 += p0+p1;
  double z1 = (1. - a1) * (1. - b1);
  double z1z2 = (1. - a1 - a2) * (1. - b1 - b2);
  double delv = m_gamma / 2. - m_gamma / 6.*(log(z1) + log(z1z2));
  m_gg2 = (p0 + p1);
  // if(z1 <= 0 || z1z2 <= 0 ) delv = 0;
  if (m_fsrmode >= 1) delv = 0;
  m_D20 += m_beta00 * m_g1 * (1); //+delv);//*wmd(a1,b1)*wmd(a2,b2); // !!! Must include delv
  m_D21 += m_beta00 * m_gg2 * (1 + delv); //*wm0(a1,b1)*wm0(a2,b2);
}

void Real_ff::D3(int i, int j, int k) {
  double A = -1;
  if (m_fsrmode >= 1) A = 1;
  Vec4D k1 = m_photons[i];
  Vec4D k2 = m_photons[j];
  Vec4D k3 = m_photons[k];
  double a1 = k1 * m_beam1 / m_p1p2;
  double b1 = k1 * m_beam2 / m_p1p2;
  double a2 = k2 * m_beam1 / m_p1p2;
  double b2 = k2 * m_beam2 / m_p1p2;
  double a3 = k3 * m_beam1 / m_p1p2;
  double b3 = k3 * m_beam2 / m_p1p2;
  // if(!m_usefulleik && m_fsrmode!=0){
  // a1 = m_y[i];
  // b1 = m_z[i];
  // a2 = m_y[j];
  // b2 = m_z[j];
  // a3 = m_y[k];
  // b3 = m_z[k];
  // }

  double ap2 = a2 / (1 - a1);
  double bp2 = b2 / (1 - b1);

  double ap3 = a3 / (1 - a1 - a2);
  double bp3 = b3 / (1 - b1 - b2);
  double p1, p2;
  if (k1.E() >= k2.E() && k1.E() >= k3.E()) {
    p1 = 0.125 * (sqr(1. - a1)) * (sqr(1. - ap2) + sqr(1. - bp2)) * (sqr(1. - ap3) + sqr(1. - bp3));
    p2 = 0.125 * (sqr(1. - b1)) * (sqr(1. - ap2) + sqr(1. - bp2)) * (sqr(1. - ap3) + sqr(1. - bp3));
  }
  else if (k2.E() >= k1.E() && k2.E() >= k3.E()) {
    p1 = 0.125 * (sqr(1. - ap2)) * (sqr(1. - a1) + sqr(1. - b1)) * (sqr(1. - ap3) + sqr(1. - bp3));
    p2 = 0.125 * (sqr(1. - bp2)) * (sqr(1. - a1) + sqr(1. - b1)) * (sqr(1. - ap3) + sqr(1. - bp3));
  }
  else if (k3.E() >= k1.E() && k3.E() >= k2.E()) {
    p1 = 0.125 * (sqr(1. - ap3)) * (sqr(1. - a1) + sqr(1. - b1)) * (sqr(1. - ap2) + sqr(1. - bp2));
    p2 = 0.125 * (sqr(1. - bp3)) * (sqr(1. - a1) + sqr(1. - b1)) * (sqr(1. - ap2) + sqr(1. - bp2));
  }
  else
  {
    p1 = 0.125 * (sqr(1. - a1)) * (sqr(1. - ap2) + sqr(1. - bp2)) * (sqr(1. - ap3) + sqr(1. - bp3));
    p2 = 0.125 * (sqr(1. - b1)) * (sqr(1. - a2) + sqr(1. - bp2)) * (sqr(1. - ap3) + sqr(1. - bp3));
  }
  m_ggg1 += p1 + p2;
  // m_ggg2 += p2;
}

void Real_ff::Beta2() {
  m_p1p2 = m_beam1 * m_beam2;
  double b20(0), b21(0);
  m_beta20 = m_beta21 = 0;
  m_g1 = m_g2  = m_gg1 = m_gg2 = m_ggg1 = 0;
  for (int i = 1; i < m_photons.size(); ++i) {
    for (int j = 0; j < i; ++j) {
      m_g1 = m_g2  = m_gg1 = m_gg2 = m_ggg1 = 0;
      Vec4D ki = m_photons[i];
      Vec4D kj = m_photons[j];
      double s1 = Eikonal(ki);
      double s2 = Eikonal(kj);
      D2(i, j);
      // D2(kj,ki);
      m_D20 *= s1 * s2;
      m_D21 *= s1 * s2;
      if (m_fsrmode == 0) {
        b20 = (m_D20 - m_beta00 * s1 * s2
               - Beta10(ki) * s2 - Beta10(kj) * s1);
        b21 = (m_D21 - m_beta01f * s1 * s2
               - Beta11(ki, i) * s2 - Beta11(kj, j) * s1);
      }
      else {
        b20 = (m_D20 - m_beta00 * s1 * s2
               - Beta10(ki) * s2 - Beta10(kj) * s1);
        b21 = (m_D21 - m_beta01f * s1 * s2
               - Beta11(ki, i) * s2 - Beta11(kj, j) * s1);

      }
      // b20 = Beta20(ki,kj);
      m_beta20 += b20 / (s1 * s2);
      m_beta21 += b21 / (s1 * s2);
      m_BETA20[std::make_pair(i, j)] = b20;
      m_BETA21[std::make_pair(i, j)] = b21;
    }

  }
  DEBUG_FUNC("Born = " << m_born << std::endl <<
             "D20 of photon k = " << m_D20 << std::endl <<
             "B_{00} = " << m_beta00 << std::endl <<
             "B_{01} = " << m_beta01 << std::endl <<
             "B_{10} = " << m_beta10 << std::endl <<
             "B_{11} = " << m_beta11 << std::endl <<
             "B_{12} = " << m_beta12 << std::endl <<
             "B_{20} = " << m_beta20 << std::endl);
  // PRINT_VAR("END SUM");
}

double Real_ff::Beta20(Vec4D k1, Vec4D k2, int i, int k) {
  double b20(0), b21(0);
  m_beta20 = m_beta21 = 0;
  m_D20 = m_D21 = 0;
  double lam = -1;
  if (m_fsrmode != 0) lam = 1;
  m_g1 = m_g2  = m_gg1 = m_gg2 = 0;
  double s1 = Eikonal(k1);
  double s2 = Eikonal(k2);
  D2(i, k);
  D2(k, i);
  m_D20 *= 0.5 * s1 * s2;
  double v01;
  if (m_fsrmode == 0) v01 = m_beta01i;
  else v01 = m_beta01f;
  double virtfac = 1;
  if (m_use_fac && m_fsrmode == 0 ) virtfac = 1 + m_gammaF / 2;
  if (m_use_fac && m_fsrmode != 0 ) virtfac = 1 + m_gammaI / 2;
  // Beta10 is divided by S so need to multiply by extra s
  b20 = (m_D20 - m_beta00 * s1 * s2
         - Beta11(k1,i) * s2 * s1/virtfac - Beta11(k2,k) * s1 * s2/virtfac);
  return b20 / (s1 * s2);
}

double Real_ff::Beta21(Vec4D k1, Vec4D k2, int i, int k) {
  double b20(0), b21(0);
  double lam = -1;
  if (m_fsrmode != 0) lam = 1;

  m_beta20 = m_beta21 = 0;
  m_g1 = m_g2  = m_gg1 = m_gg2 = 0;
  m_D20 = 0;
  m_D21 = 0;
  double s1 = Eikonal(k1);
  double s2 = Eikonal(k2);
  D2(i, k);
  D2(k, i);
  m_D20 *= 0.5 * s1 * s2; //*wmd(m_y[i],m_z[i])*wmd(m_y[k],m_z[k]);
  m_D21 *= 0.5 * s1 * s2; //*wmd(m_y[i],m_z[i])*wmd(m_y[k],m_z[k]);
  double v01;
  if (m_fsrmode == 0) v01 = m_beta01i;
  else v01 = m_beta01f;
  // Beta10 is divided bt S so need to multiply by extra s
  double virtfac = 1;
  if (m_use_fac) virtfac = 1 + m_gammaF / 2;
  if (m_fsrmode == 0) {
    m_D21 *= 1 + m_gammaF / 2.;
    b21 = (m_D21 - m_beta01 * s1 * s2
           - Beta11(k1, i) * s2 * s1 / virtfac - Beta11(k2, k) * s1 * s2 / virtfac);
  }
  else {
    // Primitaive virt corrections
    b21 = (1 + m_gammaI / 2.) * Beta20(k1, k2, i, k) * s1 * s2;
    // b21 = Beta20(k1, k2, i, k)*s1*s2;
    // (m_D21-m_beta01*s1*s2
    // -Beta11(k1,i)*s2*s1-Beta11(k2,k)*s1*s2);
  }
  return b21 / s1 / s2;
}





double Real_ff::Beta20Int(Vec4D k1, Vec4D k2, int i, int j, int mode) {
  // mode 0 IF
  // mode 1 II
  // mode 2 FF
  double y1, z1, y2, z2, y3, z3, yf1, zf1, yfp, zfp, yf2, zf2, sfac1, sfac2, d20int, bint;
  if (mode == 0) {
    y1 = m_yisr[i];
    z1 = m_zisr[i];
    yf1 = m_yfsr[j];
    zf1 = m_zfsr[j];
    y2 = yf1 / (1 + yf1 + zf1);
    z2 = zf1 / (1 + yf1 + zf1);

    d20int = sqr(1 - y1) * sqr(1 - y2)
             + sqr(1 - y1) * sqr(1 - z2)
             + sqr(1 - z1) * sqr(1 - y2)
             + sqr(1 - z1) * sqr(1 - z2);

    sfac1 = m_sisr[i];
    sfac2 = m_sfsr[j];
    d20int *= sfac1 * sfac2 * m_born / 4;
    bint = d20int - m_beta10f[j] * sfac1
           - m_beta10i[i] * sfac2
           - m_beta00 * sfac2 * sfac1;
  }
  else if (mode == 1) {
    y1 = m_yisr[i];
    z1 = m_zisr[i];
    yf1 = m_yisr[j];
    zf1 = m_zisr[j];

    y2 = y1 / (1 - yf1);
    z2 = z1 / (1 - zf1);
    y3 = yf1 / (1 - y1);
    z3 = zf1 / (1 - z1);

    d20int = chi(y1, y2, z2) + chi(z1, y3, z3); //+chi(y1,y3,z3)+chi(z1,y2,z2);
    sfac1 = m_sisr[i];
    sfac2 = m_sisr[j];

    d20int *= sfac1 * sfac2 * m_born;
    bint = d20int - m_beta10i[j] * sfac1
           - m_beta10i[i] * sfac2
           - m_beta00 * sfac2 * sfac1;
  }
  else if (mode == 2) {
    y1 = m_yfsr[i] / (1 + m_yfsr[i] + m_zfsr[i]);
    z1 = m_zfsr[i] / (1 + m_yfsr[i] + m_zfsr[i]);

    y2 = m_yfsr[j] / (1 + m_yfsr[j] + m_zfsr[j]);
    z2 = m_zfsr[j] / (1 + m_yfsr[j] + m_zfsr[j]);

    yf1 = y1 / (1 + yf2);
    zf1 = z1 / (1 + zf2);

    yf2 = y2 / (1 + yf1);
    zf2 = z2 / (1 + zf1);

    double yp1 = yf1 / (1 + yf1 + zf1);
    double zp1 = zf1 / (1 + yf1 + zf1);

    double yp2 = yf2 / (1 + yf2 + zf2);
    double zp2 = z2 / (1 + yf2 + zf2);


    d20int = chi(y1, yp2, zp2) + chi(z1, yp2, zp2);
    d20int += chi(y2, yp1, zp1) + chi(z2, yp1, zp1);
    sfac1 = m_sfsr[i];
    sfac2 = m_sfsr[j];
    d20int *= sfac1 * sfac2 * m_born;
    bint = d20int - m_beta10f[j] * sfac1
           - m_beta10f[i] * sfac2
           - m_beta00 * sfac2 * sfac1;
  }
  return bint;
}


double Real_ff::IntIF(Vec4D_Vector &k1, Vec4D_Vector &k2, std::vector<double> yi, std::vector<double> zi, std::vector<double> yf, std::vector<double> zf) {
  m_beta20Int = 0;
  m_BETA20.clear();
  m_yisr = yi;
  m_zisr = zi;
  m_yfsr = yf;
  m_zfsr = zf;
  if ((m_beta10f.size() != k2.size()) || (m_beta10i.size() != k1.size())) {
    msg_Error() << "Error in YFS interference\n";
    msg_Out() << k1.size() << " ISR photons with " << m_beta10i.size() << " Real corrections\n";
    msg_Out() << k2.size() << " FSR photons with " << m_beta10f.size() << " Real corrections\n";
  }

  for (int i = 0; i < k1.size(); i++) {
    for (int j = 0; j < k2.size(); j++) {
      m_beta20Int += Beta20Int(k1[i], k2[j], i, j, 0) / (m_sisr[i] * m_sfsr[j]);
      // m_BETA20[std::make_pair(i,j)] =  bint;///(sfac1*sfac2);
      // m_BETA20[std::make_pair(j,i)] =  bint;
    }
  }
  double beta3iif = 0;
  for (int i = 2; i < k1.size(); ++i) {
    for (int j = 0; j < i; ++j) {
      for (int k = 0; k < k2.size(); ++k) {
        if (i == j) continue;
        m_ggg1 = m_ggg2 = 0;
        Vec4D ki = k1[i];
        Vec4D kj = k1[j];
        Vec4D kk = k2[k];
        double sfac1 = m_sisr[i];
        double sfac2 = m_sisr[j];
        double sfac3 = m_sfsr[k];
        D3Int(i, j, k, 1);
        D3Int(j, i, k, 1); // permute intial photons
        double SF = sfac1 * sfac2 * sfac3;
        m_D30 = m_born * m_ggg1 * SF / 2.;
        double sub = sfac1 * Beta20Int(kj, kk, j, k, 0) +
                     sfac2 * Beta20Int(ki, kk, i, k, 0) +
                     sfac3 * Beta20Int(ki, kj, i, j, 1) +
                     sfac1 * sfac2 * m_beta10f[k] +
                     sfac3 * sfac2 * m_beta10i[i] +
                     sfac1 * sfac3 * m_beta10i[j] +
                     sfac1 * sfac2 * sfac3 * m_beta00;
        // PRINT_VAR(m_D30);
        // PRINT_VAR(sub);
        // PRINT_VAR(sfac1*sfac2*sfac3);
        beta3iif += (m_D30 - sub) / (sfac1 * sfac2 * sfac3);
      }
    }
  }
  double beta3ffi = 0;
  for (int i = 2; i < k2.size(); ++i) {
    for (int j = 0; j < i; ++j) {
      for (int k = 0; k < k1.size(); ++k) {
        // if(i==j) continue;
        m_ggg1 = m_ggg2 = 0;
        Vec4D ki = k2[i];
        Vec4D kj = k2[j];
        Vec4D kk = k1[k];
        double sfac1 = m_sfsr[i];
        double sfac2 = m_sfsr[j];
        double sfac3 = m_sisr[k];
        D3Int(k, i, j, 2);
        D3Int(k, j, i, 2); // permute final photons
        double SF = sfac1 * sfac2 * sfac3;
        m_D30 = m_born * m_ggg1 * SF / 2.;
        // mode 0 IF
        // mode 1 II
        // mode 2 FF
        double sub = sfac1 * Beta20Int(kk, kj, k, j, 0) +
                     sfac2 * Beta20Int(kk, ki, k, i, 0) +
                     sfac3 * Beta20Int(ki, kj, i, j, 2) +
                     sfac1 * sfac2 * m_beta10i[k] +
                     sfac3 * sfac2 * m_beta10f[i] +
                     sfac1 * sfac3 * m_beta10f[j] +
                     sfac1 * sfac2 * sfac3 * m_beta00;
        // PRINT_VAR(m_D30);
        // PRINT_VAR(sub);
        // PRINT_VAR(sfac1*sfac2*sfac3);
        beta3ffi += (m_D30 - sub) / (sfac1 * sfac2 * sfac3);
      }
    }
  }
  // return beta3ffi;
  return m_beta20Int;//+beta3iif+beta3ffi;
}


void Real_ff::D3Int(int i, int j, int k, int mode) {
  // mode 1 IIF, mode 2 IFF
  double A = -1;
  double ap1, bp1, ap2, bp2, ap3, bp3, p1, p2, p3, yp1, yp2;
  // if(m_fsrmode>=1) A = 1;
  Vec4D k1 = m_photons[i];
  Vec4D k2 = m_photons[j];
  Vec4D k3 = m_photons[k];
  double a1 = k1 * m_beam1 / m_p1p2;
  double b1 = k1 * m_beam2 / m_p1p2;
  double a2 = k2 * m_beam1 / m_p1p2;
  double b2 = k2 * m_beam2 / m_p1p2;
  double a3 = k3 * m_beam1 / m_p1p2;
  double b3 = k3 * m_beam2 / m_p1p2;
  if (mode == 1) {
    a1 = m_yisr[i];
    b1 = m_zisr[i];
    a2 = m_yisr[j];
    b2 = m_zisr[j];
    a3 = m_yfsr[k] / (1 + m_yfsr[k] + m_zfsr[k]);;
    b3 = m_zfsr[k] / (1 + m_yfsr[k] + m_zfsr[k]);;

    ap1 = a1 / (1 - a2);
    bp1 = b1 / (1 - b2);
    ap2 = a2 / (1 - a1);
    bp2 = b2 / (1 - b1);

    ap3 = a3 / (1 + a3 + b3);
    bp3 = b3 / (1 + a3 + b3);

    double aap3 = ap3 / (1 + bp3);
    double bbp3 = bp3 / (1 + ap3);
    // p1 = 0.25*(sqr(1.-a1))*(sqr(1.-a2)+sqr(1.-b2));
    // p2 = 0.25*(sqr(1.-b1))*(sqr(1.-a2)+sqr(1.-b2));
    p1 = chi(a1, ap2, bp2) + chi(a2, ap1, bp1);
    p2 = chi(b1, ap2, bp2) + chi(b2, ap1, ap1);
    // if(a1>a2) p1 = 0;
    // else p2 = 0;
    p3 = 0.5 * (sqr(1 - ap3) + sqr(1 - bp3));
    m_ggg1 += (p1 + p2) * p3;
    // m_p1+p2+p3;
  }
  else if (mode == 2) {
    a1 = m_yisr[i];
    b1 = m_zisr[i];


    a2 = m_yfsr[j] / (1 + m_yfsr[j] + m_zfsr[j]);
    b2 = m_zfsr[j] / (1 + m_yfsr[j] + m_zfsr[j]);

    a3 = m_yfsr[k] / (1 + m_yfsr[k] + m_zfsr[k]);
    b3 = m_zfsr[k] / (1 + m_yfsr[k] + m_zfsr[k]);

    ap1 = a1 / (1 - a2);
    bp1 = b1 / (1 - b2);

    ap2 = a2 / (1 + a3);
    bp2 = b2 / (1 + b3);

    ap3 = a3 / (1 + a2);
    bp3 = b3 / (1 + b2);

    p1 = chi(a2, ap3, bp3) + chi(b2, ap3, bp3);
    p2 = chi(a3, ap2, bp2) + chi(b3, ap2, ap2);
    // if(a3>a2) p1 = 0;
    // else p2 = 0;
    p3 = 0.5 * (sqr(1 - a1) + sqr(1 - b1));
    // m_ggg1 += (p1+p2);//*p3;//*p3;
    m_ggg1 += (p1 + p2) * p3;

  }


  // m_ggg2 += p2;
}




// double Real_ff::Beta21(Vec4D k1, Vec4D k2){
//   m_p1p2 = m_beam1*m_beam2;
//   double b20(0), b21(0);
//   m_beta20 = m_beta21 = 0;
//   m_g1 = m_g2  = m_gg1 = m_gg2 = 0;
//   double s1 = Eikonal(k1);
//   double s2 = Eikonal(k2);
//   double ai = k1*m_beam1/m_p1p2;
//   double bi = k1*m_beam2/m_p1p2;
//   double aj = k2*m_beam1/m_p1p2;
//   double bj = k2*m_beam2/m_p1p2;
//   D2(k1,k2);
//   m_D20 *= s1*s2;
//   m_D21 *= s1*s2;
//   b21 = (m_D21-m_beta00*s1*s2
//               -Beta10(k1)*s2-Beta10(k2)*s1);
//   return b21;
// }



void Real_ff::Beta3() {
  m_beta30 = 0;
  double virtfac = 1;
  if (m_use_fac && m_fsrmode == 0 ) virtfac = 1 + m_gammaF / 2;
  if (m_use_fac && m_fsrmode != 0 ) virtfac = 1 + m_gammaI / 2;

  if (m_fsrmode >= 1) return;
  double d30 = 0;
  double A = -1;
  m_p1p2 = m_beam1 * m_beam2;
  for (int k = 2; k < m_photons.size(); ++k) {
    for (int i = 1; i < k; ++i) {
      for (int j = 0; j < i; ++j) {
        m_ggg1 = m_ggg2 = 0;
        Vec4D ki = m_photons[i];
        Vec4D kj = m_photons[j];
        Vec4D kk = m_photons[k];
        double sfac1 = Eikonal(ki);//2/(m_y[i]*m_z[i])*wm0(m_y[i], m_z[i]);
        double sfac2 = Eikonal(kj);//2/(m_y[j]*m_z[j])*wm0(m_y[j], m_z[j]);
        double sfac3 = Eikonal(kk);//2/(m_y[k]*m_z[k])*wm0(m_y[k], m_z[k]);
        D3(i, j, k);
        D3(j, i, k);
        D3(i, k, j);
        D3(j, k, i);
        D3(k, i, j);
        D3(k, j, i);
        double SF = sfac1 * sfac2 * sfac3; //*wmd(m_y[i], m_z[i])*wmd(m_y[j], m_z[j])*wmd(m_y[k], m_z[k]);
        m_D30 = m_born * (m_ggg1) * SF / 6; // 6 for 3 isr diagrams
        double sub = Beta20(kj, kk, j, k) / virtfac +
                     Beta20(ki, kk, i, k) / virtfac +
                     Beta20(ki, kj, i, j) / virtfac +
                     Beta10(kk) / virtfac +
                     Beta10(ki) / virtfac +
                     Beta10(kj) / virtfac +
                     m_beta00;
        // PRINT_VAR(m_D30);
        // PRINT_VAR(sub);
        // PRINT_VAR(sfac1*sfac2*sfac3);
        sub *= sfac1 * sfac2 * sfac3;
        m_beta30 += (m_D30 - sub) / (sfac1 * sfac2 * sfac3);
      }
    }
  }
  // PRINT_VAR(m_beta30);
  DEBUG_FUNC("Born = " << m_born << std::endl <<
             "B_{00} = " << m_beta00 << std::endl <<
             "B_{10} = " << m_beta10 << std::endl <<
             "B_{11} = " << m_beta11 << std::endl <<
             "B_{12} = " << m_beta12 << std::endl <<
             "B_{20} = " << m_beta20 << std::endl <<
             "B_{21} = " << m_beta21 << std::endl <<
             "B_{30} = " << m_beta30 << std::endl);
}



double Real_ff::CalculateBeta(const Vec4D& p) {
  return Vec3D(p).Abs() / p[0];
}


void Real_ff::Sort(Vec4D_Vector &p) {
  struct {
    bool operator()(Vec4D a, Vec4D b) const
    {
      return a.PPerp() > b.PPerp();
      // return a[0] > b[0];
    }
  } cmp;

  std::sort(p.begin(), p.end(), cmp);
}
