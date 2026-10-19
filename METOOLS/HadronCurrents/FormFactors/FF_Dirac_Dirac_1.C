#include "METOOLS/HadronCurrents/FormFactors/FF_Dirac_Dirac_1.H"
#include "METOOLS/HadronCurrents/FormFactors/RChL_Functions.H"
#include "METOOLS/HadronCurrents/FormFactors/Kstar_Decays.H"
#include "METOOLS/HadronCurrents/FormFactors/A0_Decays.H"
#include "METOOLS/HadronCurrents/Tools.H"

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;


namespace {

  Complex SachsToF1(const double GE,
                    const double GM,
                    const double tau)
  {
    return Complex((GE + tau * GM) / (1.0 + tau), 0.0);
  }

}

void FF_Dirac_Dirac_1_Base::FixMode()
{
  if (m_flavs[m_pi[0]].Kfcode() == kf_p_plus &&
      m_flavs[m_pi[1]].Kfcode() == kf_p_plus)   m_mode = FF_Dirac_Dirac_1_mode::pp;

  else if (m_flavs[m_pi[0]].Kfcode() == kf_n &&
      m_flavs[m_pi[1]].Kfcode() == kf_n)   m_mode = FF_Dirac_Dirac_1_mode::nn;

}
//FixMode() identifies which particle combination the form factor is being used for, setting m_mode to pp or nn based on the flavours in m_flavs.

//Construct() uses m_mode to select the appropriate setup: Construct_pp(params) or Construct_nn(params). T
//his is where you initialise any proton-specific or neutron-specific parameters.

//operator() evaluates the form factor by selecting the requested model through m_ffmodel, such as Dipole, Kelly, BBBA05 or YAHL18.

Complex FF_Dirac_Dirac_1_Base::operator()(const ATOOLS::Vec4D_Vector& moms){
 double Q2 = (moms[m_pi[0]] + moms[m_pi[1]]).Abs2();

switch (m_ffmodel) {
case ff_model::none:       return Complex(0.0, 0.0);


case ff_model::pp_Point:   return Complex(1.0, 0.0);
case ff_model::pp_Dipole:  return F1_pp_Dipole(Q2);
case ff_model::pp_Kelly:   return F1_pp_Kelly(Q2);
case ff_model::pp_BBBA05:  return F1_pp_BBBA05(Q2);
case ff_model::pp_YAHL18:  return F1_pp_YAHL18(Q2);

case ff_model::nn_Point:   return Complex(1.0, 0.0);
case ff_model::nn_Dipole:  return F1_nn_Dipole(Q2);
case ff_model::nn_Kelly:   return F1_nn_Kelly(Q2);
case ff_model::nn_BBBA05:  return F1_nn_BBBA05(Q2);
case ff_model::nn_YAHL18:  return F1_nn_YAHL18(Q2);

case ff_model::unknown:
default:                   return Complex(0.0, 0.0);
}
}


//construct

void FF_Dirac_Dirac_1_Base::Construct(const FF_Parameters& params)
{
  switch (m_mode) {
  case FF_Dirac_Dirac_1_mode::pp: Construct_pp(params); break;
  case FF_Dirac_Dirac_1_mode::nn: Construct_nn(params); break;
  default: break;
  }
}


void FF_Dirac_Dirac_1_Base::Construct_pp(
  const FF_Parameters& params)
{
Mp = Flavour(kf_p_plus).HadMass();
mu_p = 2.7928473446;

Lambda2 = 0.71;

if (params.m_parameterMap.find("Lambda2") !=
    params.m_parameterMap.end()) {
  Lambda2 = params.m_parameterMap.at("Lambda2");
}
}


void FF_Dirac_Dirac_1_Base::Construct_nn(
  const FF_Parameters& params)
{
Mn = Flavour(kf_n).HadMass();
mu_n = -1.9130427;

Lambda2 = 0.71;

if (params.m_parameterMap.find("Lambda2") !=
    params.m_parameterMap.end()) {
  Lambda2 = params.m_parameterMap.at("Lambda2");
}
}
  


//funcs

Complex FF_Dirac_Dirac_1_Base::F1_pp_Dipole(const double& Q2) {
  if (m_mode == FF_Dirac_Dirac_1_mode::pp) {
  const double tau = Q2 / (4.0 * Mp * Mp);
  const double GD = std::pow(1.0 + Q2 / Lambda2, -2.0);

  const double GE = GD;
  const double GM = mu_p * GD;
  return SachsToF1(GE, GM, tau);
}

return Complex(0.0, 0.0);
}


Complex FF_Dirac_Dirac_1_Base::F1_nn_Dipole(const double& Q2) {
  if (m_mode == FF_Dirac_Dirac_1_mode::nn) {
  const double tau = Q2 / (4.0 * Mn * Mn);
  const double GD = std::pow(1.0 + Q2 / Lambda2, -2.0);

  const double GM = mu_n * GD;
  const double GE = -mu_n * (tau/(1+5.6*tau))* GD;
  return SachsToF1(GE, GM, tau);
  }
  
  return Complex(0.0, 0.0);
}

Complex FF_Dirac_Dirac_1_Base::F1_pp_Kelly(const double& Q2)
{
  return Complex(1.0, 0.0);
}

Complex FF_Dirac_Dirac_1_Base::F1_pp_BBBA05(const double& Q2)
{
  return Complex(1.0, 0.0);
}

Complex FF_Dirac_Dirac_1_Base::F1_pp_YAHL18(const double& Q2)
{
  return Complex(1.0, 0.0);
}

Complex FF_Dirac_Dirac_1_Base::F1_nn_Kelly(const double& Q2)
{
  return Complex(1.0, 0.0);
}

Complex FF_Dirac_Dirac_1_Base::F1_nn_BBBA05(const double& Q2)
{
  return Complex(1.0, 0.0);
}

Complex FF_Dirac_Dirac_1_Base::F1_nn_YAHL18(const double& Q2)
{
  return Complex(1.0, 0.0);
}



//getter


DECLARE_FF_GETTER( FF_Dirac_Dirac_1_Base, "FF_Dirac_Dirac_1")

FormFactor_Base * ATOOLS::Getter<FormFactor_Base,FF_Parameters,
                                 FF_Dirac_Dirac_1_Base>::
operator()(const METOOLS::FF_Parameters &params) const{
  msg_Out()<<METHOD<<":\n";

  size_t Nbaryons = 0;
  for (size_t i=0;i<params.m_pi.size();i++) {
    if (params.m_flavs[params.m_pi[i]].IsBaryon()) Nbaryons++;
  }
  if (Nbaryons!=2) return NULL;

  if ((params.m_flavs[params.m_pi[0]].Kfcode()==kf_p_plus &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_p_plus)    ||
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_n &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_n)         ) {
    return new FF_Dirac_Dirac_1_Base(params);
  }

  return NULL;
}

//CHANGE!!!!!!!!!!!!!!!!!!!!!!!!
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//
// pi^0 pi^+ form factors.
// Todos:
//  - add more theory - i.e. Gounaris-Sakurai forms and variations of RChiPT
//  - mirror parameters to run card/decay yaml files
//
// Form factors from:
// - KS, 1: (Kuehn-Santamaria model):
//   * pi pi (original version): Z.Phys.C 48 (1990) 445-452
//     (https://doi.org/10.1007/BF01572024)
//   * pi pi (TODO: KS with Gounaris-Sakurai form factors):
//     from Belle measurement & analysis:
//     (https://arxiv.org/pdf/0805.3773) 
//   * K pi: Z.Phys.C 69 (1996) 243 (vector form factor)
//     (https://doi.org/10.1007/s002880050024)
// - RChT, 2: Resonance Chiral Perturbation Theory
//   * pi pi/KK fit: Eur.Phys.J.C 79 (2019) 5, 436
//     (https://doi.org/10.1140/epjc/s10052-019-6943-9)
//           - We could also have some alternative dispersive approach, 
//             maybe as form factor model 3, based on RchT
//             (https://arxiv.org/pdf/1112.0962)
//
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////