#include "METOOLS/HadronCurrents/FormFactors/FF_N_NPi.H"
#include "METOOLS/HadronCurrents/FormFactors/FF_N_N.H"
#include "METOOLS/HadronCurrents/Tools.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

FF_N_NPi::FF_N_NPi(const FF_Parameters & params) :
  FormFactor_Base(params), p_NN(NULL)
{
  if      (m_name=="FpiNN")   m_fftype = FFtype::FpiNN;
  else if (m_name=="FPiPole") m_fftype = FFtype::FPiPole;
  else THROW(fatal_error,"FF_N_NPi has no form factor '"+m_name+"'; use "
             "FpiNN or FPiPole.");
  ValidateModel();

  GeneralModel empty;
  const GeneralModel & md = p_model ? *p_model : empty;
  m_mN    = 0.5*(m_masses[m_pi[0]]+m_masses[m_pi[1]]);
  m_gA    = md("GA",1.2754);
  m_fpi   = md("fpi",0.093);
  m_gpiNN = md("gpiNN",m_gA*m_mN/m_fpi);  // Goldberger-Treiman
  m_mrho2 = sqr(Flavour(kf_rho_770).HadMass());

  // PCAC models reuse FF_N_N's axial models: g_piNN(x) = (m_N/f_pi) F_A(x).
  if (m_fftype==FFtype::FpiNN && m_ffmodel>=20) {
    std::map<std::string,double> pmap;
    FF_Parameters fa(m_ffmodel,params.m_flavs,m_pi,pmap,"FA",p_model);
    p_NN = new FF_N_N(fa);
  }
  // F1V reuses FF_N_N's vector models: F1^p - F1^n (BOSON 24), VECTOR_FF.
  if (GetModel()==Model::F1V) {
    std::map<std::string,double> pmap;
    FF_Parameters fv(int(md("VECTOR_FF",10.)+0.5),params.m_flavs,m_pi,pmap,"F1",p_model);
    p_NN = new FF_N_N(fv);
  }

  // p_model dangles once the owning current's setup returns, as in FF_P_X.
  p_model = NULL;
  msg_Tracking()<<"### FF_N_NPi '"<<m_name<<"': model "<<m_ffmodel
                <<", F(0) = "<<FF_Assigner(0.)<<"\n";
}

FF_N_NPi::~FF_N_NPi()
{
  if (p_NN) delete p_NN;
}

void FF_N_NPi::ValidateModel() const
{
  switch (GetModel()) {
    case Model::None: case Model::Point:
      return;
    case Model::PionDipole:
    case Model::AxialDipole:
    case Model::Meyer:
    case Model::BBBA07:
    case Model::MINERvA:
    case Model::Lattice:
    case Model::DjukanovicLattice:
      if (m_fftype==FFtype::FpiNN) return;
      break;
    case Model::Rho:
    case Model::F1V:
      if (m_fftype==FFtype::FPiPole) return;
      break;
  }
  THROW(fatal_error, "Unknown FORM_FACTOR=" + std::to_string(m_ffmodel) +
        " for FF_N_NPi '" + m_name + "'");
}

////////////////////////////////////////////////////////////////////////////
// Standalone models
////////////////////////////////////////////////////////////////////////////

double FF_N_NPi::PionDipole(const double & x) const
{
  return m_gpiNN*std::pow(1.+x/0.934,-2);
}

double FF_N_NPi::Rho(const double & x) const
{
  // rho dominance, implementation_1.tex, m_rho from Sherpa
  return m_mrho2/(m_mrho2+x);
}

////////////////////////////////////////////////////////////////////////////
// Model switch: the only place a new model number is wired in.
////////////////////////////////////////////////////////////////////////////

double FF_N_NPi::FF_Assigner(const double & x) const
{
  switch (GetModel()) {
    case Model::None:       return 0.;
    case Model::Point:      return m_fftype==FFtype::FpiNN ? m_gpiNN : 1.;
    case Model::PionDipole: return PionDipole(x);
    case Model::Rho:        return Rho(x);
    case Model::F1V:        return p_NN->FF_Assigner(x);
    default:                return m_mN/m_fpi*p_NN->FF_Assigner(x); // PCAC
  }
}

Complex FF_N_NPi::operator()(const ATOOLS::Vec4D_Vector& moms)
{
  const Vec4D q(moms[m_pi[0]]-moms[m_pi[1]]);
  return Complex(FF_Assigner(-q.Abs2()),0.);
}

DECLARE_FF_GETTER(FF_N_NPi,"FF_N_NPi")

FormFactor_Base * ATOOLS::Getter<FormFactor_Base,FF_Parameters,FF_N_NPi>::
operator()(const METOOLS::FF_Parameters &params) const
{
  if (params.m_pi.size()!=2) return NULL;
  return new FF_N_NPi(params);
}
