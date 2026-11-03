#include "METOOLS/HadronCurrents/FormFactors/FF_N_Nprime.H"
#include "METOOLS/HadronCurrents/Tools.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "MODEL/Main/Model_Base.H"
#include <fstream>
#include <string>

using namespace METOOLS;
using namespace ATOOLS; 
using namespace std;

FF_N_Nprime::FF_N_Nprime(const FF_Parameters & params) :
  FormFactor_Base(params), m_vertex(FF_N_Nprime_vertex::unknown)
{
  if      (m_name=="F1prime") m_fftype = FFtype::F1prime;
  else if (m_name=="F2prime") m_fftype = FFtype::F2prime;
  else THROW(fatal_error,"FF_N_Nprime has no form factor '"+m_name+"'; use "
             "F1prime, F2prime.");
  ValidateModel();

  GeneralModel empty;
  const GeneralModel & md = p_model ? *p_model : empty;
  FixMode(md);
  m_N  = m_masses[m_pi[1]];
  m_R  = m_masses[m_pi[0]];

  // p_model dangles once the owning current's setup returns, as in FF_P_X.
  p_model = NULL;
  msg_Tracking()<<"### FF_N_Nprime '"<<m_name<<"': model "<<m_ffmodel
                <<", F(0) = "<<FF_Assigner(0.)<<", F(1 GeV^2) = "<<FF_Assigner(1.)<<"\n";
  
  CreateF1F2Histogram();
}


void FF_N_Nprime::FixMode(const GeneralModel & md)
{
  const kf_code kf0(m_flavs[m_pi[0]].Kfcode());
  const kf_code kf1(m_flavs[m_pi[1]].Kfcode());

  const bool proton_in = (kf1 == kf_p_plus);
  const bool neutron_in = (kf1 == kf_n);
  const bool N1440_out = (kf0 == kf_N_1440 || kf0 == kf_N_1440_plus);
  const bool N1535_out = (kf0 == kf_N_1535 || kf0 == kf_N_1535_plus);

  const bool proton_transition = proton_in && (kf0 == kf_N_1440_plus || kf0 == kf_N_1535_plus);
  const bool neutron_transition = neutron_in && (kf0 == kf_N_1440 || kf0 == kf_N_1535);

  const int boson(std::abs(int(md("BOSON", 22.))));
  if (boson != 22)
    THROW(fatal_error,
          "FF_N_Nprime: only electromagnetic transitions with BOSON 22 "
          "are supported.");

  if (proton_transition)
    m_vertex = FF_N_Nprime_vertex::gamma_p;
  else if (neutron_transition)
    m_vertex = FF_N_Nprime_vertex::gamma_n;
  else
    THROW(fatal_error,
          "FF_N_Nprime: unsupported nucleon-to-resonance transition "
          "for the supplied particle legs.");
}

////////////////////////////////////////////////////////////////////////////
// Vector models. 1440
////////////////////////////////////////////////////////////////////////////

void FF_N_Nprime::MAID2009_1440_p(const double & Q2,double & F1,double & F2) const
{
  // Add ref

  const double x(Q2/1); //x = Q^2/(1 GEV^2)
  const double alpha(MODEL::s_model->ScalarConstant("alpha_QED"));
  const double Qp2(sqr(m_R+m_N)+Q2);
  const double Qm2(sqr(m_R - m_N)+Q2);
  const double k((sqrt(Qm2 * Qp2)/(2*m_R)));
  const double K((sqr(m_R)-sqr(m_N))/(2*m_R));
  const double tau(Q2/sqr(m_R + m_N));
  const double R(sqrt((M_PI*alpha*Qm2)/(m_R*m_N*K)));
  const double D((m_R + m_N)*k/(sqrt(2)*Q2));

  const double A(-61.4e-3 *(1.0 + 0.871*x - 3.52*sqr(x) - 0.158*std::pow(x, 4)) * std::exp(-1.36*x));
  const double S(4.2e-3 *(1.0 + 40.0*x + 1.5*std::pow(x, 4)) * std::exp(-1.75*x));

  const double s(S/(R*D));
  const double a(A/R);
  
  F1 = (s+tau*a)/(1.+tau);
  F2 = (a-s)/(1.+tau);
}


void FF_N_Nprime::MAID2009_1440_n(const double & Q2,double & F1,double & F2) const
{
  // add ref

  const double x(Q2/1); //x = Q^2/(1 GEV^2)
  const double alpha(MODEL::s_model->ScalarConstant("alpha_QED"));
  const double Qp2(sqr(m_R+m_N)+Q2);
  const double Qm2(sqr(m_R - m_N)+Q2);
  const double k((sqrt(Qm2 * Qp2)/(2*m_R)));
  const double K((sqr(m_R)-sqr(m_N))/(2*m_R));
  const double tau(Q2/sqr(m_R + m_N));
  const double R(sqrt((M_PI*alpha*Qm2)/(m_R*m_N*K)));
  const double D((m_R + m_N)*k/(sqrt(2)*Q2));

  const double A(54.1e-3 *(1.0 + 0.95*x) * std::exp(-1.77*x));
  const double S(-41.5e-3 *(1.0 + 2.98*x) * std::exp(-1.55*x));

  const double s(S/(R*D));
  const double a(A/R);
  
  F1 = (s+tau*a)/(1.+tau);
  F2 = (a-s)/(1.+tau);
}



////////////////////////////////////////////////////////////////////////////
// Vector models. 1440
////////////////////////////////////////////////////////////////////////////

void FF_N_Nprime::MAID2009_1535_p(const double & Q2,double & F1,double & F2) const
{
  // Add ref

  const double x(Q2/1); //x = Q^2/(1 GEV^2)
  const double alpha(MODEL::s_model->ScalarConstant("alpha_QED"));
  const double e(sqrt(4*M_PI*alpha));
  const double Qp2(sqr(m_R+m_N)+Q2);
  const double Qm2(sqr(m_R-m_N)+Q2);
  const double k((sqrt(Qm2 * Qp2)/(2*m_R)));
  const double d((m_R-m_N)/(m_R+m_N));
  const double b(e*sqrt(Qp2/(8*m_N*(sqr(m_R)-sqr(m_N)))));
  const double tau(Q2/sqr(m_R + m_N));
  const double D((m_R + m_N)*sqrt(2)*k/Q2);

  const double A(66.4e-3 *(1.0 + 1.608*x) * std::exp(-0.7*x));
  const double S(-2.0e-3 *(1.0 + 23.9*x) * std::exp(-0.81*x));

  const double s(S/(b*D));
  const double a(-A/(2*b));
  
  F1 = (d*s+tau*a)/(sqr(d)+tau);
  F2 = (d*a-s)/(sqr(d)+tau);
}


void FF_N_Nprime::MAID2009_1535_n(const double & Q2,double & F1,double & F2) const
{
  // add ref

  const double x(Q2/1); //x = Q^2/(1 GEV^2)
  const double alpha(MODEL::s_model->ScalarConstant("alpha_QED"));
  const double e(sqrt(4*M_PI*alpha));
  const double Qp2(sqr(m_R+m_N)+Q2);
  const double Qm2(sqr(m_R-m_N)+Q2);
  const double k((sqrt(Qm2 * Qp2)/(2*m_R)));
  const double d((m_R-m_N)/(m_R+m_N));
  const double b(e*sqrt(Qp2/(8*m_N*(sqr(m_R)-sqr(m_N)))));
  const double tau(Q2/sqr(m_R + m_N));
  const double D((m_R + m_N)*sqrt(2)*k/Q2);

  const double A(-50.7e-3 *(1.0 + 4.75*x) * std::exp(-1.69*x));
  const double S(28.5e-3 *(1.0 + 0.36*x) * std::exp(-1.55*x));

  const double s(S/(b*D));
  const double a(-A/(2*b));
  
  F1 = (d*s+tau*a)/(sqr(d)+tau);
  F2 = (d*a-s)/(sqr(d)+tau);
}


////////////////////////////////////////////////////////////////////////////
// Model switches: the only place a new model number is wired in.
////////////////////////////////////////////////////////////////////////////
void FF_N_Nprime::Vector(const double & Q2,double & F1primep,double & F2primep,
                    double & F1primen,double & F2primen) const
{
  F1primep = F2primep = F1primen = F2primen = 0.;
  switch (GetModel()) {
    case Model::Point: F1primep = 1.; return;
    case Model::MAID2009_1440: MAID2009_1440_p(Q2,F1primep,F2primep); MAID2009_1440_n(Q2,F1primen,F2primen); return;
    case Model::MAID2009_1535: MAID2009_1535_p(Q2,F1primep,F2primep); MAID2009_1535_n(Q2,F1primen,F2primen); return;
    default: return;
  }
}



////////////////////////////////////////////////////////////////////////////
// Boson combinations: this should be the same for any form factor models
////////////////////////////////////////////////////////////////////////////
double FF_N_Nprime::FF_Assigner(const double & Q2) const
{
  // vector
  double F1primep, F2primep, F1primen, F2primen;
  Vector(Q2,F1primep, F2primep, F1primen, F2primen);
  const double Fp(m_fftype==FFtype::F1prime ? F1primep : F2primep), Fn(m_fftype==FFtype::F1prime ? F1primen : F2primen);
  switch (m_vertex) {
    case FF_N_Nprime_vertex::gamma_p: return Fp;
    case FF_N_Nprime_vertex::gamma_n: return Fn;
    default:                   return 0.;
  }
}


// histogram


void FF_N_Nprime::CreateF1F2Histogram() const
{
  const int nBins = 1000;
  const double Q2Min = 1.0e-4;
  const double Q2Max = 2.0;

  const char* modelName = "Unknown";

  switch (GetModel()) {
    case Model::Point:    modelName = "Point"; break;
    case Model::MAID2009_1440: modelName = "MAID2009_1440"; break;
    default:
      msg_Error() << "FF_N_Nprime: No histogram implementation for model "
                  << m_ffmodel << std::endl;
      return;
  }

  const int resonanceCode = static_cast<int>(std::round(m_R * 1000.0));

  const std::string prefix =
      std::string("Results/FF_N_Nprime_")
      + std::to_string(resonanceCode) + "_"
      + modelName;

  std::ofstream protonFile(prefix + "_proton.dat");
  std::ofstream neutronFile(prefix + "_neutron.dat");

  if (!protonFile || !neutronFile) {
    msg_Error() << "FF_N_Nprime: Could not create histogram files. "
                << "Check that Results exists." << std::endl;
    return;
  }

  protonFile << "# Q2_GeV2 F1primep F2primep A1/2 S1/2\n";
  neutronFile << "# Q2_GeV2 F1primen F2primen A1/2 S1/2\n";

  const double alpha = MODEL::s_model->ScalarConstant("alpha_QED");
  const double K = (sqr(m_R) - sqr(m_N)) / (2.0 * m_R);

  for (int i = 0; i <= nBins; ++i) {
    const double Q2 = Q2Min + (Q2Max - Q2Min) * i / nBins;

    double F1p, F2p, F1n, F2n;
    Vector(Q2, F1p, F2p, F1n, F2n);

    const double Qp2(sqr(m_R+m_N)+Q2);
    const double Qm2(sqr(m_R - m_N)+Q2);
    const double k((sqrt(Qm2 * Qp2)/(2*m_R)));
    const double tau(Q2/sqr(m_R + m_N));
    const double R(sqrt((M_PI*alpha*Qm2)/(m_R*m_N*K)));
    const double D((m_R + m_N)*k/(sqrt(2)*Q2));

    const double A1p = R * (F1p + F2p);
    const double S1p = R * D * (F1p - tau*F2p);

    const double A1n = R * (F1n + F2n);
    const double S1n = R * D * (F1n - tau*F2n);

    protonFile << Q2 << " " << F1p << " " << F2p << " "
               << A1p << " " << S1p << "\n";

    neutronFile << Q2 << " " << F1n << " " << F2n << " "
                << A1n << " " << S1n << "\n";
  }

  msg_Info() << "FF_N_Nprime: Form-factor and helicity-amplitude "
             << "files written for model " << modelName << std::endl;
}





Complex FF_N_Nprime::operator()(const ATOOLS::Vec4D_Vector& moms)
{
  const Vec4D q(moms[m_pi[0]]-moms[m_pi[1]]);
  return Complex(FF_Assigner(-q.Abs2()),0.);
}

DECLARE_FF_GETTER(FF_N_Nprime,"FF_N_Nprime")

FormFactor_Base * ATOOLS::Getter<FormFactor_Base,FF_Parameters,FF_N_Nprime>::
operator()(const METOOLS::FF_Parameters &params) const
{
  if (params.m_pi.size()!=2) return NULL;
  return new FF_N_Nprime(params);
}









