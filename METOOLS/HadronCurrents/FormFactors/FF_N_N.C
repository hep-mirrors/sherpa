#include "METOOLS/HadronCurrents/FormFactors/FF_N_N.H"
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

FF_N_N::FF_N_N(const FF_Parameters & params) :
  FormFactor_Base(params), m_vertex(FF_N_N_vertex::unknown)
{
  if      (m_name=="F1") m_fftype = FFtype::F1;
  else if (m_name=="F2") m_fftype = FFtype::F2;
  else if (m_name=="FA") m_fftype = FFtype::FA;
  else if (m_name=="FP") m_fftype = FFtype::FP;
  else THROW(fatal_error,"FF_N_N has no form factor '"+m_name+"'; use "
             "F1, F2, FA or FP.");
  ValidateModel();

  GeneralModel empty;
  const GeneralModel & md = p_model ? *p_model : empty;
  FixMode(md);
  m_mN  = 0.5*(m_masses[m_pi[0]]+m_masses[m_pi[1]]);
  m_mpi = Flavour(kf_pi).HadMass(); // make this available from the runcard for specific pion masses
  m_s2  = std::abs(MODEL::s_model->ComplexConstant(string("csin2_thetaW")));
  m_T3  = (m_vertex==FF_N_N_vertex::Z_p) ? 0.5 : -0.5;
  m_gA  = md("GA",1.2754);

  // p_model dangles once the owning current's setup returns, as in FF_P_X.
  p_model = NULL;
  msg_Tracking()<<"### FF_N_N '"<<m_name<<"': model "<<m_ffmodel
                <<", F(0) = "<<FF_Assigner(0.)<<", F(1 GeV^2) = "<<FF_Assigner(1.)<<"\n";
  
  CreateF1F2Histogram();
  CreateAxialRatioHistogram();
}

void FF_N_N::FixMode(const GeneralModel & md)
{
  const kf_code kf0(m_flavs[m_pi[0]].Kfcode()), kf1(m_flavs[m_pi[1]].Kfcode());
  const bool cc(kf0!=kf1), isp(kf0==kf_p_plus);
  const int  boson(std::abs(int(md("BOSON", cc ? 24. : 22.))));
  if      (boson==22 && !cc) m_vertex = isp ? FF_N_N_vertex::gamma_p : FF_N_N_vertex::gamma_n;
  else if (boson==23 && !cc) m_vertex = isp ? FF_N_N_vertex::Z_p : FF_N_N_vertex::Z_n;
  else if (boson==24 &&  cc) m_vertex = FF_N_N_vertex::W_pn;
  else THROW(fatal_error,"FF_N_N: BOSON "+std::to_string(boson)+
             " does not match the nucleon legs.");
}

////////////////////////////////////////////////////////////////////////////
// Vector models. G_E = F1 - tau F2, G_M = F1 + F2.
////////////////////////////////////////////////////////////////////////////


void FF_N_N::Dipole_p(const double & Q2,double & F1,double & F2) const
{
  // Standard dipole + Galster baseline 
  // S. Galster, H. Klein, J. Moritz, K. H. Schmidt, D. Wegener and J. Bleckwenn, 
  // Nucl. Phys. B32 (1971) 221–237, doi:10.1016/0550-3213(71)90068-X.
  const double mu_p(2.793);
  const double tau(Q2/(4.*sqr(m_mN)));
  const double GD(std::pow((1 + (Q2/0.71)), -2));
  const double GE(GD);
  const double GM(mu_p*GD);
  F1 = (GE+tau*GM)/(1.+tau);
  F2 = (GM-GE)/(1.+tau);
}


void FF_N_N::Dipole_n(const double & Q2,double & F1,double & F2) const
{
  // Standard dipole + Galster baseline 
  // S. Galster, H. Klein, J. Moritz, K. H. Schmidt, D. Wegener and J. Bleckwenn, 
  // Nucl. Phys. B32 (1971), doi:10.1016/0550-3213(71)90068-X.
  const double mu_n(-1.913);
  const double tau(Q2/(4.*sqr(m_mN)));
  const double GD(std::pow((1 + (Q2/0.71)), -2));
  const double GM(mu_n*GD);
  const double GE(-mu_n*GD*tau/(1.+5.6* tau));
  F1 = (GE+tau*GM)/(1.+tau);
  F2 = (GM-GE)/(1.+tau);
}


void FF_N_N::Kelly_p(const double & Q2,double & F1,double & F2) const
{
  // Kelly, PRC 70 (2004) 068202, doi:10.1103/PhysRevC.70.068202
  const double mu_p(2.793);
  const double tau(Q2/(4.*sqr(m_mN)));
  const double GE((1.-0.24*tau)/(1.+tau*(10.98+tau*(12.82+tau*21.97))));
  const double GM(mu_p*(1.+0.12*tau)/(1.+tau*(10.97+tau*(18.86+tau*6.55))));
  F1 = (GE+tau*GM)/(1.+tau);
  F2 = (GM-GE)/(1.+tau);
}

void FF_N_N::Kelly_n(const double & Q2,double & F1,double & F2) const
{
  // Kelly, PRC 70 (2004) 068202; G_E^n Galster-like A tau/(1+B tau) G_D
  const double mu_n(-1.913), A(1.7), B(3.3), L2(0.7174);
  const double tau(Q2/(4.*sqr(m_mN)));
  const double GE(A*tau/(1.+B*tau)/sqr(1.+Q2/L2));
  const double GM(mu_n*(1.+2.33*tau)/(1.+tau*(14.72+tau*(24.2+tau*84.1))));
  F1 = (GE+tau*GM)/(1.+tau);
  F2 = (GM-GE)/(1.+tau);
}


void FF_N_N::BBBA05_p(const double & Q2,double & F1,double & F2) const
{
  // Bradford, R. and Bodek, A. and Budd, H. and Arrington, J.
  // Nuclear Physics B 159 (2006)
  // dio: 10.1016/j.nuclphysbps.2006.08.028},
  const double mu_p(2.793);
  const double tau(Q2/(4.*sqr(m_mN)));
  const double GE((1.-0.0578*tau)/(1.+tau*(11.1+tau*(13.6+tau*33.0))));
  const double GM(mu_p*(1.+0.150*tau)/(1.+tau*(11.1+tau*(19.6+tau*7.54))));
  F1 = (GE+tau*GM)/(1.+tau);
  F2 = (GM-GE)/(1.+tau);
}

void FF_N_N::BBBA05_n(const double & Q2,double & F1,double & F2) const
{
  // Bradford, R. and Bodek, A. and Budd, H. and Arrington, J.
  // Nuclear Physics B 159 (2006)
  // dio: 10.1016/j.nuclphysbps.2006.08.028},
  const double mu_n(-1.913);
  const double tau(Q2/(4.*sqr(m_mN)));
  const double GE((tau*(1.25+tau*1.3))/(1.+tau*(-9.86+tau*(305+tau*(-758+tau*802)))));
  const double GM(mu_n*(1.+1.81*tau)/(1.+tau*(14.1+tau*(20.7+tau*68.7))));
  F1 = (GE+tau*GM)/(1.+tau);
  F2 = (GM-GE)/(1.+tau);
}


void FF_N_N::Ye2018_p(const double & Q2,double & F1,double & F2) const
{
  // Ye, Zhihong and Arrington, John and Hill, Richard J. and Lee, Gabriel
  // Physics Letters B 777 (2018)
  // dio: 10.1016/j.physletb.2017.11.023
  const double mu_p(2.793);
  const double tau(Q2/(4.*sqr(m_mN)));

  const double tc(0.0779191396), t0(-0.7);
  const std::vector<double> ak{0.239163298067,  -1.10985857441,  1.44438081306,  
    0.479569465603,  -2.28689474187,  1.12663298498,  1.25061984354,  -3.63102047159,  
    4.08221702379,  0.504097346499,  -5.08512046051,  3.96774254395,  -0.981529071103};
  const std::vector<double> bk{0.264142994136, -1.09530612212, 1.21855378178, 
    0.661136493537, -1.40567892503, -1.35641843888, 1.44702915534, 4.2356697359, 
    -5.33404565341, -2.91630052096, 8.70740306757, -5.70699994375, 1.28081437589};

    auto GEp = [&](const double & Qsq) {
    const double A(sqrt(Max(0.,tc+Qsq))), B(sqrt(tc-t0)), z((A-B)/(A+B));
    double zk(1.), sum(0.);
    for (const double & a : ak) { sum += a*zk; zk *= z; }
    return sum;
    };
    auto GMp = [&](const double & Qsq) {
    const double A(sqrt(Max(0.,tc+Qsq))), B(sqrt(tc-t0)), z((A-B)/(A+B));
    double zk(1.), sum(0.);
    for (const double & b : bk) { sum += b*zk; zk *= z; }
    return mu_p*sum;
    };

    const double GE(GEp(Q2));
    const double GM(GMp(Q2));
    
  F1 = (GE+tau*GM)/(1.+tau);
  F2 = (GM-GE)/(1.+tau);
}

void FF_N_N::Ye2018_n(const double & Q2,double & F1,double & F2) const
{
  // Ye, Zhihong and Arrington, John and Hill, Richard J. and Lee, Gabriel
  // Physics Letters B 777 (2018)
  // dio: 10.1016/j.physletb.2017.11.023
  const double mu_n(-1.913);
  const double tau(Q2/(4.*sqr(m_mN)));

  const double tc(0.0779191396), t0(-0.7);
  const std::vector<double> ak{0.048919981379,-0.064525053912,-0.240825897382,
    0.392108744873, 0.300445258602,-0.661888687179,-0.175639769687, 0.624691724461,
    -0.077684299367,-0.236003975259, 0.090401973470, 0.0, 0.0};
  const std::vector<double> bk{0.257758326959,-1.079540642058, 1.182183812195,
    0.711015085833,-1.348080936796,-1.662444025208, 2.624354426029, 1.751234494568,
    -4.922300878888, 3.197892727312,-0.712072389946, 0.0, 0.0};

    auto GEn = [&](const double & Qsq) {
    const double A(sqrt(Max(0.,tc+Qsq))), B(sqrt(tc-t0)), z((A-B)/(A+B));
    double zk(1.), sum(0.);
    for (const double & a : ak) { sum += a*zk; zk *= z; }
    return sum;
    };
    auto GMn = [&](const double & Qsq) {
    const double A(sqrt(Max(0.,tc+Qsq))), B(sqrt(tc-t0)), z((A-B)/(A+B));
    double zk(1.), sum(0.);
    for (const double & b : bk) { sum += b*zk; zk *= z; }
    return mu_n*sum;
    };

    const double GE(GEn(Q2));
    const double GM(GMn(Q2));
    
  F1 = (GE+tau*GM)/(1.+tau);
  F2 = (GM-GE)/(1.+tau);
}


////////////////////////////////////////////////////////////////////////////
// Axial models. F_A = g_A Ftilde(Q^2)/Ftilde(0), positive-F_A convention.
//        - ensure you normalize to F_A to g_A !!!! otherwise you may end up with
//          negative F_A or with differeing normalization to Delta's and N* 
////////////////////////////////////////////////////////////////////////////

double FF_N_N::AxialDipole(const double & Q2) const
{
  // https://doi.org/10.1146/annurev-nucl-102115-044720
  const double MA(1.03); 
  auto F = [MA](const double & Q2) {return std::pow(1. + Q2 / sqr(MA), -2);};
  return m_gA * F(Q2) / F(0.);
}


double FF_N_N::BBBA07(const double & Q2) const
{
  // A. Bodek, S. Avvakumov, R. Bradford and H. Budd, 
  // "Vector and Axial Nucleon Form Factors: A Duality Constrained Parameterization,” European Physical Journal C 53, 349–354 (2008).
  //dio:10.1140/epjc/s10052-007-0491-4
  const double MA(1.015);
  const std::vector<double> p{1.0000, 0.9207, 0.9795, 1.0480, 1.0516, 1.2874, 0.7707 };
  const std::vector<double> nodes{0.0, 1.0/6.0, 2.0/6.0, 3.0/6.0, 4.0/6.0, 5.0/6.0, 1.0};
  auto F = [&](const double & Qsq) {
    double xi(0.0);
    if (Qsq > 0.0) {
      const double tau(Qsq/(4.0*sqr(m_mN)));
      xi = 2.0/(1.0 + std::sqrt(1.0 + 1.0/tau));}
    double AFA(0.0);
    for (int j = 0; j < 7; ++j) {
      double Lj(1.0);
      for (int k = 0; k < 7; ++k) {
        if (k != j) {
          Lj *= (xi - nodes[k])/(nodes[j] - nodes[k]);
        }
      }
      AFA += p[j]*Lj;
    }
    const double dipole(std::pow(1.0 + Qsq/sqr(MA), -2));
    return AFA*dipole;
  };
  return m_gA*F(Q2)/F(0.);
}


double FF_N_N::Meyer(const double & Q2) const
{
  // Meyer, Betancourt, Gran, Hill, PRD 93 (2016) 113015, arXiv:1603.03048
  const double tc(9.*sqr(Flavour(kf_pi_plus).HadMass())), t0(-0.28);
  const std::vector<double> ak{0.759,-2.30,0.6,3.8,-2.3,-2.16,0.896,1.58,-0.823};
  auto F = [&](const double & Qsq) {
    const double A(sqrt(Max(0.,tc+Qsq))), B(sqrt(tc-t0)), z((A-B)/(A+B));
    double zk(1.), sum(0.);
    for (const double & a : ak) { sum += a*zk; zk *= z; }
    return sum;
  };
  return m_gA*F(Q2)/F(0.);
}


double FF_N_N::MINERvA(const double & Q2) const
{
  // T. Cai et al. (MINERvA Collaboration), Nature 614 (2023) 48–53, doi:10.1038/s41586-022-05478-3.
  // values from supplementary material
  const double tc(9.*sqr(Flavour(kf_pi_plus).HadMass())), t0(-0.75);
  const std::vector<double> ak{0.50,-1.50,1.2,0.1,-0.2,-0.46,0.4,-0.15,0.044};
  auto F = [&](const double & Qsq) {
    const double A(sqrt(Max(0.,tc+Qsq))), B(sqrt(tc-t0)), z((A-B)/(A+B));
    double zk(1.), sum(0.);
    for (const double & a : ak) { sum += a*zk; zk *= z; }
    return sum;
  };
  return m_gA*F(Q2)/F(0.);
}


double FF_N_N::Lattice(const double & Q2) const
{
  // lattice-QCD average z-expansion
  const double tc(9.*sqr(0.134)), t0(-0.50);
  const std::vector<double> ak{0.72115656,-1.70104640,0.26324902,1.53433681,
                               0.01061114,-1.49893610,0.67062898};
  auto F = [&](const double & Qsq) {
    const double A(sqrt(Max(0.,tc+Qsq))), B(sqrt(tc-t0)), z((A-B)/(A+B));
    double zk(1.), sum(0.);
    for (const double & a : ak) { sum += a*zk; zk *= z; }
    return sum;
  };
  return m_gA*F(Q2)/F(0.);
}


double FF_N_N::DjukanovicLattice(const double & Q2) const
{
  // Djukanovic, D., von Hippel, G., Koponen, J., Meyer, H. B., Ottnad, K., Schulz, T. and Wittig, H. (2022). 
  // “Isovector axial form factor of the nucleon from lattice QCD.” Physical Review D, 106, 074503.
  // DOI: 10.1103/PhysRevD.106.074503
  const double tc(sqr(3.0 * Flavour(kf_pi).HadMass())), t0(0);
  const std::vector<double> ak{1.225, -1.274, -0.379};
  auto F = [&](const double & Qsq) {
    const double A(sqrt(Max(0.,tc+Qsq))), B(sqrt(tc-t0)), z((A-B)/(A+B));
    double zk(1.), sum(0.);
    for (const double & a : ak) { sum += a*zk; zk *= z; }
    return sum;
  };
  return m_gA*F(Q2)/F(0.);
}

////////////////////////////////////////////////////////////////////////////
// Model switches: the only place a new model number is wired in.
////////////////////////////////////////////////////////////////////////////
void FF_N_N::Vector(const double & Q2,double & F1p,double & F2p,
                    double & F1n,double & F2n) const
{
  F1p = F2p = F1n = F2n = 0.;
  switch (GetModel()) {
    case Model::Point: F1p = 1.; return;
    case Model::Dipole: Dipole_p(Q2,F1p,F2p); Dipole_n(Q2,F1n,F2n); return;
    case Model::Kelly: Kelly_p(Q2,F1p,F2p); Kelly_n(Q2,F1n,F2n); return;
    case Model::BBBA05: BBBA05_p(Q2,F1p,F2p); BBBA05_n(Q2,F1n,F2n); return;
    case Model::Ye2018: Ye2018_p(Q2,F1p,F2p); Ye2018_n(Q2,F1n,F2n); return;
    default: return;
  }
}

double FF_N_N::Axial(const double & Q2) const
{
  switch (GetModel()) {
    case Model::Point:   return m_gA;
    case Model::AxialDipole:   return AxialDipole(Q2);
    case Model::Meyer:   return Meyer(Q2);
    case Model::BBBA07:   return BBBA07(Q2);
    case Model::MINERvA:   return MINERvA(Q2);
    case Model::Lattice: return Lattice(Q2);
    case Model::DjukanovicLattice: return DjukanovicLattice(Q2);
    default:             return 0.;
  }
}


//histogram


void FF_N_N::CreateF1F2Histogram() const
{
  const int nBins = 1000;
  const double Q2Max = 2.0;

  // Convert the selected model into a filename-safe string.
  const char* modelName = "Unknown";

  switch (GetModel()) {
    case Model::Point:   modelName = "Point";   break;
    case Model::Dipole:  modelName = "Dipole";  break;
    case Model::Kelly:   modelName = "Kelly";   break;
    case Model::BBBA05:   modelName = "BBBA05";   break;
    case Model::Ye2018:   modelName = "Ye2018";   break;

    case Model::AxialDipole:   modelName = "AxialDipole";   break;
    case Model::Meyer:   modelName = "Meyer";   break;
    case Model::BBBA07:   modelName = "BBBA07";   break;
    case Model::MINERvA:   modelName = "MINERvA";   break;
    case Model::Lattice: modelName = "Lattice"; break;
    case Model::DjukanovicLattice: modelName = "DjukanovicLattice"; break;
    default:             modelName = "Unknown"; break;
  }

  // Vector form factors: F1 and F2 for proton and neutron.
  const bool vectorType =
      m_fftype == FFtype::F1 || m_fftype == FFtype::F2;

  const bool vectorModel =
      GetModel() == Model::Point ||
      GetModel() == Model::Dipole ||
      GetModel() == Model::Kelly ||
      GetModel() == Model::BBBA05 ||
      GetModel() == Model::Ye2018;

      if (vectorType && vectorModel) {
        // Existing F1/F2 files
        std::ofstream protonFile(
            std::string("Results/FF_N_N_proton_") + modelName + ".dat");
        std::ofstream neutronFile(
            std::string("Results/FF_N_N_neutron_") + modelName + ".dat");
    
        // New GE/GM files
        std::ofstream protonSachsFile(
            std::string("Results/FF_N_N_proton_GE_GM_") + modelName + ".dat");
        std::ofstream neutronSachsFile(
            std::string("Results/FF_N_N_neutron_GE_GM_") + modelName + ".dat");
    
        // Check that all files opened successfully
        if (!protonFile || !neutronFile ||
            !protonSachsFile || !neutronSachsFile) {
          msg_Error() << "FF_N_N: Could not create vector form-factor files."
                      << std::endl;
          return;
        }
    
        // Column headers
        protonFile << "# Q2_GeV2 F1p F2p\n";
        neutronFile << "# Q2_GeV2 F1n F2n\n";
    
        protonSachsFile << "# Q2_GeV2 GEp GMp\n";
        neutronSachsFile << "# Q2_GeV2 GEn GMn\n";
    
        // Generate form factors over the Q2 range
        for (int i = 0; i <= nBins; ++i) {
          const double Q2 = Q2Max * i / nBins;
          const double tau = Q2 / (4. * sqr(m_mN));
    
          double F1p, F2p, F1n, F2n;
          Vector(Q2, F1p, F2p, F1n, F2n);
    
          // Write F1 and F2
          protonFile << Q2 << " " << F1p << " " << F2p << "\n";
          neutronFile << Q2 << " " << F1n << " " << F2n << "\n";
    
          // Convert to Sachs form factors
          const double GEp = F1p - tau * F2p;
          const double GMp = F1p + F2p;
    
          const double GEn = F1n - tau * F2n;
          const double GMn = F1n + F2n;
    
          // Write GE and GM
          protonSachsFile << Q2 << " " << GEp << " " << GMp << "\n";
          neutronSachsFile << Q2 << " " << GEn << " " << GMn << "\n";
        }
      }
  


  // Axial form factors: FA or FP, depending on the selected type.
  const bool axialType =
      m_fftype == FFtype::FA || m_fftype == FFtype::FP;

  const bool axialModel =
      GetModel() == Model::Point ||
      GetModel() == Model::AxialDipole ||
      GetModel() == Model::Meyer ||
      GetModel() == Model::BBBA07 ||
      GetModel() == Model::MINERvA ||
      GetModel() == Model::Lattice ||
      GetModel() == Model::DjukanovicLattice;
      

  if (axialType && axialModel) {
    const bool pseudoscalar = (m_fftype == FFtype::FP);

    const std::string filename = std::string("Results/FF_N_N_") + (pseudoscalar ? "FP_" : "FA_") + modelName + ".dat";

    std::ofstream axialFile(filename);

    if (!axialFile) {
      msg_Error() << "FF_N_N: Could not create axial form-factor file."
                  << std::endl;
      return;
    }

    axialFile << (pseudoscalar
                  ? "# Q2_GeV2 FP\n"
                  : "# Q2_GeV2 FA\n");

    for (int i = 0; i <= nBins; ++i) {
      const double Q2 = Q2Max * i / nBins;
      double value = Axial(Q2);

      if (pseudoscalar) {
        value *= 4. * sqr(m_mN) / (Q2 + sqr(m_mpi));
      }

      axialFile << Q2 << " " << value << "\n";
    }
  }

  msg_Info() << "FF_N_N: Form-factor output written for model "
             << modelName << std::endl;
}


//testing BBBA07
void FF_N_N::CreateAxialRatioHistogram() const
{
  if (m_fftype != FFtype::FA ||
      GetModel() != Model::BBBA07) return;

  const double MA(1.015);
  const double MN(m_mN);

  std::ofstream output("Results/FF_N_N_FA_over_GAD_BBBA07.dat");

  if (!output.is_open()) {
    THROW(fatal_error,"Could not open FF_N_N_FA_over_GAD_BBBA07.dat");
  }

  output << "# Q2_GeV2 xi FA_over_GAD\n";

  const int nBins(1000);
  const double Q2max(50.0);

  for (int i = 0; i <= nBins; ++i) {
    // Concentrate more bins at low Q2 to improve spacing near xi = 0.
    const double x(static_cast<double>(i) / nBins);
    const double Q2(Q2max * x * x);

    double xi(0.0);

    if (Q2 > 0.0) {
      const double tau(Q2 / (4.0 * sqr(MN)));
      xi = 2.0 / (1.0 + std::sqrt(1.0 + 1.0 / tau));
    }

    const double GAD(std::pow(1.0 + Q2 / sqr(MA), -2));
    const double FA(Axial(Q2));
    const double ratio(FA / (m_gA * GAD));

    output << Q2 << " "
           << xi << " "
           << ratio << "\n";
  }

  output.close();
}



////////////////////////////////////////////////////////////////////////////
// Boson combinations: this should be the same for any form factor models
////////////////////////////////////////////////////////////////////////////
double FF_N_N::FF_Assigner(const double & Q2) const
{
  // axial
  if (m_fftype==FFtype::FA || m_fftype==FFtype::FP) {
    double FA(Axial(Q2));
    if      (m_vertex==FF_N_N_vertex::gamma_p || m_vertex==FF_N_N_vertex::gamma_n) FA = 0.;
    else if (m_vertex==FF_N_N_vertex::Z_p     || m_vertex==FF_N_N_vertex::Z_n)     FA *= m_T3;
    return m_fftype==FFtype::FA ? FA : 4.*sqr(m_mN)/(Q2+sqr(m_mpi))*FA; // currently using dipole
  }
  // vector
  double F1p, F2p, F1n, F2n;
  Vector(Q2,F1p,F2p,F1n,F2n);
  const double Fp(m_fftype==FFtype::F1 ? F1p : F2p), Fn(m_fftype==FFtype::F1 ? F1n : F2n);
  switch (m_vertex) {
    case FF_N_N_vertex::gamma_p: return Fp;
    case FF_N_N_vertex::gamma_n: return Fn;
    case FF_N_N_vertex::Z_p:     return m_T3*(Fp-Fn)-2.*m_s2*Fp;
    case FF_N_N_vertex::Z_n:     return m_T3*(Fp-Fn)-2.*m_s2*Fn;
    case FF_N_N_vertex::W_pn:    return Fp-Fn;
    default:                   return 0.;
  }
}






Complex FF_N_N::operator()(const ATOOLS::Vec4D_Vector& moms)
{
  const Vec4D q(moms[m_pi[0]]-moms[m_pi[1]]);
  return Complex(FF_Assigner(-q.Abs2()),0.);
}

DECLARE_FF_GETTER(FF_N_N,"FF_N_N")

FormFactor_Base * ATOOLS::Getter<FormFactor_Base,FF_Parameters,FF_N_N>::
operator()(const METOOLS::FF_Parameters &params) const
{
  if (params.m_pi.size()!=2) return NULL;
  return new FF_N_N(params);
}





