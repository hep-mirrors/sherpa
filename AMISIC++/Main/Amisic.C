#include "AMISIC++/Main/Amisic.H"
#include "AMISIC++/Perturbative/MI_Processes.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace AMISIC;
using namespace ATOOLS;
using namespace std;

Amisic::Amisic() :
  p_processes(nullptr), p_soft(nullptr), m_sigmaND_norm(1.),
  m_isMinBias(false), m_ana(false)
{}

Amisic::~Amisic() {
  if (p_processes) delete p_processes;
  if (p_soft)      delete p_soft;
  if (m_ana) FinishAnalysis();
}

bool Amisic::Initialize(MODEL::Model_Base *const model,
			PDF::ISR_Handler *const isr,
			YFS::YFS_Handler *const yfs,
                        REMNANTS::Remnant_Handler * remnant_handler)
{
  InitParametersAndType(isr);
  Vec4D P = isr->GetBeam(0)->OutMomentum()+isr->GetBeam(1)->OutMomentum();
  m_S     = P.Abs2();
  m_Y     = P.Y();
  // Calculate hadronic non-diffractive cross sections, the normalization for the
  // multiple scattering probability.
  p_xsecs     = new Hadronic_XSec_Calculator(model,isr->Flav(0),isr->Flav(1));
  // Initialize the parton-level processes - currently only 2->2 scatters.  Even for
  // soft processes we need the Mass_Selector for the construction of the amplitudes
  // in the parton shower.
  p_processes = new MI_Processes(m_variable_s);
  p_processes->SetXSecCalculator(p_xsecs);
  if (mipars->GetEvtType()==evt_type::Perturbative) {
    // Initialize the hard processes (i.e. the 2->2 scatters)
    p_processes->Initialize(model,nullptr,isr,yfs);
    // Initialize the Over_Estimator - mainly fixing an effective prefactor to allow
    // for a quick'n'dirty fix to create fast estimates of the next scatter's pT^2.
    // This follows closely the algorithm in the original Sjostrand - van Zijl
    // publication.  The logic for the overestimator is based on
    // - t-channel dominance allowing to approximate differential cross sections as
    //   dsigma ~ dp_T dy f(x_1) f_(x_2) C_1 C_2 as(t)^2/(t+t_0)^2
    //   where C_1/2 are the colour factors of the incoming partons, x_1/2 are the
    //   Bjorken-x.
    // - assuming that the product of the PDFs f(x_1)f(x_2) is largest for mid-rapidity
    //   where x_1 and x_2 are identical
    m_overestimator.Initialize(p_processes);
    // Initializing the Single_Collision_Handler which creates the next scatter: it needs
    // the remnants, processes, and the overestimator
    m_singlecollision.Init(remnant_handler,p_processes,&m_overestimator);
    // Initialize everything to do with the inpact parameter dependence.
    m_impact.Initialize(remnant_handler,p_processes);
  }
  else {
    p_soft = new NonPerturbative_XSecs(remnant_handler,p_xsecs);
    // Initializing the Single_Collision_Handler which creates the next scatter: it needs
    // the remnants, processes, and the overestimator
    m_singlecollision.Init(remnant_handler,p_soft);
  }

  if (m_ana) InitAnalysis();
  m_isFirst   = true;
  m_isMinBias = false;

  return true;
}

void Amisic::InitParametersAndType(PDF::ISR_Handler *const isr) {
  mipars = new MI_Parameters();
  bool shown = false;
  // Must distinguish the hard and rescatter process.  For the latter, we take energies
  // as fixed, for the former, the energies may vary (we have to check the spectrum):
  // - if EPA is used the energies entering the ISR will vary,
  // - otherwise the energy is fixed.
  m_variable_s = isr->GetBeam(0)->On() || isr->GetBeam(1)->On();
  if (isr->Flav(0).IsHadron() && isr->Flav(1).IsHadron())
    m_type = mitype::hadron_hadron;
  else if ((isr->Flav(0).IsHadron() && isr->Flav(1).IsPhoton()) ||
           (isr->Flav(1).IsHadron() && isr->Flav(0).IsPhoton()))
    m_type = mitype::gamma_hadron;
  else if (isr->Flav(0).IsPhoton() && isr->Flav(1).IsPhoton())
    m_type = mitype::gamma_gamma;
  else
    msg_Error() << METHOD <<": unknown multiple interaction model for " <<
            isr->Flav(0) << " and " << isr->Flav(1) << "\n";
  for (size_t beam=0;beam<2;beam++) {
    if(!shown && sqr((*mipars)("pt_0"))<isr->PDF(beam)->Q2Min()) {
      msg_Error()<<"Potential error in "<<METHOD<<":\n"
		 <<"   IR cutoff of MPI model "<<(*mipars)("pt_0")
		 <<" below minimal scale of PDFs.\n"
		 <<"   Will freeze PDFs at their minimal scale: "
		 <<sqrt(isr->PDF(beam)->Q2Min())<<" GeV.\n";
      shown = true;
    }
  }
}

bool Amisic::InitMPIs(PDF::ISR_Handler *const isr, const double & scale) {
  // Initialise the MPI simulation: fixing the maximal scale for the downward evolution
  // and determining an impact parameter in SetB().
  if (m_isFirst) {
    m_isFirst = false;
    if (m_variable_s) UpdateForNewS();
    SetMassMode(1);
    SetMaxScale(scale);
    SetB();
  }
  if (!VetoEvent(scale)) return true;
  return false;
}

void Amisic::UpdateForNewS() {
  // Update if first scatter with variable c.m. energy (e.g. collisions with pomoerons
  // or resolved photons): update c.m. energy in processes, re-calculate non-diffractive
  // cross sections for normalisation, adjust prefactors for overestimators, etc..
  // TODO: will have to check if we need another longitudinal boost (hence the Y)
  Vec4D P(0.,0.,0.,0.);
  for (size_t beam=0;beam<2;beam++) {
    P += m_singlecollision.InMomentum(beam);
  }
  m_S = P.Abs2();
  m_Y = P.Y();
  m_singlecollision.UpdateSandY(m_S, m_Y);
}

void Amisic::SetB(const double & b) {
  // Generation of the next sctter in the Single_Collision_Handler depends on the
  // impact-parameter enhancement, f(b), Eq. (28), here named m_bfac.
  // It is obtained by interpolation from a look-up table in the Interaction_Probability,
  // accessed through the Impact_Parameter class.
  m_b    = (b<0.) ? m_impact.CalculateB(m_S,m_pt2) : b;
  m_bfac = ATOOLS::Max(0.,m_impact(m_S,m_b));
  m_overestimator.SetBFac(m_bfac);
}


int Amisic::InitMinBiasEvent() {
  if (m_isFirst) {
    m_isFirst      = false;
    m_producedSoft = false;
    m_isMinBias    = true;
    if (!p_soft) SetB();
    return 1;
  }
  return 0;
}

int Amisic::InitRescatterEvent() {
  if (m_isFirst) {
    m_isFirst   = false;
    m_isMinBias = true;
    if (m_variable_s) UpdateForNewS();
    SetB(m_singlecollision.B());
    m_singlecollision.SetLastPT2();
    return 1;
  }
  return 0;
}

Blob * Amisic::GenerateScatter() {
  // If a next (perturbative) scatter event has been found, the pointer to the respective
  // blob is returned.
  // TODO: we may want to think about something for rescatter events - but this is for
  //       future work.
  msg_Out()<<"**** "<<METHOD<<": "<<mipars->GetEvtType()<<"\n";
  Blob * blob = m_singlecollision.NextScatter();
  if (blob) {
    if (mipars->GetEvtType()==evt_type::Perturbative) {
      m_pt2 = m_singlecollision.LastPT2();
      blob->SetPosition(m_impact.SelectPositionForScatter(m_b));
      if (m_ana) AnalysePerturbative(false);
      m_producedSoft = false;
    }
    else {
      m_producedSoft = true;
    }
    blob->SetTypeSpec("AMISIC++ 1.1");
    return blob;
  }
  if (m_ana && mipars->GetEvtType()==evt_type::Perturbative) AnalysePerturbative(true);
  return nullptr;
}

bool Amisic::VetoEvent(const double & scale) const {
  if (scale<0.) return true;
  return false;
}

Cluster_Amplitude * Amisic::ClusterConfiguration(Blob * blob) {
  Cluster_Amplitude * ampl = Cluster_Amplitude::New();
  CreateAmplitudeLegs(ampl,blob);
  FillAmplitudeSettings(ampl);
  return ampl;
}

void Amisic::CreateAmplitudeLegs(Cluster_Amplitude * ampl,Blob * blob) {
  for (size_t i(0);i<blob->NInP()+blob->NOutP();++i) {
    size_t     id(1<<ampl->Legs().size());
    Particle * part(blob->GetParticle(i));
    ColorID    col(part->GetFlow(1),part->GetFlow(2));
    if (i<blob->NInP()) {
      ampl->CreateLeg(-part->Momentum(),part->Flav().Bar(),col.Conj(),id);
    }
    else {
      ampl->CreateLeg(part->Momentum(),part->Flav(),col,id);
    }
    ampl->Legs().back()->SetStat(0);
  }
}

void Amisic::FillAmplitudeSettings(Cluster_Amplitude * ampl) {
  double muf2 = m_singlecollision.muF2(), muq2 = muf2;
  double mur2 = m_singlecollision.muR2();
  ampl->SetNIn(2);
  ampl->SetMuR2(mur2);
  ampl->SetMuF2(muf2);
  ampl->SetMuQ2(muq2);
  ampl->SetKT2(muf2);
  ampl->SetMu2(mur2);
  ampl->SetOrderEW(0);
  ampl->SetOrderQCD(2);
  ampl->SetMS(p_processes);
}

void Amisic::CleanUpMinBias() {
  SetMaxEnergies(rpa->gen.PBeam(0)[0],rpa->gen.PBeam(1)[0]);
  SetMaxScale(rpa->gen.Ecms()/2.);
  m_isFirst   = true;
  m_isMinBias = false;
}

void Amisic::Reset() {}

void Amisic::InitAnalysis() {
  m_nscatters = 0;
  m_histos[string("N_scatters")] = new Histogram(0,0,20,20);
  m_histos[string("B")]          = new Histogram(0,0,10,100);
  m_histos[string("Bfac")]       = new Histogram(0,0,10,100);
  m_histos[string("P_T(1)")]     = new Histogram(0,0,100,100);
  m_histos[string("Y(1)")]       = new Histogram(0,-10,10,10);
  m_histos[string("Delta_Y(1)")] = new Histogram(0,0,10,10);
  m_histos[string("P_T(2)")]     = new Histogram(0,0,100,100);
  m_histos[string("Y(2)")]       = new Histogram(0,-10,10,10);
  m_histos[string("Delta_Y(2)")] = new Histogram(0,0,10,10);
}

void Amisic::FinishAnalysis() {
  Histogram * histo;
  string name;
  for (map<string,Histogram *>::iterator
	 hit=m_histos.begin();hit!=m_histos.end();hit++) {
    histo = hit->second;
    name  = string("MPI_Analysis/")+hit->first+string(".dat");
    histo->Finalize();
    histo->Output(name);
    delete histo;
  }
  m_histos.clear();
}

void Amisic::AnalysePerturbative(const bool & last) {
  if (!last) {
    if (m_nscatters==0) {
      m_histos[string("P_T(1)")]->Insert(sqrt(m_singlecollision.PT2()));
      m_histos[string("Y(1)")]->Insert(m_singlecollision.Y3());
      m_histos[string("Y(1)")]->Insert(m_singlecollision.Y4());
      m_histos[string("Delta_Y(1)")]->Insert(dabs(m_singlecollision.Y3()-
						  m_singlecollision.Y4()));
    }
    if (m_nscatters==1) {
      m_histos[string("P_T(2)")]->Insert(sqrt(m_singlecollision.PT2()));
      m_histos[string("Y(2)")]->Insert(m_singlecollision.Y3());
      m_histos[string("Y(2)")]->Insert(m_singlecollision.Y4());
      m_histos[string("Delta_Y(2)")]->Insert(dabs(m_singlecollision.Y3()-
						  m_singlecollision.Y4()));
    }
  }
  m_nscatters++;
  if (last) {
    m_histos[string("N_scatters")]->Insert(double(m_nscatters)+0.5);
    m_histos[string("B")]->Insert(m_b);
    m_histos[string("Bfac")]->Insert(m_bfac);
    m_nscatters = 0;
  }
}

void Amisic::Test() {
  double Q2start(100);
  long int n(1000000);
  m_overestimator.Test(Q2start,n);
  m_singlecollision.Test(Q2start,n);
  THROW(normal_exit,"testing complete");
}
