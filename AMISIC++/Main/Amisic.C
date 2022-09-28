#include "AMISIC++/Main/Amisic.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace AMISIC;
using namespace ATOOLS;
using namespace std;

Amisic::Amisic() : m_sigmaND_norm(1.), p_processes(NULL), m_Enorm(rpa->gen.Ecms()/2.),
                   m_isMinBias(false), m_ana(false)
{}

Amisic::~Amisic() {
  if (p_processes) delete p_processes;
  if (m_ana) FinishAnalysis();
  if (p_xsecs) delete p_xsecs;
}

bool Amisic::Initialize(MODEL::Model_Base *const model,
			PDF::ISR_Handler *const isr,
                        REMNANTS::Remnant_Handler * remnant_handler)
{
  if (!InitParameters()) return false;
  bool shown = false;
  // if EPA is used the energies entering the ISR will vary
  m_variable_s = isr->GetBeam(0)->Type() == BEAM::beamspectrum::EPA ||
                 isr->GetBeam(1)->Type() == BEAM::beamspectrum::EPA;
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
  m_pbeam0 = isr->GetBeam(0)->OutMomentum();
  m_pbeam1 = isr->GetBeam(1)->OutMomentum();

  // Calculate hadronic non-diffractive cross sections, to act as normalization for the
  // multiple scattering probability. 
  p_xsecs = new Hadronic_XSec_Calculator(m_type);
  (*p_xsecs)((m_pbeam0 + m_pbeam1).Abs2());
  // show output if the calculation is not repeated for different energies
  if (!m_variable_s)
    p_xsecs->Output();

  // Initialize the parton-level processes - currently only 2->2 scatters and use the
  // information to construct a very quick overestimator - this follows closely the
  // algorithm in the original Sjostrand - van Zijl publication.
  // The logic for the overestimator is based on
  // - t-channel dominance allowing to approximate differential cross sections as
  //   dsigma ~ dp_T dy f(x_1) f_(x_2) C_1C_2 as(t)^2/(t+t_0)^2
  //   where C_1/2 are the colour factors of the incoming partons, x_1/2 are the
  //   Bjorken-x.
  // - assuming that the product of the PDFs f(x_1)f(x_2) is largest for mid-rapidity
  //   where x_1 and x_2 are identical
  m_sigmaND_norm = (*mipars)("SigmaND_Norm");
  p_processes = new MI_Processes(m_variable_s);
  p_processes->SetSigmaND(m_sigmaND_norm * p_xsecs->XSnd());
  p_processes->SetXSecCalculator(p_xsecs);
  p_processes->Initialize(model,NULL,isr);
  
  m_overestimator.Initialize(p_processes);
  m_overestimator.SetXSnd(m_sigmaND_norm * p_xsecs->XSnd());
  
  m_singlecollision.Init(remnant_handler);
  m_singlecollision.SetMIProcesses(p_processes);
  m_singlecollision.SetOverEstimator(&m_overestimator);

  m_impact.SetProcesses(p_processes);
  m_impact.Initialize(p_processes->XShard()/(m_sigmaND_norm * p_xsecs->XSnd()));

  if (m_ana) InitAnalysis();
  CleanUp();

  return true;
}

bool Amisic::InitParameters() {
  mipars = new MI_Parameters();
  return mipars->Init();
}

bool Amisic::InitMPIs(const double & scale) {
  SetMassMode(1);
  SetMaxScale(scale);
  SetB();
  if (!VetoEvent(scale)) return true;
  //m_stop = true;
  return false;
}

void Amisic::Update(const PDF::ISR_Handler *isr) {
  if (!m_variable_s)
    return;
  // Get new energy from Beams
  m_pbeam0 = isr->GetBeam(0)->OutMomentum();
  m_pbeam1 = isr->GetBeam(1)->OutMomentum();
  double s = (m_pbeam0 + m_pbeam1).Abs2();
  double y = (m_pbeam0 + m_pbeam1).Y();
  // Calculate new cross-sections
  (*p_xsecs)(s);
  // Update the processes
  p_processes->SetSigmaND(m_sigmaND_norm * p_xsecs->XSnd());
  p_processes->Update(s);
  // Update the over-estimator (implicitly uses the processes)
  m_overestimator.Update();
  m_overestimator.SetXSnd(m_sigmaND_norm * p_xsecs->XSnd());
  // Update the single-collision handler and impact parameter
  m_singlecollision.UpdateSandY(s, y);
  m_impact.Update(p_processes->XShard()/(m_sigmaND_norm * p_xsecs->XSnd()));
}

void Amisic::SetMaxEnergies(const double & E0,const double & E1) {
  m_singlecollision.SetResidualX(E0/m_pbeam0[0],E1/m_pbeam1[0]);
}

void Amisic::SetMaxScale(const double & scale) {
  m_pt2 = sqr(scale);
  m_singlecollision.SetLastPT2(m_pt2);
}

void Amisic::SetB(const double & b) {
  // Select b and set the enhancement factor
  m_b    = (b<0.)?m_impact.SelectB(m_pt2):b;
  m_bfac = Max(0.,m_impact.Enhancement());
}
  
bool Amisic::VetoEvent(const double & scale) {
  if (scale<0.) return false;
  return false;
}

const double Amisic::ScaleMin() const {
  return (*mipars)("pt_min");
}

const double Amisic::ScaleMax() const {
  return m_pt2;
}

Blob * Amisic::GenerateScatter() {
  Blob * blob = m_singlecollision.NextScatter(m_bfac);
  if (blob) {
    m_pt2 = m_singlecollision.LastPT2();
    blob->SetPosition(m_impact.SelectPositionForScatter(m_b));
    blob->SetTypeSpec("AMISIC++ 1.1");
    if (m_ana) Analyse(false);
    return blob;
  }
  if (m_ana) Analyse(true);
  return NULL;
}

int Amisic::InitMinBiasEvent(ATOOLS::Blob_List * blobs) {
  if (m_isFirst==true) {
    m_isFirst   = false;
    m_isMinBias = true;
    m_b    = m_impact.SelectB();
    m_bfac = Max(0.,m_impact(m_b));
  }
  return 0;
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

int Amisic::ShiftMasses(ATOOLS::Cluster_Amplitude * ampl) {
  return p_processes->ShiftMasses(ampl);
}


bool Amisic::VetoScatter(Blob * blob)
{
  msg_Error()<<METHOD<<" ont implemented yet.  Will exit.\n";
  exit(1);
}

void Amisic::CleanUp() {
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

void Amisic::Analyse(const bool & last) {
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
  exit(1);
}
