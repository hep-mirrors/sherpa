#include "AMISIC++/Main/Amisic.H"
#include "AMISIC++/Tools/MI_Parameters.H"
#include "AMISIC++/Tools/Hadronic_XSec_Calculator.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Scoped_Settings.H"

using namespace AMISIC;
using namespace ATOOLS;
using namespace std;

Amisic::Amisic() :
  m_sigmaND_norm(1.), p_processes(NULL),
  m_isMinBias(false), m_ana(true)
{}

Amisic::~Amisic() {
  if (p_processes) delete p_processes;
  if (m_ana) FinishAnalysis();
}

bool Amisic::Initialize(MODEL::Model_Base *const model,
			PDF::ISR_Handler *const isr)
{
  if (!InitParameters()) return false;
  bool shown = false;
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

  // Calculate hadronic non-diffractive cross sections, to act as normalization for the
  // multiple scattering probability. 
  Hadronic_XSec_Calculator xsecs;
  xsecs();

  // Initialize the parton-level processes - currently only 2->2 scatters and use the
  // information to construct a verry quick overestimator - this follows closely the
  // algorithm in the original Sjostrand - van der Zijl publication.
  // The logic for the overestimator is based on
  // - t-channel dominance allowing to approximate differential cross sections as
  //   dsigma ~ dp_T dy f(x_1) f_(x_2) C_1C_2 as(t)^2/(t+t_0)^2
  //   where C_1/2 are the colour factors of the incoming partons, x_1/2 are the
  //   Bjorken-x.
  // - assuming that the product of the PDFs f(x_1)f(x_2) is largest for mid-rapidity
  //   where x_1 and x_2 are identical
  m_sigmaND_norm = (*mipars)("SigmaND_Norm");
  p_processes = new MI_Processes();
  p_processes->SetSigmaND(m_sigmaND_norm * xsecs.XSnd());
  p_processes->Initialize(model,NULL,isr);
  
  m_overestimator.Initialize(p_processes);
  m_overestimator.SetXSnd(m_sigmaND_norm * xsecs.XSnd());
  
  m_singlecollision.Init();
  m_singlecollision.SetMIProcesses(p_processes);
  m_singlecollision.SetOverEstimator(&m_overestimator);

  m_impact.SetProcesses(p_processes);
  m_impact.Initialize(p_processes->XShard()/(m_sigmaND_norm * xsecs.XSnd()));

  if (m_ana) InitAnalysis();
  CleanUp();
  
  return true;
}

bool Amisic::InitParameters() {
  mipars = new MI_Parameters();
  return mipars->Init();
}

void Amisic::SetMaxEnergies(const double & E1,const double & E2) {
  m_residualE1 = E1;
  m_residualE2 = E2;
  m_singlecollision.SetResidualX(2.*E1/rpa->gen.Ecms(),2.*E2/rpa->gen.Ecms());
}

void Amisic::SetMaxScale(const double & scale) {
  m_pt2 = sqr(scale);
  m_singlecollision.SetLastPT2(m_pt2);
}

void Amisic::SetB(const double & b) {
  // Select b and set the enhancement factor
  m_b    = (b<0.)?m_impact.SelectB(m_pt2):b;
  m_bfac = Max(0.,m_impact(m_b));
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
  //msg_Out()<<METHOD<<": asking Amisic to start filling Min Bias event, "
  //	   <<"pt = "<<sqrt(m_pt2)<<" for b = "<<m_b<<" and bfac = "<<m_bfac<<".\n";
  Blob * blob = m_singlecollision.NextScatter(m_bfac);
  if (blob) {
    m_pt2 = m_singlecollision.LastPT2();
    blob->SetPosition(m_impact.SelectPositionForScatter(m_b));
    blob->SetTypeSpec("AMISIC++ 1.1");
    if (m_ana) Analyse(false);
    return blob;
  }
  else {
    if (m_ana) Analyse(true);
  }
  return NULL;
}

int Amisic::InitMinBiasEvent(ATOOLS::Blob_List * blobs) {
  if (m_isFirst==true) {
    //msg_Out()<<METHOD<<": asking Amisic to provide Min Bias event.\n";
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
  msg_Out()<<METHOD<<" ont implemented yet.  Will exit.\n";
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
