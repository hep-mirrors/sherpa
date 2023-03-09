#include "AMISIC++/Perturbative/Single_Collision_Handler.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/Random.H"

using namespace AMISIC;
using namespace ATOOLS;
using namespace std;

Single_Collision_Handler::Single_Collision_Handler() :
  p_processes(NULL), p_overestimator(NULL),
  m_pt2(0.), m_pt2min(0.),  
  m_S((rpa->gen.PBeam(0)+rpa->gen.PBeam(1)).Abs2()), m_lastpt2(m_S),
  m_residualx1(1.), m_residualx2(1.), m_Ycms(0.),
  m_xt(1.), m_y3(0.), m_y4(0.), m_x1(1.), m_x2(1.),
  m_ana(false)
{}

Single_Collision_Handler::~Single_Collision_Handler() {
  if (m_ana) FinishAnalysis();
}

void Single_Collision_Handler::
Init(REMNANTS::Remnant_Handler * remnant_handler,
     MI_Processes * processes, Over_Estimator * overestimator) {
  p_remnants      = remnant_handler;
  p_processes     = processes;
  p_overestimator = overestimator;
  /////////////////////////////////////////////////////////////////////////////////
  // TODO: Will have to make pt2min initial-particle dependent
  //       (currently: photon vs proton, but treated the same)
  /////////////////////////////////////////////////////////////////////////////////
  m_pt2min        = sqr((*mipars)("pt_min"));
  if (m_ana) InitAnalysis();
}

Blob * Single_Collision_Handler::NextScatter() {
  /////////////////////////////////////////////////////////////////////////////////
  // Simple logic - bfac is provided from outside, now
  // - produce a trial kinematics (new transverse momentum smaller than the last one,
  //   supplemented with rapidities, Mandelstams ...)
  // - select a process
  // - construct its full kinematics and set the colours
  // - return a filled (at the moment 2->2 only) scattering blob
  /////////////////////////////////////////////////////////////////////////////////
  do {
    if (!SelectPT2(m_lastpt2)) return NULL;
    p_proc = p_processes->SelectProcess();
  }
  while (!p_proc || !p_proc->MakeKinematics(m_pt2,m_y3,m_y4,sqrt(m_shat)) ||
	 !p_proc->SetColours() || !TestRemnants());
  return MakeBlob();
}

bool Single_Collision_Handler::TestRemnants() const {
  /////////////////////////////////////////////////////////////////////////////////
  // Make sure there is enough energy left in the remnants
  /////////////////////////////////////////////////////////////////////////////////
  return ( p_remnants->GetRemnant(0)->TestExtract(p_proc->Flav(0), p_proc->Momentum(0)) &&
	   p_remnants->GetRemnant(1)->TestExtract(p_proc->Flav(1), p_proc->Momentum(1)) );
}

bool Single_Collision_Handler::SelectPT2(const double & pt2) {
  /////////////////////////////////////////////////////////////////////////////////
  // Generate a trial kinematics
  // - produces a trial pt2 based on a fast and crude overestimator (in TrialPT2).
  //   if it falls below the minimal value, false is returned which effectively stops
  //   the generation of secondary scatters
  // - produce rapidites for the two outgoing particles (flat distribution), reconstruct
  //   Bjorken-x for the PDFs and the Mandelstams
  // - calculate the cross section summed over all parton-level processes
  // - accept or reject the kinematics with a hit-or-miss of true over overestimated
  //   differential cross section dsigma/dpt2
  /////////////////////////////////////////////////////////////////////////////////
  if (p_processes->GetSudakov()->XSratio(m_S)<1.) return false;
  m_pt2 = pt2;
  double dsigmatrue, dsigmaapprox, weight;
  bool success(false), output(false);
  do {
    m_pt2  = p_overestimator->TrialPT2(m_pt2);
    m_muf2 = m_mur2 = m_pt2;
    msg_Out()<<"         * "<<METHOD<<"(pt^2 = "<<m_pt2<<", s = "<<m_S<<"):\n";
    if (m_pt2<m_pt2min) return false;
    if (!SelectRapidities() || !CalcXs() || !CalcMandelstams()) continue;
    dsigmatrue   = (*p_processes)(m_shat,m_that,m_uhat,m_x1,m_x2);
    dsigmaapprox = (*p_overestimator)(m_pt2, m_yvol);
    weight       = dsigmatrue/dsigmaapprox;
    msg_Out()<<"                xsec ratio = "<<weight<<" from "
	     <<dsigmatrue<<" / "<<dsigmaapprox<<"\n";
    if (m_ana) AnalyseWeight(weight);
    if (weight > ran->Get()) success = true;
  } while (!success);
  SetLastPT2(m_pt2);
  return true;
}

Blob * Single_Collision_Handler::MakeBlob() {
  //////////////////////////////////////////////////////////////////////////////////////
  // Making a hard collision blob for pertutbative scatters and fill it with scales
  // and incoming and outgoing particles.
  //////////////////////////////////////////////////////////////////////////////////////
  Blob * blob = new Blob();
  blob->SetType(btp::Hard_Collision);
  blob->SetStatus(blob_status::needs_showers);
  blob->SetId();
  blob->AddData("WeightsMap",new Blob_Data<Weights_Map>({}));
  blob->AddData("Renormalization_Scale",new Blob_Data<double>(m_pt2));
  blob->AddData("Factorization_Scale",new Blob_Data<double>(m_pt2));
  blob->AddData("Resummation_Scale",new Blob_Data<double>(m_pt2));
  for (size_t i=0;i<2;i++) blob->AddToInParticles(p_proc->GetParticle(i));
  for (size_t i=2;i<4;i++) blob->AddToOutParticles(p_proc->GetParticle(i));
  if (m_ana) Analyse(m_pt2,blob);
  return blob;
}

void Single_Collision_Handler::UpdateSandY(double s, double y) {
  //////////////////////////////////////////////////////////////////////////////////////
  // Updating everything that is needed for generation of next scatter
  // Setting the last pT^2 on s, as a default.
  //////////////////////////////////////////////////////////////////////////////////////
  m_lastpt2 = m_S = s;
  m_Ycms = y;
  p_processes->UpdateS(m_S);
  p_overestimator->UpdateS();
}

bool Single_Collision_Handler::SelectRapidities() {
  //////////////////////////////////////////////////////////////////////////////////////
  // Generate two trial rapidities and keep the volume of the rapidity range -
  // it has to be divided out of the approximated/overestimated differential cross section 
  // as the prefactor there includes the volume to allow fast pt^2 generation.
  //////////////////////////////////////////////////////////////////////////////////////
  m_xt   = sqrt(4.*m_pt2/m_S);
  if (m_xt>1.) return false;
  m_ymax = log(1./m_xt*(1.+sqrt(1.-m_xt*m_xt)));
  m_yvol = sqr(2.*m_ymax);
  m_y3   = m_ymax*(2.*ran->Get()-1.);
  m_y4   = m_ymax*(2.*ran->Get()-1.);
  return true;
}

bool Single_Collision_Handler::CalcXs() {
  //////////////////////////////////////////////////////////////////////////////////////
  // Obtain Bjorken-x from the trial rapidities and make sure they are still smaller
  // than the residual x in the remnants.  Will enter the exact PDF calculation in
  // the exact differential cross section.
  //////////////////////////////////////////////////////////////////////////////////////
  m_x1   = m_xt*(exp(m_y3)+exp(m_y4))/2.*exp(-m_Ycms);
  m_x2   = m_xt*(exp(-m_y3)+exp(-m_y4))/2.*exp(m_Ycms);
  //////////////////////////////////////////////////////////////////////////////////////
  // TODO: Misses term for incoming masses, we will hope for the best by treating them
  //       as massless - that could be improved in the future.
  //////////////////////////////////////////////////////////////////////////////////////
  msg_Debugging()<<"                x1 = "<<m_x1<<" ("<<m_residualx1<<"), "
		 <<"x2 = "<<m_x2<<" ("<<m_residualx2<<"), "
		 <<"xt = "<<m_xt<<".\n";
  if (m_x1<p_processes->PDFXmin(0) || m_x2<p_processes->PDFXmin(1) ||
      m_x1>1. || m_x2>1.) return 0.;
  return (m_x1<m_residualx1 && m_x2<m_residualx2);
}

bool Single_Collision_Handler::CalcMandelstams() {
  //////////////////////////////////////////////////////////////////////////////////////
  // Mandelstams for the exact matrix element.
  //////////////////////////////////////////////////////////////////////////////////////
  if (m_xt*m_xt>m_x1*m_x2) return false;
  double tanhy = sqrt(1.-(m_xt*m_xt)/(m_x1*m_x2)); // = tanh((y3+y4)/2)
  m_shat = m_x1*m_x2*m_S;
  m_that = -0.5*m_shat*(1.-tanhy);
  m_uhat = -0.5*m_shat*(1.+tanhy);
  return true;
}

void Single_Collision_Handler::InitAnalysis() {
  m_histos[string("weights")]     = new Histogram(0,0.,2.,200);
  m_histos[string("weights_low")] = new Histogram(0,0.,0.1,1000);
  m_histos[string("pt")]          = new Histogram(0,0.,200.,200);
  m_histos[string("flavs")]       = new Histogram(0,-6.5,6.5,13);
}

void Single_Collision_Handler::FinishAnalysis() {
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

void Single_Collision_Handler::AnalyseWeight(const double & weight) {
  m_histos[string("weights")]->Insert(weight);
  m_histos[string("weights_low")]->Insert(weight);
}
  
void Single_Collision_Handler::Analyse(const double & pt2,Blob * blob) {
  m_histos[string("pt")]->Insert(sqrt(pt2));
  Flavour flav1 = blob->OutParticle(0)->Flav();
  Flavour flav2 = blob->OutParticle(1)->Flav();
  int fl1 = size_t(flav1.Kfcode());
  if (flav1.IsAnti()) fl1 = -fl1;
  if (fl1==21) fl1=6;
  if (fl1==22) fl1=-6;
  int fl2 = size_t(flav2.Kfcode());
  if (flav2.IsAnti()) fl2 = -fl2;
  if (fl2==21) fl2=6;
  if (fl2==22) fl2=-6;
  m_histos[string("flavs")]->Insert(fl1);
  m_histos[string("flavs")]->Insert(fl2);
}
  
void Single_Collision_Handler::Test(const double & Q2,const long int & n) {
  msg_Out()<<METHOD<<" for Q^2 = "<<Q2<<", s = "<<m_S<<".\n";
  Histogram histo(0,0.0,Q2,100);
  for (long int dryrun=0;dryrun<n;dryrun++) {
    SetLastPT2(Q2);
    bool taken(false);
    while (NextScatter() && m_pt2>0) {
      if (!taken) {
	histo.Insert(m_pt2);
	taken = true;
      }
      SetLastPT2(m_pt2);
    }
  }
  histo.Finalize();
  histo.Output("True_PT2");
  msg_Out()<<METHOD<<": finished "<<n<<" dry runs.\n";
}
