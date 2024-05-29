#include "AMISIC++/Perturbative/Single_Collision_Handler.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/Random.H"

using namespace AMISIC;
using namespace ATOOLS;
using namespace std;

Single_Collision_Handler::Single_Collision_Handler() :
  m_evttype(evt_type::Perturbative),
  p_processes(NULL), p_overestimator(NULL), p_soft(NULL),
  m_pt2(0.), m_pt2min(0.),
  m_S((rpa->gen.PBeam(0)+rpa->gen.PBeam(1)).Abs2()), m_lastpt2(m_S),
  m_residualx1(1.), m_residualx2(1.), m_Ycms(0.),
  m_xt(1.), m_y3(0.), m_y4(0.), m_x1(1.), m_x2(1.),
  m_ana(false)
{}

Single_Collision_Handler::~Single_Collision_Handler() {
  if (m_ana) FinishAnalysis();
}

void Single_Collision_Handler::Init(MI_Processes * processes,
				    Over_Estimator * overestimator) {
  p_processes     = processes;
  p_integrator    = p_processes->GetIntegrator();
  p_xsecs         = p_processes->GetXSecs();
  p_overestimator = overestimator;
  for (size_t i=0;i<2;i++) p_remnants[i] = p_processes->GetRemnant(i);
  ///////////////////////////////////////////////////////////////////////////
  // TODO: Will have to make pt2min initial-particle dependent
  //       (currently: photon vs proton, but treated the same)
  ///////////////////////////////////////////////////////////////////////////
  m_pt2min        = sqr((*mipars)("pt_min"));
  if (m_ana) InitAnalysis();
}

void Single_Collision_Handler::Init(REMNANTS::Remnant_Handler * remnant_handler,
				    NonPerturbative_XSecs * soft) {
  m_evttype  = evt_type::NonPerturbative;
  p_remnants = remnant_handler;
  p_soft     = soft;
  if (m_ana) InitAnalysis();
}

void Single_Collision_Handler::PrefabricateBlob() {
  ///////////////////////////////////////////////////////////////////////////
  // Current use is in x-dependent matter overlap, which mix the generation
  // of a scatter with fixing the impact parameter.  Assuming the integrator
  // has produced a successful trial event, triggered by the "FirstB" method
  // in the Impact class, we select the specific process, create the
  // kinematics, and fill the blob.
  ///////////////////////////////////////////////////////////////////////////
  Blob * blob = nullptr;
  while (!blob) { blob = p_processes->FillHardScatterBlob(); }
  m_prefabs.push_back(blob);
}

Blob * Single_Collision_Handler::NextScatter() {
  // If AMISIC only produces (semi-inclusive) soft interactions use the
  // NonPerturbative_XSecs to generate a blob.
  if (p_soft) {
    Blob * blob = p_soft->MakeScatter();
    m_muf2 = p_soft->MuF2();
    m_mur2 = p_soft->MuR2();
    return blob;
  }
  ///////////////////////////////////////////////////////////////////////////
  // Simple logic - bfac is provided from outside, now
  // - either use (one of) the pre-fabricated blob(s), or
  // - generate the next pt^2 with the Sudakov-like hit-or-miss method
  //   and create a new blob based on the associated kinematics
  // - if successful, return a filled (at the moment 2->2 only) scattering blob
  ///////////////////////////////////////////////////////////////////////////
  Blob * blob = nullptr;
  if (!m_prefabs.empty()) {
    blob = m_prefabs.front();
    m_prefabs.pop_front();
  }
  else {
    m_pt2 = m_lastpt2;
    while (!blob) {
      if (!SelectPT2()) {
        if (force) m_pt2 = m_lastpt2;
        else return nullptr;
      }
      else {
        blob      = p_processes->FillHardScatterBlob();
        m_lastpt2 = m_pt2;
      }
    }
    if (m_ana && blob) Analyse(m_pt2,blob);
  }
  return blob;
}

bool Single_Collision_Handler::SelectPT2() {
  ///////////////////////////////////////////////////////////////////////////
  // Generate a trial kinematics, starting from m_lastpt2
  // - produces a trial pt2 based on a fast and crude overestimator (in TrialPT2).
  //   If it falls below the minimal value, false is returned which
  //   effectively stops the generation of secondary scatters
  // - produce rapidites for the two outgoing particles (flat distribution),
  //   reconstruct Bjorken-x for the PDFs and the Mandelstams
  // - calculate the cross section summed over all parton-level processes
  // - accept or reject the kinematics with a hit-or-miss of true over
  //   overestimated differential cross section dsigma/dpt2
  ///////////////////////////////////////////////////////////////////////////
  if (p_processes->XSratio(m_S)<1.) return false;
  double wt;
  do {
    m_pt2 = p_overestimator->TrialPT2(m_pt2);
    if (m_pt2<=m_pt2min) return false;
    if (!p_integrator->MakeKinematics(m_pt2,m_S)) continue;
    wt = (*p_processes)() / (*p_overestimator)(m_pt2, p_integrator->Yvol());
    if (m_ana) AnalyseWeight(wt);
  } while (wt<ran->Get());
  return true;
}

void Single_Collision_Handler::UpdateSandY(double s, double y) {
  ///////////////////////////////////////////////////////////////////////////
  // Updating everything that is needed for generation of next scatter
  // Setting the last pT^2 on s, as a default.
  m_lastpt2 = m_S = s;
  m_Ycms = y;
  if (m_evttype==evt_type::Perturbative) {
    p_processes->UpdateS(m_S);
    double xsnd_eff = p_xsecs->XSndNorm()*p_xsecs->XSnd();
    p_overestimator->UpdateS(m_S, xsnd_eff, p_processes->PT02());
  }
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
    while (NextScatter(true) && m_lastpt2>0) {
      if (!taken) {
	histo.Insert(m_lastpt2);
	taken = true;
      }
      SetLastPT2(m_pt2);
    }
  }
  histo.Finalize();
  histo.Output("True_PT2");
  msg_Out()<<METHOD<<": finished "<<n<<" dry runs.\n";
}
