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
  m_done(false), 
  m_ana(false)
{}

Single_Collision_Handler::~Single_Collision_Handler() {
  if (m_ana) FinishAnalysis();
}

void Single_Collision_Handler::
Init(MI_Processes * processes,Over_Estimator * overestimator,
     Interaction_Probability * pint, Matter_Overlap * overlap) {
  p_processes     = processes;
  p_integrator    = p_processes->GetIntegrator();
  p_xsecs         = p_processes->GetXSecs();
  m_needsTrigger  = (p_processes->TrigSize()>0);
  p_overestimator = overestimator;
  for (size_t i=0;i<2;i++) p_remnants[i] = p_processes->GetRemnant(i);
  p_pint          = pint;
  p_overlap       = overlap;
  ///////////////////////////////////////////////////////////////////////////
  // TODO: Will have to make pt2min initial-particle dependent
  //       (currently: photon vs proton, but treated the same)
  ///////////////////////////////////////////////////////////////////////////
  m_pt2min        = sqr((*mipars)("pt_min"));
  if (m_ana) InitAnalysis();
}

void Single_Collision_Handler::
Init(REMNANTS::Remnant_Handler * remnant_handler,
     NonPerturbative_XSecs * soft) {
  m_evttype  = evt_type::NonPerturbative;
  for (size_t i=0;i<2;i++) p_remnants[i] = remnant_handler->GetRemnant(i);
  p_soft     = soft;
  if (m_ana) InitAnalysis();
}

void Single_Collision_Handler::Reset() {
  for (map<double,Blob *>::iterator mit=m_prefabs.begin();
       mit!=m_prefabs.end();mit++) {
    delete mit->second;
  }
  m_prefabs.clear();
  m_nextblob = nullptr;
  m_nextpt2  = -1.;
  m_done     = false;
}


double Single_Collision_Handler::PrefabricateBlob(const int & mode) {
  ///////////////////////////////////////////////////////////////////////////
  // Current use is in producing "biased" minbias, i.e. by enforcing a
  // particular final state at parton level.
  // The method will return the relative weight for its production.
  ///////////////////////////////////////////////////////////////////////////
  double wt = 1., pt2;
  if (mode!=0) {
    Blob * blob = nullptr;
    while (!blob) {
      if (p_integrator->TrialEvent(m_S,p_overlap)) {
	wt  = p_processes->MakeTriggerBlob(blob);
      }
    }
    m_prefabs[p_integrator->PT2()] = blob;
  }
  return wt;
}


bool Single_Collision_Handler::FirstMPI(Blob * signal) {
  ///////////////////////////////////////////////////////////////////////////
  // Initialising the MPI sequence along the following steps:
  // 1. Extract the relevant kinematics variables, including in particular the
  //    hardest pt^2_veto scale, from the signal blob and adjust the effective
  //    dynamic radius in the matter overlap.
  // 2. Produce trial impact parameters m_b and trial pt^2_trial of potential 
  //    first MPI scatters, until a trial pt^2_trial<pt^2_veto for a viable
  //    scatter that is softer than the signal has been found.
  // TODO: in case the signal is independent of, i.e. different from, the 2->2
  //       scatters in the potential MPI scatters, we should not veto but
  //       just pick an impact parameter m_b and start without the veto
  //       condition.
  ///////////////////////////////////////////////////////////////////////////
  double pt2veto = sqr((*signal)["MI_Scale"]->Get<double>());
  double x1      = (*signal)["PDFInfo"]->Get<PDF_Info>().m_x1;
  double x2      = (*signal)["PDFInfo"]->Get<PDF_Info>().m_x2;
  p_overlap->FixDynamicRadius(x1,x2,pt2veto,pt2veto);
  do {
    /////////////////////////////////////////////////////////////////////////
    // Form_Factor (and by extension Matter_Overlap) has radii etc. in fm,
    // event record needs it in mm, therefore we have to divide by 10^12.
    /////////////////////////////////////////////////////////////////////////    
    do {
      m_b   = p_overlap->SelectB();
    } while (!p_overlap->SelectPositionForScatter(m_b,x1,m_pt2,x2,m_pt2,m_deltapos));
    m_pt2 = m_lastpt2 = m_S/4.;
  } while (SelectPT2()==0 && m_pt2>pt2veto);
  m_lastpt2 = m_pt2;
  if (m_pt2<=m_pt2min) m_done = true;
  return true;
}

bool Single_Collision_Handler::FirstMinBiasScatter(Blob * blob) {
  ///////////////////////////////////////////////////////////////////////////
  // If AMISIC only produces (semi-inclusive) soft interactions use the
  // NonPerturbative_XSecs to generate a blob.
  ///////////////////////////////////////////////////////////////////////////
  if (p_soft) {
    if (p_soft->MakeScatter(blob)) {
      blob->AddData("Renormalization_Scale",new Blob_Data<double>(p_soft->MuR2()));
      blob->AddData("Factorization_Scale",new Blob_Data<double>(p_soft->MuF2()));
      blob->SetTypeSpec("Soft_Process_from_AMISIC");
      m_done = true;
      return true;
    }
    THROW(fatal_error,"Couldn't fill non-perturbative scatter blob.");
  }
  ///////////////////////////////////////////////////////////////////////////
  // Start the Min bias event by repeating the two steps below, until
  // a viable first (hardest) blob is found.
  // 1. Selecting b according to db dphi b P_int(b)/sigma_ND from
  //    the Interaction_Probability class
  // 2. Produce a hard scatter, using the NextScatter method.
  ///////////////////////////////////////////////////////////////////////////
  blob->ClearAllData();
  blob->DeleteOwnedParticles();
  bool success = false;
  int  trials  = 10;
  while (!success) {
    m_b    = p_pint->SelectB(m_S);
    trials = 10;
    do {
      m_done  = false;
      m_pt2   = m_lastpt2 = m_S/4.;
      success = NextScatter(blob);
    } while (!success && (trials--)>0);
  }
  ///////////////////////////////////////////////////////////////////////////////
  // Updating the new "signal" blob of the min bias event and adding the last
  // bits of relevant information.
  // Switching Sudakov evolution back on: m_done = false
  ///////////////////////////////////////////////////////////////////////////////
  blob->AddData("Trials",new Blob_Data<size_t>(1));
  blob->AddData("Weight_Norm",new Blob_Data<double>(p_xsecs->XSnd()));
  m_done = false;
  return true;
}

bool Single_Collision_Handler::FirstRescatter(Blob * blob) {
  blob->ClearAllData();
  blob->DeleteOwnedParticles();
  m_done = false;
  m_pt2  = m_lastpt2 = m_S/4.;
  ///////////////////////////////////////////////////////////////////////////////
  // Updating the new "signal" blob of the min bias event and adding the last
  // bits of relevant information.
  // Switching Sudakov evolution back on: m_done = false
  ///////////////////////////////////////////////////////////////////////////////
  if (NextScatter(blob)) {
    blob->AddData("Trials",new Blob_Data<size_t>(1));
    blob->AddData("Weight_Norm",new Blob_Data<double>(p_xsecs->XSnd()));
    blob->AddStatus(blob_status::needs_beamRescatter);
    return true;
  }
  m_done = true;
  return false;
  ///////////////////////////////////////////////////////////////////////////
  // Todo: Add soft interactions as an option here.
  ///////////////////////////////////////////////////////////////////////////
}

bool Single_Collision_Handler::NextScatter(Blob * blob) {
  ///////////////////////////////////////////////////////////////////////////
  // Simple logic:
  // - either use (one of) the pre-fabricated blob(s), or
  // - generate the next pt^2 with the Sudakov-like hit-or-miss method
  //   and create a new blob based on the associated kinematics
  // - if successful, return a filled (at the moment 2->2 only)
  //   scattering blob
  ///////////////////////////////////////////////////////////////////////////
  if (m_done) return false;
  m_pt2 = m_lastpt2;
  DefinePossibleNext();
  int  fill = 0;
  bool done = false;
  while (!done) {
    switch (SelectPT2()) {
    case 1:
      ////////////////////////////////////////////////////////////////////////
      // pt^2 is smaller than the pt^2 of the hardest pre-fabricated blob,
      // we will copy it over and remove it from the (ordered) list.
      ////////////////////////////////////////////////////////////////////////
      CopyPrefabricatedBlob(blob);
      m_prefabs.erase(--m_prefabs.end());
      m_lastpt2 = m_pt2;
      return true;
    case 0:
      ////////////////////////////////////////////////////////////////////////
      // pt^2 is harder than the pt^2 of the hardest pre-fabricated blob,
      // or such a blob doesn't exist.  Will check if we can fill a
      // hard scatter blob, if successful we can stop.
      ////////////////////////////////////////////////////////////////////////
      fill = p_processes->FillHardScatterBlob(blob,m_lastpt2);
      if (fill==1) {
        m_lastpt2 = m_pt2;
	double x1 = 0., x2 = 0.;
	if ((*blob)["PDFInfo"]!=NULL) {
	  x1 = (*blob)["PDFInfo"]->Get<PDF_Info>().m_x1;
	  x2 = (*blob)["PDFInfo"]->Get<PDF_Info>().m_x2;
	}
	m_deltapos = Vec4D(0.,0.,0.,0.);
	p_overlap->SelectPositionForScatter(m_b,x1,m_pt2,x2,m_pt2,m_deltapos);
	blob->SetType(btp::Hard_Collision);
	blob->SetStatus(blob_status::needs_showers);
	blob->SetId();
	return true;
      }
      else {
	if (fill==-1) THROW(fatal_error,"mismatched sequence of pt^2");
      }
      break;
    case -1:
    default:
      ////////////////////////////////////////////////////////////////////////
      // pt^2 is smaller than pt2min.  
      // This means we will switch off "Sudakov evolution" at the end of the
      // method (m_done = true) - for the FirstMinBias method we will have to
      // switch it on again.
      ////////////////////////////////////////////////////////////////////////
      done = true;
      break;
    }
  }
  if (m_pt2<m_pt2min) m_done = true;
  m_lastpt2 = m_pt2;
  return false;
}

void Single_Collision_Handler::DefinePossibleNext() {
  ///////////////////////////////////////////////////////////////////////////
  // Checks the pt2-ordered list of prefabricated blobs and, if there are any 
  // left, the next pt^2 and corresponding blob are extracted from the list. 
  ///////////////////////////////////////////////////////////////////////////
  m_nextpt2  = -1.;
  m_nextblob = nullptr;
  if (!m_prefabs.empty()) {
    m_nextpt2  = (--m_prefabs.end())->first;
    m_nextblob = (--m_prefabs.end())->second;
  }
}

int Single_Collision_Handler::SelectPT2() {
  ///////////////////////////////////////////////////////////////////////////
  // Generate a trial kinematics, starting from m_lastpt2
  // - produces a trial pt2 based on a fast and crude overestimator (in
  //   TrialPT2). If it falls below the minimal value, -1 is returned,
  //    which effectively stops the generation of secondary scatters
  // - if the trial pt2 is below the (next, if present) prefabricated blob,
  //   its own pt^2 scale overwrites the trial pt2 and 1 is returned.
  // - produce rapidites for the two outgoing particles (flat distribution),
  //   reconstruct Bjorken-x for the PDFs and the Mandelstams
  // - calculate the cross section summed over all parton-level processes
  // - accept or reject the kinematics with a hit-or-miss of true over
  //   overestimated differential cross section dsigma/dpt2,
  //   if accepted 0 is returned. 
  ///////////////////////////////////////////////////////////////////////////
  while (true) {
    m_pt2 = p_overestimator->TrialPT2(m_pt2);
    if (m_pt2<=m_pt2min)  return -1;
    if (m_pt2<=m_nextpt2) { m_pt2 = m_nextpt2; return 1; }
    if (!p_integrator->MakeKinematics(m_pt2,m_S)) continue;
    p_overlap->FixDynamicRadius(p_integrator->X(0),p_integrator->X(1));
    double wt = ( p_processes->PDFnorm()*(*p_processes)()/
		  (*p_overestimator)(m_pt2,p_integrator->Yvol()) *
		  (*p_overlap)(m_b) );
    if (m_ana) AnalyseWeight(wt);
    if (wt>=ran->Get()) break;
  }
  return 0;
}

void Single_Collision_Handler::CopyPrefabricatedBlob(Blob * blob) {
  if (m_prefabs.empty() || m_nextblob==nullptr) return;
  Particle_Vector * parts = m_nextblob->InParticles();
  while (!parts->empty()) {
    Particle * part = parts->back();
    part->SetDecayBlob(nullptr);
    blob->AddToInParticles(part);
    parts->pop_back();
  }
  parts = m_nextblob->OutParticles();
  while (!parts->empty()) {
    Particle * part = parts->back();
    part->SetProductionBlob(nullptr);
    blob->AddToOutParticles(part);
    parts->pop_back();
  }
  for (String_BlobDataBase_Map::const_iterator sdit=m_nextblob->GetData().begin();
       sdit!=m_nextblob->GetData().end();sdit++) {
    if (sdit->first!="WeightsMap")
      blob->AddData(sdit->first,sdit->second->ClonePtr());
  }
  blob->SetId();
  blob->SetTypeSpec("AMISIC++ 1.1");
  blob->SetType(btp::Hard_Collision);
  blob->SetStatus(blob_status::needs_showers);
  delete m_nextblob;
}

void Single_Collision_Handler::UpdateSandY(double s, double y) {
  ///////////////////////////////////////////////////////////////////////////
  // Updating everything that is needed for generation of next scatter
  // Setting the last pT^2 on s, as a default.
  ///////////////////////////////////////////////////////////////////////////
  m_lastpt2 = m_S = s;
  m_Ycms = y;
  if (m_evttype==evt_type::Perturbative) {
    p_processes->UpdateS(m_S);
    p_overestimator->UpdateS(m_S,p_processes->PT02(),p_processes->PT2Min());
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
