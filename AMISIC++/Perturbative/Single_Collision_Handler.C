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

void Single_Collision_Handler::Init(REMNANTS::Remnant_Handler * remnant_handler,
				    NonPerturbative_XSecs * soft) {
  m_evttype  = evt_type::NonPerturbative;
  p_remnants = remnant_handler;
  p_soft     = soft;
  if (m_ana) InitAnalysis();
}

void Single_Collision_Handler::Reset() {
  while (!m_prefabs.empty()) {
    delete m_prefabs.back();
    m_prefabs.pop_back();
  }
  m_done = false;
}

bool Single_Collision_Handler::PrefabricateBlob() {
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

bool Single_Collision_Handler::FirstMPI(Blob * signal) {
  //msg_Out()<<"  * "<<METHOD<<": |"<<signal<<"|\n";
  double pt2_1  = sqr((*signal)["MI_Scale"]->Get<double>());
  double x1     = (*signal)["PDFInfo"]->Get<PDF_Info>().m_x1;
  double x2     = (*signal)["PDFInfo"]->Get<PDF_Info>().m_x2;
  p_overlap->FixDynamicRadius(x1,x2,pt2_1,pt2_1);
  //msg_Out()<<"  * "<<METHOD<<"(pt^2 = "<<pt2_1<<", x1 = "<<x1<<", x2 = "<<x2<<"): "
  //	   <<"R = "<<sqrt(p_overlap->DynamicRadius2())<<".\n";
  while (true) {
    m_b     = p_overlap->SelectB();
    m_fb    = p_pint->fb(m_S,m_b);
    m_bfac  = m_fb * p_overlap->MaxValue(m_b);
    p_overestimator->SetBFac(m_bfac);
    m_pt2   = m_S/4.;
    //msg_Out()<<"  * "<<METHOD<<" starts veto cycle: b = "<<m_b<<", "
    //<<"bfac = "<<m_bfac<<"\n";
    while (SelectPT2()) {
      Blob * second = new Blob();
      msg_Out()<<"  * "<<METHOD<<" inits |"<<second<<"|, "
	       <<Blob::Counter()<<"/"<<Particle::Counter()<<"\n";
      int fill      = p_processes->FillHardScatterBlob(second,pt2_1);
      // New blob (second) was successfull filled and fulfils requirements that
      // its pt^2 < pt^2_1, the pt^2 of the signal blob.
      // We store the new blob now to get it off the shelf in the NextScatter call.
      if (fill==1) {
	m_prefabs.push_back(second);
	m_lastpt2 = m_pt2;
	//msg_Out()<<"  * "<<METHOD<<" ends veto cycle ("<<m_done<<"): "
	//	 <<"|"<<second<<"|, "<<Blob::Counter()<<"/"<<Particle::Counter()<<"\n";
	return true;
      }
      // New blob (second) couldn't be created and needs to be deleted.
      else {
	delete second;
	//msg_Out()<<"  * "<<METHOD<<" deletes |"<<second<<"|, "
	//	 <<Blob::Counter()<<"/"<<Particle::Counter()<<"\n";
	second = nullptr;
	// Kinematics of new blob had to be vetoed, as its pt^2 > pt^2_1, which
	// means we have to select a new impact parameter b.
	if (fill==-1) break;
      }
    }
    if (m_pt2<m_pt2min) {
      m_lastpt2 = 0.;
      m_done    = true;
      //msg_Out()<<"  * "<<METHOD<<" ends veto cycle ("<<m_done<<"): |no blob|\n";
      return true;
    }
  }
  return true;
}

bool Single_Collision_Handler::FirstMinBiasScatter(Blob * first) {
  bool done = false;
  while (!done) {
    double pt2_1 = MakeHardScatterBlob(first);
    p_overlap->FixDynamicRadius(p_integrator->X(0),p_integrator->X(1));
    m_b          = p_overlap->SelectB();
    m_fb         = p_pint->fb(m_S,m_b);
    m_bfac       = m_fb * p_overlap->MaxValue(m_b);
    //msg_Out()<<"  * "<<METHOD<<"(pt^2 = "<<pt2_1<<"): "
    //	     <<"R = "<<sqrt(p_overlap->DynamicRadius2())<<".\n";
    p_overestimator->SetBFac(m_bfac);
    m_pt2        = m_S/4.;
    while (SelectPT2()) {
      Blob * second  = new Blob();
      //msg_Out()<<"  * "<<METHOD<<" inits |"<<second<<"|, "
      //       <<Blob::Counter()<<"/"<<Particle::Counter()<<"\n";
      int fill       = p_processes->FillHardScatterBlob(second,pt2_1);
      // New blob (second) was successfull filled and fulfils requirements that
      // its pt^2 < pt^2_1, the pt^2 of the signal blob.
      // We store the new blob now to get it off the shelf in the NextScatter call.
      if (fill==1) {
	m_prefabs.push_back(second);
	m_lastpt2 = m_pt2;
	m_done    = false;
	done      = true;
	break;
      }
      // New blob (second) couldn't be created and needs to be deleted.
      else {
	delete second;
	//msg_Out()<<"  * "<<METHOD<<" deletes |"<<second<<"|, "
	//	 <<Blob::Counter()<<"/"<<Particle::Counter()<<"\n";
	second = nullptr;
	// Kinematics of new blob had to be vetoed, as its pt^2 > pt^2_1, which
	// means we have to try a new starting event.
	if (fill==-1) {
	  first->ClearAllData();
	  first->DeleteOwnedParticles();
	  break;
	}
      }
      if (m_pt2<m_pt2min) {
	m_lastpt2 = 0.;
	m_done = done = true;
	break;
      }
    }
  };
  // Updating the new "signal" blob of the min bias event and adding the last
  // bits of relevant information.
  first->SetType(btp::Hard_Collision);
  first->SetStatus(blob_status::needs_showers);
  first->AddData("Trials",new Blob_Data<size_t>(1));
  first->AddData("Weight_Norm",new Blob_Data<double>(p_xsecs->XSnd()));
  msg_Out()<<METHOD<<"(first):\n"<<(*first)<<"\n";
  return true;
}

double Single_Collision_Handler::MakeHardScatterBlob(Blob * blob) {
  ///////////////////////////////////////////////////////////////////////////
  // Here we do not apply any pt2 veto as this is the first blob in a
  // min bias or rescattering sequence.
  ///////////////////////////////////////////////////////////////////////////
  do { } while (!(p_integrator->TrialEvent(m_S,p_overlap) &&
		  p_processes->FillHardScatterBlob(blob)));
  return p_integrator->PT2();
}

bool Single_Collision_Handler::NextScatter(Blob * blob) {
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
  if (m_done) return false;
  if (!m_prefabs.empty()) {
    //msg_Out()<<"  * "<<METHOD<<" copies blob: "
    //	     <<"|"<<m_prefabs.front()<<"| --> |"<<blob<<"|, "
    //	     <<Blob::Counter()<<"/"<<Particle::Counter()<<".\n";
    CopyPrefabricatedBlob(blob);
    blob->SetType(btp::Hard_Collision);
    blob->SetStatus(blob_status::needs_showers);
    msg_Out()<<METHOD<<"(prefab):\n"<<(*blob)<<"\n";
    return true;
  }
  m_pt2        = m_lastpt2;
  //msg_Out()<<"  * "<<METHOD<<" starting veto cycle (pt^2 = "<<m_pt2<<"): "
  //	   <<"b = "<<m_b<<", bfac = "<<m_bfac<<"\n";
  while (SelectPT2()) {
        if (m_pt2<m_pt2min) {
        //msg_Out()<<"  * "<<METHOD<<": no blob |"<<blob<<"|, "
      //       <<Blob::Counter()<<"/"<<Particle::Counter()<<" for "<<m_pt2<<"\n";
      m_done = true;
      return false;
        }
    int fill = p_processes->FillHardScatterBlob(blob,m_lastpt2);
    if (fill==1) {
      //msg_Out()<<"  * "<<METHOD<<": blob for "<<m_pt2<<": |"<<blob<<"|, "
      //       <<Blob::Counter()<<"/"<<Particle::Counter()<<"\n";
        m_lastpt2 = m_pt2;
      blob->SetType(btp::Hard_Collision);
      blob->SetStatus(blob_status::needs_showers);
      msg_Out()<<METHOD<<"(new one):\n"<<(*blob)<<"\n";
      return true;
    }
    else {
      if (fill==-1) THROW(fatal_error,"mismatched sequence of pt^2");
    }
  }
  //msg_Out()<<"  * "<<METHOD<<": no blob |"<<blob<<"|, "
  //	   <<Blob::Counter()<<"/"<<Particle::Counter()<<" for "<<m_pt2<<"\n";
  m_done = true;
  return false;
}

bool Single_Collision_Handler::SelectPT2() {
  ///////////////////////////////////////////////////////////////////////////
  // Generate a trial kinematics, starting from m_lastpt2
  // - produces a trial pt2 based on a fast and crude overestimator (in
  //   TrialPT2). If it falls below the minimal value, false is returned,
  //    which effectively stops the generation of secondary scatters
  // - produce rapidites for the two outgoing particles (flat distribution),
  //   reconstruct Bjorken-x for the PDFs and the Mandelstams
  // - calculate the cross section summed over all parton-level processes
  // - accept or reject the kinematics with a hit-or-miss of true over
  //   overestimated differential cross section dsigma/dpt2
  ///////////////////////////////////////////////////////////////////////////
  double wt=0.;
  do {
    m_pt2 = p_overestimator->TrialPT2(m_pt2);
    if (m_pt2<=m_pt2min) return false;
    if (!p_integrator->MakeKinematics(m_pt2,m_S)) continue;
    p_overlap->FixDynamicRadius(p_integrator->X(0),p_integrator->X(1));
    wt = ( (*p_processes)()/(*p_overestimator)(m_pt2,p_integrator->Yvol()) *
	   (*p_overlap)(m_b)/m_bfac );
    if (m_ana) AnalyseWeight(wt);
  } while (wt<ran->Get());
  return true;
}

void Single_Collision_Handler::CopyPrefabricatedBlob(ATOOLS::Blob * blob) {
  if (m_prefabs.empty()) return;
  Particle_Vector * parts = m_prefabs.front()->InParticles();
  while (!parts->empty()) {
    Particle * part = parts->back();
    part->SetDecayBlob(nullptr);
    blob->AddToInParticles(part);
    parts->pop_back();
  }
  parts = m_prefabs.front()->OutParticles();
  while (!parts->empty()) {
    Particle * part = parts->back();
    part->SetProductionBlob(nullptr);
    blob->AddToOutParticles(part);
    parts->pop_back();
  }
  for (String_BlobDataBase_Map::const_iterator sdit=
	 m_prefabs.front()->GetData().begin();
       sdit!=m_prefabs.front()->GetData().end();sdit++) {
    if (sdit->first!="WeightsMap")
      blob->AddData(sdit->first,sdit->second->ClonePtr());
  }
  blob->SetId();
  blob->SetTypeSpec("AMISIC++ 1.1");
  blob->SetType(btp::Hard_Collision);
  blob->SetStatus(blob_status::needs_showers);
  Blob * del = m_prefabs.front();
  delete m_prefabs.front();
  //msg_Out()<<"  * "<<METHOD<<" deletes |"<<del<<"|, "
  //	   <<Blob::Counter()<<"/"<<Particle::Counter()<<"\n";
  m_prefabs.pop_front();
}

void Single_Collision_Handler::UpdateSandY(double s, double y) {
  ///////////////////////////////////////////////////////////////////////////
  // Updating everything that is needed for generation of next scatter
  // Setting the last pT^2 on s, as a default.
  m_lastpt2 = m_S = s;
  m_Ycms = y;
  if (m_evttype==evt_type::Perturbative) {
    p_processes->UpdateS(m_S);

    p_overestimator->UpdateS(m_S,  p_processes->PT02());
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
  /*
  msg_Out()<<"  * "<<METHOD<<" for Q^2 = "<<Q2<<", s = "<<m_S<<".\n";
  Histogram histo(0,0.0,Q2,100);
  for (long int dryrun=0;dryrun<n;dryrun++) {
    SetLastPT2(Q2);
    bool taken(false);
    while (NextScatter() && m_lastpt2>0) {
      if (!taken) {
	histo.Insert(m_lastpt2);
	taken = true;
      }
      SetLastPT2(m_pt2);
    }
  }
  histo.Finalize();
  histo.Output("True_PT2");
  msg_Out()<<"  * "<<METHOD<<": finished "<<n<<" dry runs.\n";
  */
}
