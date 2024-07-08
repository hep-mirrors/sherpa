#include "AMISIC++/Main/Amisic.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace AMISIC;
using namespace ATOOLS;
using namespace std;

Amisic::Amisic() :
  m_processes(MI_Processes()), p_soft(nullptr),
  m_sudakovarg(Sudakov_Argument(&m_processes)),
  m_sigmaND_norm(1.),
  m_isMinBias(false), m_ana(true)
{}

Amisic::~Amisic() {
  if (mipars) delete mipars;
  if (p_soft)      delete p_soft;
  if (m_ana) FinishAnalysis();
}

bool Amisic::Initialize(MODEL::Model_Base *const model,
			PDF::ISR_Handler *const isr,
			YFS::YFS_Handler *const yfs,
                        REMNANTS::Remnant_Handler * remnant_handler)
{
  InitParametersAndType(isr,remnant_handler);
  Vec4D P = isr->GetBeam(0)->OutMomentum()+isr->GetBeam(1)->OutMomentum();
  m_S     = P.Abs2();
  m_Y     = P.Y();
  ///////////////////////////////////////////////////////////////////////////
  // Calculate hadronic non-diffractive cross sections, the normalization
  // for the multiple scattering probability.
  ///////////////////////////////////////////////////////////////////////////
  m_xsecs.Initialize(isr->Flav(0),isr->Flav(1),model);
  ///////////////////////////////////////////////////////////////////////////
  // Initialize the parton-level processes - currently only 2->2 scatters.
  ///////////////////////////////////////////////////////////////////////////
  m_processes.SetXSecCalculator(&m_xsecs);
  m_processes.SetMatterOverlap(&m_mo);
    if (mipars->GetEvtType()==evt_type::Perturbative) {
    m_processes.Initialize(model,nullptr,isr,yfs);
    ///////////////////////////////////////////////////////////////////////////
    // Calculate the ratios of hard and non-diffractive cross sections.
    ///////////////////////////////////////////////////////////////////////////
    m_xsecs.CalculateXSratios(&m_processes,p_sbins);
    ///////////////////////////////////////////////////////////////////////////
    // Prepare the integral for the "Sudakov form factor", Eq. (SZ, 37)
    ///////////////////////////////////////////////////////////////////////////
    m_sudakovarg.Initialize(p_sbins,m_mo.GetBBins());
    ///////////////////////////////////////////////////////////////////////////
    // Initialize everything to do with the inpact parameter dependence.
    // For dynamic matter overlaps this includes fixing the effective radius
    ///////////////////////////////////////////////////////////////////////////
    m_pint.Initialize(&m_mo,&m_processes,p_sbins);
    ///////////////////////////////////////////////////////////////////////////
    // Initialize the Over_Estimator - mainly fixing an effective prefactor
    // to allow for a quick'n'dirty fix to create fast estimates of the next
    // scatter's pT^2.  The logic for the overestimator is based on
    // - t-channel dominance allowing to approximate differential cross
    //   sections as
    //          dsigma ~ dp_T dy f(x_1) f_(x_2) C_1 C_2 as(t)^2/(t+t_0)^2
    //   where C_1/2 are the colour factors of the incoming partons, x_1/2
    //   are the Bjorken-x;
    // - assuming that the product of the PDFs f(x_1)f(x_2) is largest for
    //   mid-rapidity where x_1 and x_2 are identical
    ///////////////////////////////////////////////////////////////////////////
    m_overestimator.Initialize(isr,&m_processes,
                               p_sbins,m_sudakovarg.GetPT2Bins());
    ///////////////////////////////////////////////////////////////////////////
    // Initializing the Single_Collision_Handler which creates the next
    // scatter: it needs the remnants, processes, and the overestimator
    ///////////////////////////////////////////////////////////////////////////
    m_singlecollision.Init(&m_processes,&m_overestimator,&m_pint,&m_mo);
  } else {
      p_soft = new NonPerturbative_XSecs(remnant_handler,p_xsecs);
      p_soft->SetBeams(isr->GetBeam(0),isr->GetBeam(1));
      p_soft->CalculateSDependentCrossSections();
      // Initializing the Single_Collision_Handler which creates the next scatter: it needs
      // the remnants, processes, and the overestimator
      m_singlecollision.Init(remnant_handler,p_soft);
      remnant_handler->SetType(REMNANTS::strat::simple);
  }
  if (m_ana) InitAnalysis();
  m_isFirst   = true;
  m_isMinBias = false;
  return true;
}

void Amisic::InitParametersAndType(PDF::ISR_Handler *const isr,
				   REMNANTS::Remnant_Handler * remnants) {
  mipars = new MI_Parameters();
  bool shown = false;
  // Must distinguish the hard and rescatter process.  For the latter, we
  // take energies as fixed, for the former, the energies may vary (we have
  // to check the spectrum):
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
  else {
    msg_Error() <<METHOD<<": unknown multiple interaction model for "
		<<isr->Flav(0) << " and " << isr->Flav(1) << "\n"
		<<"   Will terminate the run.\n";
    THROW(normal_exit,"Inconsistent set-up in AMISIC.");
  }
  Vec4D P = isr->GetBeam(0)->OutMomentum()+isr->GetBeam(1)->OutMomentum();
  m_S     = P.Abs2();
  m_Y     = P.Y();
  size_t nbins = size_t((*mipars)("nS_bins"));
  p_sbins      = (m_variable_s ?
		  new axis(nbins, 4.*Max(1., sqr(2.*(*mipars)("pt_min"))), m_S,
			   axis_mode::linear) :
		  new axis(1, m_S, (1.+1.e-6)*m_S, axis_mode::linear));

  size_t nspace = 15;
  if (m_variable_s) {
    msg_Info()<<"   | Variable centre-of mass energies for "; nspace-=3;
  }
  else { msg_Info()<<"   | Fixed centre-of mass energies for "; }
  switch (m_type) {
  case mitype::hadron_hadron: msg_Info()<<"hadron-hadron collisions."; break;
  case mitype::gamma_hadron:  msg_Info()<<"photon-hadron collisions."; break;
  case mitype::gamma_gamma:   msg_Info()<<"photon-photon collisions."; break;
  default:                    msg_Info()<<"unknown particles.       "; break;
  }
  msg_Info()<<std::string(nspace,' ')<<"|\n"<<"   "
	    <<std::string(77,'=')<<"\n\n";
  ///////////////////////////////////////////////////////////////////////////
  // Initializing the matter overlap, specific for the UE/MB model, cf.
  // Eqs. (SZ, 21-22) and (CS, 11), respectively.
  ///////////////////////////////////////////////////////////////////////////
  m_mo.Initialize(remnants,isr);
  m_variable_b = m_mo.IsDynamic();
if (m_variable_s && m_variable_b)
    THROW(fatal_error,
	  "B-dependent form factors not implemented yet for variable s.");

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

bool Amisic::InitMPIs(const double & ptmax,const double & x1,const double & x2,
		      const double & scale) {
  ///////////////////////////////////////////////////////////////////////////
  // Initialise the MPI simulation: fixing the maximal scale, take from the
  // hard event, for the downward evolution and determining a
  // kinematics-dependent impact parameter in SetB().
  ///////////////////////////////////////////////////////////////////////////
  if (m_isFirst) {
    m_isFirst = false;
    if (m_variable_s) UpdateForNewS();
    m_mo.SetKRadius(m_pint.K(m_S));
  }
  if (!VetoEvent(scale)) return true;
  return false;
}

bool Amisic::InitMinBiasEvent(Blob_List * blobs) {
  ///////////////////////////////////////////////////////////////////////////
  // Initialise the MinBias simulation:
  // - update for new centre-of-mass energy and rapidity
  //   TODO: will have to use rapidity, postponed until pp works (again)
  // - create a first scatter (and store a second one, if necessary) and
  //   add the blob to the blob list.
  ///////////////////////////////////////////////////////////////////////////
  if (m_isFirst) {
    m_isFirst   = false;
    m_producedSoft = false;
    m_isMinBias = true;
    if (m_variable_s) UpdateForNewS();
    SetMaxScale(rpa->gen.Ecms()/2.);
    if (!p_soft) SetFirstB();
    if (m_mo.IsDynamic()) m_singlecollision.PrefabricateBlob();
    m_mo.SetKRadius(m_pint.K(m_S));
    return 1;
  }
  return 0;
}

bool Amisic::InitRescatterEvent() {
  ///////////////////////////////////////////////////////////////////////////
  // Initialise the MinBias simulation: fixing the maximal scale = S/4, the
  // kinematic limit, for the downward evolution and take the impact parameter
  // from the already existing scatter.
  // TODO: we may have to check the logical flow here.
  ///////////////////////////////////////////////////////////////////////////
  if (m_isFirst) {
    m_isFirst   = false;
    m_isMinBias = true;
    if (m_variable_s) UpdateForNewS();
    // SetFirstB(m_singlecollision.B());
    SetMaxScale2(sqr(rpa->gen.Ecms()/2.));
    m_singlecollision.SetLastPT2();
    return true;
  }
  return false;
}

void Amisic::UpdateForNewS() {
  ///////////////////////////////////////////////////////////////////////////
  // Update if first scatter with variable c.m. energy (e.g. collisions with
  // pomerons or resolved photons): update c.m. energy in processes,
  // re-calculate non-diffractive cross sections for normalisation, adjust
  // prefactors for overestimators, etc..
  // TODO: will have to check if we need another longitudinal boost (hence
  //       the Y)
  ///////////////////////////////////////////////////////////////////////////
  Vec4D P(0.,0.,0.,0.);
  for (size_t beam=0;beam<2;beam++) {
    P += m_singlecollision.InMomentum(beam);
  }
  m_S = P.Abs2();
  m_Y = P.Y();
  m_singlecollision.UpdateSandY(m_S, m_Y);
}

void Amisic::SetB(const double & x1,const double & x2,const double & scale) {
  ///////////////////////////////////////////////////////////////////////////
  // It is obtained by interpolation from a look-up table in the
  // Interaction_Probability.
  ///////////////////////////////////////////////////////////////////////////
  msg_Out()<<METHOD<<"\n";
  exit(1);
}

bool Amisic::FirstMinBias(Blob * blob) {
  //msg_Out()<<METHOD<<"(s = "<<m_S<<"): should produce only "
  //	   <<m_xsecs.XSratio(m_S)<<" hard events.\n";
  if (m_xsecs.XSratio(m_S)<1.) return false;
  if (m_singlecollision.FirstMinBiasScatter(blob)) {
    m_b    = m_singlecollision.B();
    m_bfac = m_singlecollision.BFac();
    blob->UnsetStatus(blob_status::needs_minBias);
    return true;
  }
  return false;
}

bool Amisic::FirstMPI(Blob * blob) {
  if (m_singlecollision.FirstMPI(blob)) {
    m_b    = m_singlecollision.B();
    m_bfac = m_singlecollision.BFac();
    return true;
  }
  return false;
}

bool Amisic::GenerateScatter(const size_t & type,Blob * blob) {
  // TODO_MERGE check this function implementation
  ///////////////////////////////////////////////////////////////////////////
  // If a next (perturbative) scatter event has been found, the pointer to
  // the respective blob is returned.
  // TODO: we may want to think about something for rescatter events -
  //       but this is for future work.
  ///////////////////////////////////////////////////////////////////////////
  Blob * blob = nullptr;
  if (force) m_singlecollision.UpdateSandY(m_S, m_Y);
  blob = m_singlecollision.NextScatter(force);
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
    blob->UnsetStatus(blob_status::needs_minBias);
    return blob;
  }
  if (m_ana && mipars->GetEvtType()==evt_type::Perturbative) AnalysePerturbative(true);
  return nullptr;
  // old version below
  /*if (!m_singlecollision.Done()) {
    //msg_Out()<<"*** "<<METHOD<<"("<<type<<") for "<<blob->Id()<<" |"<<blob<<"| ,"
    //	     <<Blob::Counter()<<"/"<<Particle::Counter()<<"\n";
    //     <<(*blob)<<"\n";
    bool outcome = false;
    switch (type) {
    case 3:
      outcome = ( (blob->Type()==btp::Signal_Process) ?
		  FirstMPI(blob) :
		  m_singlecollision.NextScatter(blob) );
      break;
    case 2:
      msg_Out()<<METHOD<<": rescatter not yet implemented.\n";
      exit(1);
      break;
    case 1:
      if (!m_isMinBias) THROW(fatal_error,"Conflict for "+to_string(int(type))+
			      " vs "+to_string(int(m_isMinBias)));
      outcome = FirstMinBias(blob);
      break;
    default: THROW(fatal_error,"Unknown type: "+to_string(type));
    }
    if (outcome) {
      AddInformationToBlob(blob);
      if (m_singlecollision.Done() && m_ana) {
	msg_Out()<<"*** "<<METHOD<<"(1) with "<<m_nscatters<<"\n";
	Analyse(true);
      }
      //msg_Out()<<"*** "<<METHOD<<" successful for blob |"<<blob<<"|, "
      //       <<Blob::Counter()<<"/"<<Particle::Counter()<<"\n";
      return true;
    }
    //msg_Out()<<"*** "<<METHOD<<" couldn't add to the blob |"<<blob<<"|, "
    //	     <<Blob::Counter()<<"/"<<Particle::Counter()<<"\n";
    msg_Out()<<"*** "<<METHOD<<"(2) with "<<m_nscatters<<"\n";
    if (m_ana) Analyse(true);
  }
  //else {
  //msg_Out()<<"*** "<<METHOD<<" already done, obsolete blob |"<<blob<<"|, "
  //	     <<Blob::Counter()<<"/"<<Particle::Counter()<<"\n";
  //}
  return false;
}

void Amisic::AddInformationToBlob(ATOOLS::Blob * blob) {
  m_pt2     = m_singlecollision.LastPT2();
  double x1 = 0., x2 = 0.;
  if ((*blob)["PDFInfo"]!=NULL) {
    x1 = (*blob)["PDFInfo"]->Get<PDF_Info>().m_x1;
    x2 = (*blob)["PDFInfo"]->Get<PDF_Info>().m_x2;
  }
  blob->SetPosition(m_mo.SelectPositionForScatter(m_b,x1,m_pt2,x2,m_pt2));
  if (m_ana && mipars->GetEvtType()==evt_type::Perturbative) AnalysePerturbative(false,blob);
   */
}

bool Amisic::VetoEvent(const double & scale) const {
  ///////////////////////////////////////////////////////////////////////////
  // So far this has not been properly filled.
  // TODO: maybe we need to add another weight in the spirit of
  //       the (1-x1-y1-...) seen in the double-parton scattering formalism.
  ///////////////////////////////////////////////////////////////////////////
  if (scale<0.) return true;
  return false;
}

Cluster_Amplitude * Amisic::ClusterConfiguration(Blob * blob) {
  ///////////////////////////////////////////////////////////////////////////
  // Simplified version of cluster amplitude generation (only 2->2 scatter
  // so far) to connect with the parton shower.
  ///////////////////////////////////////////////////////////////////////////
  Cluster_Amplitude * ampl = Cluster_Amplitude::New();
  CreateAmplitudeLegs(ampl,blob);
  double muf2 = (*blob)["Factorization_Scale"]->Get<double>(), muq2 = muf2;
  double mur2 = (*blob)["Renormalization_Scale"]->Get<double>();
  ampl->SetNIn(2);
  ampl->SetMuR2(mur2);
  ampl->SetMuF2(muf2);
  ampl->SetMuQ2(muq2);
  ampl->SetKT2(muf2);
  ampl->SetMu2(mur2);
  ampl->SetOrderEW(0);
  ampl->SetOrderQCD(2);
  ampl->SetMS(&m_processes);
  return ampl;
}

void Amisic::CreateAmplitudeLegs(Cluster_Amplitude * ampl,Blob * blob) {
  for (size_t i(0);i<blob->NInP()+blob->NOutP();++i) {
    size_t     id(1<<ampl->Legs().size());
    Particle * part(blob->GetParticle(i));
    ColorID    col(part->GetFlow(1),part->GetFlow(2));
    if (i<blob->NInP()) ampl->CreateLeg(-part->Momentum(),part->Flav().Bar(),
					col.Conj(),id);
    else ampl->CreateLeg(part->Momentum(),part->Flav(),col,id);
    ampl->Legs().back()->SetStat(0);
  }
}

void Amisic::CleanUpMinBias() {
  SetMaxEnergies(rpa->gen.PBeam(0)[0],rpa->gen.PBeam(1)[0]);
  SetMaxScale2(sqr(rpa->gen.Ecms()/2.));
  m_isFirst   = true;
  m_isMinBias = false;
  m_singlecollision.Reset();
}

void Amisic::Reset() {
  m_singlecollision.Reset();
}

void Amisic::InitAnalysis() {
  m_Nev = m_nscatters = m_Nscat = 0;
  m_histos[string("N_scatters")] = new Histogram(0,0,50,50);
  m_histos[string("N_B0_1")]     = new Histogram(0,0,50,50);
  m_histos[string("N_B1_2")]     = new Histogram(0,0,50,50);
  m_histos[string("N_B2_4")]     = new Histogram(0,0,50,50);
  m_histos[string("N_B4_8")]     = new Histogram(0,0,50,50);
  m_histos[string("B")]          = new Histogram(0,0,10,100);
  m_histos[string("Bfac")]       = new Histogram(0,0,10,100);
  m_histos[string("P_T1")]       = new Histogram(0,0,100,100);
  m_histos[string("Y1")]         = new Histogram(0,-10,10,10);
  m_histos[string("Delta_Y1")]   = new Histogram(0,0,10,10);
  m_histos[string("P_T2")]       = new Histogram(0,0,100,100);
  m_histos[string("Y2")]         = new Histogram(0,-10,10,10);
  m_histos[string("Delta_Y2")]   = new Histogram(0,0,10,10);
}

void Amisic::FinishAnalysis() {
  msg_Info()<<METHOD<<": <nscatters> = "
	    <<double(m_Nscat)/double(m_Nev)<<" vs. "
	    <<m_xsecs.XSratio()<<".\n";
  Histogram * histo;
  string name;
  for (map<string,Histogram *>::iterator hit=m_histos.begin();
       hit!=m_histos.end();hit++) {
    histo = hit->second;
    name  = string("MPI_Analysis/")+hit->first+string(".dat");
    histo->Finalize();
    histo->Output(name);
    delete histo;
  }
  m_histos.clear();
}

void Amisic::AnalysePerturbative(const bool & last,Blob * blob) {
  if (!last) {
    if (!blob) return;
    Vec4D mom[2];
    for (size_t i=0;i<2;i++) mom[i] = blob->OutParticle(i)->Momentum();
    if (m_nscatters==0) {
      m_histos[string("P_T1")]->Insert(sqrt(mom[0].PPerp2()));
      m_histos[string("Y1")]->Insert(mom[0].Y());
      m_histos[string("Y1")]->Insert(mom[1].Y());
      m_histos[string("Delta_Y1")]->Insert(dabs(mom[0].Y()-mom[1].Y()));
    }
    if (m_nscatters==1) {
      m_histos[string("P_T2")]->Insert(sqrt(mom[0].PPerp2()));
      m_histos[string("Y2")]->Insert(mom[0].Y());
      m_histos[string("Y2")]->Insert(mom[1].Y());
      m_histos[string("Delta_Y2")]->Insert(dabs(mom[0].Y()-mom[1].Y()));
    }
    m_nscatters++;
  }
  if (last) {
    m_Nscat += int(m_nscatters);  m_Nev++;
    m_histos[string("N_scatters")]->Insert(double(m_nscatters)+0.5);
    m_histos[string("B")]->Insert(m_b*rpa->hBar()*rpa->c()*1.e12);
    m_histos[string("Bfac")]->Insert(m_bfac);
    if (m_b<=1.)
      m_histos[string("N_B0_1")]->Insert(double(m_nscatters)+0.5);
    else if (m_b<=2.)
      m_histos[string("N_B1_2")]->Insert(double(m_nscatters)+0.5);
    else if (m_b<=4.)
      m_histos[string("N_B2_4")]->Insert(double(m_nscatters)+0.5);
    else if (m_b<=8.)
      m_histos[string("N_B4_8")]->Insert(double(m_nscatters)+0.5);
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


/*
  void Multiple_Interactions::ExtractMPIStartingConditions(ATOOLS::Blob * blob) {
  ////////////////////////////////////////////////////////////////////////////
  // Extracting relevant information to start MPIs from the *** first ***
  // signal blob in the event (I assume we will be able to produce hard DPIs):
  // m_ptmax  = maximal pt allowed for the MPIs (possibly modified by external
  //            factor): Trivially the factorization scale - plus the
  //	        renormalization scale if the two do not coincide.
  // m_x[0,1] = Bjorken-x of the two incoming particles.
  // m_scale  = factorisation scale of the process.
  // The latter two are used if we have x (and scale)-dependent hadron radii.
  ////////////////////////////////////////////////////////////////////////////
  Blob_Data_Base * facinfo = (*blob)["Factorization_Scale"];
  Blob_Data_Base * reninfo = (*blob)["Renormalization_Scale"];
  if (facinfo==NULL || reninfo==NULL)
    THROW(fatal_error,"No starting scale info in signal blob");
  double ptfac = sqrt(facinfo->Get<double>());
  double ptren = sqrt(reninfo->Get<double>());
  m_ptmax = ptfac+(!IsZero(ptfac-ptren)?ptren:0.);
  m_scale = ptfac;
  for (size_t i=0;i<2;i++) {
    m_x[i] = ( blob->InParticle(i)->Momentum()[0]/
	       p_activeMI->Remnants()->GetRemnant(i)->GetBeam()
	       ->InMomentum()[0] );
  }
}


*/
