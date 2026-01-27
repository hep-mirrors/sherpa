#include "AMISIC++/Main/Amisic.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include <algorithm>
#include <cmath>
#include <fstream> // output
#include <sstream> // output

using namespace AMISIC;
using namespace ATOOLS;
using namespace std;

Amisic::Amisic() :
  m_processes(MI_Processes()), p_soft(nullptr),
  m_sigmaND_norm(1.),
  m_Nscatters(0), m_producedSoft(false), m_isMinBias(false),
  m_ana(false),
  m_n_variations(1), 
  m_mpi_scatter_count(0) // output
{}

Amisic::~Amisic() {
  if (mipars) delete mipars;
  if (p_soft) delete p_soft;
  if (m_ana) FinishAnalysis();
  // output
  if (m_b_weight_file.is_open()) {
    m_b_weight_file.close();
  }
  if (m_total_weight_file.is_open()) {
    m_total_weight_file.close();
  }
  // output
}

bool Amisic::Initialize(MODEL::Model_Base *const model,
			PDF::ISR_Handler *const isr,
			YFS::YFS_Handler *const yfs,
                        REMNANTS::Remnant_Handler * remnant_handler)
{
  msg_Info()<<"   "<<std::string(77,'=')<<"\n"
	    <<"   | "<<METHOD<<std::string(56,' ')<<"|\n";
  InitParametersAndType(isr,remnant_handler);
  InitParameterVariations();
  ///////////////////////////////////////////////////////////////////////////
  // Calculate hadronic non-diffractive cross sections, the normalization
  // for the multiple scattering probability.
  ///////////////////////////////////////////////////////////////////////////
  m_xsecs.Initialize(isr->Flav(0),isr->Flav(1),model,m_evttype);
  if (m_evttype==evt_type::Perturbative) {
    ///////////////////////////////////////////////////////////////////////////
    // Initialize the parton-level processes - currently only 2->2 scatters.
    ///////////////////////////////////////////////////////////////////////////
    m_processes.SetXSecCalculator(&m_xsecs);
    m_processes.Initialize(model,nullptr,isr,yfs);
    ///////////////////////////////////////////////////////////////////////////
    // Calculate the ratios of hard and non-diffractive cross sections.
    ///////////////////////////////////////////////////////////////////////////
    m_xsecs.CalculateXSratios(&m_processes,p_sbins);
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
    m_overestimator.Initialize(isr,&m_processes,p_sbins,&m_pint,&m_mo);
    ///////////////////////////////////////////////////////////////////////////
    // Initializing the Single_Collision_Handler which creates the next
    // scatter: it needs the processes, overestimator, interaction probability
    // and matter overlap
    ///////////////////////////////////////////////////////////////////////////
    m_singlecollision.Init(&m_processes,&m_overestimator,&m_pint,&m_mo);
    m_singlecollision.SetReweightCallback(
      [this](bool accepted, double prob) { this->AcceptRejectReweighting(accepted, prob); }
    );
  }
  else {
    p_soft = new NonPerturbative_XSecs(remnant_handler,&m_xsecs);
    p_soft->SetBeams(isr->GetBeam(0),isr->GetBeam(1));
    p_soft->CalculateSDependentCrossSections();
    ///////////////////////////////////////////////////////////////////////////
    // Initializing the Single_Collision_Handler which creates the next
    // scatter: it needs the remnants and soft processes
    ///////////////////////////////////////////////////////////////////////////
    m_singlecollision.Init(remnant_handler,p_soft);
    remnant_handler->SetType(REMNANTS::strat::simple);
  }
  /////////////////////////////////////////////////////////////////////////////
  // TODO: make sure we can inclusively generate min bias events mixing
  //       perturbative and non-perturbative (i.e. diffractive) channels.
  /////////////////////////////////////////////////////////////////////////////
  if (m_ana) InitAnalysis();
  m_isFirst   = true;
  m_isMinBias = false;
  return true;
}

void Amisic::InitParametersAndType(PDF::ISR_Handler *const isr,
				   REMNANTS::Remnant_Handler * remnants) {
  mipars = new MI_Parameters();
  bool shown = false;
  ///////////////////////////////////////////////////////////////////////////
  // Must distinguish the hard and rescatter process.  For the latter, we
  // take energies as fixed, for the former, the energies may vary (we have
  // to check the spectrum):
  // - if EPA is used the energies entering the ISR will vary,
  // - otherwise the energy is fixed.
  //
  // TODO: fix/check things up for pomerons - another interesting case
  ///////////////////////////////////////////////////////////////////////////
  m_evttype    = mipars->GetEvtType();
  m_variable_s = ( isr->Id()!=PDF::isr::bunch_rescatter &&
		   ( isr->GetBeam(0)->Type() == BEAM::beamspectrum::EPA ||
		     isr->GetBeam(1)->Type() == BEAM::beamspectrum::EPA ) );
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
  Vec4D P  = isr->GetBeam(0)->OutMomentum()+isr->GetBeam(1)->OutMomentum();
  m_S      = P.Abs2();
  m_Y      = P.Y();
  m_weight = 1.;
  size_t nbins = size_t((*mipars)["nS_bins"]);
  p_sbins      = (m_variable_s ?
		  new axis(nbins, 4.*Max(1., sqr((*mipars)("pt_min"))),
			   (1.+1.e-6)*m_S, axis_mode::log) :
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
      msg_Info()<<"   "<<string(77,'-')<<"\n"
		<<"   |  Potential error in "<<METHOD<<":"<<string(24,' ')<<"|\n"
		<<"   |  IR cutoff of MPI model "
		<<setw(5)<<setprecision(3)<<(*mipars)("pt_0")
		<<" GeV below minimal scale of PDFs."<<string(12,' ')<<"|\n"
		<<"   |  Will freeze PDFs at their minimal scale: "
		<<setw(5)<<setprecision(3)<<sqrt(isr->PDF(beam)->Q2Min())
		<<" GeV."<<string(22,' ')<<"|\n"
		<<"   "<<std::string(77,'-')<<"\n\n";
      shown = true;
    }
  }
  m_maxNscatters = size_t((*mipars)["nMaxScatters"]);
}

void Amisic::InitParameterVariations() {
  ///////////////////////////////////////////////////////////////////////////
  // Initialize the parameter variations for MPI reweighting.
  ///////////////////////////////////////////////////////////////////////////
  if (!mipars) return;
  m_sigma_nd_variations = (*mipars).GetVariationVector("SigmaND_Norm");
  m_pt0_variations      = (*mipars).GetVariationVector("pt_0");
  m_ptmin_variations    = (*mipars).GetVariationVector("pt_min");
  m_eta_variations      = (*mipars).GetVariationVector("eta");
  m_n_variations = std::max({m_sigma_nd_variations.size(),
                              m_pt0_variations.size(),
                              m_ptmin_variations.size(),
                              m_eta_variations.size(),
                              m_mo.MatterFormVariationSize(),
                              size_t(1)});
  m_sigma_nd_variations.resize(m_n_variations, m_sigma_nd_variations[0]);
  m_pt0_variations.resize(m_n_variations, m_pt0_variations[0]);
  m_ptmin_variations.resize(m_n_variations, m_ptmin_variations[0]);
  m_eta_variations.resize(m_n_variations, m_eta_variations[0]);
  
  ResetVariationWeights();

  m_xsecs.SetVariations(m_n_variations, m_sigma_nd_variations);
  m_pint.SetVariations(m_n_variations, m_sigma_nd_variations);
  m_overestimator.SetVariations(m_n_variations);

  // output
  const bool print_files = false;
  if (!print_files) return;
  std::string filename;
  if (m_n_variations == 3) {
    filename = "b_weights_rew.dat";
  } else {
    double nominal_value = m_sigma_nd_variations[0];
    std::ostringstream oss;
    oss << "b_weights_" << nominal_value << ".dat";
    filename = oss.str();
  }
  m_b_weight_file.open(filename);
  if (m_b_weight_file.is_open()) {
    if (m_n_variations == 3) {
      m_b_weight_file << "# b_value w_b_var0 w_b_var1 w_b_var2\n";
    } else {
      m_b_weight_file << "# b_value\n";
    }
    m_b_weight_file << std::scientific << std::setprecision(10);
  }
  std::string total_filename;
  if (m_n_variations == 3) {
    total_filename = "total_weights_rew.dat";
  } else {
    double nominal_value = m_sigma_nd_variations[0];
    std::ostringstream oss;
    oss << "total_weights_" << nominal_value << ".dat";
    total_filename = oss.str();
  }
  m_total_weight_file.open(total_filename);
  if (m_total_weight_file.is_open()) {
    if (m_n_variations == 3) {
      m_total_weight_file << "# n_mpi w_total_var0 w_total_var1 w_total_var2\n";
    } else {
      m_total_weight_file << "# n_mpi\n";
    }
    m_total_weight_file << std::scientific << std::setprecision(10);
  }
  // output
}

bool Amisic::InitMPIs(const double & ptmax,const double & x1,const double & x2,
		      const double & scale) {
  ///////////////////////////////////////////////////////////////////////////
  // Initialise the MPI simulation: fixing the maximal scale, take from the
  // hard event, for the downward evolution and determining a
  // kinematics-dependent impact parameter in SetB().
  ///////////////////////////////////////////////////////////////////////////
  if (!m_isFirst) return false;
  m_producedSoft = false;
  m_isFirst      = false;
  m_Nscatters    = 0;
  if (m_variable_s) UpdateForNewS();
  if (VetoEvent(scale)) return false;
  return true;
}

bool Amisic::InitMinBiasEvent(Blob_List * blobs) {
  ///////////////////////////////////////////////////////////////////////////
  // Initialise the MinBias simulation:
  // - update for new centre-of-mass energy and rapidity
  // - create a first scatter (and store a second one, if necessary) and
  //   add the blob to the blob list.
  ///////////////////////////////////////////////////////////////////////////
  if (!m_isFirst) return false;
  m_isFirst      = false;
  m_producedSoft = false;
  m_isMinBias    = true;
  m_Nscatters    = 0;
  if (m_variable_s) UpdateForNewS();
  if (!p_soft && m_singlecollision.NeedsTrig()) {
    m_weight *= m_singlecollision.PrefabricateBlob(1);
  }
  return true;
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
    m_singlecollision.SetB();
    m_singlecollision.SetLastPT2();
  }
  return true;
}

void Amisic::SetB(const double & b) { m_singlecollision.SetB(b); }

bool Amisic::FirstRescatter(Blob * blob) {
  ///////////////////////////////////////////////////////////////////////////
  // Currently only *** perturbative *** rescattering.
  // TODO: add diffractive/elastic rescattering
  ///////////////////////////////////////////////////////////////////////////
  return m_singlecollision.FirstRescatter(blob);
}

bool Amisic::FirstMinBias(Blob * blob) {
  if (m_evttype==evt_type::Perturbative) {
    if (m_xsecs.XSratio(m_S)<1.) return false;
  }
  if (m_singlecollision.FirstMinBiasScatter(blob)) {
    m_b            = m_singlecollision.B();
    m_producedSoft = m_singlecollision.Done();
    blob->UnsetStatus(blob_status::needs_minBias);
    return true;
  }
  return false;
}

bool Amisic::FirstMPI(Blob * blob) {
  ///////////////////////////////////////////////////////////////////////////
  // Refactored FirstMPI with loop at Amisic level.
  // This allows lambda ratios to be computed with the actual impact parameter
  // BEFORE Sudakov evolution, ensuring correct reweighting.
  ///////////////////////////////////////////////////////////////////////////
  UpdateForNewS();
  double pt2veto = sqr((*blob)["MI_Scale"]->Get<double>());
  if (!m_singlecollision.InitFirstMPI(blob)) return false;
  int sudakov_result;
  do {
    m_singlecollision.SelectNewB();
    m_b = m_singlecollision.B();
    ImpactParameterReweighting(m_S);
  } while (m_singlecollision.RunFirstMPISudakov(pt2veto) == 1);
  return true;
}


bool Amisic::GenerateScatter(const size_t & type,Blob * blob) {
  ///////////////////////////////////////////////////////////////////////////
  // If a next (perturbative) scatter event has been found, the pointer to
  // the respective blob is returned.
  ///////////////////////////////////////////////////////////////////////////
  if (!m_singlecollision.Done() && m_Nscatters<m_maxNscatters) {
    bool outcome = false;
    switch (type) {
    case 3:
      outcome = ((blob->Type()==btp::Signal_Process) ?
		 FirstMPI(blob) : m_singlecollision.NextScatter(blob) );
      break;
    case 2:
      if (!m_isMinBias) THROW(fatal_error,"Conflict for "+to_string(int(type))+
			      " vs "+to_string(int(m_isMinBias)));
      outcome = FirstRescatter(blob);
      break;
    case 1:
      if (!m_isMinBias) THROW(fatal_error,"Conflict for "+to_string(int(type))+
			      " vs "+to_string(int(m_isMinBias)));
      outcome = FirstMinBias(blob);
      break;
    default: THROW(fatal_error,"Unknown type: "+to_string(type));
    }
    if (outcome) {
      if (type==3 && blob->Type()!=btp::Signal_Process) {
        ++m_mpi_scatter_count; // output
      }
      AddInformationToBlob(blob);
      m_Nscatters++;
      if (m_singlecollision.Done()) {
	m_Nscatters = 0;
	if (m_ana) AnalysePerturbative(true);
      }
      return true;
    }
    if (m_ana) AnalysePerturbative(true);
  }
  return false;
}

void Amisic::UpdateForNewS() {
  ///////////////////////////////////////////////////////////////////////////
  // Update if first scatter with variable c.m. energy (e.g. collisions with
  // pomerons or resolved photons): update c.m. energy in processes,
  // re-calculate non-diffractive cross sections for normalisation, adjust
  // prefactors for radii in matter overlaps etc., etc..
  // TODO: will have to check if we need another longitudinal boost (hence
  //       the Y)
  ///////////////////////////////////////////////////////////////////////////
  Vec4D P(0.,0.,0.,0.);
  for (size_t beam=0;beam<2;beam++) P += m_singlecollision.InMomentum(beam);
  m_S = P.Abs2();
  m_Y = P.Y();
  if (m_evttype==evt_type::Perturbative) m_mo.SetKRadius(m_pint.K(m_S));
  m_singlecollision.UpdateSandY(m_S, m_Y);
}

void Amisic::AddInformationToBlob(ATOOLS::Blob * blob) {
  if (m_evttype==evt_type::Perturbative) {
    m_pt2     = m_singlecollision.LastPT2();
    double x1 = 0., x2 = 0.;
    if ((*blob)["PDFInfo"]!=NULL) {
      x1 = (*blob)["PDFInfo"]->Get<PDF_Info>().m_x1;
      x2 = (*blob)["PDFInfo"]->Get<PDF_Info>().m_x2;
    }
    blob->SetPosition(m_mo.SelectPositionForScatter(m_b,x1,m_pt2,x2,m_pt2));
    if (m_ana) AnalysePerturbative(false,blob);
  }
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
  m_singlecollision.Reset();
  m_isFirst   = true;
  m_isMinBias = false;
  m_Nscatters = 0;
}

void Amisic::Reset() {
  m_weight    = 1.;
  m_Nscatters = 0;
  m_singlecollision.Reset();
}

void Amisic::InitAnalysis() {
  m_nev = m_nscatters = m_nscat = 0;
  m_histos[string("N_scatters")] = new Histogram(0,0,50,50);
  m_histos[string("N_B0_1")]     = new Histogram(0,0,50,50);
  m_histos[string("N_B1_2")]     = new Histogram(0,0,50,50);
  m_histos[string("N_B2_4")]     = new Histogram(0,0,50,50);
  m_histos[string("N_B4_8")]     = new Histogram(0,0,50,50);
  m_histos[string("B")]          = new Histogram(0,0,10,100);
  m_histos[string("P_T1")]       = new Histogram(0,0,100,100);
  m_histos[string("Y1")]         = new Histogram(0,-10,10,10);
  m_histos[string("Delta_Y1")]   = new Histogram(0,0,10,10);
  m_histos[string("P_T2")]       = new Histogram(0,0,100,100);
  m_histos[string("Y2")]         = new Histogram(0,-10,10,10);
  m_histos[string("Delta_Y2")]   = new Histogram(0,0,10,10);
}

void Amisic::FinishAnalysis() {
  if (m_evttype==evt_type::Perturbative) {
    msg_Info()<<METHOD<<": <nscatters> = "
	      <<double(m_nscat)/double(m_nev)<<" vs. "
	      <<m_xsecs.XSratio()<<".\n";
  }
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
    m_nscat += int(m_nscatters);  m_nev++;
    m_histos[string("N_scatters")]->Insert(double(m_nscatters)+0.5);
    m_histos[string("B")]->Insert(m_b*rpa->hBar()*rpa->c()*1.e12);
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

void Amisic::ResetVariationWeights() {
  ///////////////////////////////////////////////////////////////////////////
  // Reset variation weights to 1.0 for a new event.
  ///////////////////////////////////////////////////////////////////////////
  m_variation_weights.resize(m_n_variations);
  std::fill(m_variation_weights.begin(), m_variation_weights.end(), 1.);
  m_b_weights.assign(m_n_variations, 1.);
  m_lambda_ratios.assign(m_n_variations, 1.);
  m_sudakov_weights.assign(m_n_variations, 1.);
  m_mpi_scatter_count = 0; // output
}

void Amisic::ImpactParameterReweighting(const double & s) {
  ///////////////////////////////////////////////////////////////////////////
  // Cache per-event information needed for Sudakov reweighting
  // and compute impact-parameter reweighting factors.
  ///////////////////////////////////////////////////////////////////////////
  if (m_variation_weights.empty()) return;

  // b is selected according to d^2b O(b)
  // w_b = O_var(b) / O_nom(b) = lambda_var / lambda_nom = lambda_ratio
  m_mo.SetMatterFormVariationIndex(0);
  const double K_nom = m_pint.K(s);
  const double overlap_nom = m_mo.EvaluateAt(m_b, K_nom);

  // Compute lambda(b) ratios for each variation
  for (size_t ivar=0; ivar<m_n_variations; ++ivar) {
    m_mo.SetMatterFormVariationIndex(ivar);
    const double K_var = m_pint.KVariation(s, ivar);
    const double overlap_var = m_mo.EvaluateAt(m_b, K_var);
    double lambda_ratio = overlap_var / overlap_nom;
    if (!std::isfinite(lambda_ratio) || lambda_ratio<=0.) lambda_ratio = 1.;
    m_lambda_ratios[ivar] = lambda_ratio;
    
    // This accounts for all b values tried in FirstMPI loop
    m_b_weights[ivar] *= lambda_ratio;
  }
  m_mo.SetMatterFormVariationIndex(0);
}

void Amisic::AcceptRejectReweighting(const bool accepted, const double prob_nom) {
  ///////////////////////////////////////////////////////////////////////////
  // Callback from Single_Collision_Handler for each accept/reject event
  // during Sudakov evolution of MPI scatters.
  // For accepted scatter: multiply weight by p_var/p_nom
  // For rejected scatter: multiply weight by (1-p_var)/(1-p_nom)
  ///////////////////////////////////////////////////////////////////////////
  if (m_sudakov_weights.empty()) return;

  const double xs_nom = (m_processes)();

  for (size_t ivar=0; ivar<m_n_variations; ++ivar) {
    // p_nom = dsigma/dpt2 * overlap(b,K_var) / overestimator
    // The ratio is: p_var/p_nom = overlap_var/overlap_nom * (dsigma/dpt2)_nom / (dsigma/dpt2)_var
    //                           = lambda_ratio * xs_ratio
    // p_var = 0 for pt2 <= (ptmin2)_var
    double prob_var = 0.;
    if (m_singlecollision.PT2() > sqr(m_ptmin_variations[ivar])) {
      (&m_processes)->SetPT02(sqr(m_pt0_variations[ivar]));
      const double xs_var = (m_processes)();
      (&m_processes)->SetPT02(sqr(m_pt0_variations[0]));
      prob_var = prob_nom * m_lambda_ratios[ivar] * xs_var / xs_nom;
    }
    if (accepted) {
      // Accepted: weight *= p_var / p_nom
      const double ratio = prob_var / prob_nom;
      if (std::isfinite(ratio) && ratio >= 0.) {
        m_sudakov_weights[ivar] *= ratio;
      }
    } else {
      // Rejected: weight *= (1 - p_var) / (1 - p_nom)
      const double ratio = (1. - prob_var) / (1. - prob_nom);
      if (std::isfinite(ratio) && ratio >= 0.) {
        m_sudakov_weights[ivar] *= ratio;
      }
    }
  }
}

void Amisic::ApplyVariationWeights(ATOOLS::Blob * blob) {
  ///////////////////////////////////////////////////////////////////////////
  // Compute and apply variation weights to the event.
  // The total weight is:
  // w_total = w_b * w_(n|b) = w_b * w_sudakov
  // where w_sudakov accounts for the changed MPI multiplicity distribution
  // and w_b accounts for the changed impact parameter sampling distribution.
  ///////////////////////////////////////////////////////////////////////////
  if (m_variation_weights.empty()) return;

  // output
  std::vector<double> w_b_vec(m_variation_weights.size());
  std::vector<double> w_sudakov_vec(m_variation_weights.size());
  std::vector<double> w_total_vec(m_variation_weights.size());

  for (size_t ivar=0; ivar<m_variation_weights.size(); ++ivar) {
    // The impact-parameter reweighting factor w_b is accumulated in
    // ImpactParameterReweighting() calls, stored in m_b_weights
    const double w_b = m_b_weights[ivar];
    
    // During the MPI evolution, the Sudakov weights are accumulated in
    // AcceptRejectReweighting() calls, stored in m_sudakov_weights
    const double w_sudakov = m_sudakov_weights[ivar];

    // w_total = w_b * w_(n|b) = w_b * w_sudakov
    const double w_total = w_b * w_sudakov;
    m_variation_weights[ivar] *= w_total;

    // output
    w_b_vec[ivar]       = w_b;
    w_total_vec[ivar]   = w_total;
    // output
  }
  // output
  if (m_b_weight_file.is_open()) {
    m_b_weight_file << m_b;
    if (m_variation_weights.size() == 3) {
      m_b_weight_file << " " << w_b_vec[0] << " " << w_b_vec[1] << " " << w_b_vec[2];
    }
    m_b_weight_file << "\n";
  }
  if (m_total_weight_file.is_open()) {
    m_total_weight_file << m_mpi_scatter_count;
    if (m_variation_weights.size() == 3) {
      m_total_weight_file << " " << w_total_vec[0] << " " << w_total_vec[1] << " " << w_total_vec[2];
    }
    m_total_weight_file << "\n";
  }
  // output
  
  if (blob == NULL) {
    return;
  }
  auto &wgt_map = (*blob)["WeightsMap"]->Get<Weights_Map>();
  
  for(size_t ivar=0; ivar<m_n_variations; ++ivar) {
    const std::string name {"v" + std::to_string(ivar)};
    const double wgt = m_variation_weights[ivar];
    wgt_map["MPI"][name] = wgt;
  }
  ResetVariationWeights();
}

