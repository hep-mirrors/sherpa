#include "CFPSHOWER++/Shower/Shower.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "CFPSHOWER++/Calculators/SF_Base.H"
#include "CFPSHOWER++/Shower/Cluster_Definitions.H"
#include "CFPSHOWER++/Tools/CFP_Parameters.H"
#include "CFPSHOWER++/Tools/Kernel_Constructor.H"
#include "CFPSHOWER++/Tools/Kernel_Info.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Single_Vertex.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Default_Reader.H"
#include "ATOOLS/Org/My_Limits.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

using namespace CFPSHOWER;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;

Shower::Shower():
  p_cluster(new Cluster_Definitions(this))
{}

Shower::~Shower()
{
  Reset();
  for (size_t i=0;i<5;i++) {
    while (!p_kernels[i]->empty()) {
      delete p_kernels[i]->begin()->second;
      p_kernels[i]->erase(p_kernels[i]->begin());
    }
  }
  delete p_cluster;
}

void Shower::Reset() {
  while (!m_splittings.empty()) { delete m_splittings.back(); m_splittings.pop_back(); }
}

bool Shower::Evolve(Configuration * config) {
  p_config = config;
  // Initialisation - will have to be shifted to a new method (or some Reset()).
  m_weight = 1.; m_nem = 0;
  if (p_config->T0()<0.) p_config->SetT0(m_t0min);
  // Logic here:
  // 1. As long as we find possible splittings (p_winner) and we have not exhausted
  //    the number of possible emissions (m_nmax_em>m_nem), we continue with the
  //    evolution of the configuration, in step 2.
  // 2. We iterate over all partons in the configuration that have not yet been
  //    switched off and check if they can evolve in Evolve(Parton *):
  //    - If successful this yields a winning splitting, p_winner.  Then we add
  //      the combined acceptance weight for the splitter and rejection weights for
  //      all splittings that failed for this and other partons in the configuration,
  //      in AddWeight().  The winner sets the maximal t for the next evolution step,
  //      and we have to reset the minimal t, t0.  Finally we perform the actual
  //      splitting, in PerformSplitting(), by adding the splitting products and
  //      updating kinematics and colours.  
  //    - If we are not successful we go to step 3.
  // 3. If the evolution of the configuration is over we add a final weight,
  //    consisting of all rejection weights of failed splittings.
  do {
    p_winner = NULL;
    size_t i=0;
    for (Parton_List::reverse_iterator pit=p_config->rbegin();
	 pit!=p_config->rend();pit++) {
      if ((*pit)->On() && !Evolve(*pit)) continue;
    }
    if (p_winner) {
      AddWeight(p_winner->t(0));
      p_config->SetT(p_winner->t(0));
      p_config->SetT0(m_t0min);
      if (PerformSplitting()) {
	m_splittings.push_back(p_winner);
	m_nem++;
      }
      else {
	for (Parton_List::iterator pit=p_config->begin();pit!=p_config->end();pit++) {
	  (*pit)->ClearWeights();
	}
	delete p_winner;
      }
    }
  } while (p_winner && m_nem<m_nmax_em);
  AddWeight(p_config->T0());
  //msg_Out()<<METHOD<<" finished, weight = "<<m_weight<<", "<<m_nem<<" emissions.\n"
  // 	   <<"#########################################################\n"
  //	   <<"#########################################################\n"
  //	   <<"#########################################################\n";
  //if (dabs(m_weight-1.)>1.e-6) {
    //msg_Out()<<METHOD<<" fails with weight = "<<(dabs(m_weight-1)*1.e6)<<" * 10^-6. "
    //	     <<"Will continue the run & hope for the best.\n";
    //exit(1);
  //}
  return true;
}

void Shower::AddWeight(const double & t) {
  double stepweight = 1., partweight;
  for (Parton_List::iterator pit=p_config->begin();pit!=p_config->end();pit++) {
    Parton * part = (*pit);
    stepweight *= partweight = part->GetWeight(t);
    part->ClearWeights();
  }
  m_weight *= stepweight;
}

bool Shower::Evolve(Parton * splitter) {
  // Initialise the integrated splitting kernel, summed over all spectators (at
  // least one in QCD, and maximally two).  If splittings are kinematically allowed,
  // the sum will be larger than one and we will generate a test splitting through
  // GenerateTestSplitting().  This test includes splitting parameters t and z,
  // plus corresponding kinematics information (y, z, as well as all masses
  // involved), from which we also construct the four-momenta of the decay products
  // and the spectator.  The kinematics will subsequently be fully realised
  // in PerformSplitting(), by adding the decay products to the configuration, by
  // switching off the splitter, and by updating the spectator kinematics and the
  // colour connections.
  
  if (splitter->GetSpectators()->size()<=0) return false;
  bool success = false;
  double sum  = InitialiseIntegrals(splitter);
  if (sum>0.) {
    Splitting * split = GenerateTestSplitting(splitter,sum);
    //msg_Out()<<METHOD<<" yields: "<<split<<".\n";
    if (split) {
      success = true;
      // we found a splitting - if we have not found a possible splitting yet,
      // this splitting is the winner.  If we already have a winner, there are
      // two possibilities:
      // 1. the t of the current splitting is larger than the winner to date.
      //    In this case, we delete the winner and replace it with the current trial
      //    splitting (this should become obsolete now).
      // 2. the t of the current splitter is below the current winner's t.  Then
      //    we just delete the current splitting.  
      if (p_winner==NULL) p_winner = split;
      else if (split->t(0) > p_winner->t(0)) {
	delete p_winner; p_winner = split;
      }
      else delete split;
      // If we have a winning splitter, all contenders must have a larger t, so
      // the winning t is the lower limit t0 of the next round of attempts.
      if (p_winner) { p_config->SetT0(p_winner->t(0)); }
    }
  }
  return success;
}

double Shower::InitialiseIntegrals(Parton * splitter) {
  const Parton_List * spectators = splitter->GetSpectators();
  while (!m_integrals.empty()) {
    m_integrals.back().clear(); m_integrals.pop_back();
  }
  m_integrals.resize(spectators->size());
  m_splitkernels.resize(spectators->size());
  double sum = 0.;
  size_t i = 0;
  for (Parton_List::const_iterator spit=spectators->begin();
       spit!=spectators->end();spit++) {
    Parton * spectator = (*spit);
    // Initialise a container (Splitting) split holding the information defining the
    // potential splitting.
    Splitting split(splitter,spectator,p_config->T(),p_config->T0());
    // Select the list of applicable kernels depending on the kinematic configuration
    // and the flavour of the splitter.  
    kernel_type::code type = GetCode((splitter->Beam()>0),(spectator->Beam()>0));
    map<Flavour, Kernels *>::iterator kit =
      p_kernels[int(type)]->find(splitter->Flav());
    if (kit==p_kernels[int(type)]->end()) {
      // didn't find any meaningful splitting for the type (FF, FI, IF, or II):
      // add a zero weight.
      vector<double> help; help.resize(1); help.push_back(0.);
      m_integrals[i] = help;
    }
    else {
      // found a meaningful splitting for the type (FF, FI, IF, or II):
      // calculate the integrals and store the vector of results for
      // later selextion of winning kernel
      Kernels * kernels = kit->second;
      if (kernels) sum += kernels->CalcIntegrals(split,p_msel);
      m_integrals[i]    = kernels->GetIntegrals();
      m_splitkernels[i] = kernels; 
    }
    i++;
  }
  return sum;
}
  
Splitting * Shower::GenerateTestSplitting(Parton * splitter,const double & sum) {
  const Parton_List * spectators = splitter->GetSpectators();
  double t = p_config->T(), tstart = t, t0 = p_config->T0();
  while (t>t0) {
    t *= exp(log(ran->Get())/sum);
    if (t<t0) return NULL;
    double disc = sum * ran->Get(), specsum;
    size_t i = 0;
    for (Parton_List::const_iterator spit=spectators->begin();
	 spit!=spectators->end();spit++) {
      // Select the list of applicable kernels depending on the kinematic configuration
      // and the flavour of the splitter.  Initialise a container "Splitting" for the
      // information defining the potential splitting.
      bool active = true;
      disc -= specsum = m_integrals[i].back();
      if (disc<=0. && active) {
	disc = specsum * ran->Get();
	for (size_t j=0;j<m_integrals[i].size()-1;j++) {
	  disc -= m_integrals[i][j];
	  if (disc<=0. && active) {
	    Parton * spectator = (*spit);
	    Kernel * kernel    = (*m_splitkernels[i])[j];
	    Splitting * split  = new Splitting(splitter,spectator,t,t0);
	    split->Set_tstart(tstart);
	    split->SetKinScheme(m_kinscheme);
	    // Veto and adding of weights is embedded here.
	    // Generate z, phi, run the veto algorithm, and construct the kinematics.
	    // Veto (accept/reject) with exact operator divided by over estimator.
	    // If everything works out, the splitting is allowed and we keep it.
	    // Rejected splittings add to the overall rejection weight related to
	    // the splitter parton.
	    if (kernel->Generate(*split,*p_config,p_msel,m_weightover)) return split;
	    delete split;
	    active = false;
	    break;
	  }
	}
      }
      i++;
      if (!active) break;
    }
  }
  return NULL;
}

bool Shower::PerformSplitting() {
  if (p_winner->IsEndPoint()) {
    if (!p_winner->GetKernel()->UpdateSystem(*p_winner,*p_config)) {
      msg_Error()<<METHOD<<" failed to update to endpoint kinematics for:\n"<<(*p_winner)
		 <<"   Return false and hope for the best.\n";
      return false;
    }
  }
  if (!p_winner->GetKernel()->FillOffsprings((*p_winner))) {
    msg_Error()<<METHOD<<" failed to set colours for splitting:\n"<<(*p_winner)
	       <<"   Return false and hope for the best.\n";
    return false;
  }
  EstablishSpectators(p_winner->GetSplitter());
  //EstablishSoftPartners();
  return p_winner->GetKernel()->UpdateSystem(*p_winner,*p_config);
}

void Shower::EstablishSoftPartners() {
  if (p_winner->NPartons()>2) return;
  Parton * offspring[2], * spectator = p_winner->GetSpectator();
  bool establishsoftpartners = false;
  for (size_t i=0;i<2;i++) {
    offspring[i] = p_winner->GetParton(i);
    if (offspring[i]->Flav().IsGluon()) establishsoftpartners = true;
  }
  if (!establishsoftpartners) return;
  if (offspring[0]->Flav().IsGluon() && !offspring[0]->Flav().IsGluon())
    swap(offspring[0],offspring[1]);
  offspring[0]->SetSoftPartner(0,offspring[1]); 
  offspring[0]->SetSoftPartner(1,NULL);         
  offspring[0]->SetSoftPartner(2,spectator);    
  offspring[1]->SetSoftPartner(0,offspring[0]); 
  offspring[1]->SetSoftPartner(1,spectator);    
  offspring[1]->SetSoftPartner(2,NULL);         
  spectator->SetSoftPartner(0,offspring[1]);
  spectator->SetSoftPartner(1,NULL);
  spectator->SetSoftPartner(2,offspring[0]);
}
  
void Shower::EstablishSpectators(Parton * splitter) {
  const Parton_List * spectators = splitter->GetSpectators();
  vector<unsigned int> cols;  cols.resize(2);
  for (size_t i=0;i<3;i++) {
    Parton * offspring = p_winner->GetParton(i);
    if (offspring==NULL) continue;
    for (size_t col=0;col<2;col++) { cols[col] = offspring->GetColor()[col]; }
    for (size_t j=i+1;j<3;j++) {
      Parton * compare = p_winner->GetParton(j);
      if (compare==NULL) continue;
      for (size_t col=0;col<2;col++) {
	if (cols[col]!=0 && cols[col]==compare->GetColor()[1-col]) {
	  if (!offspring->FindSpectator(compare)) {
	    offspring->AddSpectator(compare);
	    compare->AddSpectator(offspring);
	  }
	}
      }
    }
    for (Parton_List::const_iterator spec=spectators->begin();
	 spec!=spectators->end();spec++) {
      Parton * spectator = (*spec);
      for (size_t col=0;col<2;col++) {
	if (cols[col]!=0 && cols[col]==spectator->GetColor()[1-col]) {
	  if (!offspring->FindSpectator(spectator)) {
	    offspring->AddSpectator(spectator);
	  }
	  if (!spectator->FindSpectator(offspring)) spectator->AddSpectator(offspring);
	}
      }
    }
    p_config->push_back(offspring);
  }
  for (Parton_List::const_iterator spec=spectators->begin();
       spec!=spectators->end();spec++) {
    (*spec)->RemoveSpectator(splitter);
  }
  splitter->SwitchOff();
}

bool Shower::Init(MODEL::Model_Base * const model,
		  PDF::ISR_Handler * const isr)
{
  // Shower parameters and switches are pulled from the map in the CFP_Parameter
  // - parton shower cutoffs for FS and IS showering
  // - order of the splitting function and details of terms included
  // - scale setting schemes
  // - eventually recoil schemes
  m_t0[0]       = (*cfp_pars)("pt2min(FS)");
  m_t0[1]       = (*cfp_pars)("pt2min(IS)");
  m_t0min       = Min(m_t0[0], m_t0[1]);
  m_kinscheme   = (*cfp_pars)["kinematics"];
  m_nmax_em     = (*cfp_pars)["max_emissions"];
  m_weightover  = 3.;

  Kernel_Constructor kernelconstructor(this);
  kernelconstructor.Init(model,isr);
  for (size_t i=0;i<5;i++) p_kernels[i] = kernelconstructor(i);
  //msg_Out()<<"Will exit in "<<METHOD<<"\n";
  //exit(1);
  return true;
}

