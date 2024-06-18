#include "SHERPA/Single_Events/Jet_Evolution.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Phys/NLO_Types.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "REMNANTS/Main/Remnant_Handler.H"
#include "SHERPA/PerturbativePhysics/Perturbative_Interface.H"
#include <map>

using namespace SHERPA;
using namespace ATOOLS;
using namespace PDF;
using namespace std;

Jet_Evolution::Jet_Evolution(Matrix_Element_Handler * me,
                             Hard_Decay_Handler * harddecs,
                             Decay_Handler_Base * decs,
                             const MI_Handler_Map *mihandlers,
                             const Soft_Collision_Handler_Map *schandlers,
                             const Shower_Handler_Map &showers,
                             REMNANTS::Remnant_Handler_Map &remnanthandlers) {
  Shower_Handler_Map::const_iterator shIter = showers.find(isr::hard_process);
  m_name = string("Jet_Evolution:") + shIter->second->ShowerGenerator();
  m_type = eph::Perturbative;
  FillPerturbativeInterfaces(me,harddecs,decs,mihandlers,schandlers,showers,remnanthandlers);
}

void Jet_Evolution::FillPerturbativeInterfaces(Matrix_Element_Handler * me,
					       Hard_Decay_Handler* harddecs,
					       Decay_Handler_Base * decs,
					       const MI_Handler_Map * mis,
					       const Soft_Collision_Handler_Map * scs,
					       const Shower_Handler_Map & showers,
					       REMNANTS::Remnant_Handler_Map & rhs) {
  REMNANTS::Remnant_Handler * remnants = NULL; 
  if (rhs.find(isr::hard_process)!=rhs.end()) remnants = rhs[isr::hard_process];
  else msg_Error()<<"Error in "<<METHOD<<":\n"
		  <<"  No remnant handling found for hard part of the process.\n"
		  <<"  Continue and hope for the best.\n";
  
  Shower_Handler_Map::const_iterator shower = showers.find(isr::hard_process);
  if (shower!=showers.end() && me) {
    m_pertinterfaces["SignalMEs"] = new Perturbative_Interface(me, harddecs, shower->second);
    m_pertinterfaces["SignalMEs"]->SetRemnantHandler(remnants);
  }

  shower = showers.find(isr::hard_subprocess);
  if (shower!=showers.end()) {
    m_pertinterfaces["HadronDecays"] = new Perturbative_Interface(decs, shower->second);
    MI_Handler_Map::const_iterator mihandler = mis->find(isr::hard_subprocess);
    if (mihandler!=mis->end()) {
      m_pertinterfaces["MPIs"] = new Perturbative_Interface(mihandler->second, shower->second);
      m_pertinterfaces["MPIs"]->SetRemnantHandler(remnants);
    }
    Soft_Collision_Handler_Map::const_iterator schandler = scs->find(isr::hard_subprocess);
    if (schandler!=scs->end()) {
      m_pertinterfaces["SoftCollisions"] = new Perturbative_Interface(schandler->second, shower->second);
      m_pertinterfaces["SoftCollisions"]->SetRemnantHandler(remnants);
    }
  }

  shower = showers.find(isr::bunch_rescatter);
  if (shower!=showers.end()) {
    if (rhs.find(isr::bunch_rescatter)!=rhs.end()) remnants = rhs[isr::bunch_rescatter];
    else msg_Error()<<"Error in "<<METHOD<<":\n"
		    <<"  No remnant handling found for bunch rescattering.\n"
		    <<"  Continue and hope for the best.\n";
    MI_Handler_Map::const_iterator mihandler = mis->find(isr::bunch_rescatter);
    if (mihandler!=mis->end()) {
      m_pertinterfaces["BR_MPIs"] = new Perturbative_Interface(mihandler->second, shower->second);
      m_pertinterfaces["BR_MPIs"]->SetRemnantHandler(remnants);
    }
    Soft_Collision_Handler_Map::const_iterator schandler = scs->find(isr::bunch_rescatter);
    if (schandler!=scs->end()) {
      m_pertinterfaces["BR_SoftCollisions"] = new Perturbative_Interface(schandler->second, shower->second);
      m_pertinterfaces["BR_SoftCollisions"]->SetRemnantHandler(remnants);
    }
  }
}

Jet_Evolution::~Jet_Evolution() {
  while (m_pertinterfaces.size() > 0) {
    if (m_pertinterfaces.begin()->second != NULL)
      delete m_pertinterfaces.begin()->second;
    m_pertinterfaces.erase(m_pertinterfaces.begin());
  }
}

Return_Value::code Jet_Evolution::Treat(Blob_List *bloblist) {
  if (bloblist->empty()) {
    msg_Error() << "Potential error in Jet_Evolution::Treat." << endl
                << "   Incoming blob list contains " << bloblist->size()
                << " entries." << endl
                << "   Continue and hope for the best." << endl;
    return Return_Value::Error;
  }
  PertInterfaceIter piIter;
  bool hit(false), found(true);
  while (found) {
    found = false;
    for (size_t i = 0; i < bloblist->size(); ++i) {
      Blob *meblob = (*bloblist)[i];
      if (meblob->Has(blob_status::needs_showers) &&	  
          meblob->Type() != btp::Hard_Decay) {
	piIter = SelectInterface(meblob);
        switch (AttachShowers(meblob, bloblist, piIter->second)) {
        case Return_Value::Success:
          found = hit = true;
          break;
        case Return_Value::New_Event:
          return Return_Value::New_Event;
        case Return_Value::Retry_Event:
          return Return_Value::Retry_Event;
        case Return_Value::Nothing:
          return Return_Value::Nothing;
        case Return_Value::Error:
          return Return_Value::Error;
        default:
          msg_Error() << "ERROR in " << METHOD
                      << ": Unexpected status of AttachShowers for\n"
                      << (*meblob)
                      << "   Return 'Error' and hope for the best.\n";
          return Return_Value::Error;
        }
      }
    }
    if (found)
      hit = true;
    Reset();
  }
  if (hit) {
    // enable shower generator independent FS QED correction to ME
    // TODO: check first, whether shower did FS QED
    if (!bloblist->FourMomentumConservation()) {
      msg_Tracking() << METHOD << " found four momentum conservation error.\n";
      return Return_Value::New_Event;
    }
    Blob * showerblob = bloblist->FindLast(btp::Shower);
    showerblob->AddStatus(blob_status::needs_extraQED);
    return Return_Value::Success;
  }
  // Capture potential problem with empty remnants here.
  // This should only happen after retrying an event has been called.  In this
  // case we find the last (and hopefully only) shower blob and extract its
  // initiators.
  Blob *showerblob = bloblist->FindLast(btp::Shower);
  if (showerblob!=NULL && showerblob->Has(blob_status::needs_beams)) {
    Blob * meblob = showerblob->InParticle(0)->ProductionBlob();
    REMNANTS::Remnant_Handler * remnants =
      SelectInterface(meblob)->second->RemnantHandler();
    if (meblob->Type()!=btp::Hadron_Decay &&
	!remnants->ExtractShowerInitiators(showerblob))
      return Return_Value::New_Event;
  }
  return Return_Value::Nothing;
}


PertInterfaceIter Jet_Evolution::SelectInterface(Blob * blob) {
  string tag("");
  switch (int(blob->Type())) {
  case (int(btp::Signal_Process)):
  case (int(btp::Hard_Decay)):
    tag = string("SignalMEs");
    MODEL::as->SetActiveAs(PDF::isr::hard_process);
    break;
  case (int(btp::Hard_Collision)):
    if (blob->Has(blob_status::needs_beamRescatter)) tag = string("BR_");
    tag += string("MPIs");
    if (blob->TypeSpec() == "MinBias" || blob->TypeSpec()=="Shrimps")
      tag += string("SoftCollisions");
    MODEL::as->SetActiveAs(PDF::isr::hard_subprocess);
    break;
  case (int(btp::Hadron_Decay)):
    tag = string("HadronDecays");
    MODEL::as->SetActiveAs(PDF::isr::hard_subprocess);
    break;
  default:
    msg_Error() << "ERROR in " << METHOD << ": "
		<< "Do not have an interface for this type of blob.\n"
		<< (*blob) << "\n   Will abort.\n";
    THROW(fatal_error, "No perturbative interface found.");
  }
  PertInterfaceIter piIter = m_pertinterfaces.find(tag);
  if (piIter == m_pertinterfaces.end()) {
    msg_Error() << "Error in Jet_Evolution::Treat: "
		<< "No Perturbative_Interface found for type " << tag
		<< "\n"
		<< "   Abort the run.\n";
    THROW(fatal_error, "No perturbative interface found.");
  }
  return piIter;
}

Return_Value::code
Jet_Evolution::AttachShowers(Blob *blob, Blob_List *bloblist,
                             Perturbative_Interface *pertinterface) {
  p_remnants = pertinterface->RemnantHandler();
  if (!pertinterface->Shower()->On() ||
      (pertinterface->MEHandler() &&
       pertinterface->MEHandler()->Process()->Info().m_nlomode ==
           nlo_mode::fixedorder)) {
    AftermathOfNoShower(blob, bloblist);
    Blob * noshowerblob = bloblist->FindLast(btp::Shower);
    noshowerblob->AddStatus(blob_status::needs_extraQED);
    return Return_Value::Nothing;
  }
  int shower(0);
  Return_Value::code stat(pertinterface->DefineInitialConditions(blob, bloblist));
  if (stat == Return_Value::New_Event || stat == Return_Value::Retry_Event) {
    pertinterface->CleanUp();
    return stat;
  }
  if (blob->Type() != ::btp::Hadron_Decay) {
    msg_Debugging() << METHOD << "(): Setting scale for MI {\n";
    double scale(0.0);
    Cluster_Amplitude *ampl(pertinterface->Amplitude());
    while (ampl->Next())
      ampl = ampl->Next();
    msg_Debugging() << *ampl << "\n";
    scale = sqrt(ampl->MuQ2());
    blob->AddData("MI_Scale", new Blob_Data<double>(scale));
    msg_Debugging() << "} -> p_T = " << scale << "\n";
  }
  switch (stat) {
  case Return_Value::Success:
    if (blob->Type() != ::btp::Hadron_Decay)
      DefineInitialConditions(blob, bloblist, pertinterface);
    if (blob->NInP() == 1) shower = pertinterface->PerformDecayShowers();
    else if (blob->NInP() == 2) {
      shower = pertinterface->PerformShowers();
      blob->UnsetStatus(blob_status::needs_beamRescatter);
    }
    switch (shower) {
    case 1:
      // No Sudakov rejection
      Reset();
      if (AftermathOfSuccessfulShower(blob, bloblist, pertinterface)) {
        pertinterface->CleanUp();
        return Return_Value::Success;
      }
      blob->SetStatus(blob_status::inactive);
      CleanUp();
      return Return_Value::New_Event;
    case 0:
      // Sudakov rejection
      Reset();
      CleanUp();
      return Return_Value::New_Event;
    default:
      THROW(fatal_error, "Invalid return value from shower");
    }
  case Return_Value::Nothing:
    if (AftermathOfNoShower(blob, bloblist)) {
      pertinterface->CleanUp();
      return Return_Value::Success;
    }
    blob->SetStatus(blob_status::inactive);
    CleanUp();
    return Return_Value::New_Event;
  case Return_Value::Error:
    msg_Error() << "ERROR in " << METHOD << ":" << std::endl
                << "   DefineInitialConditions yields an error for "
                << std::endl
                << (*blob) << "   Return 'Error' and hope for the best."
                << std::endl;
    blob->SetStatus(blob_status::inactive);
    CleanUp();
    return Return_Value::Error;
  default:
    msg_Error() << "ERROR in " << METHOD << ":" << std::endl
                << "   Unexpected status of DefineInitialConditions for "
                << std::endl
                << (*blob) << "   Return 'Error' and hope for the best."
                << std::endl;
    blob->SetStatus(blob_status::inactive);
    CleanUp();
    return Return_Value::Error;
  }
  return Return_Value::Error;
}

bool Jet_Evolution::AftermathOfNoShower(Blob *blob, Blob_List *bloblist) {
  Blob *noshowerblob = new Blob();
  noshowerblob->SetType(btp::Shower);
  for (size_t i = 0; i < blob->GetInParticles().size(); ++i) {
    noshowerblob->AddToOutParticles(blob->InParticle(i));
    noshowerblob->AddToInParticles(new Particle(*blob->InParticle(i)));
    noshowerblob->InParticle(i)->SetBeam(i);
    blob->InParticle(i)->SetStatus(part_status::decayed);
  }
  for (size_t i = 0; i < blob->GetOutParticles().size(); ++i) {
    if (blob->OutParticle(i)->DecayBlob())
      continue;
    noshowerblob->AddToInParticles(blob->OutParticle(i));
    noshowerblob->AddToOutParticles(new Particle(*blob->OutParticle(i)));
    blob->OutParticle(i)->SetStatus(part_status::decayed);
  }
  noshowerblob->SetStatus(blob_status::needs_beams |
                          blob_status::needs_hadronization);
  if (blob->Type() != btp::Hadron_Decay) {
    noshowerblob->AddStatus(blob_status::needs_reconnections);
  }
  noshowerblob->SetId();
  noshowerblob->SetTypeSpec("No_Shower");
  bloblist->push_back(noshowerblob);
  blob->SetStatus(blob_status::inactive);
  return p_remnants->ExtractShowerInitiators(noshowerblob);
}

bool Jet_Evolution::AftermathOfSuccessfulShower(Blob *blob, Blob_List *bloblist,
						Perturbative_Interface *pertinterface) {
  if (blob->NInP() == 1 && blob->Type() != btp::Hadron_Decay)
    blob->InParticle(0)->SetInfo('h');
  pertinterface->FillBlobs();
  blob->UnsetStatus(blob_status::needs_showers);
  Blob *showerblob =
      (!pertinterface->Shower()->On() ? CreateMockShowerBlobs(blob, bloblist)
                                  : bloblist->FindLast(btp::Shower));
  if (showerblob==NULL || blob->Type()== btp::Hadron_Decay) return true;
  showerblob->AddStatus(blob_status::needs_reconnections);
  return p_remnants->ExtractShowerInitiators(showerblob);
}

ATOOLS::Blob *Jet_Evolution::CreateMockShowerBlobs(Blob *const meblob,
                                                   Blob_List *const bloblist) {
  Blob *ISRblob = NULL;
  if (meblob->NInP() != 1) {
    for (int i = 0; i < 2; i++) {
      // new ISR Blob
      ISRblob = new Blob();
      ISRblob->SetType(btp::Shower);
      ISRblob->SetStatus(blob_status::needs_beams);
      Particle *part = new Particle(*meblob->InParticle(i));
      part->SetStatus(part_status::decayed);
      part->SetBeam(int(meblob->InParticle(1 - i)->Momentum()[3] >
                        meblob->InParticle(i)->Momentum()[3]));
      ISRblob->AddToInParticles(part);
      ISRblob->AddToOutParticles(meblob->InParticle(i));
      meblob->InParticle(i)->SetStatus(part_status::decayed);
      ISRblob->SetId();
      bloblist->insert(bloblist->begin(), ISRblob);
    }
  }
  for (int i = 0; i < meblob->NOutP(); i++) {
    Blob *FSRblob = new Blob();
    FSRblob->SetType(btp::Shower);
    FSRblob->SetStatus(blob_status::needs_hadronization);
    if (meblob->Type() != btp::Hadron_Decay) {
      FSRblob->AddStatus(blob_status::needs_reconnections);
    }
    Particle *part = new Particle(*meblob->OutParticle(i));
    if (meblob->OutParticle(i)->DecayBlob()) {
      Blob *dec = meblob->OutParticle(i)->DecayBlob();
      if (dec->Type() == btp::Hard_Decay) {
        dec->RemoveInParticle(meblob->OutParticle(i));
        dec->AddToInParticles(part);
      }
    }
    FSRblob->AddToInParticles(meblob->OutParticle(i));
    meblob->OutParticle(i)->SetStatus(part_status::decayed);
    FSRblob->AddToOutParticles(part);
    FSRblob->SetId();
    bloblist->push_back(FSRblob);
  }
  return ISRblob;
}

void Jet_Evolution::CleanUp(const size_t &mode) {
  for (PertInterfaceIter piIter = m_pertinterfaces.begin();
       piIter != m_pertinterfaces.end(); ++piIter) {
    piIter->second->CleanUp();
  }
}

void Jet_Evolution::Reset() {
  for (PertInterfaceIter piIter = m_pertinterfaces.begin();
       piIter != m_pertinterfaces.end(); ++piIter) {
    piIter->second->Shower()->GetISRHandler()->Reset(0);
    piIter->second->Shower()->GetISRHandler()->Reset(1);
  }
}

bool Jet_Evolution::DefineInitialConditions(const Blob *blob,
                                            const Blob_List *bloblist,
                                            Perturbative_Interface *pertinterface) {
  Reset();
  msg_Debugging() << METHOD << "(): {\n";
  for (::Blob_List::const_iterator blit = bloblist->begin();
       blit != bloblist->end(); ++blit) {
    if ((*blit)->Type() == ::btp::Shower) {
      // Update(*blit,0, pertinterface);
      // Update(*blit,1, pertinterface);
    }
  }
  msg_Debugging() << "}\n";
  return true;
}

void Jet_Evolution::Update(Blob *blob, const size_t beam,
                           Perturbative_Interface *pertinterface) {
  size_t cbeam = 0;
  for (int i = 0; i < blob->NInP(); ++i) {
    Particle *cur = blob->InParticle(i);
    if (!cur->Flav().Strong() || cur->ProductionBlob())
      continue;
    if (cbeam == beam) {
      msg_Debugging() << "  " << *cur << ", beam = " << beam << "\n";
      if (!p_remnants->Extract(cur, beam))
        msg_Error() << METHOD << ": Cannot extract particle:\n"
                    << (*cur) << "\n  from:\n"
                    << p_remnants[beam].GetRemnant(beam)->GetBeam()->Bunch()
                    << "\n";
      return;
    }
    ++cbeam;
  }
}

void Jet_Evolution::Finish(const string &) {}
