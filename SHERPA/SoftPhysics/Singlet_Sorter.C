#include "SHERPA/SoftPhysics/Singlet_Sorter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Phys/Particle.H"

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Singlet_Sorter::Singlet_Sorter() {}
Singlet_Sorter::~Singlet_Sorter() {
  ResetPartLists();
}

void Singlet_Sorter::ResetPartLists() {
  while (!m_partlists.empty()) {
    m_partlists.back()->clear();
    delete m_partlists.back();
    m_partlists.pop_back();
  }  
  m_partlists.clear();
  m_hadrons.clear();
}

Return_Value::code Singlet_Sorter::operator()(Blob_List * bloblist) {
  // Logic:
  // 1. Harvest particles in HarvestParticles
  //    fill particles into lists, one for each hadron or tau decay (to keep track of
  //    decay vertices) and one for everything else.  These lists are the accumulated in
  //    plists, with the default list = plists[0] for all coloured particles not
  //    coming from hadron/tau decays.
  // 2. For the hadrons there is a separate list, hadrons.  Its particles will be filled
  //    into a separate blob, in DealWithHadrons.
  //    TODO: I have to check for partonic hadron decays if I get the connections right.
  // 3. The plists will be decomposed into singlets and filled into one blob

  ResetPartLists();
  if (!HarvestParticles(bloblist)) return Return_Value::New_Event;
  if (m_partlists.size()==0 || (*m_partlists.begin())->empty()) return Return_Value::Nothing;
  while (!m_partlists.empty()) {
    p_partlist = m_partlists.front();  
    if (DecomposeIntoSinglets()) {
      Blob * blob = MakeBlob();
      if (blob) bloblist->push_back(blob);
    }
    else {
      msg_Error()<<"Error in "<<METHOD<<" failed to decompose particle list into singlet.\n"
		 <<"   Reset list, return Error and hope for the best.\n";
      ResetPartLists();
      return Return_Value::Error;
    }
    delete p_partlist; /// does this help?
    m_partlists.pop_front();
  }
  return Return_Value::Success;
}

bool Singlet_Sorter::HarvestParticles(Blob_List * bloblist) {
  for (Blob_List::iterator blit=bloblist->begin();
       blit!=bloblist->end();++blit) {
    if ((*blit)->Has(blob_status::needs_reconnections) ||
	(*blit)->Has(blob_status::needs_hadronization)) {
      Blob* upstream_blob=(*blit)->UpstreamBlob();
      if (upstream_blob && upstream_blob->Type()==btp::Hadron_Decay) {
	p_partlist = new Part_List;
        m_partlists.push_back(p_partlist);
      }
      else {
	m_partlists.push_back(new Part_List);
	p_partlist = m_partlists.front();
      }
      if (!FillParticleLists(*blit)) return false;
      (*blit)->UnsetStatus(blob_status::needs_reconnections |
			   blob_status::needs_hadronization);
    }
  }
  if (m_partlists.size()==1 && m_partlists.front()->empty()) m_partlists.pop_front();
  DealWithHadrons(bloblist);
  return true;
}

bool Singlet_Sorter::FillParticleLists(Blob * blob) {
  for (int i=0;i<blob->NOutP();i++) {
    Particle * part = blob->OutParticle(i); 
    if (part->Status()==part_status::active && 
	part->Info()!='G' && part->Info()!='I') {
      if (part->GetFlow(1)!=0 || part->GetFlow(2)!=0) {
	if (part->GetFlow(1)==part->GetFlow(2)) {
	  msg_Error()<<"Error in "<<METHOD<<": blob with funny colour assignements:\n"
		     <<"   "<<(*part)<<"\n"
		     <<"   Will demand new event and hope for the best.\n";
	  return false;
	}
	p_partlist->push_back(part);
	part->SetStatus(part_status::fragmented);
      }
      else if (part->Flav().Kfcode()==kf_tau || part->Flav().IsHadron())
	m_hadrons.push_back(part);
    }
  }
  return true;
}

void Singlet_Sorter::DealWithHadrons(Blob_List * bloblist) {
  if (m_hadrons.size()>0) {
    Blob * blob = new Blob();
    blob->SetId();
    blob->SetType(btp::Fragmentation);
    blob->SetStatus(blob_status::needs_hadrondecays);
    while (!m_hadrons.empty()) {
      Particle * part = m_hadrons.back();
      blob->AddToInParticles(part);
      part->SetStatus(part_status::decayed);
      blob->AddToOutParticles(new Particle((*part)));
      blob->GetOutParticles().back()->SetStatus(part_status::active);
      m_hadrons.pop_back();
    }
    bloblist->push_back(blob);
  }
}

bool Singlet_Sorter::DecomposeIntoSinglets() {
  Part_List sorted;
  while (!p_partlist->empty()) {
    if (!NextSinglet(sorted,true) &&
	!NextSinglet(sorted,false)) {
      msg_Error()<<"Error in "<<METHOD<<" particles left in list.\n";
      for (Part_List::iterator pit=p_partlist->begin();pit!=p_partlist->end();pit++)
	msg_Error()<<"  "<<(**pit)<<"\n";
      return false;
    }
  }
  p_partlist->splice(p_partlist->begin(),sorted);
  return true;
}

bool Singlet_Sorter::NextSinglet(Part_List & sorted,const bool triplet) {
  Particle  * part = NULL;
  for (Part_Iterator pit=p_partlist->begin();pit!=p_partlist->end();++pit) {
    if ((*pit)->GetFlow(1)!=0 && ((triplet && (*pit)->GetFlow(2)==0) ||
				  (!triplet && (*pit)->GetFlow(2)!=0))) {
      part = (*pit);
      p_partlist->erase(pit);
      break;
    }
  }
  if (part) {
    sorted.push_back(part);
    size_t col = triplet?0:part->GetFlow(1);
    do {
      part = FindNext(part->GetFlow(1));
      if (part) sorted.push_back(part);
    } while (part && part->GetFlow(1)!=col);
    return true;
  }
  return false;
}

Particle * Singlet_Sorter::FindNext(const size_t col) {
  for (Part_Iterator pit=p_partlist->begin();pit!=p_partlist->end();++pit) {
    if ((*pit)->GetFlow(2)==col) {
      Particle * part = (*pit);
      p_partlist->erase(pit);
      return part;
    }
  }
  return NULL;
}

Blob * Singlet_Sorter::MakeBlob() {
  if (p_partlist->empty()) return NULL;
  Blob * blob = new Blob();
  blob->SetId();
  blob->SetType(btp::Fragmentation);
  blob->SetStatus(blob_status::needs_hadronization);
  Particle * part = p_partlist->front();
  Blob     * prod = part->ProductionBlob(), * up(prod?prod->UpstreamBlob():NULL);
  if (!up || (up && up->Type()!=btp::Hadron_Decay)) {
    blob->AddStatus(blob_status::needs_reconnections);
  }
  while (!p_partlist->empty()) {
    blob->AddToInParticles(p_partlist->front());
    p_partlist->pop_front();
  }
  return blob;
}
