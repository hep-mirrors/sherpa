#include "SHERPA/Single_Events/Hadronization.H"
#include "ATOOLS/Org/Message.H"

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Hadronization::Hadronization(Colour_Reconnection_Handler * reconnections,
			     Fragmentation_Handler * fragmentation) :
  p_reconnectionhandler(reconnections), p_fragmentationhandler(fragmentation)
{
  m_name = std::string("Hadronization:")+p_fragmentationhandler->FragmentationModel();
  m_type = eph::Hadronization;
}

Hadronization::~Hadronization() {}

Return_Value::code Hadronization::Treat(ATOOLS::Blob_List *bloblist,double &weight) 
{
  if (bloblist->empty()) {
    msg_Error()<<"Hadronization::Treat("<<bloblist<<","<<weight<<"): "<<endl
	       <<"   Blob list contains "<<bloblist->size()<<" entries."<<endl
	       <<"   Continue and hope for the best."<<endl;
    return Return_Value::Error;
  }
  switch (int(ExtractSinglets(bloblist))) {
  case int(Return_Value::Success) : break;
  case int(Return_Value::Nothing) : return Return_Value::Nothing;
  case int(Return_Value::Error)   : return Return_Value::Error;
  default :
    msg_Error()<<"ERROR in "<<METHOD<<":"<<std::endl
	       <<"   ExtractSinglets yields unknown return value."<<std::endl
	       <<"   Return 'Retry_Event' and hope for the best."<<std::endl;
    return Return_Value::Retry_Event;
  }
  (*p_reconnectionhandler)(bloblist);
  return p_fragmentationhandler->PerformFragmentation(bloblist);
}


void Hadronization::CleanUp(const size_t & mode) {
  p_reconnectionhandler->CleanUp(mode);
}

void Hadronization::Finish(const std::string &) {}


Return_Value::code Hadronization::ExtractSinglets(Blob_List * bloblist)
{
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
  std::vector<SP(Part_List)> plists;
  plists.push_back(new Part_List);
  if (!HarvestParticles(bloblist,&plists)) return Return_Value::New_Event;
  if (plists[0]->empty() && plists.size()<2) {
    msg_Debugging()<<"WARNING in Lund_Interface::PrepareFragmentationBlob:\n"
		   <<"   No coloured particle found leaving shower blobs.\n";
    return Return_Value::Nothing;
  }
  Return_Value::code ret(Return_Value::Success);
  for (size_t i=0; i<plists.size(); ++i) {
    if (plists[i]->empty()) continue;
    SP(Part_List) plist=plists[i];
    vector<SP(Part_List)> partlists; 
    if (DecomposeIntoSinglets(plist,partlists)) {
      bloblist->push_back(MakeBlob(partlists));
      ret=Return_Value::Success;
    }
    else {
      ret=Return_Value::Error;
      break;
    }
  }
  return ret;
}

bool Hadronization::HarvestParticles(Blob_List * bloblist, 
				     vector<SP(Part_List)> * const plists) {
  vector<Particle*> hadrons;
  for (Blob_List::iterator blit=bloblist->begin();
       blit!=bloblist->end();++blit) {
    if ((*blit)->Has(blob_status::needs_hadronization)) {
      SP(Part_List) plist((*plists)[0]);
      Blob* upstream_blob=(*blit)->UpstreamBlob();
      if (upstream_blob && upstream_blob->Type()==btp::Hadron_Decay) {
        plist=new Part_List;
        plists->push_back(plist);
      }
      if (!FillParticleList((*blit),plist,hadrons)) {
	msg_Out()<<(*bloblist)<<"\n";
	exit(1);
	return false;
      }
      (*blit)->UnsetStatus(blob_status::needs_reconnections |
			   blob_status::needs_hadronization);
    }
  }
  DealWithHadrons(bloblist,hadrons);
  return true;
}

bool Hadronization::FillParticleList(Blob * blob, SP(Part_List) & plist,
				     vector<Particle*> hadrons) {
  for (int i=0;i<blob->NOutP();i++) {
    Particle * part = blob->OutParticle(i); 
    if (part->Status()==part_status::active && 
	part->Info()!='G' && part->Info()!='I') {
      if (part->GetFlow(1)!=0 || part->GetFlow(2)!=0) {
	if (part->GetFlow(1)==part->GetFlow(2)) {
	  msg_Error()<<"Error in "<<METHOD<<":\n"
		     <<"   Blob with funny colour assignements.\n"
		     <<"   Will demand new event and hope for the best.\n";
	  return false;
	}
	plist->push_back(part);
	part->SetStatus(part_status::fragmented);
      }
      else if (part->Flav().Kfcode()==kf_tau || part->Flav().IsHadron()) {
	hadrons.push_back(part);
      }
    }
  }
}

void Hadronization::DealWithHadrons(Blob_List * bloblist, 
				    vector<Particle*> hadrons) {
  if (hadrons.size()>0) {
    Blob * blob = new Blob();
    blob->SetId();
    blob->SetType(btp::Fragmentation);
    blob->SetStatus(blob_status::needs_hadrondecays);
    for (size_t i=0;i<hadrons.size();++i) {
      blob->AddToInParticles(hadrons[i]);
      hadrons[i]->SetStatus(part_status::decayed);
      blob->AddToOutParticles(new Particle((*hadrons[i])));
      blob->GetOutParticles().back()->SetStatus(part_status::active);
    }
    bloblist->push_back(blob);
  }
}

bool Hadronization::
DecomposeIntoSinglets(SP(Part_List) plist,vector<SP(Part_List)> & partlists) {
  Part_List * pli(NULL);
  //for (Part_List::iterator pit=plist->begin();pit!=plist->end();pit++)
  //  msg_Out()<<"  "<<(**pit)<<"\n";
  while (!plist->empty()) {
    pli = NextSinglet(plist,true);
    if (!pli) pli = NextSinglet(plist,false);
    if (pli) partlists.push_back(pli);
    else {
      if (!plist->empty()) {
	msg_Error()<<"Error in "<<METHOD<<" particles left in list.\n";
	for (Part_List::iterator pit=plist->begin();pit!=plist->end();pit++)
	  msg_Error()<<"  "<<(**pit)<<"\n";
	exit(1);
      }
    }
  }
  return true;
}

Part_List * Hadronization::NextSinglet(SP(Part_List) plist,bool triplet) {
  Part_List * pli = NULL;
  Particle * part = NULL;
  for (Part_Iterator pit=plist->begin();pit!=plist->end();++pit) {
    if ((*pit)->GetFlow(1)!=0 && ((triplet && (*pit)->GetFlow(2)==0) ||
				  (!triplet && (*pit)->GetFlow(2)!=0))) {
      part = (*pit);
      //msg_Out()<<"Start singlet with "<<part->Number()<<"("<<triplet<<").\n";
      plist->erase(pit);
      break;
    }
  }
  if (part) {
    pli  = new Part_List;  
    pli->push_back(part);
    size_t col = triplet?0:part->GetFlow(1);
    do {
      part = FindNext(plist,part->GetFlow(1));
      if (part) {
	pli->push_back(part);
	//msg_Out()<<"               add "<<part->Number()<<"("<<triplet<<").\n";
      }
    } while (part && part->GetFlow(1)!=col);
  }
  return pli;
}

Particle * Hadronization::FindNext(SP(Part_List) plist,const size_t col) {
  for (Part_Iterator pit=plist->begin();pit!=plist->end();++pit) {
    if ((*pit)->GetFlow(2)==col) {
      Particle * part = (*pit);
      plist->erase(pit);
      return part;
    }
  }
  return NULL;
}

Blob * Hadronization::MakeBlob(vector<SP(Part_List)> & partlists) {
  Blob * blob = new Blob();
  blob->SetId();
  blob->SetType(btp::Fragmentation);
  Particle * part = (*partlists.begin())->front();
  Blob     * prod = part->ProductionBlob(), * up(prod->UpstreamBlob());
  blob->SetStatus(blob_status::needs_hadronization);
  if (!up || (up && up->Type()!=btp::Hadron_Decay)) {
    blob->AddStatus(blob_status::needs_reconnections);
  }
  for (vector<SP(Part_List)>::iterator pliter=partlists.begin();
       pliter!=partlists.end();pliter++) {
    while (!(*pliter)->empty()) {
      blob->AddToInParticles((*pliter)->front());
      (*pliter)->pop_front();
    }
  }
  return blob;
}
