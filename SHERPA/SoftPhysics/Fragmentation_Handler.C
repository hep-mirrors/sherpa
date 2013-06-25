#include "SHERPA/SoftPhysics/Fragmentation_Handler.H"

#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Smart_Pointer.H"
#include "ATOOLS/Org/Return_Value.H"
#ifdef USING__Ahadic
#include "AHADIC++/Main/Ahadic.H"
#endif

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;
#ifdef USING__Ahadic
using namespace AHADIC;
#include "AHADIC++/Tools/Hadron_Init.H"
#else
#error Ahadic must be enabled for filling the hadron tables.
#endif

Fragmentation_Handler::Fragmentation_Handler(string _dir,string _file):
  m_dir(_dir), m_file(_file), m_mode(0)
#ifdef USING__Ahadic
  ,p_ahadic(NULL)
#endif
#ifdef USING__PYTHIA
  ,p_lund(NULL)
#endif
{
  Data_Reader dr(" ",";","!","=");
  dr.AddComment("#");
  dr.AddWordSeparator("\t");
  dr.AddIgnore("[");
  dr.AddIgnore("]");
  dr.SetInputPath(m_dir);
  dr.SetInputFile(m_file);
  m_fragmentationmodel=dr.GetValue<string>("FRAGMENTATION",string("Ahadic"));
  m_shrink=dr.GetValue<int>("COMPRESS_PARTONIC_DECAYS",1);
  m_flagpartonics=dr.GetValue<int>("FLAG_PARTONIC_DECAYS",0);
  if (m_fragmentationmodel==string("Lund")) {
#ifndef USING__PYTHIA
    THROW(fatal_error, "Fragmentation/decay interface to Pythia has not been "+
          string("enabled during compilation (./configure --enable-pythia)."));
#else
    m_sfile=dr.GetValue<string>("LUND_FILE",string("Lund.dat"));
    Hadron_Init init;
    init.Init();
    init.OverWriteProperties(dr);
    ATOOLS::OutputHadrons(msg->Tracking());
    p_lund = new Lund_Interface(m_dir,m_sfile);
    m_mode=1;
    exh->AddTerminatorObject(this);
    // hack for particle initialization, because we don't want to replicate
    // this method in the obsolete Lund Interface.
    return;
#endif
  }
#ifdef USING__Ahadic
  else if (m_fragmentationmodel==string("Ahadic")) {
    m_sfile=dr.GetValue<string>("AHADIC_FILE",m_file);
    Hadron_Init init;
    init.Init();
    init.OverWriteProperties(dr);
    ATOOLS::OutputHadrons(msg->Tracking());
    p_ahadic = new AHADIC::Ahadic(m_dir,m_sfile,true);
    m_mode=2;
    exh->AddTerminatorObject(this);
    return;
  }
#endif
  else if (m_fragmentationmodel==string("Off") ||
           m_fragmentationmodel==string("None") ||
           m_fragmentationmodel==string("0")) return;
  THROW(critical_error,"Fragmentation model not implemented.");
}
   
Fragmentation_Handler::~Fragmentation_Handler() 
{
#ifdef USING__PYTHIA
  if (p_lund!=NULL)   { delete p_lund;   p_lund   = NULL;   }
#endif
#ifdef USING__Ahadic
  if (p_ahadic!=NULL) { delete p_ahadic; p_ahadic = NULL;   }
#endif
  exh->RemoveTerminatorObject(this);
}

void Fragmentation_Handler::PrepareTerminate() 
{
  std::string path(rpa->gen.Variable("SHERPA_STATUS_PATH"));
  if (path=="") return;
  Copy(m_dir+"/"+m_sfile,path+"/"+m_sfile);
}

Return_Value::code 
Fragmentation_Handler::PerformFragmentation(Blob_List *bloblist,
					    Particle_List *particlelist) 
{
  if (m_mode==0 || bloblist->size()==0) return Return_Value::Nothing;

  Return_Value::code success;
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
  switch (m_mode) {
#ifdef USING__PYTHIA
  case 1  : 
    success = p_lund->Hadronize(bloblist);
    if (m_shrink>0 && success==Return_Value::Success) Shrink(bloblist);
    return success;
#endif
#ifdef USING__Ahadic
  case 2  : 
    success = p_ahadic->Hadronize(bloblist);
    if (success!=Return_Value::Success &&
	success!=Return_Value::Nothing) {
      msg_Tracking()<<"Potential problem in "<<METHOD<<":\n"<<(*bloblist)<<"\n";
    }
    if (m_shrink>0 && success==Return_Value::Success) Shrink(bloblist);
    return success;
#endif
  default : 
    msg_Error()<<"ERROR in "<<METHOD<<":\n"
	       <<"   Unknown hadronization model in mode = "<<m_mode<<".\n"
	       <<"   Abort the run.\n";
    THROW(critical_error,"Fragmentation model not implemented.");
  }
}

void Fragmentation_Handler::Shrink(Blob_List * bloblist) {
  list<Blob *> deleteblobs;
  Particle_Vector * parts;
  for (Blob_List::reverse_iterator blit=bloblist->rbegin();
       blit!=bloblist->rend();++blit) {
    Blob * blob = (*blit);
    if (blob->Type()==btp::Fragmentation) {
      if (blob->Status()!=blob_status::needs_hadrondecays) return;
      Blob * showerblob(blob->InParticle(0)->ProductionBlob());
      Blob * decblob(showerblob->InParticle(0)->ProductionBlob());
      showerblob->DeleteInParticles(0);
      showerblob->DeleteOutParticles(0);
      deleteblobs.push_back(blob);
      deleteblobs.push_back(showerblob);
      while (!blob->GetOutParticles().empty()) {
	Particle * part = 
	  blob->RemoveOutParticle(blob->GetOutParticles().front());
	decblob->AddToOutParticles(part);
      }
      if (m_flagpartonics>0) decblob->AddData("Partonic",new Blob_Data<int>(1));
    }
  }
  for (list<Blob *>::iterator blit=deleteblobs.begin();
       blit!=deleteblobs.end();blit++) bloblist->Delete((*blit));
}



Return_Value::code Fragmentation_Handler::ExtractSinglets(Blob_List * bloblist)
{
  Particle  * part(NULL);
  Blob      * blob(NULL);
  std::vector<SP(Part_List)> plists;
  plists.push_back(new Part_List);
  for (Blob_List::iterator blit=bloblist->begin();
       blit!=bloblist->end();++blit) {
    if ((*blit)->Has(blob_status::needs_hadronization)) {
      // If not coming from hadron decays, fill default plists[0].
      // If from hadron decays, create separate plist for each blob
      // such that there is a one-to-one correspondence between
      // fragmentation outcome and hadron decay. This is needed for setting
      // the correct vertex position, and to reject exclusive final states
      SP(Part_List) plist(plists[0]);
      Blob* upstream_blob=(*blit)->UpstreamBlob();
      if (upstream_blob && upstream_blob->Type()==btp::Hadron_Decay) {
        plist=new Part_List;
        plists.push_back(plist);
      }
      
      std::vector<Particle*> taus;
      for (int i=0;i<(*blit)->NOutP();i++) {
	part = (*blit)->OutParticle(i); 
	if (part->Status()==part_status::active && 
	    part->Info()!='G' && part->Info()!='I') {
	  if (part->GetFlow(1)!=0 || part->GetFlow(2)!=0) {
	    if (part->GetFlow(1)==part->GetFlow(2)) {
	      msg_Error()<<"Error in "<<METHOD<<":\n"
			 <<"   Blob with funny colour assignements.\n"
			 <<"   Will demand new event and hope for the best.\n";
	      return Return_Value::New_Event;
	    }
	    plist->push_back(part);
	    part->SetStatus(part_status::fragmented);
	  }
	  else if (part->Flav().Kfcode()==kf_tau || part->Flav().IsHadron()) {
            taus.push_back(part);
	  }
	}
      }
      if (taus.size()>0) {
        blob = new Blob();
        blob->SetId();
        blob->SetType(btp::Fragmentation);
        blob->SetStatus(blob_status::needs_hadrondecays);
        for (size_t i=0;i<taus.size();++i) {
          blob->AddToInParticles(taus[i]);
          taus[i]->SetStatus(part_status::decayed);
          blob->AddToOutParticles(new Particle((*taus[i])));
          blob->GetOutParticles().back()->SetStatus(part_status::active);
        }
        bloblist->push_back(blob);
      }
      (*blit)->UnsetStatus(blob_status::needs_hadronization);
    }
  }
  if (plists[0]->empty() && plists.size()<2) {
    msg_Debugging()<<"WARNING in Lund_Interface::PrepareFragmentationBlob:\n"
		   <<"   No coloured particle found leaving shower blobs.\n";
    return Return_Value::Nothing;
  }
  
  
  Return_Value::code ret(Return_Value::Success);
  for (size_t i=0; i<plists.size(); ++i) {
    if (plists[i]->empty()) continue;
    SP(Part_List) plist=plists[i];
    int  col1, col2;
    bool hit1, hit2;
    Part_List * pli(NULL);
    vector<SP(Part_List)> partlists; 
    int plsize;
    do {
      plsize=plist->size();
      hit1    = false;
      for (Part_Iterator pit=plist->begin();pit!=plist->end();++pit) {
	col1 = (*pit)->GetFlow(1);
	col2 = (*pit)->GetFlow(2);
	if (col1!=0 && col1==col2) {
	  msg_Error()<<"Error in "<<METHOD<<":\n"
		     <<"   Blob with funny colour assignements.\n"
		     <<"   Will demand new event and hope for the best.\n";
	  return Return_Value::New_Event;
	}
	if (col1!=0 && col2==0) {
	  hit1 = true;
	  pli  = new Part_List;
	  pli->push_back((*pit));
	  pit  = plist->erase(pit);
	  partlists.push_back(pli);
	  do {
	    hit2 = false;
	    for (Part_Iterator pit1=plist->begin();pit1!=plist->end();++pit1) {
	      if ((int)((*pit1)->GetFlow(2))==col1) {
		col1 = (*pit1)->GetFlow(1);
		pli->push_back((*pit1));
		pit1 = plist->erase(pit1);
		hit2 = true;
		break;
	      }
	    }
	  } while (hit2 && col1!=0);
	}
	if (hit1) break;
      }
      if (!hit1) {
	for (Part_Iterator pit=plist->begin();pit!=plist->end();++pit) {
	  col1 = (*pit)->GetFlow(1);
	  col2 = (*pit)->GetFlow(2);
	  if (col1!=0 && col2!=0) {
	    hit1 = true;
	    pli  = new Part_List;
	    pli->push_back((*pit));
	    pit  = plist->erase(pit);
	    partlists.push_back(pli);
	    do {
	      hit2 = false;
	      for (Part_Iterator pit1=plist->begin();
		   pit1!=plist->end();++pit1) {
		if ((int)((*pit1)->GetFlow(2))==col1) {
		  col1 = (*pit1)->GetFlow(1);
		  pli->push_back((*pit1));
		pit1 = plist->erase(pit1);
		hit2 = true;
		break;
		}
	      }
	    } while (hit2 && col1!=col2);
	  }
	  if (hit1) break;
	}
      }
      if (!hit1 && plist->size()==plsize) {
	msg_Error()<<"Error in "<<METHOD<<":\n"
		   <<"   Will throw new event.\n";
	return Return_Value::New_Event;
      }
    } while(plist->size()>0);
    
    if (plist->empty()) {
      blob = new Blob();
      blob->SetId();
      blob->SetType(btp::Fragmentation);
      blob->SetStatus(blob_status::needs_hadronization);
      bloblist->push_back(blob);
      for (vector<SP(Part_List)>::iterator pliter=partlists.begin();
	   pliter!=partlists.end();pliter++) {
	while (!(*pliter)->empty()) {
	  blob->AddToInParticles((*pliter)->front());
	  (*pliter)->pop_front();
	}
      }
      
      ret=Return_Value::Success;
    }
    else {
      ret=Return_Value::Error;
      break;
    }
  }
  return ret;
}

