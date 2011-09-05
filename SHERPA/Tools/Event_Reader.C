#include "SHERPA/Tools/Event_Reader.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include <iostream>

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;



Event_Reader::Event_Reader(const std::string & path,const std::string & file) :
  Event_Reader_Base(path,file), 
  m_generator(std::string("Unknown generator")),
  f_gz(false), m_inputmode(0), m_eventmode(0), m_phasemode(-1)
{
  std::string filename=m_path+m_file;
  msg_Out()<<" -> "<<filename<<"\n";

  size_t pst=filename.find(".gz");
  size_t lng=filename.length();

  if(pst!=std::string::npos) {
    if(pst==lng-3) {
      f_gz=true;
      int stat=0;
      stat=system((std::string("gunzip ")+filename).c_str());
      filename.resize(pst);
      m_file.resize(m_file.length()-3);
    } else {
      msg_Error()<<"Error in "<<__PRETTY_FUNCTION__<<":\n"
		 <<"   Cannot handle event file "<<filename<<".\n"
		 <<"   Will abort the run. Please check."<<std::endl;
      abort();
    }
  }

  pst=m_add.find("/");
  if(pst!=std::string::npos) {
    m_add.resize(pst+1);
    m_file=m_file.substr(pst+1);
  }

  msg_Out()<<" -> "<<filename<<" , "<<m_add<<" , "<<m_file<<"\n";
  p_instream=new std::ifstream(filename.c_str());

  if(!p_instream->good()) {
    msg_Error()<<"ERROR: Event file "<<filename<<" not found."<<std::endl
	       <<"   Will abort the run."<<std::endl;
    abort();
  }

  InitialSettings();
  msg_Out()<<"Generator  = "<<m_generator<<"\n"
	   <<"Input mode = "<<m_inputmode<<"\n"
	   <<"Event mode = "<<m_eventmode<<"\n"
	   <<"Phase mode = "<<m_phasemode<<std::endl;
  if(f_gz) msg_Out()<<"G(un)zip handling is switched on for the whole run!\n......."<<std::endl;
}





void Event_Reader::InitialSettings() 
{
  std::string  buffer;
  unsigned int pos;
  for (;;) {
    if (!p_instream->eof()) {
      getline(*p_instream,buffer);
      buffer += std::string(" ");
      if (buffer.find("Generated events start here.")!=std::string::npos) break;
      if (buffer[0] != '%' && buffer[1] != '%' && buffer.length()>0) {
	if (buffer.find("Generator   = ")!=std::string::npos) {
	  pos = buffer.find("Generator   = ");
	  buffer = buffer.substr(pos+14);
	  while(buffer.length()>0) {
	    if (buffer[0]==' ') buffer = buffer.substr(1);
	    else {
	      pos = buffer.find(string(" "));
	      if (pos>0) m_generator = buffer.substr(0,pos);
	      break;
	    }
	  }
	}
	if (buffer.find("Output mode = ")!=std::string::npos) {
	  pos = buffer.find("Output mode = ");
	  buffer = buffer.substr(pos+14);
	  while(buffer.length()>0) {
	    if (buffer[0]==' ') buffer = buffer.substr(1);
	    else {
	      pos    = buffer.find(string(" "));
	      buffer = buffer.substr(0,pos);
	      MyStrStream str;
	      str<<buffer;
	      str>>m_inputmode;
	      break;
	    }
	  }
	}
	if (buffer.find("Event mode  = ")!=std::string::npos) {
	  pos = buffer.find("Event mode  = ");
	  buffer = buffer.substr(pos+14);
	  while(buffer.length()>0) {
	    if (buffer[0]==' ') buffer = buffer.substr(1);
	    else {
	      pos    = buffer.find(string(" "));
	      buffer = buffer.substr(0,pos);
	      MyStrStream str;
	      str<<buffer;
	      str>>m_eventmode;
	      break;
	    }
	  }
	}
	if (buffer.find("Phase mode  = ")!=std::string::npos) {
	  pos = buffer.find("Phase mode  = ");
	  buffer = buffer.substr(pos+14);
	  while(buffer.length()>0) {
	    if (buffer[0]==' ') buffer = buffer.substr(1);
	    else {
	      pos    = buffer.find(string(" "));
	      buffer = buffer.substr(0,pos);
	      MyStrStream str;
	      str<<buffer;
	      str>>m_phasemode;
	      break;
	    }
	  }
	}
      }
    }
  }
}





void Event_Reader::CloseFile() {
  if(p_instream) {
    p_instream->close();
    delete p_instream;
    p_instream=NULL;
    int stat=0;
    if(f_gz) stat=system((std::string("gzip ")+m_path+m_add+m_file).c_str());
  } else {
    msg_Error()<<__PRETTY_FUNCTION__<<":\n   Warning: File already closed."<<std::endl;
  }
}





bool Event_Reader::FillBlobs(Blob_List * blobs) 
{
  bool result;
  long nev=rpa->gen.NumberOfEvents();
  switch (m_phasemode) {
  case -1:
    result=ReadInEvent(blobs);
    if(nev==rpa->gen.NumberOfGeneratedEvents()) CloseFile();
    return result;
  }
  msg_Error()<<"Error in Event_Reader::FillBlobs."<<std::endl
	     <<"   Phasemode = "<<m_phasemode<<" is not specified so far."<<std::endl
	     <<"   Don't know what to do with file : "<<m_file<<", will abort the run."
	     <<std::endl;
  abort();
  return false;
}





bool Event_Reader::ReadInEvent(Blob_List * blobs) 
{
  std::string    buffer,tmp;
  unsigned int   pos;
  bool           found_event, new_file;
  m_weight = 1.;
  for (;;) {
    if (!p_instream->eof()) {
      new_file = found_event = false;
      getline(*p_instream,buffer);
      buffer += std::string(" ");
      if (buffer[0] != '%' && buffer[1] != '%' && buffer.length()>0) {
	while(buffer.length()>0) {
	  if (buffer[0]==' ') buffer = buffer.substr(1);
	  else break;
	}
	if (buffer.find("Next file  =")!=std::string::npos) {
	  pos = buffer.find("Next file  =");
	  tmp = buffer.substr(pos+12);
	  msg_Info()<<"End of file. Go to next file |"<<tmp<<"|"<<std::endl;
	  while(tmp.length()>0) {
	    if (tmp[0]==' ') tmp = tmp.substr(1);
	    else {
	      CloseFile();
	      pos      = tmp.find(string(" "));
	      tmp      = tmp.substr(0,pos);
	      m_file   = tmp;
	      new_file = true;
	      break;
	    }
	  }
	  if (new_file) {
	    std::string filename=m_path+m_add+m_file;
	    int stat=0;
	    if(f_gz) stat=system((std::string("gunzip ")+filename+std::string(".gz")).c_str());
	    msg_Info()<<" => "<<filename<<"\n";
	    p_instream=new std::ifstream(filename.c_str());
	    if(!p_instream->good()) {
	      msg_Error()<<"ERROR: Event file "<<filename<<" not found."<<std::endl
			 <<"   Will abort the run."<<std::endl;
	      abort();
	    }
	    continue;
	  }
	}
	if (buffer.find("Event :")!=std::string::npos) {
	  pos = buffer.find("Event :");
	  tmp = buffer.substr(pos+8);
	  while(tmp.length()>0) {
	    if (tmp[0]==' ') tmp = tmp.substr(1);
	    else {
	      pos = tmp.find(string(" "));
	      tmp = tmp.substr(0,pos);
	      MyStrStream str;
	      str<<tmp;
	      str>>m_evtnumber;
	      found_event = true;
	      break;
	    }
	  }
	}
	if (buffer.find("Weight :")!=std::string::npos) {
	  pos = buffer.find("Weight :");
	  tmp = buffer.substr(pos+8);
	  while(tmp.length()>0) {
	    if (tmp[0]==' ') tmp = tmp.substr(1);
	    else {
	      pos = tmp.find(string(" "));
	      tmp = tmp.substr(0,pos);
	      MyStrStream str;
	      str<<tmp;
	      str>>m_weight;
	      break;
	    }
	  }
	}
	if (found_event) break;
      }
    }
  }
  switch (m_inputmode) {
  case 0: 
    // HepEvt Format
    switch (m_eventmode) {
    case 0: 
      // Simple Events
      return ReadInSimpleHepEvtEvent(blobs);
    }
  }
  msg_Error()<<"Error in Event_Reader::ReadInEvent."<<std::endl
	     <<"   Input/Event mode = "<<m_inputmode<<" / "<<m_eventmode
	     <<" is not specified so far."<<std::endl
	     <<"   Don't know what to do with file : "<<m_file
	     <<", will abort the run."<<std::endl;
  abort();
  return false;
}





bool Event_Reader::ReadInSimpleHepEvtEvent(Blob_List * blobs) 
{
  std::string    buffer,tmp,IS,FS;
  unsigned int   pos;
  int            kfc, mother, minhadron=0;
  vector<int>    ISc;
  vector<int>    FSc;
  Particle*      part=NULL;
  Blob         * hardblob, * showerblob, * hadronblob;
  bool           fhard, fshower, fhadron;

  hardblob         = new Blob();
  hardblob->SetType(btp::Signal_Process);
  hardblob->SetTypeSpec(m_generator);
  hardblob->SetId();
  hardblob->SetPosition(Vec4D(0.,0.,0.,0.));
  hardblob->SetStatus(blob_status::code(30));
  hardblob->SetWeight(m_weight);
  blobs->push_back(hardblob);
  hardblob->AddData("ME_Weight",new Blob_Data<double>(m_weight));

  showerblob         = new Blob();
//   showerblob->SetType(btp::FS_Shower);
  showerblob->SetTypeSpec(m_generator);
  showerblob->SetId();
  showerblob->SetPosition(Vec4D(0.,0.,0.,0.));
  showerblob->SetStatus(blob_status::code(28));
  blobs->push_back(showerblob);

  hadronblob         = new Blob();
  hadronblob->SetType(btp::Fragmentation);
  hadronblob->SetTypeSpec(m_generator);
  hadronblob->SetId();
  hadronblob->SetPosition(Vec4D(0.,0.,0.,0.));
  hadronblob->SetStatus(blob_status::needs_hadrondecays);
  blobs->push_back(hadronblob);

  for (;;) {
    if(!p_instream->eof()) {
      getline(*p_instream,buffer);
      buffer  += std::string(" ");
      if(buffer[0] != '%' && buffer[1]!='%' && buffer.length()>0) {
	while(buffer.length()>0) {
	  if (buffer[0]==' ') buffer = buffer.substr(1);
	  else break;
	}
	if (buffer.find("End event")!=std::string::npos) break;
	if (buffer.find("Event :")!=std::string::npos) {
	  msg_Error()<<"Error in Event_Reader::ReadInSimpleHepEvtEvent\n"
		     <<"   Found Event start for event end.\n"
		     <<"   Consider the current event finished and continue."<<std::endl;
	  break;
	}

	part=NULL;
	fhard=fshower=fhadron=true;

	if (buffer.find("->")!=std::string::npos) {    // Fill signal blob
	  pos = buffer.find("->");
	  IS  = buffer.substr(0,pos);
	  FS  = buffer.substr(pos+2)+std::string(" ");
	  ISc.clear();
	  FSc.clear();
	  while(IS.length()>0) {
	    if (IS[0]==' ') IS = IS.substr(1);
	    else {
	      pos = IS.find(string(" "));
	      tmp = IS.substr(0,pos);
	      MyStrStream str;
	      str<<tmp;
	      str>>kfc;
	      ISc.push_back(kfc);
	      IS  = IS.substr(pos);
	    }
	  }
	  while(FS.length()>0) {
	    if (FS[0]==' ') FS = FS.substr(1);
	    else {
	      pos = FS.find(string(" "));
	      tmp = FS.substr(0,pos);
	      MyStrStream str;
	      str<<tmp;
	      str>>kfc;
	      if (kfc!=0) FSc.push_back(kfc);
	      FS  = FS.substr(pos);
	    }
	  }
	  size_t hardsize = ISc.size()+FSc.size();
	  int control;
	  Vec4D cms = Vec4D(0.,0.,0.,0.);
	  for (size_t i=0;i<hardsize;i++) {
	    getline(*p_instream,buffer);
	    if (i<ISc.size()) control = ISc[i];
	                 else control = FSc[i-ISc.size()]; 
	    part=TranslateFromInput(buffer,mother,control);
	    if(part==NULL) break;
	    if (i<ISc.size()) {
	      hardblob->AddToInParticles(part);
	      cms += part->Momentum();
	      part->SetStatus(part_status::decayed);
	      part->SetInfo('G');
	    }
	    else {
	      hardblob->AddToOutParticles(part);
	      showerblob->AddToInParticles(part);
	      part->SetStatus(part_status::active);
	      part->SetInfo('H');
	      minhadron = part->Number()+1;
	    }
	  }
	  if(part==NULL) fhard=false;
	  else hardblob->SetCMS(cms);
	}	// End of fill signal blob
	else {    // Fill other blobs
	  part=TranslateFromInput(buffer,mother,0);
	  if(!part) { fshower=fhadron=false;}
	  else {
	    if(part->Status()!=part_status::active &&
	       part->Status()!=part_status::decayed &&
	       part->Status()!=part_status::fragmented) delete part;
	    else {
	      if(part->Flav().IsHadron()) {
		part->SetInfo('P');
		hadronblob->AddToOutParticles(part);
	      }
	      else {
		part->SetInfo('F');
		if(part->Status()!=part_status::active) part->SetInfo('f');
		if((!part->Flav().IsPhoton()) || (part->Flav().IsPhoton() &&  mother<minhadron)) {
		  showerblob->AddToOutParticles(part);
		  hadronblob->AddToInParticles(part);
		  if(!part->Flav().Strong() && !part->Flav().IsDiQuark())
		    hadronblob->AddToOutParticles(new Particle(*part));
		  minhadron=part->Number();
		}
		else hadronblob->AddToOutParticles(part);
	      }
	    }
	  }
	}    //End of fill other blobs
	if(fhard && fshower && fhadron);
	else {
	  msg_Info()<<__PRETTY_FUNCTION__<<":\n "
		    <<"Warning: Event could not be reconstructed..."
		    <<"fhard="<<fhard<<",  fshower="<<fshower<<",  fhadron="<<fhadron<<std::endl;
	  return false;
	}
      }
    }
  }
  return true;
}





Particle* Event_Reader::TranslateFromInput(std::string buffer, int& mother, const int control)
{
  unsigned int     pos;
  std::string      tmp;
  int              flag;
  double           number;
  vector<int>      flags;
  vector<double>   numbers;
  buffer += std::string(" ");
  while(buffer.length()>0) { 
    if (buffer[0]==' ') buffer = buffer.substr(1); 
    else {
      pos = buffer.find(string(" "));
      tmp = buffer.substr(0,pos);
      MyStrStream str;
      str<<tmp;
      if (tmp.find(string("."))!=std::string::npos) {
	str>>number;
	numbers.push_back(number);
      }
      else {
	str>>flag;
	flags.push_back(flag);
      }
      buffer  = buffer.substr(pos);
    }
  }
  if (control != 0) {
    if (flags[1]!=control) {
      msg_Error()<<"Error in Event_Reader::TranslateFromInput."<<std::endl
		 <<"   Particle ID and control number do not coincide : "<<flags[1]<<" vs. "
		 <<control<<std::endl
		 <<"   in event number "<<m_evtnumber<<" of file "<<m_file<<"."<<std::endl
		 <<"   Abort the run and check."<<std::endl;
      abort();
    }
  }
  Particle* part=NULL;
  if (flags.size()<4 || numbers.size()<4) {
    msg_Error()<<"Error in "<<__PRETTY_FUNCTION__<<":\n"
	       <<"   Not enough information provided for particle construction: "
	       <<"   ("<<flags.size()<<" "<<numbers.size()<<")\n"
	       <<"   in event number "<<m_evtnumber<<" of file "<<m_file<<".\n"
	       <<"   Will return (nil)-particle."<<std::endl;
    return part;    //abort();
  }
  Flavour flav;
  flav.FromHepEvt(flags[1]);
  Vec4D momentum=Vec4D(numbers[0],numbers[1],numbers[2],numbers[3]);
  part=new Particle(flags[0],flav,momentum);
  part->SetStatus(part_status::code(flags[2]));
  part->SetFinalMass(numbers[4]);
  mother=flags[3];
  return part;
}
