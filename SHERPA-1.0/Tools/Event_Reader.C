#include "Event_Reader.H"
#include "MyStrStream.H"
#include "Message.H"
#include <iostream>

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;


Event_Reader::Event_Reader(const std::string & path,const std::string & file) :
  m_path(path), m_file(file), m_generator(std::string("Unknown generator")),
  m_inputmode(0), m_eventmode(0), m_phasemode(-1)
{
  std::string filename = m_path + m_file; 
  p_instream = new std::ifstream(filename.c_str()); 

  if (!p_instream->good()) {
    msg.Error()<<"ERROR: Event file "<<filename<<" not found."<<std::endl
	       <<"   Will abort the run."<<std::endl;
    abort();
  }
  InitialSettings();
  msg.Events()<<"Generator  = "<<m_generator<<std::endl
	      <<"Input mode = "<<m_inputmode<<std::endl
	      <<"Event mode = "<<m_eventmode<<std::endl
	      <<"Phase mode = "<<m_phasemode<<std::endl;
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
	

bool Event_Reader::FillBlobs(Blob_List * blobs) 
{
  switch (m_phasemode) {
  case -1: 
    return ReadInEvent(blobs);
  }
  msg.Error()<<"Error in Event_Reader::FillBlobs."<<std::endl
	     <<"   Phasemode = "<<m_phasemode<<" is not specified so far."<<std::endl
	     <<"   Don't know what to do with file : "<<m_file<<", will abort the run."<<std::endl;
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
	  msg.Events()<<"End of file. Go to next file |"<<tmp<<"|"<<std::endl;
	  while(tmp.length()>0) {
	    if (tmp[0]==' ') tmp = tmp.substr(1);
	    else {
	      pos      = tmp.find(string(" "));
	      tmp      = tmp.substr(0,pos);
	      m_file   = tmp;
	      new_file = true;
	      break;
	    }
	  }
	  if (new_file) {
	    std::string filename = m_path + m_file; 
	    p_instream->close();
	    delete p_instream;
	    p_instream = new std::ifstream(filename.c_str()); 
	    if (!p_instream->good()) {
	      msg.Error()<<"ERROR in  Event_Reader::ReadInEvent."<<std::endl
		       <<"   New event file "<<filename<<" not found."<<std::endl
		       <<"   Will abort the run."<<std::endl;
	      abort();
	    }
	    break;
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
  msg.Error()<<"Error in Event_Reader::ReadInEvent."<<std::endl
	     <<"   Input/Event mode = "<<m_inputmode<<" / "<<m_eventmode<<" is not specified so far."<<std::endl
	     <<"   Don't know what to do with file : "<<m_file<<", will abort the run."<<std::endl;
  abort();
  return false;
}

bool Event_Reader::ReadInSimpleHepEvtEvent(Blob_List * blobs) 
{
  std::string    buffer,tmp,IS,FS;
  unsigned int   pos;
  int            kfc, minhadron,mother;
  vector<int>    ISc;
  vector<int>    FSc;
  Particle     * part;
  Blob         * hardblob, * showerblob, * hadronblob;
  bool           endevent,hard,shower,hadron;
  endevent  = hard = shower = hadron = false;
  minhadron = 0;

  hardblob         = new Blob();
  hardblob->SetType(btp::Signal_Process);
  hardblob->SetTypeSpec(m_generator);
  hardblob->SetId(0);
  hardblob->SetPosition(Vec4D(0.,0.,0.,0.));
  hardblob->SetStatus(1);
  hardblob->SetBeam(-1);
  hardblob->SetWeight(m_weight);
  blobs->push_back(hardblob);
  hardblob->AddData("ME_Weight",new Blob_Data<double>(m_weight));

  showerblob         = new Blob();
  showerblob->SetType(btp::FS_Shower);
  showerblob->SetTypeSpec(m_generator);
  showerblob->SetId(1);
  showerblob->SetPosition(Vec4D(0.,0.,0.,0.));
  showerblob->SetStatus(1);
  showerblob->SetBeam(-1);
  blobs->push_back(showerblob);

  hadronblob         = new Blob();
  hadronblob->SetType(btp::Fragmentation);
  hadronblob->SetTypeSpec(m_generator);
  hadronblob->SetId(2);
  hadronblob->SetPosition(Vec4D(0.,0.,0.,0.));
  hadronblob->SetStatus(1);
  hadronblob->SetBeam(-1);
  blobs->push_back(hadronblob);

  for (;;) {
    if (!p_instream->eof()) {
      endevent = false;
      getline(*p_instream,buffer);
      buffer  += std::string(" ");
      if (buffer[0] != '%' && buffer[1]!='%' && buffer.length()>0) {
	while(buffer.length()>0) {
	  if (buffer[0]==' ') buffer = buffer.substr(1);
	  else break;
	}
	if (buffer.find("End event")!=std::string::npos) {
	  endevent = true;
	  break;
	}
	if (buffer.find("Event :")!=std::string::npos) {
	  msg.Error()<<"Error in Event_Reader::ReadInSimpleHepEvtEvent"<<std::endl
		     <<"   Found Event start for event end."<<std::endl
		     <<"   Consider the current event finished and continue."<<std::endl;
	  endevent = true;
	  break;
	}

	// Fill signal blob
	if (buffer.find("->")!=std::string::npos) {
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
	      FSc.push_back(kfc);
	      FS  = FS.substr(pos);
	    }
	  }
	  int hardsize = ISc.size()+FSc.size();
	  int control;
	  Vec4D cms = Vec4D(0.,0.,0.,0.);
	  for (int i=0;i<hardsize;i++) {
	    getline(*p_instream,buffer);
	    if (i<ISc.size()) control = ISc[i];
	                 else control = FSc[i-ISc.size()]; 
	    part    = TranslateFromInput(buffer,mother,control);
	    if (i<ISc.size()) {
	      hardblob->AddToInParticles(part);
	      cms += part->Momentum();
	      part->SetStatus(2);
	      part->SetInfo('G');
	    }
	    else {
	      hardblob->AddToOutParticles(part);
	      showerblob->AddToInParticles(part);
	      part->SetStatus(1);
	      part->SetInfo('H');
	      minhadron = part->Number()+1;
	    }
	  }
	  hardblob->SetCMS(cms);
	  hard = true;
	}	// End of fill signal blob
	else {
	  // Fill other blobs
	  part    = TranslateFromInput(buffer,mother,0);
	  if (part->Status()>2) delete part;
	  else {
	    if (part->Flav().IsHadron()) {
	      part->SetInfo('P');
	      hadronblob->AddToOutParticles(part);
	      hadron = true;
	    }
	    else {
	      part->SetInfo('F');
	      if (part->Status()>1) part->SetInfo('f');
	      if ((!part->Flav().IsPhoton()) || 
		  (part->Flav().IsPhoton() &&  mother<minhadron)) {
		showerblob->AddToOutParticles(part);
		hadronblob->AddToInParticles(part);
		if (!part->Flav().Strong() &&
		    !part->Flav().IsDiQuark()) hadronblob->AddToOutParticles(new Particle(part));
		minhadron = part->Number();
		shower = true;
	      }
	      else {
		hadronblob->AddToOutParticles(part);
	      }
	    }
	  }
	}
      }
    }
    if (endevent) {
      if (hard && shower && hadron) return true;
      return false;
    }
  }
  return true;
}



Particle * Event_Reader::TranslateFromInput(std::string buffer,int & mother,const int control)
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
      msg.Error()<<"Error in Event_Reader::TranslateFromInput."<<std::endl
		 <<"   Particle ID and control number do not coincide : "<<flags[1]<<" vs. "<<control<<std::endl
		 <<"   in event number "<<m_evtnumber<<" of file "<<m_file<<"."<<std::endl
		 <<"   Abort the run and check."<<std::endl;
      abort();
    }
  }
  if (flags.size()<4 || numbers.size()<4) {
    msg.Error()<<"Error in Event_Reader::TranslateFromInput."<<std::endl
	       <<"   Not enough information provided for particle construction :"
	       <<"   ("<<flags.size()<<" "<<numbers.size()<<")"<<std::endl
	       <<"   in event number "<<m_evtnumber<<" of file "<<m_file<<"."<<std::endl
	       <<"   Abort the run and check."<<std::endl;
    abort();
  }
  Flavour flav; flav.FromHepEvt(flags[1]);
  Vec4D   momentum = Vec4D(numbers[0],numbers[1],numbers[2],numbers[3]);
  Particle * part  = new Particle(flags[0],flav,momentum);
  part->SetStatus(flags[2]);
  mother           = flags[3];
  return part;
}

