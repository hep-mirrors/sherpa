#include "HepEvt_Interface.H"

#include "Run_Parameter.H"
#include "Message.H"
#include "MyStrStream.H"

#include <iomanip>
#include <stdio.h>

using namespace ATOOLS;
using namespace std;

HepEvt_Interface::HepEvt_Interface(int _generator) : 
  p_instream(NULL), p_outstream(NULL), 
  m_evtnumber(0), m_nhep(-1), m_filesize(0), m_evtcount(0), 
  m_generator(gtp::code(_generator)), p_phep(NULL), p_vhep(NULL),
  p_jmohep(NULL), p_jdahep(NULL), p_isthep(NULL), p_idhep(NULL)
{ }

HepEvt_Interface::HepEvt_Interface(gtp::code _generator) : 
  p_instream(NULL), p_outstream(NULL), 
  m_evtnumber(0), m_nhep(-1), m_filesize(0), m_evtcount(0), 
  m_generator(_generator), p_phep(NULL), p_vhep(NULL),
  p_jmohep(NULL), p_jdahep(NULL), p_isthep(NULL), p_idhep(NULL)
{ }

HepEvt_Interface::HepEvt_Interface(bool _io,int _mode,
				   const std::string & _path, 
				   const std::string & _file,
				   const int _filesize) : 
  m_io(_io), m_mode(_mode), m_path(_path), m_file(_file), 
  p_instream(NULL), p_outstream(NULL), 
  m_evtnumber(0), m_nhep(-1), m_filesize(_filesize), m_evtcount(0), 
  m_generator(gtp::Unspecified)
{
  // io = true : Output mode, Sherpa2HepEvt
  std::string filename = m_path+std::string("/")+m_file;
  if (m_io) {
    if (m_mode>=10) {
      m_mode      -= 10;
    }
    else {
      filename += std::string(".0.evts"); 
      p_outstream = new std::ofstream(filename.c_str(),std::ios::out);
      if (!p_outstream->good()) { 
	msg.Error()<<"ERROR in HepEvt_Interface."<<std::endl
		   <<"   Could not open event file "<<filename<<"."<<std::endl
		   <<"   Will abort the run."<<std::endl;
	abort();
      }
      p_outstream->precision(10);
    }
  }
  else {
    p_instream = new std::ifstream(filename.c_str()); 
    if (!p_instream->good()) {
      msg.Error()<<"ERROR in HepEvt_Interface."<<std::endl
		 <<"   Event file "<<filename<<" not found."<<std::endl
		 <<"   Will abort the run."<<std::endl;
      abort();
    }
    std::string gentype;
    (*p_instream)>>gentype>>m_filesize;
    if (gentype==std::string("Sherpa")) m_generator = gtp::Sherpa;
    if (gentype==std::string("Herwig")) m_generator = gtp::Herwig;
    if (gentype==std::string("Pythia")) m_generator = gtp::Pythia;
    m_evtcount=0;
  }
  p_phep     = new double[5*s_maxentries];
  p_vhep     = new double[4*s_maxentries];
  p_jmohep   = new int[2*s_maxentries];
  p_jdahep   = new int[2*s_maxentries];
  p_isthep   = new int[s_maxentries];
  p_idhep    = new int[s_maxentries];
}

HepEvt_Interface::~HepEvt_Interface() 
{
  if (p_outstream) {
    p_outstream->close();
    delete p_outstream; p_outstream=NULL;
  }
  if (p_instream) {
    p_instream->close();
    delete p_instream; p_instream=NULL;
  }
  if (p_jmohep) { delete [] p_jmohep; p_jmohep=NULL; }
  if (p_jdahep) { delete [] p_jdahep; p_jdahep=NULL; }
  if (p_isthep) { delete [] p_isthep; p_isthep=NULL; }
  if (p_idhep)  { delete [] p_idhep;  p_idhep=NULL;  }
  if (p_phep)   { delete [] p_phep;   p_phep=NULL;   }
  if (p_vhep)   { delete [] p_vhep;   p_vhep=NULL;   }
}

/*------------------------------------------------------------------------------------
  Sherpa to HepEvt methods
------------------------------------------------------------------------------------*/


void HepEvt_Interface::ChangeOutStream(std::string & filename, long int evtsperfile)
{
  if (p_outstream->is_open()) p_outstream->close();
  p_outstream->open(filename.c_str(),std::ios::out);
  if (!p_outstream->good()) { 
    msg.Error()<<"ERROR in HepEvt_Interface::ChangeOutStream"<<std::endl
	       <<"   Could not change to event file "<<filename<<"."<<std::endl
	       <<"   Will abort the run."<<std::endl;
    abort();
  }
  p_outstream->precision(10);
  (*p_outstream)<<"Pythia "<<evtsperfile<<std::endl;
}

void HepEvt_Interface::ChangeOutStream()
{
  if (p_outstream->is_open()) p_outstream->close();
  std::string filename = m_path+"/"+m_file+"."+ToString(int(m_evtnumber/m_filesize))+".evts";
  p_outstream->open(filename.c_str(),std::ios::out);
  if (!p_outstream->good()) { 
    msg.Error()<<"ERROR in HepEvt_Interface::ChangeOutStream"<<std::endl
	       <<"   Could not change to event file "<<filename<<"."<<std::endl
	       <<"   Will abort the run."<<std::endl;
    abort();
  }
  p_outstream->precision(10);
}




bool HepEvt_Interface::Sherpa2HepEvt(Blob_List * const _blobs) {
  m_convertS2H.clear();

  m_evtnumber++;
  if ((m_evtnumber-1)%m_filesize==0 && p_outstream!=NULL) ChangeOutStream();

  int nhep = 0;

  ISBlobs2HepEvt(_blobs,nhep);
  HardBlob2HepEvt(_blobs,nhep);
  FSBlobs2HepEvt(_blobs,nhep);
  FragmentationBlob2HepEvt(_blobs,nhep);
  HadronDecayBlobs2HepEvt(_blobs,nhep);
  
  m_nhep=nhep;

  switch (m_mode) {
  case 0 :  break;
  case 2 :  WriteReducedHepEvt(nhep); break;
  case 3 :  WriteFormatedHepEvt(nhep); break;
  default:  WriteFullHepEvt(nhep);
  }
  return true;
}

void HepEvt_Interface::PrintHepEvtEvent(int nhep) 
{
  for (int i=0;i<nhep;++i) {
    msg.Out()<<i+1<<"  "<<p_isthep[i]<<" "<<p_idhep[i]<<" "<<p_jmohep[2*i]<<" "<<p_jmohep[2*i+1]
	     <<" "<<p_jdahep[2*i]<<" "<<p_jdahep[2*i+1]<<" \n ";
    for (int j=0;j<5;++j) msg.Out()<<p_phep[5*i+j]<<" ";
    msg.Out()<<"\n ";
    for (int j=0;j<4;++j) msg.Out()<<p_vhep[4*i+j]<<" ";
    msg.Out()<<"\n";
  }  
}

void HepEvt_Interface::WriteFullHepEvt(int nhep)
{
  (*p_outstream)<<"  "<<m_evtnumber<<" "<<nhep<<" \n";
  for (int i=0;i<nhep;++i) {
    (*p_outstream)<<i+1<<"  "<<p_isthep[i]<<" "<<p_idhep[i]<<" "<<p_jmohep[2*i]<<" "<<p_jmohep[2*i+1]
		  <<" "<<p_jdahep[2*i]<<" "<<p_jdahep[2*i+1]<<" \n ";
    for (int j=0;j<5;++j) (*p_outstream)<<p_phep[5*i+j]<<" ";
    (*p_outstream)<<"\n ";
    for (int j=0;j<4;++j) (*p_outstream)<<p_vhep[4*i+j]<<" ";
    (*p_outstream)<<"\n";
  }
}

void HepEvt_Interface::WriteReducedHepEvt(int nhep)
{
  (*p_outstream)<<"  "<<m_evtnumber<<" "<<nhep<<" \n";
  for (int i=0;i<nhep;++i) {
    (*p_outstream)<<i+1<<"  "<<p_isthep[i]<<" "<<p_idhep[i]<<" "<<p_jmohep[2*i]<<" "<<p_jmohep[2*i+1]
	       <<" "<<p_jdahep[2*i]<<" "<<p_jdahep[2*i+1]<<" \n ";
    for (int j=0;j<5;++j) (*p_outstream)<<p_phep[5*i+j]<<" ";
    (*p_outstream)<<"\n ";
    for (int j=0;j<4;++j) (*p_outstream)<<p_vhep[4*i+j]<<" ";
    (*p_outstream)<<"\n";
  }
}

void HepEvt_Interface::WriteFormatedHepEvt(int nhep)
{
  (*p_outstream)<<" "<<std::setw(4)<<nhep<<" \n";
  for (int i=0;i<nhep;++i) {
    (*p_outstream)<<" "<<std::setw(8)<<p_isthep[i]<<" "<<std::setw(8)<<p_idhep[i]
		  <<" "<<std::setw(4)<<p_jmohep[2*i]<<" "<<std::setw(4)<<p_jmohep[2*i+1]
		  <<" "<<std::setw(4)<<p_jdahep[2*i]<<" "<<std::setw(4)<<p_jdahep[2*i+1]<<" \n ";
    (*p_outstream)<<std::setprecision(10);
    (*p_outstream)<<std::setiosflags(std::ios::fixed);
    for (int j=0;j<5;++j) (*p_outstream)<<std::setw(16)<<p_phep[5*i+j]<<" ";
    (*p_outstream)<<"\n ";
    for (int j=0;j<4;++j) (*p_outstream)<<std::setw(16)<<p_vhep[4*i+j]<<" ";
    (*p_outstream)<<"\n";
    (*p_outstream)<<std::resetiosflags(std::ios::fixed);
  }
}

void HepEvt_Interface::ISBlobs2HepEvt(Blob_List * const _blobs,int & _nhep) {
  for (int beam=0;beam<2;beam++) {
    for (Blob_List::const_iterator bit=_blobs->begin(); bit!=_blobs->end();++bit) {
      if ((*bit)->Type()==btp::Bunch && (*bit)->Beam()==beam) {
	if ((*bit)->NInP()!=1) {
 	  msg.Error()<<"Error in HepEvt_Interface::ISBlobs2HepEvt."<<endl
		     <<"   Bunch blob with more than one incoming particle !"<<endl
		     <<(*bit)<<endl;
	  abort();
	}
	if ((*bit)->NOutP()>1) {
	  Particle2HepEvt((*bit)->InParticle(0),_nhep);
	  for (int j=0;j<(*bit)->NOutP();j++) Particle2HepEvt((*bit)->OutParticle(j),_nhep);
	  EstablishRelations((*bit));
	}
      }
      if ((*bit)->Type()==btp::Beam && (*bit)->Beam()==beam) {
	if ((*bit)->NInP()!=1) {
	  msg.Error()<<"Error in HepEvt_Interface::ISBlobs2HepEvt."<<endl
		     <<"   Beam Remnant blob with more than one incoming particle !"<<endl
		     <<(*bit)<<endl;
	  abort();
	}
	if ((*bit)->NOutP()>1) {
	  Particle2HepEvt((*bit)->InParticle(0),_nhep);
	  for (int j=0;j<(*bit)->NOutP();j++) Particle2HepEvt((*bit)->OutParticle(j),_nhep);
	  EstablishRelations((*bit));
	}
      }
      if ((*bit)->Type()==btp::IS_Shower && (*bit)->Beam()==beam) {
	if ((*bit)->NInP()!=1) {
	  msg.Error()<<"Error in HepEvt_Interface::ISBlobs2HepEvt."<<endl
		     <<"   IS blob with more than one incoming particle !"<<endl
		     <<(*bit)<<endl;
	  abort();
	}
	Particle2HepEvt((*bit)->InParticle(0),_nhep);
	for (int j=0;j<(*bit)->NOutP();j++) Particle2HepEvt((*bit)->OutParticle(j),_nhep);
	EstablishRelations((*bit));
      }
    }
  }
}


void HepEvt_Interface::HardBlob2HepEvt(Blob_List * const _blobs,int & _nhep) {
  int mo,da;
  for (Blob_List::const_iterator bit=_blobs->begin(); bit!=_blobs->end();++bit) {
    if ((*bit)->Type()==btp::ME_PS_Interface_IS) {
      if ((*bit)->NInP()!=2 || (*bit)->NOutP()!=2) {
	msg.Error()<<"Error in HepEvt_Interface::HardBlob2HepEvt."<<endl
		   <<"   ME_PS_Interface_IS blob with other than 2->2 particles !"<<endl
		   <<(*bit)<<endl;
	abort();
      }
      else {
	for (int i=0;i<2;i++) {
	  Particle2HepEvt((*bit)->InParticle(i),_nhep);
	  Particle2HepEvt((*bit)->OutParticle(i),_nhep);
	  mo = m_convertS2H[(*bit)->InParticle(i)];
	  da = m_convertS2H[(*bit)->OutParticle(i)];
	  for (int j=0;j<2;j++) {
	    p_jmohep[2*da+j] = mo+1; p_jdahep[2*mo+j] = da+1;
	  } 
	}
      }
    }
    if ((*bit)->Type()==btp::Signal_Process) {
      if ((*bit)->NInP()!=2) {
	msg.Error()<<"Error in HepEvt_Interface::HardBlob2HepEvt."<<endl
		   <<"   Hard ME blob with other than 2 incoming particles !"<<endl
		   <<(*bit)<<endl;
	abort();
      }
      if ((*bit)->NOutP()>=2) {
	Particle2HepEvt((*bit)->InParticle(0),_nhep);
	Particle2HepEvt((*bit)->InParticle(1),_nhep);
	for (int j=0;j<(*bit)->NOutP();j++) Particle2HepEvt((*bit)->OutParticle(j),_nhep);
	EstablishRelations((*bit));
      }
    }
    if ((*bit)->Type()==btp::ME_PS_Interface_FS) {
      if ((*bit)->NInP()<2 || (*bit)->NOutP()!=(*bit)->NInP()) {
	msg.Error()<<"Error in HepEvt_Interface::HardBlob2HepEvt."<<endl
		   <<"   ME_PS_Interface_IS blob with other than 2->2 particles !"<<endl
		   <<(*bit)<<endl;
	abort();
      }
      else {
	for (int i=0;i<(*bit)->NOutP();i++) {
	  Particle2HepEvt((*bit)->InParticle(i),_nhep);
	  Particle2HepEvt((*bit)->OutParticle(i),_nhep);
	  mo = m_convertS2H[(*bit)->InParticle(i)];
	  da = m_convertS2H[(*bit)->OutParticle(i)];
	  for (int j=0;j<2;j++) {
	    p_jmohep[2*da+j] = mo+1; p_jdahep[2*mo+j] = da+1;
	  } 
	}
      }
    }
  }
}

void HepEvt_Interface::FSBlobs2HepEvt(Blob_List * const _blobs,int & _nhep) {
  for (Blob_List::const_iterator bit=_blobs->begin(); bit!=_blobs->end();++bit) {
    if ((*bit)->Type()==btp::FS_Shower && (*bit)->NInP()==1) {
      for (int j=0;j<(*bit)->NOutP();j++) Particle2HepEvt((*bit)->OutParticle(j),_nhep);
      EstablishRelations((*bit));
    } 
  }
}

void HepEvt_Interface::FragmentationBlob2HepEvt(Blob_List * const _blobs,int & _nhep) {
  for (Blob_List::const_iterator bit=_blobs->begin(); bit!=_blobs->end();++bit) {
    if ((*bit)->Type()==btp::Fragmentation) {
      String2HepEvt((*bit),_nhep);;
    }
  }
}

void HepEvt_Interface::HadronDecayBlobs2HepEvt(Blob_List * const _blobs,int & _nhep) {
  for (Blob_List::const_iterator bit=_blobs->begin(); bit!=_blobs->end();++bit) {
    if ((*bit)->Type()==btp::Hadron_Decay) {
      if ((*bit)->NInP()!=1) {
	msg.Error()<<"Error in HepEvt_Interface::HadronDecays2HepEvt."<<endl
		   <<"   Decay blob with other than 1 incoming particles !"<<endl
		   <<(*bit)<<endl;
	abort();
      }
      if ((*bit)->NOutP()>=2) {
	Particle2HepEvt((*bit)->InParticle(0),_nhep);
	for (int j=0;j<(*bit)->NOutP();j++) Particle2HepEvt((*bit)->OutParticle(j),_nhep);
	EstablishRelations((*bit));
      }      
      else if ((*bit)->NOutP()==1 &&
	       ((*bit)->InParticle(0)->Flav().Kfcode()==311 ||
		(*bit)->InParticle(0)->Flav().Kfcode()==511) ) {
	// KK or BB mixing !!!!
	Particle2HepEvt((*bit)->InParticle(0),_nhep);
	for (int j=0;j<(*bit)->NOutP();j++) Particle2HepEvt((*bit)->OutParticle(j),_nhep);
	EstablishRelations((*bit));
      }
      else {
	msg.Error()<<"Warning : Potential error in HepEvt_Interface::HadronDecays2HepEvt."<<endl
		   <<"   Decay blob for 1 -> 1 process with no identified mxing !"<<std::endl;
      }
    }
  }
}

void HepEvt_Interface::Particle2HepEvt(Particle * const _part,int & _nhep)
{
  int number = m_convertS2H.count(_part);
  if (number>0) return;
  p_idhep[_nhep]    = _part->Flav().HepEvt();
  p_jmohep[2*_nhep] = p_jmohep[2*_nhep+1] = 0; 
  p_jdahep[2*_nhep] = p_jdahep[2*_nhep+1] = 0; 
        
  for (short int j=1; j<4; ++j) p_phep[(j-1)+_nhep*5] = _part->Momentum()[j];
  p_phep[3+_nhep*5] = _part->Momentum()[0];
  double pabs = (_part->Momentum()).Abs2();
  if (pabs<0) p_phep[4+_nhep*5] = 0.;
         else p_phep[4+_nhep*5] = sqrt(pabs);
  if (_part->ProductionBlob()!=NULL) {
    for (short int j=1; j<4; ++j) p_vhep[(j-1)+_nhep*4] = _part->XProd()[j];
    p_vhep[3+_nhep*4]     = _part->XProd()[0];
  }
  else {
    for (short int j=1; j<4; ++j) p_vhep[(j-1)+_nhep*4] = 0.;
    p_vhep[3+_nhep*4]     = 0.;
  }
  if (_part->DecayBlob()!=NULL) p_isthep[_nhep] = 2;
                           else p_isthep[_nhep] = 1;

  m_convertS2H.insert(std::make_pair(_part,_nhep));
  _nhep++;
}

void HepEvt_Interface::String2HepEvt(Blob * const _string,int & _nhep)
{
  p_idhep[_nhep]    = 92;
  p_jmohep[2*_nhep] = p_jmohep[2*_nhep+1] = 0; 
  p_jdahep[2*_nhep] = p_jdahep[2*_nhep+1] = 0; 
  for (short int j=1; j<4; ++j) p_phep[(j-1)+_nhep*5] = _string->CMS()[j];
  p_phep[3+_nhep*5] = _string->CMS()[0];
  double pabs = (_string->CMS()).Abs2();
  if (pabs<0) p_phep[4+_nhep*5] = 0.;
         else p_phep[4+_nhep*5] = sqrt(pabs);
  for (short int j=1;j<4;j++) p_vhep[(j-1)+_nhep*4] = _string->Position()[j];
  p_vhep[3+_nhep*4] = _string->Position()[0];
  p_isthep[_nhep]   = 2;

  Particle * incoming, * outgoing;
  int number, stringnumber = _nhep;
  for (int i=0;i<_string->NInP();i++) {
    incoming = _string->InParticle(i);
    number   = m_convertS2H[incoming];
    if (i==0)                 p_jmohep[2*_nhep]   = number+1;
    if (i==_string->NInP()-1) p_jmohep[2*_nhep+1] = number+1;
    p_jdahep[2*number] = p_jdahep[2*number+1]       = _nhep+1; 
  }
  
  _nhep++;
  for (int i=0;i<_string->NOutP();i++) {
    outgoing = _string->OutParticle(i);
    Particle2HepEvt(outgoing,_nhep);
    if (i==0)                  p_jdahep[2*stringnumber]   = _nhep;
    if (i==_string->NOutP()-1) p_jdahep[2*stringnumber+1] = _nhep;
    p_jmohep[2*m_convertS2H[outgoing]] = p_jmohep[2*m_convertS2H[outgoing]+1] = stringnumber+1;
  }
}


void HepEvt_Interface::EstablishRelations(Blob * const _blob) {
  int mothers[2];
  int daughters[2];
  mothers[0]  = mothers[1]   = 0;
  for (int i=0;i<_blob->NInP();i++) {
    mothers[i] = m_convertS2H[_blob->InParticle(i)];
  }
  if (_blob->NOutP()>0) {
    daughters[0] = m_convertS2H[_blob->OutParticle(0)];
    if (_blob->NOutP()>1) daughters[1] = m_convertS2H[_blob->OutParticle(_blob->NOutP()-1)];
    else daughters[1] = m_convertS2H[_blob->OutParticle(0)];
  }
  else daughters[0] = daughters[1] = 0;
  
  if (_blob->NInP()>0) {
    for (int i=0;i<_blob->NInP();i++) {
      p_jdahep[2*mothers[i]]   = daughters[0]+1;
      p_jdahep[2*mothers[i]+1] = daughters[1]+1;
    }
  }
  if (_blob->NOutP()>0) {
    for (int i=0;i<_blob->NOutP();i++) {
      for (int j=0;j<_blob->NInP();j++){ 
	p_jmohep[2*m_convertS2H[_blob->OutParticle(i)]+j] = mothers[j]+1; 
      }
    }
  }
}

void HepEvt_Interface::EstablishRelations(Particle * const _mother,
					  Blob * const _blob) {
  int mother;
  int daughters[2];
  mother = m_convertS2H[_mother];
  if (_blob->NOutP()>0) {
    daughters[0] = m_convertS2H[_blob->OutParticle(0)];
    if (_blob->NOutP()>1) daughters[1] = m_convertS2H[_blob->OutParticle(_blob->NOutP()-1)];
    else daughters[1] = m_convertS2H[_blob->OutParticle(0)];
  }
  else daughters[0]   = daughters[1] = 0;
  
  p_jdahep[2*mother]    = daughters[0]+1;
  p_jdahep[2*mother+1]  = daughters[1]+1;
  if (_blob->NOutP()>0) {
    for (int i=0;i<_blob->NOutP();i++) {
      p_jmohep[2*m_convertS2H[_blob->OutParticle(i)]]   = mother+1;
      p_jmohep[2*m_convertS2H[_blob->OutParticle(i)]+1] = 0;
    }
  }
}



/*------------------------------------------------------------------------------------
  HepEvt to Sherpa methods
------------------------------------------------------------------------------------*/

bool HepEvt_Interface::HepEvt2Sherpa(Blob_List * const blobs) {
  bool okay;
  if (!m_convertH2S.empty()) {
    for (Translation_Map::iterator piter=m_convertH2S.begin();
	 piter!=m_convertH2S.end();piter++) {
      if (piter->second.second) {
	delete (piter->second.first); piter->second.first=NULL; 
      }
    }
    m_convertH2S.clear();
  }

  if (p_instream) ReadHepEvt(blobs);
  else { for (int i=0;i<m_nhep;++i) HepEvt2Particle(i); }
  switch (m_generator)  {
    case gtp::Herwig:  okay = ConstructBlobsFromHerwig(blobs); break;
    case gtp::Pythia:  okay = ConstructBlobsFromPythia(blobs); break;
    case gtp::Sherpa:  okay = ConstructBlobs(blobs); break;
    default:
      msg.Error()<<"Error in HepEvt_Interface::ReadHepEvt."<<std::endl
		 <<"   Generator type unspecified : "<<m_generator<<std::endl
		 <<"   Will abort the run."<<std::endl;
      abort();
      
  }
  okay = okay && IdentifyBlobs(blobs);
  m_evtcount++;
  if (p_instream && m_evtcount%m_filesize==0) OpenNewHepEvtFile();
  return okay; 
}

void HepEvt_Interface::ReadHepEvt(Blob_List * const blobs) 
{
  int number = 0;
  if (p_instream) {
    *p_instream>>m_evtnumber>>m_nhep;
    for (int i=0;i<m_nhep;++i) {
      (*p_instream)>>number>>p_isthep[i]>>p_idhep[i]>>p_jmohep[2*i]>>p_jmohep[2*i+1]
		   >>p_jdahep[2*i]>>p_jdahep[2*i+1];
      (*p_instream)>>p_phep[5*i+0]>>p_phep[5*i+1]>>p_phep[5*i+2]>>p_phep[5*i+3]>>p_phep[5*i+4];
      (*p_instream)>>p_vhep[4*i+0]>>p_vhep[4*i+1]>>p_vhep[4*i+2]>>p_vhep[4*i+3];
      HepEvt2Particle(i);
    }
  }
}

void HepEvt_Interface::OpenNewHepEvtFile() 
{
  std::string file, filename;
  (*p_instream)>>file;
  filename =  m_path+std::string("/")+file; 
  p_instream->close();
  delete p_instream;
  p_instream = new std::ifstream(filename.c_str()); 
  if (!p_instream->good()) {
    msg.Error()<<"ERROR in HepEvt_Interface."<<std::endl
	       <<"   Event file "<<filename<<" not found."<<std::endl
	       <<"   Will abort the run."<<std::endl;
    abort();
  }
  std::string gentype;
  (*p_instream)>>gentype>>m_filesize;
  if (gentype==std::string("Sherpa")) m_generator = gtp::Sherpa;
  if (gentype==std::string("Herwig")) m_generator = gtp::Herwig;
  if (gentype==std::string("Pythia")) m_generator = gtp::Pythia;
  m_evtcount=0;
}

//-----------------------------------------------------------------------------

void HepEvt_Interface::HepEvt2Particle(const int pos)
{
  /*
    std::cout<<pos<<": stat,id "<<p_isthep[pos]<<","<<p_idhep[pos]
    <<"; mos : "<<p_jmohep[2*pos]<<","<<p_jmohep[2*pos+1]
    <<"; das : "<<p_jdahep[2*pos]<<","<<p_jdahep[2*pos+1]<<std::endl
    <<"    mom: "<<p_phep[5*pos+3]<<" "<<p_phep[5*pos+0]<<" "<<p_phep[5*pos+1]
    <<" "<<p_phep[5*pos+2]<<" "<<p_phep[5*pos+4]<<std::endl;
  */
  Flavour flav;
  if ((m_generator==gtp::Herwig) &&
      (p_idhep[pos]==94 || p_idhep[pos]==0)) flav=Flavour(kf::none);
  else flav.FromHepEvt(p_idhep[pos]);
  Vec4D momentum     = Vec4D(p_phep[3+pos*5],p_phep[0+pos*5],p_phep[1+pos*5],p_phep[2+pos*5]);
  Particle * newpart = new Particle(pos+1,flav,momentum);
  newpart->SetStatus(p_isthep[pos]);
  m_convertH2S.insert(std::make_pair(pos,std::make_pair(newpart,true)));
}

bool HepEvt_Interface::ConstructBlobsFromHerwig(ATOOLS::Blob_List * const blobs)
{
  bool signalwarner=true,breakit=false,first_fsr=false;
  int helper, helper1;
  ATOOLS::Particle * part, * mother;
  ATOOLS::Blob * blob, * help, *fsr,
    * signal = new ATOOLS::Blob(), 
    * beam1  = new ATOOLS::Blob(), 
    * beam2  = new ATOOLS::Blob(), 
    * isr1   = new ATOOLS::Blob(), 
    * isr2   = new ATOOLS::Blob(), 
    * fsr1   = new ATOOLS::Blob(), 
    * fsr2   = new ATOOLS::Blob(),
    * cf     = new ATOOLS::Blob();
  signal->SetType(btp::Signal_Process);
  beam1->SetType(btp::Beam);
  beam2->SetType(btp::Beam);
  isr1->SetType(btp::IS_Shower);
  isr2->SetType(btp::IS_Shower);
  fsr1->SetType(btp::FS_Shower);
  fsr2->SetType(btp::FS_Shower);
  cf->SetType(btp::Cluster_Formation);

  blobs->push_back(beam1);  beam1->SetId();
  blobs->push_back(beam2);  beam2->SetId();
  blobs->push_back(isr1);   isr1->SetId();
  blobs->push_back(isr2);   isr2->SetId();
  blobs->push_back(signal); signal->SetId();
  blobs->push_back(fsr1);   fsr1->SetId();
  blobs->push_back(fsr2);   fsr2->SetId();
  blobs->push_back(cf);     cf->SetId();

  Translation_Map::iterator piter, miter;
  for (int i=0;i<m_nhep;i++) {
    piter = m_convertH2S.find(i);
    if (piter==m_convertH2S.end() || !piter->second.second) continue;
    part = piter->second.first;
    if (part->Status()!=2 && part->Status()!=1) part->SetStatus(3);
    switch (p_isthep[i]) {
    case 101:
      beam1->AddToInParticles(part);
      part->SetStatus(2);
      piter->second.second = false;
      break;
    case 102:
      beam2->AddToInParticles(part);
      part->SetStatus(2);
      piter->second.second = false;
      break;
    case 121:
      isr1->AddToOutParticles(part);
      signal->AddToInParticles(part);
      part->SetStatus(2);
      part->SetInfo('G');
      piter->second.second = false;
      break;
    case 122:
      isr2->AddToOutParticles(part);
      signal->AddToInParticles(part);
      part->SetStatus(2);
      part->SetInfo('G');
      piter->second.second = false;
      break;
    case 123:
      first_fsr = true;
    case 124:
      if (p_isthep[p_jmohep[2*i]-1]==125 || 
	  p_isthep[p_jmohep[2*i]-1]==195 || p_isthep[p_jmohep[2*i]-1]==155) break;
      if (part->Flav()==Flavour(kf::none) && 
	  p_idhep[p_jdahep[2*i]-1]==94 && p_isthep[p_jdahep[2*i]-1]==144) {
	helper = p_jdahep[2*(p_jdahep[2*(p_jdahep[2*i]-1)]-1)]-1;
	if (p_idhep[helper]==0 && p_isthep[helper]==155) {
	  signalwarner=false;
	  //std::cout<<"Gotcha : "<<p_jdahep[2*i]-1<<" : "<<p_idhep[p_jdahep[2*i]-1]
	  //	   <<" -> "<<helper<<", "<<p_jdahep[2*helper]-1
	  //	   <<" -> "<<p_jdahep[2*(p_jdahep[2*helper]-1)]-1<<std::endl;
	  for (int k=0;k<2;k++) {
	    helper1 = p_jdahep[2*helper+k]-1;
	    piter = m_convertH2S.find(helper1);
	    if (piter==m_convertH2S.end() || !piter->second.second) continue;
	    part = piter->second.first;
	    piter->second.second = false;
	    part->SetStatus(2);
	    part->SetInfo('H');
	    signal->AddToOutParticles(part);
	    if (k==0) { 
	      //std::cout<<"Add "<<helper1+1<<" to "<<first_fsr<<std::endl;
	      if (first_fsr) fsr1->AddToInParticles(part);
	                else fsr2->AddToInParticles(part);
	    }
	    else {
	      fsr = new ATOOLS::Blob();
	      blobs->push_back(fsr);
	      fsr->SetType(btp::FS_Shower);
	      fsr->SetId();
	      fsr->AddToInParticles(part);
	    }
	    if (p_jdahep[2*helper1]>p_jdahep[2*helper1+1]) {
	      helper1 = p_jdahep[2*helper1]-1;
	      piter   = m_convertH2S.find(helper1);
	      if (piter==m_convertH2S.end() || !piter->second.second) continue;
	      part = piter->second.first;
	      piter->second.second = false;
	      part->SetStatus(2);
	      part->SetInfo('f');
	      if (k==0) { 
		if (first_fsr) fsr1->AddToOutParticles(part);
		          else fsr2->AddToOutParticles(part);
	      }
	      else fsr->AddToOutParticles(part);	    
	      help = new ATOOLS::Blob();
	      blobs->push_back(help);
	      help->SetType(btp::Hard_Decay);
	      help->SetId();
	      help->AddToInParticles(part);
	      for (int j=p_jdahep[2*helper1]-1;j<=p_jdahep[2*helper1+1]-1;j++) {
		//std::cout<<j+1<<" in "<<p_jdahep[2*helper1]<<"..."<<p_jdahep[2*helper1+1]<<std::endl;
		piter = m_convertH2S.find(j);
		if (piter==m_convertH2S.end() || !piter->second.second) continue;
		piter->second.second = false;
		part = piter->second.first;
		help->AddToOutParticles(part);
		part->SetInfo('h');
		if (p_isthep[j]==1) part->SetStatus(1);
		else {
		  blob = new ATOOLS::Blob();
		  blobs->push_back(blob);
		  blob->SetType(btp::FS_Shower);
		  blob->SetId();
		  blob->AddToInParticles(part);
		  part->SetStatus(2);
		  //std::cout<<"---------------"<<j+1<<" -> "<<p_jdahep[2*j]<<" ->"
		  //	   <<p_jdahep[2*(p_jdahep[2*j]-1)]<<","<<p_jdahep[2*(p_jdahep[2*j]-1)+1]
		  //	   <<"---------------------------"<<std::endl;
		  for (int l=p_jdahep[2*(p_jdahep[2*j]-1)];l<=p_jdahep[2*(p_jdahep[2*j]-1)+1];l++) {
		    //std::cout<<"Find this : "<<p_jdahep[2*j]<<" -> "<<l<<std::endl;
		    breakit=true;
		    piter = m_convertH2S.find(l-1);
		    if (piter==m_convertH2S.end() || !piter->second.second) continue;
		    part = piter->second.first;
		    piter->second.second = false;
		    blob->AddToOutParticles(part);
		    cf->AddToInParticles(part);
		    part->SetStatus(2);
		    part->SetInfo('f');
		  }
		}
	      }
	    }
	    //if (k==0) std::cout<<(*signal)<<(*fsr2)<<(*help)<<std::endl;
	    //if (k==1) std::cout<<(*signal)<<(*fsr)<<(*help)<<std::endl;
	  }
	  //if (breakit) {
	  // std::cout<<(*blobs)<<std::endl;
	  // abort();
	  //}
	  first_fsr = false;
	  break;
	}
	else{
	  msg.Error()<<"Error in HepEvt_Interface::ConstructBlobsFromHerwig : "<<std::endl
		     <<"   Unexpected feature in HepEvt, will continue & hope for the best."<<std::endl;
	  first_fsr = false;
	  break;
	}
      }
      //std::cout<<"Check124 "<<i+1<<" -> "<<p_jdahep[2*i]<<" ("<<p_isthep[p_jdahep[2*i]-1]<<")"
      //	       <<" <- "<<p_jmohep[2*i]<<" ("<<p_isthep[p_jmohep[2*i]-1]<<")"<<std::endl;
      signal->AddToOutParticles(part);
      if (first_fsr) fsr1->AddToInParticles(part);
	        else fsr2->AddToInParticles(part);
      part->SetStatus(2);
      piter->second.second = false;
      part->SetInfo('H');
      if (p_isthep[p_jdahep[2*i]-1]==195 || p_isthep[p_jdahep[2*i]-1]==3) {
	//helper = i;
	helper = p_jdahep[2*i]-1;
	//std::cout<<"Try 1 : "<<i+1<<"-> "<<helper+1<<std::endl;
	if (p_isthep[p_jdahep[2*i]-1]==3) helper = (p_jdahep[2*helper]-1);
	//std::cout<<"Try 2 : "<<i+1<<"-> "<<helper+1<<std::endl;
	piter = m_convertH2S.find(helper);
	if (piter==m_convertH2S.end() || !piter->second.second) continue;
	piter->second.second = false;
	part = piter->second.first;
	part->SetStatus(2);
	part->SetInfo('f');
	if (first_fsr) fsr1->AddToOutParticles(part);
	          else fsr2->AddToOutParticles(part);
	help = new ATOOLS::Blob();
	blobs->push_back(help);
	help->SetType(btp::Hard_Decay);
	help->SetId();
	help->AddToInParticles(part);
	for (int j=p_jdahep[2*helper]-1;j<=p_jdahep[2*helper+1]-1;j++) {
	  //std::cout<<helper+1<<" -> "<<j+1<<" into decay blob."<<std::endl;
	  piter = m_convertH2S.find(j);
	  if (piter==m_convertH2S.end() || !piter->second.second) continue;
	  part = piter->second.first;
	  help->AddToOutParticles(part);
	  piter->second.second = false;
	  part->SetInfo('h');
	  if (p_isthep[j]==1) part->SetStatus(1);
	  else {
	    blob = new ATOOLS::Blob();
	    blobs->push_back(blob);
	    blob->SetType(btp::FS_Shower);
	    blob->SetId();
	    blob->AddToInParticles(part);
	    part->SetStatus(2);
	    //std::cout<<"---------------"<<j+1<<" -> "<<p_jdahep[2*j]<<" ->"
	    //	     <<p_jdahep[2*(p_jdahep[2*j]-1)]<<","<<p_jdahep[2*(p_jdahep[2*j]-1)+1]
	    //	     <<"---------------------------"<<std::endl;
	    for (int l=p_jdahep[2*(p_jdahep[2*j]-1)];l<=p_jdahep[2*(p_jdahep[2*j]-1)+1];l++) {
	      //std::cout<<"Find this : "<<p_jdahep[2*j]<<" -> "<<l<<std::endl;
	      breakit=true;
	      piter = m_convertH2S.find(l-1);
	      if (piter==m_convertH2S.end() || !piter->second.second) continue;
	      part = piter->second.first;
	      piter->second.second = false;
	      blob->AddToOutParticles(part);
	      cf->AddToInParticles(part);
	      part->SetStatus(2);
	      part->SetInfo('f');
	    }
	  }
	}
	//std::cout<<(*blobs)<<std::endl<<std::endl<<std::endl<<std::endl;
	//abort();
      }
      first_fsr = false;
      break;
    case 125:
    case 155:
      break;
    case 141:
      beam1->AddToOutParticles(part);
      isr1->AddToInParticles(part);
      part->SetStatus(2);
      part->SetInfo('I');
      piter->second.second = false;
      break;
    case 142:
      beam2->AddToOutParticles(part);
      isr2->AddToInParticles(part);
      part->SetStatus(2);
      part->SetInfo('I');
      piter->second.second = false;
      break;
    case 143:
      // Outgoing jet 
      break;
    case 144:
      // Outgoing jet 
      break;
    case 160:
      break;
    case 158:
    case 159:
    case 161:
    case 162:
      break;
    case 181:
    case 182:
    case 183:
    case 184:
    case 185:
    case 186:
      part->SetStatus(2);
      part->SetInfo('C');
      cf->AddToOutParticles(part);
      if (p_jdahep[2*i]!=0) {
	blob = new ATOOLS::Blob();
	blobs->push_back(blob);
	blob->SetType(btp::Cluster_Decay);
	blob->SetId();
	blob->AddToInParticles(part);
      }
      piter->second.second = false;
      break;
    case 195:
    case 196:
    case 197:
      part->SetStatus(1);
      part->SetInfo('P');
      miter = m_convertH2S.find(p_jmohep[2*i]-1);
      if (miter==m_convertH2S.end()) continue;
      mother = miter->second.first;
      blob   = mother->DecayBlob();
      if (!blob) {
	blob = new ATOOLS::Blob();
	blobs->push_back(blob);
	blob->SetType(btp::Cluster_Decay);
	blob->SetId();
	blob->AddToInParticles(mother);
      }
      blob->AddToOutParticles(part);
      piter->second.second = false;
      blob->SetPosition(Vec4D(p_vhep[4*i+3],p_vhep[4*i+0],p_vhep[4*i+1],p_vhep[4*i+2]));
      if (p_jdahep[2*i]!=0) {
	part->SetStatus(2);
	part->SetInfo('p');
	blob = new ATOOLS::Blob();
	blobs->push_back(blob);
	blob->SetType(btp::Hadron_Decay);
	blob->SetId();
	blob->AddToInParticles(part);
      }
      break;
    case 198:
      part->SetStatus(1);
      part->SetInfo('D');
      miter = m_convertH2S.find(p_jmohep[2*i]-1);
      if (miter==m_convertH2S.end()) continue;
      mother = miter->second.first;
      blob   = mother->DecayBlob();
      if (!blob) {
	blob = new ATOOLS::Blob();
	blobs->push_back(blob);
	blob->SetType(btp::Hadron_Decay);
	blob->SetId();
	blob->AddToInParticles(mother);
      }
      blob->AddToOutParticles(part);
      piter->second.second = false;
      blob->SetPosition(Vec4D(p_vhep[4*i+3],p_vhep[4*i+0],p_vhep[4*i+1],p_vhep[4*i+2]));
      if (p_jdahep[2*i]!=0) {
	part->SetStatus(2);
	part->SetInfo('d');
	blob = new ATOOLS::Blob();
	blobs->push_back(blob);
	blob->SetType(btp::Hadron_Decay);
	blob->SetId();
	blob->AddToInParticles(part);
      }
      break;
    case 199:
      part->SetStatus(1);
      part->SetInfo('D');
      miter = m_convertH2S.find(p_jmohep[2*i]-1);
      if (miter==m_convertH2S.end()) continue;
      mother = miter->second.first;
      blob   = mother->DecayBlob();
      if (!blob) {
	blob = new ATOOLS::Blob();
	blobs->push_back(blob);
	blob->SetType(btp::Hadron_Decay);
	blob->SetId();
	blob->AddToInParticles(mother);
      }
      blob->AddToOutParticles(part);
      piter->second.second = false;
      blob->SetPosition(Vec4D(p_vhep[4*i+3],p_vhep[4*i+0],p_vhep[4*i+1],p_vhep[4*i+2]));
      if (p_jdahep[2*i]!=0) {
	part->SetStatus(2);
	part->SetInfo('d');
	blob = new ATOOLS::Blob();
	blobs->push_back(blob);
	blob->SetType(btp::Hadron_To_Parton);
	blob->SetId();
	blob->AddToInParticles(part);
	help = new ATOOLS::Blob();
	blobs->push_back(help);
	help->SetType(btp::Cluster_Formation);
	help->SetId();
	FollowDaughters(i,blob,help);
      }
      else {
	msg.Error()<<"Error in HepEvt_Interface::ConstructBlobsFromHerwig."<<std::endl
		   <<"   Mother of heavy hadron has no decay blob yet."<<std::endl
		   <<"   Will abort."<<Particle::Counter()<<" / "<<Blob::Counter()<<std::endl;
      }
      break;      
    case 200:
      part->SetStatus(1);
      part->SetInfo('D');
      miter = m_convertH2S.find(p_jmohep[2*i]-1);
      if (miter==m_convertH2S.end()) continue;
      mother = miter->second.first;
      blob   = mother->DecayBlob();
      if (!blob) {
	msg.Error()<<"Error in HepEvt_Interface::ConstructBlobsFromHerwig."<<std::endl
		   <<"   Mother of heavy hadron flavour has no decay blob yet."<<std::endl
		   <<"   Will create a new blob and hope for the best."
		   <<Particle::Counter()<<" / "<<Blob::Counter()<<std::endl;
	blob = new ATOOLS::Blob();
	blobs->push_back(blob);
	blob->SetType(btp::Hadron_Decay);
	blob->SetId();
	blob->AddToInParticles(mother);
      }
      blob->AddToOutParticles(part);
      piter->second.second = false;
      blob->SetPosition(Vec4D(p_vhep[4*i+3],p_vhep[4*i+0],p_vhep[4*i+1],p_vhep[4*i+2]));
      if (p_jdahep[2*i]!=0) {
	part->SetStatus(2);
	part->SetInfo('d');
	blob = new ATOOLS::Blob();
	blobs->push_back(blob);
	blob->SetType(btp::Hadron_Mixing);
	blob->SetId();
	blob->AddToInParticles(part);
      }
      else {
	msg.Error()<<"Error in HepEvt_Interface::ConstructBlobsFromHerwig."<<std::endl
		   <<"   Mother of heavy hadron has no decay blob yet."<<std::endl
		   <<"   Will abort."<<Particle::Counter()<<" / "<<Blob::Counter()<<std::endl;
      }
      break;      
    case 1:
      if (!piter->second.second) continue;
      miter = m_convertH2S.find(p_jmohep[2*i]-1);
      if (miter==m_convertH2S.end()) continue;
      mother = miter->second.first;
      blob   = mother->DecayBlob();
      if (blob) {
	mother->SetStatus(2);
	if (mother->Info()=='P') mother->SetInfo('p');
	else if (mother->Info()=='D') mother->SetInfo('d');
	blob->AddToOutParticles(part);
	piter->second.second = false;
	blob->SetPosition(Vec4D(p_vhep[4*i+3],p_vhep[4*i+0],p_vhep[4*i+1],p_vhep[4*i+2]));
      }
      part->SetStatus(1);
      part->SetInfo('D');
      break;
    case 2:
      if (p_isthep[p_jmohep[2*i]-1]==141) isr1->AddToOutParticles(part);
      if (p_isthep[p_jmohep[2*i]-1]==142) isr2->AddToOutParticles(part);
      if (p_isthep[p_jmohep[2*i]-1]==143) fsr1->AddToOutParticles(part);
      if (p_isthep[p_jmohep[2*i]-1]==144) fsr2->AddToOutParticles(part);
      cf->AddToInParticles(part); 
      part->SetStatus(2); 
      piter->second.second = false; 
    default : break;
    }
  }
  if (signalwarner) {
    if (signal->NOutP()!=2 || signal->NInP()<2) {
      msg.Error()<<"Error in HepEvt_Interface::ConstructBlobsFromHerwig"<<std::endl
		 <<"   Signal is funny: "<<signal->NInP()<<" -> "<<signal->NOutP()<<std::endl
		 <<"   ====================================="<<std::endl
		 <<(*signal)<<std::endl
		 <<"   ====================================="<<std::endl
		 <<"   Clear blobs, return false and hope for the best."<<Particle::Counter()<<" / "<<Blob::Counter()<<std::endl
		 <<(*blobs)<<std::endl;

      if (!blobs->empty()) {
	for (Blob_List::iterator blit=blobs->begin();blit!=blobs->end();++blit) delete (*blit);
	blobs->clear();
      }
      DeleteObsolete(1);
      return false;
    }
  }
  DeleteObsolete(-1);
  return true;
}

void HepEvt_Interface::FollowDaughters(const int & i,ATOOLS::Blob * hadron,ATOOLS::Blob * clusters) 
{
  Translation_Map::iterator piter;
  Particle * part;
  int end = Max(p_jdahep[2*i]-1,p_jdahep[2*i+1]-1);
  for (int j=p_jdahep[2*i]-1;j<=end;j++) {
    if (j>m_nhep || j<0) {
      msg.Error()<<"Error in HepEvt_Interface::FollowDaughters("<<i<<")"<<std::endl
		 <<"   counter "<<j<<" out of bounds ("<<m_nhep<<")."<<std::endl
		 <<"   Will abort the run."<<std::endl;
      abort();
    }
    if (p_isthep[j]==183) {
      // Final cluster
      piter = m_convertH2S.find(j);
      if (piter==m_convertH2S.end()) {
	msg.Error()<<"Error in HepEvt_Interface::FollowDaughters."<<std::endl
		   <<"   Mix up of secondary clusters."<<std::endl
		   <<"   Will abort."<<std::endl;
	abort();
      }
      if (piter->second.second) {
	piter->second.second = false;
	part = piter->second.first;
	part->SetStatus(2); 
	part->SetInfo('C');
	clusters->AddToOutParticles(part);
      } 
      piter = m_convertH2S.find(i);
      if (piter==m_convertH2S.end()) {
	msg.Error()<<"Error in HepEvt_Interface::FollowDaughters."<<std::endl
		   <<"   Mix up of secondary clusters."<<std::endl
		   <<"   Will abort."<<std::endl;
	abort();
      }
      if (piter->second.second) {
	piter->second.second = false;
	part = piter->second.first;
	part->SetStatus(2); 
	part->SetInfo('f');
	hadron->AddToOutParticles(part);
	clusters->AddToInParticles(part);
      } 
    }
    else if (p_isthep[j]==1) {
      // Final lepton
      piter = m_convertH2S.find(j);
      if (piter==m_convertH2S.end()) {
	msg.Error()<<"Error in HepEvt_Interface::FollowDaughters."<<std::endl
		   <<"   Mix up of secondary clusters."<<std::endl
		   <<"   Will abort."<<std::endl;
	abort();
      }
      if (piter->second.second) {
	piter->second.second = false;
	part = piter->second.first;
	part->SetStatus(1); 
	part->SetInfo('F');
	hadron->AddToOutParticles(part);
      } 
    } 
    else if (p_isthep[j]==160) {
      // Final spectator
      piter = m_convertH2S.find(j);
      if (piter==m_convertH2S.end()) {
	msg.Error()<<"Error in HepEvt_Interface::FollowDaughters."<<std::endl
		   <<"   Mix up of secondary clusters."<<std::endl
		   <<"   Will abort."<<std::endl;
	abort();
      }
      if (piter->second.second) {
	piter->second.second = false;
	part = piter->second.first;
	part->SetStatus(2); 
	part->SetInfo('S');
	hadron->AddToOutParticles(part);
	clusters->AddToInParticles(part);
      } 
    } 
    else {
      piter = m_convertH2S.find(j);
      if (piter==m_convertH2S.end()) {
	msg.Error()<<"Error in HepEvt_Interface::FollowDaughters."<<std::endl
		   <<"   Mix up of secondary clusters."<<std::endl
		   <<"   Will abort."<<std::endl;
	abort();
      }
      p_isthep[j]=155;
      FollowDaughters(j,hadron,clusters);
    }
  }
}

bool HepEvt_Interface::ConstructBlobsFromPythia(ATOOLS::Blob_List * const blobs)
{
  ATOOLS::Blob     * signal, * blob, * productionblob, * decayblob;
  ATOOLS::Particle * part, * mother, * daughter, * partner;

  std::vector<Blob *>     _blobs;
  std::set<int>           signalints, intermeds, outs, ins;
  std::set<int>::iterator sit;
  std::vector<int>        ueints;

  Translation_Map::iterator piter, miter, partiter, diter;
  for (int i=0;i<m_nhep;++i) {
    piter = m_convertH2S.find(i);
    if (piter==m_convertH2S.end()) continue;
    part = piter->second.first;
    //std::cout<<"Particle : "<<i<<" : "<<part->Flav()<<" "<<part->Status()<<std::endl;
    if (part->Status()==3) signalints.insert(i);
    if (part->Status()!=3 && (part->Flav().IsGluon() || part->Flav().IsQuark()) &&
	(p_jmohep[2*i+0]==0 && p_jmohep[2*i+1]==0)) {
      ueints.push_back(i);
    }
  }
  bool create,test,inout;
  int  intermed, pint, testint, lsize;
  for (sit=signalints.begin();sit!=signalints.end();sit++) {
    piter = m_convertH2S.find((*sit));
    if (piter==m_convertH2S.end()) continue;
    part = piter->second.first;
    if (IsZero(part->Momentum()[0])) {
      msg.Error()<<"WARNING : Error in HepEvt_Interface::ConstructBlobsFromPythia."<<std::endl
		 <<"    Signal particles with zero energy: Looks like a nonsensical event."<<std::endl
		 <<"    Will return .false. and hope that event is discarded."<<std::endl;
      return false;
    }
  }
  
  if (signalints.size()>0) {
    signal = new Blob;
    signal->ClearAllData();
    signal->AddData("ME_Weight",new Blob_Data<double>(1.));
    signal->SetType(btp::Signal_Process);
    _blobs.push_back(signal);
    for (sit=signalints.begin();sit!=signalints.end();sit++) {
      //piter = m_convertH2S.find((*sit));
      /*std::cout<<"Analyse : "<<(*sit)<<" : "
	       <<p_jmohep[2*(*sit)]<<" "<<p_jmohep[2*(*sit)+1]<<" "
	       <<p_jdahep[2*(*sit)]<<" "<<p_jdahep[2*(*sit)+1]<<std::endl;*/
      if (p_jmohep[2*(*sit)]!=0 && p_jmohep[2*(*sit)+1]!=0) { 
	intermeds.insert((*sit)); 
	ins.insert(p_jmohep[2*(*sit)]-1);
	ins.insert(p_jmohep[2*(*sit)+1]-1);
      }
      else {
	if (intermeds.find(p_jmohep[2*(*sit)]-1)!=intermeds.end()) {
	  outs.insert((*sit)); 
	}
      }
    }
    for (sit=intermeds.begin();sit!=intermeds.end();sit++) signalints.erase((*sit));
    for (sit=ins.begin();sit!=ins.end();sit++)             signalints.erase((*sit));
    for (sit=outs.begin();sit!=outs.end();sit++)           signalints.erase((*sit));

    /*
      std::cout<<"Signalins : ";
      for (sit=signalints.begin();sit!=signalints.end();sit++) std::cout<<(*sit)<<" ";
      std::cout<<std::endl;
      std::cout<<"Ins : ";
      for (sit=ins.begin();sit!=ins.end();sit++) std::cout<<(*sit)<<" ";
      std::cout<<std::endl;
      std::cout<<"Testints : ";
      for (sit=intermeds.begin();sit!=intermeds.end();sit++) std::cout<<(*sit)<<" ";
      std::cout<<std::endl;
      std::cout<<"Outs : ";
      for (sit=outs.begin();sit!=outs.end();sit++) std::cout<<(*sit)<<" ";
      std::cout<<std::endl;
    */
    for (sit=outs.begin();sit!=outs.end();sit++) {
      piter = m_convertH2S.find((*sit));
      if (piter==m_convertH2S.end()) continue;
      signal->AddToOutParticles(piter->second.first);
      piter->second.second = false;
      piter = m_convertH2S.find(p_jmohep[2*(*sit)]-1);
      if (piter==m_convertH2S.end()) continue;
      piter->second.second = false;
    }
    std::vector<int> transfer;
    for (sit=intermeds.begin();sit!=intermeds.end();sit++) {
      piter = m_convertH2S.find((*sit));
      if (piter==m_convertH2S.end()) continue;
      if (piter->second.second) {
	signal->AddToOutParticles(piter->second.first);
	transfer.push_back((*sit));
      }
      piter->second.second = false;
    }
    for (int i=0;i<transfer.size();i++) {
      intermeds.erase(transfer[i]);
      outs.insert(transfer[i]);
    }
    for (sit=ins.begin();sit!=ins.end();sit++) {
      piter = m_convertH2S.find((*sit));
      if (piter==m_convertH2S.end()) continue;
      signal->AddToInParticles(piter->second.first);
      piter->second.second = false;
      pint = (*sit);
      do {
	test = false;
	blob = new Blob;
	_blobs.push_back(blob);
	blob->AddToOutParticles(piter->second.first);
	//std::cout<<"Test this : "<<pint<<" <- "<<p_jmohep[2*pint]-1<<std::endl;
	pint = p_jmohep[2*pint]-1;
	if (pint>-1) {
	  if (signalints.find(pint)!=signalints.end()) {
	    pint  = (*(signalints.find(pint)));
	    //std::cout<<pint<<" <- "<<pint<<std::endl;
	    piter = m_convertH2S.find(pint);
	    if (piter==m_convertH2S.end()) break;
	    blob->AddToInParticles(piter->second.first);
	    piter->second.second = false;
	    test = true;
	  }
	}
      } while(test);
    }

    for (;;) {
      blobs->push_back(_blobs.back());
      _blobs.back()->SetId();
      _blobs.pop_back();
      if (_blobs.empty()) break;
    }

    //std::cout<<(*blobs)<<std::endl;
  }    


  Blob * beam1=NULL, * beam2=NULL;
  for (Blob_List::iterator biter=blobs->begin();biter!=blobs->end();biter++) {
    if ((*biter)->NInP()==0 && (*biter)->NOutP()==1) {
      Particle * part = (*biter)->OutParticle(0);
      if (part!=NULL && part->Flav().IsHadron()) {
	if (part->Momentum()[3]>0) {
	  if (beam1==NULL) beam1 = part->DecayBlob();
	  else {
	    msg.Error()<<"WARNING : Error in HepEvt_Interface::ConstructBlobsFromPythia."<<std::endl
		       <<"    too many transition blobs (1) found for hadron->partons for u.e.."<<std::endl
		       <<"    Will return .false. and hope that event is discarded."<<std::endl;
	    return false;
	  }
	}
	else {
	  if (beam2==NULL) beam2 = part->DecayBlob();
	  else {
	    msg.Error()<<"WARNING : Error in HepEvt_Interface::ConstructBlobsFromPythia."<<std::endl
		       <<"    too many transition blobs (2) found for hadron->partons for u.e.."<<std::endl
		       <<"    Will return .false. and hope that event is discarded."<<std::endl;
	    return false;
	  }
	}
      }
      /*
      if (part!=NULL && (part->Flav()==rpa.gen.Bunch(0)||part->Flav()==rpa.gen.Bunch(1))) {
	if (part->Momentum()[3]>0) {
	  if (beam1==NULL) beam1 = part->DecayBlob();
	  else {
	    msg.Error()<<"WARNING : Error in HepEvt_Interface::ConstructBlobsFromPythia."<<std::endl
		       <<"    too many transition blobs (1) found for incoming beam1."<<std::endl
		       <<(**biter)<<std::endl<<(*part)<<std::endl<<(*beam1)
		       <<"    Will return .false. and hope that event is discarded."<<std::endl;
	    abort();
	    return false;
	  }
	}
	else {
	  if (beam2==NULL) beam2 = part->DecayBlob();
	  else {
	    msg.Error()<<"WARNING : Error in HepEvt_Interface::ConstructBlobsFromPythia."<<std::endl
		       <<"    too many transition blobs (2) found for incoming beam2."<<std::endl
		       <<"    Will return .false. and hope that event is discarded."<<std::endl;
	    return false;
	  }
	}
      }
      */
    }
  }
  if (beam1==NULL || beam2==NULL) {
    msg.Error()<<"WARNING : Error in HepEvt_Interface::ConstructBlobsFromPythia."<<std::endl
	       <<"    too little transition blobs (1/2) found for hadron->partons for u.e. "
	       <<"("<<(beam1==NULL)<<"/"<<(beam2==NULL)<<")"<<std::endl
	       <<"    Will return .false. and hope that event is discarded."<<std::endl;
    return false;
  }
  int ueadd = m_nhep;
  double pt2, E, pl, omega1,omega2;
  lsize = ueints.size();
  while(lsize>0) {
    test = false;
    blob = new ATOOLS::Blob();
    blob->SetType(btp::Hard_Collision);
    // Find pairs of particles.
    pint = ueints[lsize-1];
    piter = m_convertH2S.find(pint);
    if (piter==m_convertH2S.end()) continue;
    part = piter->second.first;
    pt2  = part->Momentum().PPerp2();
    ueints.pop_back();
    for (int i=0;i<lsize;i++) {
      partiter = m_convertH2S.find(ueints[i]);
      if (partiter==m_convertH2S.end()) continue;
      partner = partiter->second.first;
      if (IsEqual(partner->Momentum().PPerp2(),pt2)) {
	blob->AddToOutParticles(part);
	blob->AddToOutParticles(partner);
	piter->second.second = false;
	partiter->second.second = false;
	blobs->push_back(blob);
	E      = part->Momentum()[0]+partner->Momentum()[0];
	pl     = part->Momentum()[3]+partner->Momentum()[3];
	omega1 = (E+pl)/2.;
	omega2 = (E-pl)/2.;
	mother = new Particle(ueadd++,Flavour(kf::gluon),Vec4D(omega1,0.,0.,omega1));
	mother->SetStatus(2);
	blob->AddToInParticles(mother);
	beam1->AddToOutParticles(mother);
	mother = new Particle(ueadd++,Flavour(kf::gluon),Vec4D(omega2,0.,0.,-omega2));
	mother->SetStatus(2);
	blob->AddToInParticles(mother);
	beam2->AddToOutParticles(mother);
	for (int j=i;j<lsize-1;j++) {
	  ueints[j] = ueints[j+1];
	  ueints.pop_back();
	}
	lsize-=2;
	test = true;
	break;
      }  
    }
    if (!test) {
      delete blob;
      msg.Error()<<"WARNING : Error in HepEvt_Interface::ConstructBlobsFromPythia."<<std::endl
		 <<"    no partner with identical pt found for particle from u.e.."<<std::endl
		 <<"    Will return .false. and hope that event is discarded."<<std::endl;
      return false;
    }
  }

  //std::cout<<"Before decays, showers, etc."<<std::endl;
  for (int i=0;i<m_nhep;++i) {
    piter = m_convertH2S.find(i);
    if (piter==m_convertH2S.end()) continue;
    part = piter->second.first;
    if (part->Status()!=3) {
      productionblob = part->ProductionBlob();
      pint           = p_jmohep[2*i]-1;
      miter          = m_convertH2S.find(pint);
      if (miter==m_convertH2S.end()) continue;
      mother = miter->second.first;
      //std::cout<<"Check for : "<<i+1<<"("<<part->Flav()<<") <- "<<p_jmohep[2*i]<<std::endl;
      if (intermeds.find(p_jmohep[2*i]-1)!=intermeds.end() && miter->second.second==false) {
	//std::cout<<"Continue."<<std::endl;
	continue;
      }
      if (productionblob==NULL) {
	create=(mother==NULL);
	if (!create) create=create||(mother->DecayBlob()==NULL); 
	if (create) {
	  productionblob = new ATOOLS::Blob();
	  productionblob->SetPosition(ATOOLS::Vec4D(p_vhep[4*i+3],p_vhep[4*i+0],
						    p_vhep[4*i+1],p_vhep[4*i+2]));
	  productionblob->AddToOutParticles(part);
	  piter->second.second = false;
	  if (mother!=NULL) {
	    productionblob->AddToInParticles(mother);
	    miter->second.second = false;
	  }
	  productionblob->SetId();
	  blobs->push_back(productionblob);
	}
	else {
	  mother->DecayBlob()->AddToOutParticles(part);
	  piter->second.second = false;
	  productionblob = part->ProductionBlob();
	}
      }
    }
    else part->SetStatus(2);
 
    decayblob = part->DecayBlob();
    test = false;
    for (unsigned int j=0;j<2;++j) {
      if (decayblob==NULL) {
	pint     = p_jdahep[2*i+j]-1;
	diter = m_convertH2S.find(pint);
	if (diter==m_convertH2S.end()) continue;
	daughter = diter->second.first;
	if (daughter!=NULL) {
	  if (daughter->ProductionBlob()!=NULL) {
	    if (!test) {
	      daughter->ProductionBlob()->AddToInParticles(part);
	      piter->second.second = false;
	      test        = true;
	    }
	  }
	  else {
	    decayblob = new ATOOLS::Blob();
	    decayblob->SetPosition(ATOOLS::Vec4D(p_vhep[4*pint+3],p_vhep[4*pint+0],
						 p_vhep[4*pint+1],p_vhep[4*pint+2]));
	    decayblob->AddToOutParticles(daughter);
	    decayblob->AddToInParticles(part);
	    diter->second.second = false;
	    piter->second.second = false;
	    decayblob->SetId();
	    blobs->push_back(decayblob);
	  }
	}
      }
    }
  }

  for (int i=0;i<m_nhep;++i) {
    piter = m_convertH2S.find(i);
    if (piter==m_convertH2S.end()) continue;
    daughter = piter->second.first;
    miter = m_convertH2S.find(p_jmohep[2*i]-1);
    if (miter==m_convertH2S.end()) continue;
    mother = miter->second.first;
    if (mother!=NULL && miter->second.second) {
      if (mother->DecayBlob()==NULL) {
	daughter->ProductionBlob()->AddToInParticles(mother);
	miter->second.second = false;
     }
    }
  }

  while(!intermeds.empty()) {
    pint = (*intermeds.begin());
    piter = m_convertH2S.find(pint);
    if (piter==m_convertH2S.end()) continue;
    delete piter->second.first;
    intermeds.erase(pint);
  };
  return true;
}

bool HepEvt_Interface::ConstructBlobs(ATOOLS::Blob_List * const blobs)
{
  ATOOLS::Blob     * productionblob;
  ATOOLS::Particle * part, * mother, * daughter;

  Translation_Map::iterator piter, miter, diter;
  for (int i=0;i<m_nhep;++i) {
    piter          = m_convertH2S.find(i);
    if (piter==m_convertH2S.end()) continue;
    part           = piter->second.first;
    productionblob = part->ProductionBlob();
    miter          = m_convertH2S.find(i);
    if (miter==m_convertH2S.end()) mother = NULL;
    else mother    = miter->second.first;
    if (productionblob==NULL) {
      bool createblob=(mother==NULL);
      if (!createblob) createblob=createblob||(mother->DecayBlob()==NULL); 
      if (createblob) {
	productionblob = new ATOOLS::Blob();
	productionblob->SetPosition(ATOOLS::Vec4D(p_vhep[4*i+3],p_vhep[4*i+0],
						  p_vhep[4*i+1],p_vhep[4*i+2]));
	productionblob->AddToOutParticles(part);
	piter->second.second = false;
	if (mother!=NULL) {
	  productionblob->AddToInParticles(mother);
	  miter->second.second = false;
	}
	productionblob->SetId();
	blobs->push_back(productionblob);
      }
      else {
	mother->DecayBlob()->AddToOutParticles(part);
	piter->second.second = false;
	productionblob=part->ProductionBlob();
      }
    }
    ATOOLS::Blob *decayblob=part->DecayBlob();
    for (unsigned int j=0;j<2;++j) {
      if (decayblob==NULL) {
	diter    = m_convertH2S.find(p_jdahep[2*i+j]-1);
	if (diter==m_convertH2S.end()) daughter = NULL;
	else daughter = diter->second.first;
	if (daughter!=NULL) {
	  if (daughter->ProductionBlob()!=NULL) {
	    daughter->ProductionBlob()->AddToInParticles(part);
	    piter->second.second = false;
	  }
	}
      }
    }
  }
  for (int i=0;i<m_nhep;++i) {
    piter    = m_convertH2S.find(i);
    if (piter==m_convertH2S.end()) continue;
    part = piter->second.first;
    miter    = m_convertH2S.find(p_jmohep[2*i]-1);
    if (miter==m_convertH2S.end()) mother = NULL;
    else mother = miter->second.first;
    if (mother!=NULL) {
      if (mother->DecayBlob()==NULL) {
	part->ProductionBlob()->AddToInParticles(mother);
	miter->second.second = false;
      }
    }
  }
  return true;
}

bool HepEvt_Interface::IdentifyBlobs(ATOOLS::Blob_List * const blobs) 
{
  int                       counter;
  bool                      test;
  Blob_List::iterator       biter, biter2;
  Particle                * search, * dummy;
  Blob                    * prod, * meps, * dec;
  std::vector<Blob *>       obsoletes;
  std::vector<Particle *>   incomings;
  incomings.clear();

  //std::cout<<(*blobs)<<std::endl<<"---------------------------------------------------------"<<std::endl;

  // Beams and bunches
  for (biter=blobs->begin();biter!=blobs->end();biter++) {
    if ((*biter)->NInP()==0 && (*biter)->NOutP()==1) {
      ATOOLS::Particle *incoming=(*biter)->OutParticle(0);
      (*biter)->RemoveOutParticle(incoming);
      delete *biter;
      blobs->erase(biter);
      ATOOLS::Blob * beam=incoming->DecayBlob();
      if (incoming->Flav().IsHadron()) {
	if (beam->OutParticle(0)->DecayBlob()->OutParticle(0)->Flav().IsHadron()) {
	  beam->SetType(btp::Bunch);
	  beam->OutParticle(0)->DecayBlob()->SetType(btp::Beam);
	}
	else {
	  beam->SetType(btp::Beam);
	}
      }
      else {
	ATOOLS::kf::code in=incoming->Flav().Kfcode();
	for (int i=0;i<beam->NOutP();++i) {
	  ATOOLS::kf::code out=beam->OutParticle(i)->Flav().Kfcode();
	  if (in==kf::e && out==kf::photon) {
	    beam->SetType(btp::Bunch);
	    beam->OutParticle(0)->DecayBlob()->SetType(btp::Beam);
	  }
	  /*
	  if (in==kf::e && out==kf::e) {
	    beam->SetType(btp::Bunch);
	    beam->OutParticle(0)->DecayBlob()->SetType(btp::Beam);
	  }
	  */
	}
      }
    }
  }
  //std::cout<<"Beams and bunches"<<std::endl;

  // (IS) Shower Blobs
  counter = 0;
  for (biter=blobs->begin();biter!=blobs->end();biter++) {
    if ((*biter)->NInP()==1 && (*biter)->Type()==btp::Unspecified) {
      search = (*biter)->InParticle(0);
      if (search->ProductionBlob()->Type()==btp::Beam) {
	(*biter)->SetType(btp::IS_Shower); 
	counter++; 
      }
    }
    if (counter==2) { incomings.clear(); break; }
  }
  //std::cout<<"IS"<<std::endl;

  
  // ME Blob
  for (biter=blobs->begin();biter!=blobs->end();biter++) {
    if ((*biter)->NInP()==2 && (*biter)->Type()==btp::Unspecified) {
      //std::cout<<"Check for signal : "<<std::endl<<((**biter))<<std::endl;
      test = false;
      if (incomings.size()==2) {
	if (((*biter)->InParticle(0)==incomings[0] && 
	     (*biter)->InParticle(1)==incomings[1] ) || 
	    ((*biter)->InParticle(0)==incomings[1] && 
	     (*biter)->InParticle(1)==incomings[0] ) ) {
	  test = true;
	  incomings.clear();
	}
	//std::cout<<"Unspecified 2 incomings: Test = true "<<test<<std::endl;
      }
      else {
	counter = 0;
	meps    = NULL;
	for (int i=0;i<2;i++) {
	  search = (*biter)->InParticle(i);
	  prod   = search->ProductionBlob();
	  if (prod->NInP()==1 && prod->NOutP()==1 &&
	      prod->Type()==btp::Unspecified) {
	    if (prod->InParticle(0)->ProductionBlob()->Type()==btp::IS_Shower) {
	      if (meps==NULL) {
		meps = prod;
		meps->SetType(btp::ME_PS_Interface_IS);
		counter++;
	      }
	      else {
		meps->AddToInParticles(prod->RemoveInParticle(0));
		meps->AddToOutParticles(prod->RemoveOutParticle(0));
		counter++;
		for (biter2=blobs->begin();biter2!=blobs->end();biter2++) {
		  if ((*biter2)==prod) { 
		    delete (*biter2);
		    blobs->erase(biter2); 
		    break; 
		  }
		}
	      }
	    }
	  }
	}
	if (counter==2) test = true;
	//std::cout<<"Unspecified mot 2 incomings: Test = true "<<test<<std::endl;
      }
      if (test) {
	(*biter)->SetType(btp::Signal_Process);
	meps     = NULL;
	for (int i=0;i<(*biter)->NOutP();i++) {
	  search = (*biter)->OutParticle(i);
	  dec    = search->DecayBlob();
	  if (dec->NInP()==1 && dec->NOutP()==1 &&
	      dec->Type()==btp::Unspecified) {
	    dummy = dec->OutParticle(0); 
	    if (dummy->DecayBlob()->Type()==btp::Unspecified) {
	      if (meps==NULL) {
		meps = dec;
		meps->SetType(btp::ME_PS_Interface_FS);
		dummy->DecayBlob()->SetType(btp::FS_Shower);
	      }
	      else {
		meps->AddToInParticles(dec->RemoveInParticle(0));
		dec->RemoveOutParticle(0);
		meps->AddToOutParticles(dummy);
		dummy->DecayBlob()->SetType(btp::FS_Shower);
		for (biter2=blobs->begin();biter2!=blobs->end();biter2++) {
		  if ((*biter2)==dec) { 
		    delete (*biter2);
		    blobs->erase(biter2); 
		    break; 
		  }
		}
	      }
	    }
	  }
	  else if (dec->NInP()==1 && dec->NOutP()>1 &&
		   dec->Type()==btp::Unspecified) dec->SetType(btp::FS_Shower);
	}
	break;
      }
    }
    else if ((*biter)->Type()==btp::Signal_Process) {
      //std::cout<<"Is signal : "<<std::endl<<((**biter))<<std::endl;
      for (int i=0;i<(*biter)->NOutP();i++) {
	search = (*biter)->OutParticle(i);
	dec    = search->DecayBlob();
	if (dec->NInP()==1 && dec->Type()==btp::Unspecified) {
	  if (dec->NOutP()!=1) 
	    dec->SetType(btp::FS_Shower);
	  else {
	    meps = dec->OutParticle(0)->DecayBlob();
	    if (meps==NULL) 
	      dec->SetType(btp::FS_Shower);	    
	    else if (meps->OutParticle(0)->Flav()==Flavour(kf::string))
	      dec->SetType(btp::FS_Shower);	    
	  }
	}
      }
    }
  }
  //std::cout<<"ME"<<std::endl;

  // Fragmentation blob
  Flavour cluster; cluster.FromHepEvt(91);
  Flavour string;  string.FromHepEvt(92);
  for (biter=blobs->begin();biter!=blobs->end();biter++) {
    if ((*biter)->NOutP()==1 && (*biter)->Type()==btp::Unspecified &&
	(*biter)->OutParticle(0)->Flav()==string) {
      search = (*biter)->OutParticle(0);
      for (biter2=blobs->begin();biter2!=blobs->end();biter2++) {
	if ((*biter2)->NInP()==1 &&
	    (*biter2)->Type()==btp::Unspecified &&
	    (*biter2)->InParticle(0)==search &&
	    (*biter2)!=(*biter)) {
	  (*biter)->RemoveOutParticle(0);
	  (*biter2)->RemoveInParticle(0);
	  delete search; search=NULL;
	  for (int i=(*biter2)->NOutP()-1;i>=0;i--) {
	    dummy = (*biter2)->RemoveOutParticle(i);
	    (*biter)->AddToOutParticles(dummy);
	  }
	  (*biter)->SetType(btp::Fragmentation);
	  delete (*biter2);
	  blobs->erase(biter2);
	  break;
	}
      }
    }
  }
  //std::cout<<"Frag"<<std::endl;

  // Hadron blobs
  for (biter=blobs->begin();biter!=blobs->end();biter++) {
    if ((*biter)->NInP()==1 && 
	(*biter)->Type()==btp::Unspecified &&
	(*biter)->InParticle(0)->Flav().IsHadron()) {
      (*biter)->SetType(btp::Hadron_Decay);
    }
  }
  //std::cout<<"Hadron"<<std::endl;

  Blob *nirwana = new Blob();
  // Nirwana particles
  for (biter=blobs->begin();biter!=blobs->end();biter++) {
    for (size_t i=0;i<(size_t)(*biter)->NOutP();++i) {
      Particle *part=(*biter)->OutParticle(i);
      if (part->Status()!=1 && part->DecayBlob()==NULL) {
	nirwana->AddToInParticles(part);
      }
    }
  }
  nirwana->SetTypeSpec("Nirwana");
  blobs->push_back(nirwana);
  //std::cout<<"Nirwana"<<std::endl;


  int blobid = 0;
  for (biter=blobs->begin();biter!=blobs->end();biter++) {
    (*biter)->SetId(-blobid);
    blobid++;
  }

  return true;
}

void HepEvt_Interface::DeleteObsolete(const int mode)
{
  //std::cout<<"Delete them !! "<<mode<<std::endl;
  if (!m_convertH2S.empty()) {
    for (Translation_Map::iterator piter=m_convertH2S.begin();
	 piter!=m_convertH2S.end();piter++) {
      switch(mode) {
      case 0:
	delete (piter->second.first); piter->second.first=NULL; 
	break;  
      case 1:
	if (!piter->second.first->ProductionBlob() &&
	    !piter->second.first->DecayBlob() &&
	    piter->second.first) {
	  delete (piter->second.first); piter->second.first=NULL;
	} 
	break;  
      default:
	if (piter->second.second) {
	  delete (piter->second.first); piter->second.first=NULL; 
	}
	break;
      }
    }
    m_convertH2S.clear();
  }
}

