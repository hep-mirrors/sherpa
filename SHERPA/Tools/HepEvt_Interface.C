#include "SHERPA/Tools/HepEvt_Interface.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Org/My_MPI.H"

#include <iomanip>
#include <stdio.h>
#include <cassert>

using namespace ATOOLS;
using namespace SHERPA;
using namespace std;

HepEvt_Interface::HepEvt_Interface(int _generator) : 
  p_instream(NULL), p_outstream(NULL),
  m_evtnumber(0), m_nhep(-1), m_filesize(0), m_evtcount(0), 
  m_generator(gtp::code(_generator)), p_phep(NULL), p_vhep(NULL),
  p_jmohep(NULL), p_jdahep(NULL), p_isthep(NULL), p_idhep(NULL),
  p_pythiatranslator(NULL)
{ 
  if (m_generator==gtp::Pythia)      
    p_pythiatranslator = new Pythia_HepEvt_Translator(this);
}

HepEvt_Interface::HepEvt_Interface(gtp::code _generator) : 
  p_instream(NULL), p_outstream(NULL),
  m_evtnumber(0), m_nhep(-1), m_filesize(0), m_evtcount(0), 
  m_generator(_generator), p_phep(NULL), p_vhep(NULL),
  p_jmohep(NULL), p_jdahep(NULL), p_isthep(NULL), p_idhep(NULL),
  p_pythiatranslator(NULL)
{ 
  if (m_generator==gtp::Pythia)      
    p_pythiatranslator = new Pythia_HepEvt_Translator(this);
}

HepEvt_Interface::HepEvt_Interface() :
  p_instream(NULL), p_outstream(NULL),
  m_evtnumber(0), m_nhep(-1), m_evtcount(0), m_generator(gtp::Unspecified)
{
  // io = true : Output mode, Sherpa2HepEvt
  //   std::string filename = m_path+std::string("/")+m_file;
  //   if (m_io) {
  //     if (m_mode>=10) {
  //       m_mode      -= 10;
  //     }
  //     else {
  //       filename += std::string(".0.evts"); 
  //       p_outstream = new std::ofstream(filename.c_str(),std::ios::out);
  //       if (!p_outstream->good()) { 
  // 	msg_Error()<<"ERROR in HepEvt_Interface."<<std::endl
  // 		   <<"   Could not open event file "<<filename<<"."<<std::endl
  // 		   <<"   Will abort the run."<<std::endl;
  // 	Abort();
  //       }
  //       p_outstream->precision(10);
  //     }
  //   }
  //   else {
  //     p_instream = new std::ifstream(filename.c_str()); 
  //     if (!p_instream->good()) {
  //       msg_Error()<<"ERROR in HepEvt_Interface."<<std::endl
  // 		 <<"   Event file "<<filename<<" not found."<<std::endl
  // 		 <<"   Will abort the run."<<std::endl;
  //       Abort();
  //     }
  //     std::string gentype;
  //     (*p_instream)>>gentype>>m_filesize;
  //     if (gentype==std::string("Sherpa")) m_generator = gtp::Sherpa;
  //     if (gentype==std::string("Pythia")) m_generator = gtp::Pythia;
  //     m_evtcount=0;
  //   }
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
    msg_Error()<<"ERROR in HepEvt_Interface::ChangeOutStream"<<std::endl
	       <<"   Could not change to event file "<<filename<<"."<<std::endl
	       <<"   Will abort the run."<<std::endl;
    Abort();
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
    msg_Error()<<"ERROR in HepEvt_Interface::ChangeOutStream"<<std::endl
	       <<"   Could not change to event file "<<filename<<"."<<std::endl
	       <<"   Will abort the run."<<std::endl;
    Abort();
  }
  p_outstream->precision(10);
}




bool HepEvt_Interface::Sherpa2HepEvt(Blob_List * const _blobs) {
  m_convertS2H.clear();

  m_evtnumber++;
  //if ((m_evtnumber-1)%m_filesize==0 && p_outstream!=NULL) ChangeOutStream();

  for(size_t w=0; w<_blobs->size(); ++w) {
    Blob* blob=(*_blobs)[w]; Blob* fsblob(NULL), * fragblob(NULL);
    if(blob->Type()==btp::Hadron_Decay) {
      for(size_t u=0; u<(size_t)blob->NOutP(); ++u)
	if(blob->ConstOutParticle(u)->RefFlav().IsQuark() ||
	   blob->ConstOutParticle(u)->RefFlav().IsDiQuark() ||
	   blob->ConstOutParticle(u)->RefFlav().IsGluon()) {
	  if(u!=0) {
	    //cout<<*_blobs<<endl;
	    PRINT_INFO("Warning(1): Unusual hadron-decay blob, "<<
		       "cannot yet correct for it."); break;}
	  fsblob=blob->ConstOutParticle(0)->DecayBlob();
	  assert(blob->Id()<fsblob->Id()); bool flg(false);
	  for(size_t v=1; v<(size_t)blob->NOutP(); ++v) {
	    if(blob->ConstOutParticle(v)->DecayBlob()!=fsblob) { flg=1; break;}
	    if(blob->ConstOutParticle(v)->RefFlav().Strong() ||
	       blob->ConstOutParticle(u)->RefFlav().IsDiQuark() ||
	       blob->ConstOutParticle(u)->RefFlav().IsGluon());
	    else { flg=1; break;}
	  }
	  if(flg) {
	    PRINT_INFO("Warning(2): Unusual hadron-decay blob, "<<
		       "cannot yet correct for it."); break;}
	  if(fsblob->NInP()!=blob->NOutP()) {
	    PRINT_INFO("Warning(3): Unusual hadron-decay blob, "<<
		       "cannot yet correct for it."); break;}
	  fragblob=fsblob->ConstOutParticle(0)->DecayBlob();
	  assert(fsblob->Id()<fragblob->Id()); flg=false;
	  for(size_t v=1; v<(size_t)fsblob->NOutP(); ++v)
	    if(fsblob->ConstOutParticle(v)->DecayBlob()!=fragblob) {
	      flg=1; break;}
	  if(flg) {
	    PRINT_INFO("Warning(4): Unusual hadron-decay blob, "<<
		       "cannot yet correct for it."); break;}
	  if(fragblob->NInP()!=fsblob->NOutP()) {
	    PRINT_INFO("Warning(5): Unusual hadron-decay blob, "<<
		       "cannot yet correct for it."); break;}
	  //cout<<" ..testing.. "<<_blobs->size()<<"\n";
	  //cout<<"\n"<<*blob<<"\n"<<*fsblob<<"\n"<<*fragblob<<"\n";
	  blob->RemoveOutParticles(); fragblob->RemoveInParticles();
	  size_t Nfrag(fragblob->NOutP());
	  for(size_t v=0; v<Nfrag; ++v)
	    blob->AddToOutParticles(fragblob->RemoveOutParticle
				    (fragblob->NOutP()-1));
	  //cout<<"\n\n"<<*blob<<"\n"<<*fsblob<<"\n"<<*fragblob<<endl;
	  assert(_blobs->Delete(fsblob)); assert(_blobs->Delete(fragblob));
	  //cout<<*_blobs<<endl;
	}
    }
  }

  int nhep = 0;

  ISBlobs2HepEvt(_blobs,nhep);
  HardBlob2HepEvt(_blobs,nhep);
  ShowerBlobs2HepEvt(_blobs,nhep);
  //FSBlobs2HepEvt(_blobs,nhep);
  QEDBlobs2HepEvt(_blobs,nhep);
  FragmentationBlob2HepEvt(_blobs,nhep);
  HadronDecayBlobs2HepEvt(_blobs,nhep);
  
  m_nhep=nhep;

  msg_Debugging()<<"HEI::OTF: proc weight: "
		 <<_blobs->Weight()<<"\n";
  SetWeight(_blobs->Weight());
  Blob *signal(_blobs->FindFirst(btp::Signal_Process));
  if (signal) {
    Blob_Data_Base *facscale((*signal)["Factorisation_Scale"/*"MI_Scale"*/]);//!
    if (facscale) {
      SetQ2(facscale->Get<double>());
    }
    else THROW(fatal_error,"No factorisation scale information.");
    
    Particle_Vector inparts = signal->GetInParticles();
    if(inparts.size()==2) {
      SetFl1(inparts[0]->Flav().Kfcode());
      double E1 = inparts[0]->Momentum()[0];
      double Ebeam1 = rpa->gen.PBeam(0)[0];
      Setx1(E1/Ebeam1);
      SetFl2(inparts[1]->Flav().Kfcode());
      double E2 = inparts[1]->Momentum()[0];
      double Ebeam2 = rpa->gen.PBeam(1)[0];
      Setx2(E2/Ebeam2);
    }
    else THROW(fatal_error,"Not two signal particles.");
    
  }
  else THROW(fatal_error,"No signal process.");

  //   switch (m_mode) {
  //   case 0 :  break;
  //   case 2 :  WriteReducedHepEvt(nhep); break;
  //   case 3 :  WriteFormatedHepEvt(nhep); break;
  //   case 4 :  WriteD0HepEvt(nhep); break;
  //   default:  WriteFullHepEvt(nhep);
  //   }
  return true;
}

void HepEvt_Interface::WriteFullHepEvt(std::ostream& ostr, int nhep)
{
  ostr<<"  "<<m_evtnumber<<" "<<nhep<<"\n";
  for (int i=0;i<nhep;++i) {
    ostr<<i+1<<"  "<<p_isthep[i]<<" "<<p_idhep[i]<<" "<<p_jmohep[2*i]<<" "<<p_jmohep[2*i+1]
		  <<" "<<p_jdahep[2*i]<<" "<<p_jdahep[2*i+1]<<" \n ";
    for (int j=0;j<5;++j) ostr<<p_phep[5*i+j]<<" ";
    ostr<<"\n ";
    for (int j=0;j<4;++j) ostr<<p_vhep[4*i+j]<<" ";
    ostr<<"\n";
  }
}

void HepEvt_Interface::WriteD0HepEvt(std::ostream& ostr, int nhep)
{
  ostr<<"  "<<m_evtnumber<<" "<<nhep<<" "<<"\n";
  ostr<<"    "<<m_weight<<" "<<m_Q2<<" "<<m_x1<<" "<<m_x2<<" "<<m_fl1<<" "<<m_fl2<<"\n";
  for (int i=0;i<nhep;++i) {
    ostr<<i+1<<"  "<<p_isthep[i]<<" "<<p_idhep[i]<<" "<<p_jmohep[2*i]<<" "<<p_jmohep[2*i+1]
        <<" "<<p_jdahep[2*i]<<" "<<p_jdahep[2*i+1]<<" \n ";
    for (int j=0;j<5;++j) ostr<<p_phep[5*i+j]<<" ";
    ostr<<"\n ";
    for (int j=0;j<4;++j) ostr<<p_vhep[4*i+j]<<" ";
    ostr<<"\n";
  }
}

void HepEvt_Interface::WriteFormatedHepEvt(std::ostream& ostr, int nhep)
{
  ostr<<" "<<std::setw(4)<<nhep<<" \n";
  for (int i=0;i<nhep;++i) {
    ostr<<" "<<std::setw(8)<<p_isthep[i]<<" "<<std::setw(8)<<p_idhep[i]
		  <<" "<<std::setw(4)<<p_jmohep[2*i]<<" "<<std::setw(4)<<p_jmohep[2*i+1]
		  <<" "<<std::setw(4)<<p_jdahep[2*i]<<" "<<std::setw(4)<<p_jdahep[2*i+1]<<" \n ";
    ostr<<std::setprecision(10);
    ostr<<std::setiosflags(std::ios::fixed);
    for (int j=0;j<5;++j) ostr<<std::setw(16)<<p_phep[5*i+j]<<" ";
    ostr<<"\n ";
    for (int j=0;j<4;++j) ostr<<std::setw(16)<<p_vhep[4*i+j]<<" ";
    ostr<<"\n";
    ostr<<std::resetiosflags(std::ios::fixed);
  }
}

void HepEvt_Interface::ISBlobs2HepEvt(Blob_List * const _blobs,int & _nhep) {
  for (int beam=0;beam<2;beam++) {
    for (Blob_List::const_iterator bit=_blobs->begin(); bit!=_blobs->end();++bit) {
      if ((*bit)->Type()==btp::Bunch && (*bit)->InParticle(0)->Beam()==beam) {
	if ((*bit)->NInP()!=1) {
 	  msg_Error()<<"Error in HepEvt_Interface::ISBlobs2HepEvt."<<endl
		     <<"   Bunch blob with more than one incoming particle !"<<endl
		     <<(*bit)<<endl;
	  Abort();
	}
	if ((*bit)->NOutP()>1) {
	  Particle2HepEvt((*bit)->InParticle(0),_nhep);
	  for (int j=0;j<(*bit)->NOutP();j++) Particle2HepEvt((*bit)->OutParticle(j),_nhep);
	  EstablishRelations((*bit));
	}
      }
      if ((*bit)->Type()==btp::Beam && (*bit)->InParticle(0)->Beam()==beam) {
	if ((*bit)->NInP()!=1) {
	  msg_Error()<<"Error in HepEvt_Interface::ISBlobs2HepEvt."<<endl
		     <<"   Beam Remnant blob with more than one incoming particle !"<<endl
		     <<(*bit)<<endl;
	  Abort();
	}
	if ((*bit)->NOutP()>1) {
	  Particle2HepEvt((*bit)->InParticle(0),_nhep);
	  for (int j=0;j<(*bit)->NOutP();j++) Particle2HepEvt((*bit)->OutParticle(j),_nhep);
	  EstablishRelations((*bit));
	}
      }
//       if ((*bit)->Type()==btp::IS_Shower && (*bit)->Beam()==beam) {
// 	if ((*bit)->NInP()!=1) {
// 	  msg_Error()<<"Error in HepEvt_Interface::ISBlobs2HepEvt."<<endl
// 		     <<"   IS blob with more than one incoming particle !"<<endl
// 		     <<(*bit)<<endl;
// 	  Abort();
// 	}
// 	Particle2HepEvt((*bit)->InParticle(0),_nhep);
// 	for (int j=0;j<(*bit)->NOutP();j++) Particle2HepEvt((*bit)->OutParticle(j),_nhep);
// 	EstablishRelations((*bit));
//       }
    }
  }
}


void HepEvt_Interface::HardBlob2HepEvt(Blob_List * const _blobs,int & _nhep) {
//   int mo,da;
  for(Blob_List::const_iterator bit=_blobs->begin(); bit!=_blobs->end();
      ++bit) {
//     if ((*bit)->Type()==btp::ME_PS_Interface_IS) {
//       if ((*bit)->NInP()!=2 || (*bit)->NOutP()!=2) {
// 	msg_Error()<<"Error in HepEvt_Interface::HardBlob2HepEvt."<<endl
// 		   <<"   ME_PS_Interface_IS blob with other than 2->2 particles !"<<endl
// 		   <<(*bit)<<endl;
// 	Abort();
//       }
//       else {
// 	for (int i=0;i<2;i++) {
// 	  Particle2HepEvt((*bit)->InParticle(i),_nhep);
// 	  Particle2HepEvt((*bit)->OutParticle(i),_nhep);
// 	  mo = m_convertS2H[(*bit)->InParticle(i)];
// 	  da = m_convertS2H[(*bit)->OutParticle(i)];
// 	  for (int j=0;j<2;j++) {
// 	    p_jmohep[2*da+j] = mo+1; p_jdahep[2*mo+j] = da+1;
// 	  } 
// 	}
//       }
//     }
    if ((*bit)->Type()==btp::Signal_Process ||
	(*bit)->Type()==btp::Hard_Collision) {
      if ((*bit)->NInP()!=2) {
	msg_Error()<<"Error in HepEvt_Interface::HardBlob2HepEvt."<<endl
		   <<"   Hard ME blob with other than 2 incoming particles !\n"
		   <<(*bit)<<endl;
	Abort();
      }
      if ((*bit)->NOutP()>=2) {
	Particle2HepEvt((*bit)->InParticle(0),_nhep);
	Particle2HepEvt((*bit)->InParticle(1),_nhep);
	for(int j=0;j<(*bit)->NOutP(); ++j)
	  Particle2HepEvt((*bit)->OutParticle(j),_nhep);
	EstablishRelations((*bit));
      }
    }
//     if ((*bit)->Type()==btp::ME_PS_Interface_FS) {
//       if ((*bit)->NInP()<2 || (*bit)->NOutP()!=(*bit)->NInP()) {
// 	msg_Error()<<"Error in HepEvt_Interface::HardBlob2HepEvt."<<endl
// 		   <<"   ME_PS_Interface_IS blob with other than n->n particles !"<<endl
// 		   <<(*bit)<<endl;
// 	Abort();
//       }
//       else {
// 	for (int i=0;i<(*bit)->NOutP();i++) {
// 	  Particle2HepEvt((*bit)->InParticle(i),_nhep);
// 	  Particle2HepEvt((*bit)->OutParticle(i),_nhep);
// 	  mo = m_convertS2H[(*bit)->InParticle(i)];
// 	  da = m_convertS2H[(*bit)->OutParticle(i)];
// 	  for (int j=0;j<2;j++) {
// 	    p_jmohep[2*da+j] = mo+1; p_jdahep[2*mo+j] = da+1;
// 	  } 
// 	}
//       }
//     }
  }
}

void HepEvt_Interface::ShowerBlobs2HepEvt(ATOOLS::Blob_List * const _blobs,
					  int & _nhep) {
  for(Blob_List::const_iterator bit=_blobs->begin();
      bit!=_blobs->end(); ++bit) {
    if((*bit)->Type()==btp::Shower
       //&& ((*bit)->NInP()==1 || (*bit)->NInP()==2)
       //&& (*bit)->Beam()==beam
       ) {
      //for(int j=0; j<(*bit)->NInP(); ++j)
      //  Particle2HepEvt((*bit)->InParticle(j),_nhep);
      for(int j=0; j<(*bit)->NOutP(); ++j) {
	//cout<<(*bit)->OutParticle(j)->Info()<<endl;
	if((*bit)->OutParticle(j)->Info()!='G')
	  Particle2HepEvt((*bit)->OutParticle(j),_nhep);
      }
      EstablishRelationsModified((*bit));
    }
  }
}

void HepEvt_Interface::QEDBlobs2HepEvt(ATOOLS::Blob_List * const _blobs,
				       int & _nhep) {
  for(Blob_List::const_iterator bit=_blobs->begin();
      bit!=_blobs->end(); ++bit) {
    if((*bit)->Type()==btp::QED_Radiation) {
      /*if((*bit)->NInP()<=2)*/ {
	for(int j=0; j<(*bit)->NOutP(); ++j)
	  Particle2HepEvt((*bit)->OutParticle(j),_nhep);
	EstablishRelations((*bit));
	//if((*bit)->NInP()>2) cout<<*_blobs<<endl;
      }
      //else
      //PRINT_INFO("Warning: QED Radiation blob with more than two incoming "<<
      //	   "particles. Translation skipped.");
    }
  }
}

void HepEvt_Interface::FSBlobs2HepEvt(Blob_List * const _blobs,int & _nhep) {
  for (Blob_List::const_iterator bit=_blobs->begin(); bit!=_blobs->end();++bit) {
    //if ((*bit)->Type()==btp::FS_Shower &&
    //((*bit)->NInP()==1 || (*bit)->NInP()==2)) {
    //for (int j=0;j<(*bit)->NOutP();j++) Particle2HepEvt((*bit)->OutParticle(j),_nhep);
    //EstablishRelations((*bit));
    //}
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
	msg_Error()<<"Error in HepEvt_Interface::HadronDecays2HepEvt."<<endl
		   <<"   Decay blob with other than 1 incoming particles !"<<endl
		   <<(*bit)<<endl;
	Abort();
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
	msg_Error()<<"Warning : Potential error in HepEvt_Interface::HadronDecays2HepEvt."<<endl
		   <<"   Decay blob for 1 -> 1 process with no identified mxing !"<<std::endl;
      }
    }
  }
}

void HepEvt_Interface::Particle2HepEvt(Particle * const _part,int & _nhep)
{
  int number = m_convertS2H.count(_part);
  if (number>0) return;
  p_idhep[_nhep]    = (long int)_part->Flav();
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


void HepEvt_Interface::EstablishRelations(Blob* const _blob) {
  int mothers[2];
  int daughters[2];
  int inum(_blob->NInP());
  mothers[0] = mothers[1] = 0;
  switch(inum) {
  case  0: break;
  case  2: mothers[1] = m_convertS2H[_blob->InParticle(1)];
  case  1: mothers[0] = m_convertS2H[_blob->InParticle(0)]; break;
  default: //PRINT_INFO("Note, evt#"<<m_evtnumber
                      //<<": blob (type="<<_blob->Type()
                      //<<") with "<<inum<<" incoming particles.");
    mothers[0] = m_convertS2H[_blob->InParticle(0)];
    mothers[1] = m_convertS2H[_blob->InParticle(inum-1)];
  }
  if(_blob->NOutP()>0) {
    daughters[0] = m_convertS2H[_blob->OutParticle(0)];
    if(_blob->NOutP()>1)
      daughters[1] = m_convertS2H[_blob->OutParticle(_blob->NOutP()-1)];
    else daughters[1] = m_convertS2H[_blob->OutParticle(0)];
  }
  else daughters[0] = daughters[1] = 0;

  if(inum>2) {
    for(int i=1; i<inum-1; ++i) {
      int tmp=m_convertS2H[_blob->InParticle(i)];
      p_jdahep[2*tmp]   = daughters[0]+1;
      p_jdahep[2*tmp+1] = daughters[1]+1;
    }
    inum=2;
  }
  for(int i=0; i<inum; ++i) {
    p_jdahep[2*mothers[i]]   = daughters[0]+1;
    p_jdahep[2*mothers[i]+1] = daughters[1]+1;
  }
  for(int i=0; i<_blob->NOutP(); ++i) for (int j=0; j<inum; ++j) { 
    p_jmohep[2*m_convertS2H[_blob->OutParticle(i)]+j] = mothers[j]+1; 
  }
}


void HepEvt_Interface::EstablishRelationsModified(Blob * const _blob) {
  int mothers[2];
  int daughters[2];
  int incoming[2];
  mothers[0]=mothers[1]=0;
  for(int i=0, ii=0; i<_blob->NInP(); ++i) {
    if(_blob->InParticle(i)->Info()=='I') {
      assert(ii<2); mothers[ii]=m_convertS2H[_blob->InParticle(i)]; ++ii;}
  }
  incoming[0]=incoming[1]=0;
  if(_blob->NOutP()>0) {
    int outn(0);
    while(outn<_blob->NOutP() && _blob->OutParticle(outn)->Info()=='G') {
      if(outn==0) incoming[0]=m_convertS2H[_blob->OutParticle(0)];
      ++outn;
    }
    if(outn<_blob->NOutP()) {
      if(outn>0) incoming[1]=m_convertS2H[_blob->OutParticle(outn-1)];
      daughters[0]=m_convertS2H[_blob->OutParticle(outn)];
      if(_blob->NOutP()>outn+1)
	daughters[1]=m_convertS2H[_blob->OutParticle(_blob->NOutP()-1)];
      else daughters[1]=m_convertS2H[_blob->OutParticle(outn)];
    }
    else daughters[0]=daughters[1]=0;
  }
  else daughters[0]=daughters[1]=0;

  if(mothers[0]) {
    p_jdahep[2*mothers[0]]   = incoming[0]+1;  //daughters[0]+1;
    p_jdahep[2*mothers[0]+1] = incoming[1]+1;} //daughters[1]+1;}
  if(mothers[1]) {
    p_jdahep[2*mothers[1]]   = incoming[0]+1;  //daughters[0]+1;
    p_jdahep[2*mothers[1]+1] = incoming[1]+1;} //daughters[1]+1;}
  if(_blob->NInP()>0) {
    for(int i=0; i<_blob->NInP(); ++i) {
      if(_blob->InParticle(i)->Info()!='I') {
	p_jdahep[2*m_convertS2H[_blob->InParticle(i)]] = daughters[0]+1;
	p_jdahep[2*m_convertS2H[_blob->InParticle(i)]+1] = daughters[1]+1;
      }
    }
  }
  if(_blob->NOutP()>0) {
    for(int i=0; i<_blob->NOutP(); ++i) {
      //if(_blob->OutParticle(i)->Info()!='G') {
        if(mothers[0])
	  p_jmohep[2*m_convertS2H[_blob->OutParticle(i)]] = mothers[0]+1;
	if(mothers[1])
	  p_jmohep[2*m_convertS2H[_blob->OutParticle(i)]+1] = mothers[1]+1;
      //}
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
  bool okay(false);
  //std::cout<<METHOD<<" : "<<p_instream<<" "<<m_nhep<<std::endl;
  switch (m_generator)  {
    case gtp::Pythia:  okay = p_pythiatranslator->ConstructBlobs(blobs); break;
  case gtp::Sherpa:  //okay = ConstructBlobs(blobs); break;
    default:
      msg_Error()<<"Error in HepEvt_Interface::ReadHepEvt."<<std::endl
		 <<"   Generator type unspecified : "<<m_generator<<std::endl
		 <<"   Will abort the run."<<std::endl;
      Abort();
      
  }
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
    msg_Error()<<"ERROR in "<<METHOD<<" : "<<std::endl
	       <<"   Event file "<<filename<<" not found."<<std::endl
	       <<"   Will abort the run."<<std::endl;
    Abort();
  }
  std::string gentype;
  (*p_instream)>>gentype>>m_filesize;
  if ((gentype==std::string("Sherpa") && m_generator!=gtp::Sherpa) ||
      (gentype==std::string("Pythia") && m_generator!=gtp::Pythia)) {
    msg_Error()<<"ERROR in "<<METHOD<<" : "<<std::endl
	       <<"   Types do not match : "<<gentype<<" vs. "<<int(m_generator)<<std::endl
	       <<"   Abort the run."<<std::endl;
    Abort();
  }

  m_evtcount=0;
}

//-----------------------------------------------------------------------------



