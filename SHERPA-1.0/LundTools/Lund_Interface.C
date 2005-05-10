#include "Lund_Interface.H"

#include "Lund_Wrapper.H"
#include "ISR_Handler.H"
#include "Data_Reader.H"
#include "Particle.H"
#include "Run_Parameter.H"
#include "Random.H"
#include "Message.H"
#include "Exception.H"
#include "Running_AlphaS.H"
#include "MyStrStream.H"
#include <list>
#include "Message.H"

using namespace SHERPA;
using namespace std;

bool Lund_Interface::s_exportas=false;
bool Lund_Interface::s_exportpdf=false;

size_t Lund_Interface::s_errors=0;
size_t Lund_Interface::s_maxerrors=0;

ATOOLS::Blob_List *Lund_Interface::s_bloblist=NULL; 
PDF::ISR_Handler *Lund_Interface::s_isrhandler=NULL; 

Lund_Interface::Lund_Interface(std::string _m_path,std::string _m_file,bool sherpa):
  m_path(_m_path),m_file(_m_file), m_maxtrials(1),
  p_hepevt(NULL), 
  m_compress(true),m_writeout(false),
  p_phep(new double[5*4000]),
  p_vhep(new double[4*4000]),
  p_jmohep(new int[2*4000]),
  p_jdahep(new int[2*4000])
{
  double win;
  std::string beam[2], frame("CMS");
  ATOOLS::Flavour flav[2];
  for (size_t i=0;i<2;++i) flav[i]=ATOOLS::rpa.gen.Bunch(i);
  if (flav[0]==ATOOLS::kf::e && flav[1]==ATOOLS::kf::p_plus) {
    beam[0]="e-";
    beam[1]="p+";
  }
  else if (flav[0]==ATOOLS::kf::p_plus && flav[1]==ATOOLS::kf::e) {
    beam[0]="p+";
    beam[1]="e-";
  }
  else if (flav[0]==ATOOLS::kf::e && flav[1]==ATOOLS::kf::photon) {
    if (flav[0].IsAnti()) beam[0]="e+"; else beam[0]="e-";
    beam[1]="gamma";
    pysubs.msub[33]=1;    
  }
  else if (flav[0]==ATOOLS::kf::photon && flav[1]==ATOOLS::kf::e) {
    beam[0]="gamma";
    if (flav[1].IsAnti()) beam[1]="e+"; else beam[1]="e-";
    pysubs.msub[33]=1;    
  }
  else if (flav[0]==ATOOLS::kf::photon && flav[1]==ATOOLS::kf::photon) {
    for (size_t i=0;i<2;++i) beam[i]="gamma";
    pysubs.msub[57]=1;    
  }
  else if (flav[0].Kfcode()==ATOOLS::kf::e && flav[1].Kfcode()==ATOOLS::kf::e) {
    for (size_t i=0;i<2;++i) if (flav[i].IsAnti()) beam[i]="e+"; else beam[i]="e-";
    pysubs.msub[0]=1;    
    pypars.mstp[47]=1;
    pydat1.mstj[100]=5;
  }
  else {
    for (size_t i=0;i<2;++i) if (flav[i].IsAnti()) beam[i]="p-"; else beam[i]="p+";
    pysubs.msub[0]=1;    
    pypars.mstp[47]=1;
    pydat1.mstj[100]=5;
  }
  win=ATOOLS::rpa.gen.Ecms();
  s_maxerrors=ATOOLS::rpa.gen.NumberOfEvents();
  std::vector<std::vector<double> > help;
  ATOOLS::Data_Reader *reader = new ATOOLS::Data_Reader("=",";","!");
  reader->SetMatrixType(reader->MTransposed);
  reader->SetInputPath(m_path);
  reader->SetInputFile(m_file);
  reader->AddIgnore("(");
  reader->AddIgnore(")");
  reader->AddIgnore(",");
  if (!sherpa) {
    if (!reader->ReadFromFile(pysubs.msel,"MSEL")) pysubs.msel=1;
  }
  else {
    pysubs.msel=0;
  }
  reader->MatrixFromFile(help,"MSUB");
  for (size_t i=0;i<help.size();++i) {
    if (help[i].size()>1) if ((int)help[i][0]>0) pysubs.msub[(int)help[i][0]-1]=(int)help[i][1];
  }
  reader->MatrixFromFile(help,"KFIN");
  for (size_t i=0;i<help.size();++i) {
    if (help[i].size()>2) {
      if (((int)help[i][0]>0)&&((int)help[i][1]>-41)) {
	pysubs.kfin[(int)help[i][1]+40][(int)help[i][0]-1]=(int)help[i][2];
      }
    }
  }
  reader->MatrixFromFile(help,"CKIN");
  for (size_t i=0;i<help.size();++i) {
    if (help[i].size()>1) if ((int)help[i][0]>0) pysubs.ckin[(int)help[i][0]-1]=help[i][1];
  }
  reader->MatrixFromFile(help,"MSTJ");
  for (size_t i=0;i<help.size();++i) {
    if (help[i].size()>1) if ((int)help[i][0]>0) pydat1.mstj[(int)help[i][0]-1]=(int)help[i][1];
  }
  reader->MatrixFromFile(help,"MSTP");
  for (size_t i=0;i<help.size();++i) {
    if (help[i].size()>1) if ((int)help[i][0]>0) pypars.mstp[(int)help[i][0]-1]=(int)help[i][1];
  }
  reader->MatrixFromFile(help,"MSTU");
  for (size_t i=0;i<help.size();++i) {
    if (help[i].size()>1) if ((int)help[i][0]>0) pydat1.mstu[(int)help[i][0]-1]=(int)help[i][1];
  }
  reader->MatrixFromFile(help,"MFUDGE");
  for (size_t i=0;i<help.size();++i) {
    if (help[i].size()>1) if ((int)help[i][0]>0) cfudge.mfudge[(int)help[i][0]-1]=(int)help[i][1];
  }
  reader->MatrixFromFile(help,"PARP");
  for (size_t i=0;i<help.size();++i) {
    if (help[i].size()>1) if ((int)help[i][0]>0) pypars.parp[(int)help[i][0]-1]=help[i][1];
  }
  reader->MatrixFromFile(help,"PARJ");
  for (size_t i=0;i<help.size();++i) {
    if (help[i].size()>1) if ((int)help[i][0]>0) pydat1.parj[(int)help[i][0]-1]=help[i][1];
  }
  reader->MatrixFromFile(help,"PARU");
  for (size_t i=0;i<help.size();++i) {
    if (help[i].size()>1) if ((int)help[i][0]>0) pydat1.paru[(int)help[i][0]-1]=help[i][1];
  }
  reader->MatrixFromFile(help,"PFUDGE");
  for (size_t i=0;i<help.size();++i) {
    if (help[i].size()>1) if ((int)help[i][0]>0) cfudge.pfudge[(int)help[i][0]-1]=help[i][1];
  }
  reader->MatrixFromFile(help,"MDME");
  for (size_t i=0;i<help.size();++i) {
    if (help[i].size()>2) {
      if ((int)help[i][0]>0 && abs((int)help[i][1]<2)) {
		pydat3.mdme[(int)help[i][1]-1][(int)help[i][0]-1]=(int)help[i][2];
      }
    }
  }
  reader->MatrixFromFile(help,"MDCYKF");
  for (size_t i=0;i<help.size();++i) {
    if (help[i].size()>2) {
      if ((int)help[i][0]>0 && abs((int)help[i][1]<2)) {
		msg_Tracking()<<"Lund_Interface::Lund_Interface(..): "
					  <<"Set MDCY("<<pycomp((int)help[i][0])<<","<<(int)help[i][1]
					  <<") ( from KF code "<<(int)help[i][0]<<" ) to "<<(int)help[i][2]<<std::endl;
		pydat3.mdcy[(int)help[i][1]-1][pycomp((int)help[i][0])-1]=(int)help[i][2];
      }
    }
  }
  // the next lines replace the apyinit_ call
  if (sherpa) {
    hepevt.nhep=100;
    for (int i=pydat3.mdcy[23-1][2-1];i<pydat3.mdcy[23-1][2-1]+pydat3.mdcy[23-1][3-1];++i) {
      if (abs(pydat3.kfdp[i-1][1-1])>=2) pydat3.mdme[i-1][1-1]=ATOOLS::Min(0,pydat3.mdme[i-1][1-1]);
    }
    pyinit(frame.c_str(),beam[0].c_str(),beam[1].c_str(),win);
  }
  // replacement ends here
  if (ATOOLS::msg.LevelIsDebugging()) ListLundParameters();
  if (!sherpa) {
    int helpi;
    if (reader->ReadFromFile(helpi,"EXPORT_ALPHAS")) s_exportas=(bool)helpi;
    if (reader->ReadFromFile(helpi,"EXPORT_PDF")) s_exportpdf=(bool)helpi;
    int orderas;
    double asmz, asdef, mz;  
    reader->SetInputFile("Model.dat");
    if (!reader->ReadFromFile(orderas,"ORDER_ALPHAS")) orderas=0;
    if (!reader->ReadFromFile(asmz,"ALPHAS(MZ)")) asmz=0.1188;
    if (!reader->ReadFromFile(asdef,"ALPHAS(default)")) asdef=asmz;
    mz=91.188;
    reader->SetInputFile("Particle.dat");
    std::vector<std::vector<double> > helpdvv;
    reader->MatrixFromFile(helpdvv,"");
    for (size_t i=0;i<helpdvv.size();++i) {
      if (helpdvv[i][0]==24. && helpdvv.size()>1) mz=helpdvv[i][1];
    }
    MODEL::as = new MODEL::Running_AlphaS(asmz,mz*mz,orderas);
    MODEL::as->SetDefault(asdef);
    p_hepevt = new ATOOLS::HepEvt_Interface(ATOOLS::gtp::Pythia);
    if (pypars.mstp[105-1]==0) p_hepevt->SetHadronized(false);
    pyinit(frame.c_str(),beam[0].c_str(),beam[1].c_str(),win);
    if (reader->ReadFromFile(m_outfile,"OUTPUT_FILE")) {
      m_writeout = true;
      if (!reader->ReadFromFile(m_evtsperfile,"EVENTS_PER_FILE")) m_evtsperfile=1000;
      NextFile(true);
      int helper;
      if (!reader->ReadFromFile(helper,"COMPRESS")) helper=1;
      m_compress=(bool)helper;
    }
  }
  delete reader;
  // print decay channel list
  //   pylist(12);
  // for test purposes - leave only one tau decay channel open
//  LeaveOpenOnlyDC(93);
}

void Lund_Interface::SwitchOfDecays(ATOOLS::kf::code kfc)
{
  pydat3.mdcy[1-1][pycomp(int(kfc))-1]=0;
  msg_Tracking()<<"Lund_Interface::SwitchOfDecays: "<<kfc<<std::endl;
}

void Lund_Interface::LeaveOpenOnlyDC( int idc )
{
  for( int i=89; i<=142; i++ ) pydat3.mdme[1-1][i-1] = 0; 
  pydat3.mdme[1-1][idc-1] = 1;
}

Lund_Interface::~Lund_Interface()
{
  NextFile(false);
  if (p_hepevt) { 
    p_hepevt->SetNhep(0);
    p_hepevt->SetIsthep(NULL);
    p_hepevt->SetIdhep(NULL);
    p_hepevt->SetJmohep(NULL);
    p_hepevt->SetJdahep(NULL);
    p_hepevt->SetPhep(NULL);
    p_hepevt->SetVhep(NULL);
    delete p_hepevt; p_hepevt = NULL; 
  }
  if (p_jmohep) { delete p_jmohep; p_jmohep = NULL; }
  if (p_jdahep) { delete p_jdahep; p_jdahep = NULL; }
  if (p_phep)   { delete p_phep;   p_phep   = NULL; }
  if (p_vhep)   { delete p_vhep;   p_vhep   = NULL; }
}

bool Lund_Interface::Hadronize(ATOOLS::Blob_List *bloblist,
			       ATOOLS::Particle_List *particlelist) 
{
  if (ExtractSinglets(bloblist,particlelist)) {
//     for (std::list<ATOOLS::Particle_List *>::iterator lit=m_partlists.begin();
// 	 lit!=m_partlists.end();++lit) PRINT_INFO(**lit);
    if (m_partlists.size()==0) {
	  return true; // no colours in event
	}
    ATOOLS::Blob * blob = new ATOOLS::Blob();
    blob->SetId();
    blob->SetType(ATOOLS::btp::Fragmentation);
    blob->SetTypeSpec("Pythia_v6.214");
    bloblist->push_back(blob);
//     PRINT_INFO("before prepareblob");
    
    int nhep = PrepareFragmentationBlob(blob);
    for (size_t trials=0;trials<m_maxtrials;++trials) {
//       PRINT_INFO("before frag");
      if (StringFragmentation(blob,bloblist,particlelist,nhep)) return true;
      if (m_maxtrials>1) 
		ATOOLS::msg.Error()<<"Error in Lund_Interface::Hadronize ."<<std::endl
						   <<"   Hadronization failed. Retry event."<<std::endl;
    }
  }
//     PRINT_INFO("before end");
  
  ATOOLS::msg.Error()<<"Error in Lund_Interface::Hadronize ."<<std::endl
					 <<"   Hadronization failed. Retry event."<<std::endl;
  while (bloblist->size()>0) {
    delete *bloblist->begin();
    bloblist->erase(bloblist->begin());
  }
  return false;
}

void Lund_Interface::Reset()
{
  if (m_partlists.size()>0) {
    for (std::list<ATOOLS::Particle_List *>::iterator plit=m_partlists.begin();
	 plit!=m_partlists.end();plit++) {
      do { (*plit)->pop_back(); } while (!(*plit)->empty());
      delete (*plit); 
    }
  }
  m_partlists.clear();
}


bool Lund_Interface::ExtractSinglets(ATOOLS::Blob_List *bloblist,ATOOLS::Particle_List *pl) 
{
  Reset();
  std::list<ATOOLS::Particle *> plist;
  ATOOLS::Particle * part;
  for (ATOOLS::Blob_List::iterator blit=bloblist->begin();blit!=bloblist->end();++blit) {
    if ((*blit)->Type()==ATOOLS::btp::FS_Shower || 
		(*blit)->Type()==ATOOLS::btp::IS_Shower ||
		(*blit)->Type()==ATOOLS::btp::Hard_Collision ||
		(*blit)->Type()==ATOOLS::btp::Beam) {
 	      for (int i=0;i<(*blit)->NOutP();i++) {
			part = (*blit)->OutParticle(i); 
			if (part->Status()==1 && part->DecayBlob()==NULL && 
	    		(part->GetFlow(1)!=0 || part->GetFlow(2)!=0)) {
	  		  plist.push_back(part);
			}
			if (part->Status()==1 && part->DecayBlob()==NULL) {
			  if (part->Flav()==ATOOLS::Flavour(ATOOLS::kf::tau) ||
			      part->Flav()==ATOOLS::Flavour(ATOOLS::kf::tau).Bar()) {
			    plist.push_back(part);
			  }
			}
	      }
    }
  }
  if (plist.empty()) {
    ATOOLS::msg.Debugging()<<"WARNING in Lund_Interface::PrepareFragmentationBlob:"<<std::endl
			   <<"   No coloured particle found leaving shower blobs."<<std::endl;
    return true;
  }

  int  col1, col2;
  bool hit1, hit2;
  size_t actsize;
  ATOOLS::Particle_List * pli;
  do {
    actsize = plist.size();
    hit1    = false;
    //std::cout<<"1st round : "<<plist.size()<<std::endl; 
    for (std::list<ATOOLS::Particle *>::iterator pit1=plist.begin();pit1!=plist.end();++pit1) {
      col1 = (*pit1)->GetFlow(1);
      col2 = (*pit1)->GetFlow(2);
      //std::cout<<"Test : "<<col1<<" "<<col2<<std::endl;
      if (col1!=0 && col2==0) {
		hit1 = true;
		pli  = new ATOOLS::Particle_List;
		pli->push_back((*pit1));
		pit1 = plist.erase(pit1);
		//std::cout<<"Part_List "<<m_partlists.size()+1<<" : "<<pli->size()<<std::endl;
		m_partlists.push_back(pli);
		do {
		  hit2 = false;
		  for (std::list<ATOOLS::Particle *>::iterator pit2=plist.begin();pit2!=plist.end();++pit2) {
			if ((int)((*pit2)->GetFlow(2))==col1) {
			  col1 = (*pit2)->GetFlow(1);
			  pli->push_back((*pit2));
			  //std::cout<<"Part_List "<<m_partlists.size()<<" : "<<pli->size()<<std::endl;
			  plist.erase(pit2);
			  hit2 = true;
			  break;
			}
		  }
		} while (hit2);
      }
      if (hit1) break;
    } 
  } while (actsize!=plist.size());
  

  do {
    actsize = plist.size();
    hit1    = false;
    //std::cout<<"2nd round : "<<plist.size()<<std::endl; 
    for (std::list<ATOOLS::Particle*>::iterator pit1=plist.begin();pit1!=plist.end();++pit1) {
      col1 = (*pit1)->GetFlow(1);
      col2 = (*pit1)->GetFlow(2);
      if (col1!=0 && col2!=0) {
		hit1 = true;
		pli  = new ATOOLS::Particle_List;
		pli->push_back((*pit1));
		pit1 = plist.erase(pit1);
		//std::cout<<"New Part_List "<<m_partlists.size()+1<<" : "<<pli->size()<<" : "<<col1<<" "<<col2<<std::endl;
		m_partlists.push_back(pli);
		do {
		  hit2 = false;
		  for (std::list<ATOOLS::Particle *>::iterator pit2=plist.begin();pit2!=plist.end();++pit2) {
			if ((int)((*pit2)->GetFlow(2))==col1) {
			  col1 = (*pit2)->GetFlow(1);
			  col2 = (*pit2)->GetFlow(2);
			  pli->push_back((*pit2));
			  //std::cout<<"New particle "<<col1<<" "<<col2<<std::endl;
			  //std::cout<<"Part_List "<<m_partlists.size()<<" : "<<pli->size()<<std::endl;
			  pit2 = plist.erase(pit2);
			  hit2 = true;
			  break;
			}
		  }
		} while (hit2);
      }
      if (hit1) break;
    }
  } while (actsize!=plist.size());
  
  
  for (std::list<ATOOLS::Particle*>::iterator pit1=plist.begin(); pit1!=plist.end();) {
    //     std::cout<<" Test start "<<(*pit1)<<std::endl;
    //     std::cout<<" Test end "<<std::endl;
	if ((*pit1)->Flav()==ATOOLS::Flavour(ATOOLS::kf::tau) ||
		(*pit1)->Flav()==ATOOLS::Flavour(ATOOLS::kf::tau).Bar()) {
	  pli  = new ATOOLS::Particle_List;
	  pli->push_back((*pit1));
	  m_partlists.push_back(pli);
	  pit1 = plist.erase(pit1);
	}
	else pit1++;
  }	  
  if (!plist.empty()) {
    ATOOLS::msg.Error()<<"ERROR in Lund_Interface::ExtractSinglets : "<<std::endl
		       <<"   Failed to arrange all particles leaving shower blobs in colour singlets."<<std::endl
		       <<"   Remaining particles are : "<<std::endl;
    for (std::list<ATOOLS::Particle*>::iterator pit1=plist.begin();pit1!=plist.end();++pit1) 
      ATOOLS::msg.Error()<<"   "<<(*pit1)<<std::endl;
    ATOOLS::msg.Error()<<"   Reset, return false and hope for the best."<<std::endl;
    Reset();
    return false;
  }

  /*  
  std::cout<<"===================="<<m_partlists.size()<<" Lists : "<<std::endl;
  for (std::list<ATOOLS::Particle_List *>::iterator plit=m_partlists.begin();
       plit!=m_partlists.end();plit++) {
    std::cout<<(**plit)<<std::endl;
    std::cout<<"-----------------------------------------------------"<<std::endl;
  }
  */

  return true;
}

int Lund_Interface::PrepareFragmentationBlob(ATOOLS::Blob *blob) 
{
  int nhep = 0;
  if (nhep==0) {
    hepevt.idhep[nhep]=ATOOLS::Flavour(ATOOLS::kf::photon).HepEvt();
    for (short int j=1;j<4;++j) hepevt.phep[nhep][j-1]=blob->CMS()[j];
    hepevt.phep[nhep][3]=blob->CMS()[0];
    double pabs=(blob->CMS()).Abs2();
    if (pabs<0) hepevt.phep[nhep][4]=0.0;
    else hepevt.phep[nhep][4]=sqrt(pabs);
    for (short int j=0;j<4;++j) hepevt.vhep[nhep][j]=0.0;
    hepevt.isthep[nhep]=1;
    hepevt.jmohep[nhep][0]=0;
    hepevt.jmohep[nhep][1]=0;
    hepevt.jdahep[nhep][0]=0;
    hepevt.jdahep[nhep][1]=0;
  }
  
  ATOOLS::Particle * part;

  for (std::list<ATOOLS::Particle_List *>::iterator plit=m_partlists.begin();
       plit!=m_partlists.end();) {
    part = *(*plit)->begin();
    if (part->Flav()==ATOOLS::Flavour(ATOOLS::kf::tau) ||
		part->Flav()==ATOOLS::Flavour(ATOOLS::kf::tau).Bar()) {
      while ((*plit)->size()>0) {
		part = *(*plit)->begin();
		blob->AddToInParticles(part);
		AddPartonToString(part,nhep);
		(*plit)->pop_front();
      }
      delete (*plit);
      plit = m_partlists.erase(plit);
    }
    else plit++;
  }

  // Colourful stuff
  if (m_partlists.size()==0) return nhep;
  for (std::list<ATOOLS::Particle_List *>::iterator plit=m_partlists.begin();
       plit!=m_partlists.end();) {
    part = *(*plit)->begin();
    if (part->GetFlow(1)!=0 && part->GetFlow(2)==0) {
      while ((*plit)->size()>0) {
		part = *(*plit)->begin();
		blob->AddToInParticles(part);
		AddPartonToString(part,nhep);
		(*plit)->pop_front();
      }
      delete (*plit);
      plit = m_partlists.erase(plit);
    }
    else plit++;
  }
  
  if (m_partlists.size()==0) return nhep;
  ATOOLS::Particle * help,* help1;
  for (std::list<ATOOLS::Particle_List *>::iterator plit=m_partlists.begin();
       plit!=m_partlists.end();) {
    part = *(*plit)->begin();
    if (part->GetFlow(1)!=0 && part->GetFlow(2)!=0) {
      blob->AddToInParticles(part);
      ATOOLS::Flavour            flav = ATOOLS::Flavour(ATOOLS::kf::d);
      if (ATOOLS::ran.Get()<0.5) flav = ATOOLS::Flavour(ATOOLS::kf::u);
      help = new ATOOLS::Particle(-1,flav,0.5*part->Momentum());
      help->SetProductionBlob(part->ProductionBlob());
      AddPartonToString(help,nhep);
      delete help;
      help1 = part;
      (*plit)->pop_front();
	  
      while ((*plit)->size()>0) {
		part = *(*plit)->begin();
		blob->AddToInParticles(part);
		AddPartonToString(part,nhep);
		(*plit)->pop_front();
      }
      help = new ATOOLS::Particle(-1,flav.Bar(),0.5*help1->Momentum());
      help->SetProductionBlob(help1->ProductionBlob());
      AddPartonToString(help,nhep);
      delete help;
	  
      delete (*plit);
      plit = m_partlists.erase(plit);
    }
    else {
      ATOOLS::msg.Error()<<"ERROR in Lund_Interface::PrepareFragmentationBlob : "<<std::endl
						 <<"   No octet found although expected, particle is : "<<std::endl
						 <<(*plit)<<std::endl
						 <<"   will abort."<<std::endl;
      abort();
    }
  }
  
  Reset();
  return nhep;
}



bool Lund_Interface::StringFragmentation(ATOOLS::Blob *blob,ATOOLS::Blob_List *bloblist,
					 ATOOLS::Particle_List *pl,int nhep) 
{
  hepevt.nevhep=0;
  hepevt.nhep=nhep;
  pyhepc(2);
  pydat1.mstu[70-1]=1;
  pydat1.mstu[71-1]=hepevt.nhep;
  size_t errors=s_errors;
  //std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<std::endl;
  //pylist(1);
  pyexec();
  //std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<std::endl;
  //pylist(2);
  //std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<std::endl;
  pyhepc(1);
  pydat1.mstu[70-1]=2;
  pydat1.mstu[72-1]=hepevt.nhep;
  if (ATOOLS::msg.LevelIsDebugging()) {
    ATOOLS::msg.Out()<<std::endl<<std::endl;
  }
  FillPrimaryHadronsInBlob(blob,bloblist,pl);
  return s_errors==errors;
}

void Lund_Interface::AddPartonToString(ATOOLS::Particle *parton,int &nhep)
{
  hepevt.idhep[nhep]=parton->Flav().HepEvt();
  for (short int j=1; j<4; ++j) hepevt.phep[nhep][j-1]=parton->Momentum()[j];
  hepevt.phep[nhep][3]=parton->Momentum()[0];
  double pabs=(parton->Momentum()).Abs2();
  if (pabs<0) hepevt.phep[nhep][4]=0.0;
  else hepevt.phep[nhep][4]=sqrt(pabs);
  for (short int j=1;j<4;++j) hepevt.vhep[nhep][j-1]=parton->XProd()[j];
  hepevt.vhep[nhep][3]=parton->XProd()[0];
  hepevt.isthep[nhep]=1;
  hepevt.jmohep[nhep][0]=0;
  hepevt.jmohep[nhep][1]=0;
  hepevt.jdahep[nhep][0]=0;
  hepevt.jdahep[nhep][1]=0;
  nhep++;
}

void Lund_Interface::FillPrimaryHadronsInBlob(ATOOLS::Blob *blob,ATOOLS::Blob_List *bloblist,
					      ATOOLS::Particle_List *pl)
{
  m_secondarymap.clear();
//  ATOOLS::Blob *decay;
  ATOOLS::Particle *particle;
  ATOOLS::Flavour flav;
  ATOOLS::Vec4D momentum, position;
  for (int i=0;i<hepevt.nhep;++i) {
    if ((hepevt.isthep[i]!=2)&&(hepevt.isthep[i]!=1)&&(hepevt.isthep[i]!=149)) continue;
    if (hepevt.idhep[i]==93) flav=ATOOLS::Flavour(ATOOLS::kf::cluster);
    else flav.FromHepEvt(hepevt.idhep[i]);
	if ((flav==ATOOLS::Flavour(ATOOLS::kf::tau) ||
		 flav==ATOOLS::Flavour(ATOOLS::kf::tau).Bar()) &&
		hepevt.jmohep[i][0]==0 && hepevt.jmohep[i][1]==0) {
	  /*
		std::cout<<"HepEvt : "<<i<<" "<<hepevt.isthep[i]<<" "<<hepevt.idhep[i]
		<<" "<<hepevt.jmohep[i][0]<<" "<<hepevt.jmohep[i][1]
		<<" "<<hepevt.jdahep[i][0]<<" "<<hepevt.jdahep[i][1]
		<<" "<<flav<<" -> new method."<<std::endl;
	  */
	  FillPrimaryTauInBlob(i,blob,bloblist,pl);
	  continue;
	}
    if (flav==ATOOLS::Flavour(ATOOLS::kf::string) || 
		flav==ATOOLS::Flavour(ATOOLS::kf::cluster)) {
      for (int j=hepevt.jdahep[i][0]-1;j<hepevt.jdahep[i][1];j++) {
		flav.FromHepEvt(hepevt.idhep[j]);
		if (!flav.IsHadron()) continue;
		momentum=ATOOLS::Vec4D(hepevt.phep[j][3],hepevt.phep[j][0],
							   hepevt.phep[j][1],hepevt.phep[j][2]);
		position=ATOOLS::Vec4D(hepevt.vhep[j][3],hepevt.vhep[j][0],
							   hepevt.vhep[j][1],hepevt.vhep[j][2]);
		particle = new ATOOLS::Particle(-1,flav,momentum);
		if (pl) particle->SetNumber(pl->size());
		else particle->SetNumber(0);
		particle->SetStatus(1);
		particle->SetInfo('P');
		blob->SetPosition(position);
		if (pl) pl->push_back(particle);
		blob->AddToOutParticles(particle);
		if (hepevt.jdahep[j][0]!=0 && hepevt.jdahep[j][1]!=0) {
		  m_secondarymap[particle]=j;
		}
      }
    }
  }
  blob->SetStatus(0);
}

void Lund_Interface::FillPrimaryTauInBlob(int pos,ATOOLS::Blob *blob,
										  ATOOLS::Blob_List *bloblist,
										  ATOOLS::Particle_List *pl)
{
  ATOOLS::Particle *particle;
  ATOOLS::Flavour flav;
  ATOOLS::Vec4D momentum, position;

  // Reconstruct tau and push it through the blob
  flav.FromHepEvt(hepevt.idhep[pos]);
  momentum=ATOOLS::Vec4D(hepevt.phep[pos][3],hepevt.phep[pos][0],
						 hepevt.phep[pos][1],hepevt.phep[pos][2]);
  position=ATOOLS::Vec4D(hepevt.vhep[pos][3],hepevt.vhep[pos][0],
							   hepevt.vhep[pos][1],hepevt.vhep[pos][2]);
  particle = new ATOOLS::Particle(-1,flav,momentum);
  if (pl) particle->SetNumber(pl->size());
     else particle->SetNumber(0);
  particle->SetStatus(1);
  particle->SetInfo('P');
  blob->SetPosition(position);
  if (pl) pl->push_back(particle);
  blob->AddToOutParticles(particle);
  if (hepevt.jdahep[pos][0]!=0 && hepevt.jdahep[pos][1]!=0) {
	m_secondarymap[particle]=pos;
  }
  blob->SetStatus(0);
}

bool Lund_Interface::FindDecay( ATOOLS::Particle * part )
{
  if (m_secondarymap.find(part)==m_secondarymap.end()) return false; 
  return true;
}  

void Lund_Interface::PerformDecay( ATOOLS::Particle * part, 
	                               ATOOLS::Blob_List * blob_list, 
								   ATOOLS::Particle_List * part_list )
{
  // choose decay channel and kinematics
  int pos = m_secondarymap[part];					// find particle in secondary map
  int daughter1 = hepevt.jdahep[pos][0]-1,
      daughter2 = hepevt.jdahep[pos][1];
  const int n = daughter2-daughter1+1;

  // create blob
  ATOOLS::Blob * blob;											// decay blob
  blob = new ATOOLS::Blob();
  blob->SetStatus(0);
  blob->SetType( ATOOLS::btp::Hadron_Decay );
  blob->SetTypeSpec( "Pythia_v6.214" );
  blob->SetId();
  blob->AddToInParticles( part );
  if( part->Info() == 'P' ) part->SetInfo('p');
  if( part->Info() == 'D' ) part->SetInfo('d');
  blob_list->push_back( blob );
  ATOOLS::Particle * particle;									// daughter part.
  ATOOLS::Vec4D momentum;										// daughter mom.
  ATOOLS::Vec4D position;										// daughter pos.
  ATOOLS::Flavour flav;											// daughter flav.

  // treat every daughter
  for (int i=daughter1;i<daughter2;++i) {
    flav.FromHepEvt( hepevt.idhep[i] );
    momentum=ATOOLS::Vec4D( hepevt.phep[i][3],
		                    hepevt.phep[i][0],
			                hepevt.phep[i][1],
							hepevt.phep[i][2] );
    position=ATOOLS::Vec4D( hepevt.vhep[i][3],
		                    hepevt.vhep[i][0],
			                hepevt.vhep[i][1],
							hepevt.vhep[i][2] );
	particle = new ATOOLS::Particle( -1, flav, momentum );
	if( part_list ) particle->SetNumber( part_list->size() );
	else particle->SetNumber( 0 );
	particle->SetStatus(1);
	particle->SetInfo('D');
	blob->SetPosition( position );
	if( part_list ) part_list->push_back( particle ); 
	blob->AddToOutParticles( particle );
	// check if daughter can be treated as well
	if( FindDecay(particle) ) {
	  blob->SetStatus(1);
	  PerformDecay( particle, blob_list, part_list );
	}
  }
}
  
bool Lund_Interface::FillDecay(ATOOLS::Particle * part,ATOOLS::Blob_List *bloblist,
							   ATOOLS::Particle_List *pl)
{
  msg_Tracking()<<"Lund_Interface::FillDecay()"<<endl;
  ATOOLS::Blob *decay;
  ATOOLS::Particle *particle;
  if (m_secondarymap.find(part)==m_secondarymap.end()) return; 
  int pos = m_secondarymap[part];
  decay = new ATOOLS::Blob();
  decay->SetStatus(1);
  decay->SetType(ATOOLS::btp::Hadron_Decay);
  decay->SetTypeSpec("Pythia_v6.214");
  decay->SetId();
  decay->AddToInParticles(particle);
  if (particle->Info()=='P') particle->SetInfo('p');
  if (particle->Info()=='D') particle->SetInfo('d');
  particle->SetStatus(2);
  bloblist->push_back(decay);
  FillSecondaryHadronsInBlob(decay,bloblist,hepevt.jdahep[pos][0]-1,hepevt.jdahep[pos][1],pl);
}


void Lund_Interface::FillSecondaryHadronsInBlob(ATOOLS::Blob *blob,ATOOLS::Blob_List *bloblist,
												int daughter1,int daughter2,ATOOLS::Particle_List *pl) 
{
  ATOOLS::Blob *decay;
  ATOOLS::Particle *particle;
  ATOOLS::Flavour flav;
  ATOOLS::Vec4D momentum, position;
  for (int i=daughter1;i<daughter2;++i) {
    flav.FromHepEvt(hepevt.idhep[i]);
    momentum=ATOOLS::Vec4D(hepevt.phep[i][3],hepevt.phep[i][0],
			   hepevt.phep[i][1],hepevt.phep[i][2]);
    position=ATOOLS::Vec4D(hepevt.vhep[i][3],hepevt.vhep[i][0],
			   hepevt.vhep[i][1],hepevt.vhep[i][2]);
    particle = new ATOOLS::Particle(-1,flav,momentum);
    if (pl) particle->SetNumber(pl->size());
       else particle->SetNumber(0);
    particle->SetStatus(1);
    particle->SetInfo('D');
    blob->SetPosition(position);
    if (pl) pl->push_back(particle);
    blob->AddToOutParticles(particle);
    if (hepevt.jdahep[i][0]!=0 && hepevt.jdahep[i][1]!=0) {
      decay = new ATOOLS::Blob();
      decay->SetStatus(1);
      decay->SetType(ATOOLS::btp::Hadron_Decay);
      decay->SetTypeSpec("Pythia_v6.214");
      decay->SetId();
      decay->AddToInParticles(particle);
      if (particle->Info()=='P') particle->SetInfo('p');
      if (particle->Info()=='D') particle->SetInfo('d');
      particle->SetStatus(2);
      bloblist->push_back(decay);
	  FillSecondaryHadronsInBlob(decay,bloblist,hepevt.jdahep[i][0]-1,hepevt.jdahep[i][1],pl);
    }
  }
}


void Lund_Interface::Error(const int error)
{
  ++s_errors;
  if (s_errors>s_maxerrors) {
    THROW(critical_error,"Pythia calls PYERRM("+
	  ATOOLS::ToString(error)+")");
  }
  else {
    ATOOLS::msg.Error()<<"Lund_Interface::Error("<<error<<") "<<ATOOLS::om::red
		       <<"Pythia calls PYERRM("<<error<<") in event "
		       <<ATOOLS::rpa.gen.NumberOfDicedEvents()<<"."
		       <<ATOOLS::om::reset<<std::endl;
    if (ATOOLS::msg.LevelIsDebugging()) {
      msg_Tracking()<<*s_bloblist<<std::endl;
      pylist(2);
    }
  }
}

void Lund_Interface::NextFile(const bool newfile) 
{
  if (!m_writeout) return; 
  std::string oldfile;
  bool oldfileexists=false;
  std::ofstream *outfile=p_hepevt->GetOutStream();
  if (outfile!=NULL) {
    oldfileexists=true;
    oldfile=m_outfile+ATOOLS::ToString(m_curfile)+std::string(".evts");
    if (newfile) 
      (*outfile)<<(m_outfile+ATOOLS::ToString(++m_curfile)+std::string(".evts"))<<std::endl;
    if (m_compress) {
      system((std::string("gzip ")+oldfile+std::string(".gz ")+oldfile).c_str());
      system((std::string("rm ")+oldfile).c_str());
    }
  }
  if (!newfile) {
    if (p_hepevt) { 
      p_hepevt->SetNhep(0);
      p_hepevt->SetIsthep(NULL);
      p_hepevt->SetIdhep(NULL);
      p_hepevt->SetJmohep(NULL);
      p_hepevt->SetJdahep(NULL);
      p_hepevt->SetPhep(NULL);
      p_hepevt->SetVhep(NULL);
      delete p_hepevt; 
      p_hepevt = NULL; 
    }
    return;
  }
  std::string file = std::string(m_outfile+ATOOLS::ToString(m_curfile)+std::string(".evts"));
  p_hepevt->ChangeOutStream(file,m_evtsperfile);
}

bool Lund_Interface::OneEvent(ATOOLS::Blob_List * const blobs,double &weight)
{
  bool okay = false;
  for (int i=0;i<200;i++) {
    pyevnt();
    pyhepc(1);
    weight=1.;  //*=pypars.pari[10];
    for (int i=0;i<hepevt.nhep;i++) {
      for (int j=0;j<2;j++) {
	p_jmohep[2*i+j] = hepevt.jmohep[i][j]; 
	p_jdahep[2*i+j] = hepevt.jdahep[i][j];
      } 
      for (int j=0;j<5;j++) p_phep[5*i+j] = hepevt.phep[i][j];
      for (int j=0;j<4;j++) p_vhep[4*i+j] = hepevt.vhep[i][j];
    }
    p_hepevt->SetNhep(hepevt.nhep);
    p_hepevt->SetIsthep(hepevt.isthep);
    p_hepevt->SetIdhep(hepevt.idhep);
    p_hepevt->SetJmohep(p_jmohep);
    p_hepevt->SetJdahep(p_jdahep);
    p_hepevt->SetPhep(p_phep);
    p_hepevt->SetVhep(p_vhep);
    if (ATOOLS::msg.LevelIsDebugging()) pylist(3);
    if (p_hepevt->HepEvt2Sherpa(blobs)) { 
      okay = true; 
      break; 
    }
  }
  return okay;
} 
