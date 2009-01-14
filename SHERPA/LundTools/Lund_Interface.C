#include "Lund_Interface.H"

#include "Lund_Wrapper.H"
#include "ISR_Handler.H"
#include "Data_Reader.H"
#include "Particle.H"
#include "Blob.H"
#include "Blob_List.H"
#include "Run_Parameter.H"
#include "Random.H"
#include "Message.H"
#include "Exception.H"
#include "Running_AlphaS.H"
#include "MyStrStream.H"
#include <list>
#include <cassert>
#include "Message.H"
#include "HepEvt_Interface.H"
#include "Mass_Handler.H"

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

bool Lund_Interface::s_exportas=false;
bool Lund_Interface::s_exportpdf=false;

size_t Lund_Interface::s_errors=0;
size_t Lund_Interface::s_maxerrors=0;

int * Lund_Interface::s_saved_mrpy=new int[6];
double * Lund_Interface::s_saved_rrpy=new double[100];

ATOOLS::Blob_List *Lund_Interface::s_bloblist=NULL; 
PDF::ISR_Handler *Lund_Interface::s_isrhandler=NULL; 

Lund_Interface::Lund_Interface(string _m_path,string _m_file,bool sherpa):
  m_path(_m_path),m_file(_m_file), m_maxtrials(2),
  p_hepevt(NULL), 
  m_compress(true),m_writeout(false),
  p_phep(new double[5*10000]),
  p_vhep(new double[4*10000]),
  p_jmohep(new int[2*10000]),
  p_jdahep(new int[2*10000])
{
  exh->AddTerminatorObject(this);
  double win;
  string beam[2], frame("CMS");
  Flavour flav[2];
  for (size_t i=0;i<2;++i) flav[i]=rpa.gen.Bunch(i);
  if (flav[0]==kf_e && flav[1]==kf_p_plus) {
    beam[0]="e-";
    beam[1]="p+";
  }
  else if (flav[0]==kf_p_plus && flav[1]==kf_e) {
    beam[0]="p+";
    beam[1]="e-";
  }
  else if (flav[0]==kf_e && flav[1]==kf_photon) {
    if (flav[0].IsAnti()) beam[0]="e+"; else beam[0]="e-";
    beam[1]="gamma";
    spsubs.msub[33]=1;    
  }
  else if (flav[0]==kf_photon && flav[1]==kf_e) {
    beam[0]="gamma";
    if (flav[1].IsAnti()) beam[1]="e+"; else beam[1]="e-";
    spsubs.msub[33]=1;    
  }
  else if (flav[0]==kf_photon && flav[1]==kf_photon) {
    for (size_t i=0;i<2;++i) beam[i]="gamma";
    spsubs.msub[57]=1;    
  }
  else if (flav[0].Kfcode()==kf_e && flav[1].Kfcode()==kf_e) {
    for (size_t i=0;i<2;++i) if (flav[i].IsAnti()) beam[i]="e+"; else beam[i]="e-";
    spsubs.msub[0]=1;    
    sppars.mstp[47]=1;
    spdat1.mstj[100]=5;
  }
  else {
    for (size_t i=0;i<2;++i) if (flav[i].IsAnti()) beam[i]="p-"; else beam[i]="p+";
    spsubs.msub[0]=1;    
    sppars.mstp[47]=1;
    spdat1.mstj[100]=5;
  }
  win=rpa.gen.Ecms();
  s_maxerrors=rpa.gen.NumberOfEvents();
  vector<vector<double> > help;
  Data_Reader *reader = new Data_Reader(" ",";","!","=");
  reader->AddWordSeparator("\t");
  reader->SetInputPath(m_path);
  reader->SetInputFile(m_file);
  reader->AddIgnore("(");
  reader->AddIgnore(")");
  reader->AddIgnore(",");
  if (!sherpa) {
    if (!reader->ReadFromFile(spsubs.msel,"MSEL")) spsubs.msel=1;
  }
  else {
    spsubs.msel=0;
  }
  reader->MatrixFromFile(help,"MSUB");
  for (size_t i=0;i<help.size();++i) {
    if (help[i].size()>1) if ((int)help[i][0]>0) spsubs.msub[(int)help[i][0]-1]=(int)help[i][1];
  }
  reader->MatrixFromFile(help,"KFIN");
  for (size_t i=0;i<help.size();++i) {
    if (help[i].size()>2) {
      if (((int)help[i][0]>0)&&((int)help[i][1]>-41)) {
	spsubs.kfin[(int)help[i][1]+40][(int)help[i][0]-1]=(int)help[i][2];
      }
    }
  }
  reader->MatrixFromFile(help,"CKIN");
  for (size_t i=0;i<help.size();++i) {
    if (help[i].size()>1) if ((int)help[i][0]>0) spsubs.ckin[(int)help[i][0]-1]=help[i][1];
  }
  reader->MatrixFromFile(help,"MSTJ");
  for (size_t i=0;i<help.size();++i) {
    if (help[i].size()>1) if ((int)help[i][0]>0) spdat1.mstj[(int)help[i][0]-1]=(int)help[i][1];
  }
  reader->MatrixFromFile(help,"MSTP");
  for (size_t i=0;i<help.size();++i) {
    if (help[i].size()>1) if ((int)help[i][0]>0) sppars.mstp[(int)help[i][0]-1]=(int)help[i][1];
  }
  reader->MatrixFromFile(help,"MSTU");
  for (size_t i=0;i<help.size();++i) {
    if (help[i].size()>1) if ((int)help[i][0]>0) spdat1.mstu[(int)help[i][0]-1]=(int)help[i][1];
  }
  reader->MatrixFromFile(help,"MFUDGE");
  for (size_t i=0;i<help.size();++i) {
    if (help[i].size()>1) if ((int)help[i][0]>0) cfudge.mfudge[(int)help[i][0]-1]=(int)help[i][1];
  }
  reader->MatrixFromFile(help,"PARP");
  for (size_t i=0;i<help.size();++i) {
    if (help[i].size()>1) if ((int)help[i][0]>0) sppars.parp[(int)help[i][0]-1]=help[i][1];
  }
  reader->MatrixFromFile(help,"PARJ");
  bool found21(false), found41(false), found42(false);
  for (size_t i=0;i<help.size();++i) {
    if (help[i].size()>1) 
      if ((int)help[i][0]>0) {
        spdat1.parj[(int)help[i][0]-1]=help[i][1];
        if((int)help[i][0]==21) found21=true;
        if((int)help[i][0]==41) found41=true;
        if((int)help[i][0]==42) found42=true;
      }
  }
  if(!found21) spdat1.parj[21-1]=0.36; // sigma
  if(!found41) spdat1.parj[41-1]=0.25; // a
  if(!found42) spdat1.parj[42-1]=0.75; // b
  reader->MatrixFromFile(help,"PARU");
  for (size_t i=0;i<help.size();++i) {
    if (help[i].size()>1) if ((int)help[i][0]>0) spdat1.paru[(int)help[i][0]-1]=help[i][1];
  }
  reader->MatrixFromFile(help,"PFUDGE");
  for (size_t i=0;i<help.size();++i) {
    if (help[i].size()>1) if ((int)help[i][0]>0) cfudge.pfudge[(int)help[i][0]-1]=help[i][1];
  }
  reader->MatrixFromFile(help,"MDME");
  for (size_t i=0;i<help.size();++i) {
    if (help[i].size()>2) {
      if ((int)help[i][0]>0 && abs((int)help[i][1]<2)) {
	spdat3.mdme[(int)help[i][1]-1][(int)help[i][0]-1]=(int)help[i][2];
      }
    }
  }
  reader->MatrixFromFile(help,"MDCYKF");
  for (size_t i=0;i<help.size();++i) {
    if (help[i].size()>2) {
      if ((int)help[i][0]>0 && abs((int)help[i][1]<2)) {
	msg_Tracking()<<"Lund_Interface::Lund_Interface(..): "
		      <<"Set MDCY("<<spcomp((int)help[i][0])<<","<<(int)help[i][1]
		      <<") ( from KF code "<<(int)help[i][0]<<" ) to "<<(int)help[i][2]<<endl;
	spdat3.mdcy[(int)help[i][1]-1][spcomp((int)help[i][0])-1]=(int)help[i][2];
      }
    }
  }
  // the next lines replace the aspinit_ call
  if (sherpa) {
    hepevt.nhep=100;
    for (int i=spdat3.mdcy[2-1][23-1];i<spdat3.mdcy[2-1][23-1]+spdat3.mdcy[3-1][23-1];++i) {
      if (abs(spdat3.kfdp[1-1][i-1])>=2) spdat3.mdme[1-1][i-1]=Min(0,spdat3.mdme[1-1][i-1]);
    }
    spinit(frame.c_str(),beam[0].c_str(),beam[1].c_str(),100.0);
  }
  // replacement ends here
  if (msg_LevelIsDebugging()) ListLundParameters();
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
    MODEL::as = new MODEL::Running_AlphaS(asmz,mz*mz,orderas);
    MODEL::as->SetDefault(asdef);
    p_hepevt = new HepEvt_Interface(gtp::Pythia);
    if (sppars.mstp[105-1]==0) p_hepevt->SetHadronized(false);
    spinit(frame.c_str(),beam[0].c_str(),beam[1].c_str(),win);
    if (reader->ReadFromFile(m_outfile,"OUTPUT_FILE")) {
      m_writeout = true;
      if (!reader->ReadFromFile(m_evtsperfile,"EVENTS_PER_FILE")) m_evtsperfile=1000;
      NextFile(true);
      int helper;
      if (!reader->ReadFromFile(helper,"COMPRESS")) helper=1;
      m_compress=(bool)helper;
    }
  }
  splist(0);
  delete reader;
}

bool Lund_Interface::IsAllowedDecay(kf_code can)
{
  if (can==kf_tau) return false;
  if (spcomp(int(can))<501 && spdat3.mdcy[1-1][spcomp(int(can))-1]==1) return true;
  return false;
}

void Lund_Interface::SwitchOffDecays(kf_code kfc)
{
  int kc = spcomp(int(kfc));
  if(kc>500) return;
  spdat3.mdcy[1-1][kc-1]=0;
}

void Lund_Interface::AdjustProperties(Flavour flav)
{
  int kc = spcomp(int(flav.Kfcode()));
  if(kc>500) return;
  // adjust mass
  double pythiamass = spdat2.pmas[1-1][kc-1];
  double sherpamass = flav.PSMass();
  flav.SetMass(pythiamass);
  if( !(abs(sherpamass-pythiamass)/sherpamass < 1.e-2) ) {
    msg_Info()<<METHOD<<" Adjusted mass of "<<flav<<" ("<<flav.Kfcode()
        <<") from "<<sherpamass<<" to "<<pythiamass<<" to allow Pythia decays."<<endl;
  }
}

void Lund_Interface::SwitchOffMassSmearing()
{
  spdat1.mstj[24-1]=0;
}

double Lund_Interface::DiceMass(Flavour flav, double min, double max)
{
  int kc = spcomp(flav.Kfcode())-1;
  double peak = spdat2.pmas[1-1][kc];
  double w_cut = spdat2.pmas[3-1][kc];
  if(w_cut == 0.0) {
    if(peak<min-Accu() || peak>max+Accu()) return -1.0;
    else return peak;
  }
  else {
    Mass_Handler masshandler(flav);
    double finalmin = min>(peak-w_cut)? min : peak-w_cut;
    double finalmax = max<(peak+w_cut)? max : peak+w_cut;
    if(finalmin>finalmax) return -1.0;
    else return masshandler.GetMass(finalmin,finalmax);
  }
}

Lund_Interface::~Lund_Interface()
{
  exh->RemoveTerminatorObject(this);
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
  if (p_jmohep) { delete [] p_jmohep; p_jmohep = NULL; }
  if (p_jdahep) { delete [] p_jdahep; p_jdahep = NULL; }
  if (p_phep)   { delete [] p_phep;   p_phep   = NULL; }
  if (p_vhep)   { delete [] p_vhep;   p_vhep   = NULL; }
}

Return_Value::code Lund_Interface::Hadronize(Blob_List *bloblist) 
{
  int nhep(0);
  for (Blob_List::iterator blit=bloblist->begin();blit!=bloblist->end();++blit) {
    if ((*blit)->Has(blob_status::needs_hadronization)) {
      if ((*blit)->Type()!=btp::Fragmentation) {
	msg_Error()<<"Error in "<<METHOD<<":"<<std::endl
		   <<"   Blob with status=needs_hadronizazion and "<<std::endl
		   <<"   type different from 'Fragmentation' illegally found."<<std::endl
		   <<"   Will retry the event and hope for the best."<<std::endl;
	return Return_Value::Retry_Event;
      }
      (*blit)->SetTypeSpec("Pythia_v6.214");

      Return_Value::code result=Return_Value::Nothing;
      for(size_t trials=0; trials<m_maxtrials; ++trials) {
	nhep = PrepareFragmentationBlob(*blit);
	int errs=spdat1.mstu[23-1];
	result=StringFragmentation(*blit,bloblist,nhep);
	if(result==Return_Value::Success) break;
	assert(result==Return_Value::Retry_Phase);
	if(trials+1<m_maxtrials) {
	  msg_Error()<<"Error in "<<METHOD<<"."<<endl
		     <<"   Hadronization failed. Retry it."<<endl;
	  spdat1.mstu[23-1]=errs;   //New try, set back error sum.
	  continue;
	}
	msg_Error()<<"Error in "<<METHOD<<"."<<endl
		   <<"   Hadronization failed. Retry the event."<<endl;
	return Return_Value::Retry_Event;
      }
      //Up to this point fragmentation blob remained unchanged.
      //Therefore, could be used as the safe initial state.
      FillFragmentationBlob(*blit);
    }
  }
  return Return_Value::Success;
}

Return_Value::code Lund_Interface::PerformDecay(Blob * blob)
{
  if (!blob->Has(blob_status::needs_hadrondecays) ||
      blob->NInP()!=1 ||
      blob->InParticle(0)->Status()!=part_status::active)
  {
    msg_Error()<<METHOD<<" returns Error."<<endl;
    msg_Error()<<" blob->Status()="<<blob->Status()<<endl;
    msg_Error()<<" blob->NInP()="<<blob->NInP()<<endl;
    msg_Error()<<" part->Status()="<<blob->InParticle(0)->Status()<<endl;
    return Return_Value::Error;
  }

  Particle * part = blob->InParticle(0);
  Flavour fl = part->Flav();
  int kc = spcomp(int(fl.Kfcode()))-1;
  double peak = spdat2.pmas[1-1][kc];
  double w_cut = spdat2.pmas[3-1][kc];
  if( part->FinalMass()+Accu() < peak-w_cut || part->FinalMass()-Accu() > peak+w_cut) {
    return Return_Value::Retry_Method;
  }

  int nhep(0);
  int idhep = fl.HepEvt();
  switch (idhep) {
  case  9000111: idhep =  10111;   break;
  case  9000211: idhep =  10211;   break;
  case -9000211: idhep = -10211;   break;
  case  9010221: idhep =  10221;   break;
  case    10111: idhep =  9000111; break;
  case    10211: idhep =  9000211; break;
  case   -10211: idhep = -9000211; break;
  case    10221: idhep =  10331; break;
  case    10331: idhep =  9010221; break;
  }
  hepevt.idhep[nhep] = idhep;
  for (short int j=1;j<4;++j) hepevt.phep[nhep][j-1]=part->Momentum()[j];
  hepevt.phep[nhep][3] = part->Momentum()[0];
  hepevt.phep[nhep][4] = part->FinalMass();
  for (short int j=0;j<4;++j) hepevt.vhep[nhep][j]=0.0;
  hepevt.isthep[nhep]=1;
  hepevt.jmohep[nhep][0]=0;
  hepevt.jmohep[nhep][1]=0;
  hepevt.jdahep[nhep][0]=0;
  hepevt.jdahep[nhep][1]=0;
  nhep++;
  
  hepevt.nevhep=0;
  hepevt.nhep=nhep;
  sphepc(2);
  spdat1.mstu[70-1]=1;
  spdat1.mstu[71-1]=hepevt.nhep;
  int ip(1);
  spdecy(ip);
  if (spdat1.mstu[24-1]!=0) {
    msg_Error()<<"ERROR in "<<METHOD<<" : "<<std::endl
	       <<"   SPDECY call results in error code : "<<spdat1.mstu[24-1]<<std::endl
	       <<"   for decay of "<<fl<<" ("<<fl.HepEvt()<<" -> "<<idhep<<")"<<std::endl;
    if (spdat1.mstu[23-1]<int(rpa.gen.NumberOfDicedEvents()/100) ||
	rpa.gen.NumberOfDicedEvents()<200) {
      msg_Error()<<"   Up to now: "<<spdat1.mstu[23-1]<<" errors, try new event."<<std::endl;
      return Return_Value::Retry_Method;
    }
    msg_Error()<<"   Up to now: "<<spdat1.mstu[23-1]<<" errors, abort the run."<<std::endl;
    THROW(critical_error,"Too many errors in lund decay.");
  }
  part->SetStatus(part_status::decayed);
  sphepc(1);
  FillOutgoingParticlesInBlob(blob);
  return Return_Value::Success;
} 

int Lund_Interface::PrepareFragmentationBlob(Blob * blob) 
{
  int nhep = 0;
  hepevt.idhep[nhep]=Flavour(kf_photon).HepEvt();
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
  
  // gluon splittings
  for (int i(0);i<blob->NInP();++i) {
    Particle * part = blob->InParticle(i);
  if (part->GetFlow(1)!=0 && part->GetFlow(2)!=0) {
    Flavour            flav = Flavour(kf_d);
    if (ran.Get()<0.5) flav = Flavour(kf_u);
      Particle *help1(new Particle(-1,flav,0.5*part->Momentum()));
      Particle *help2(new Particle(-1,flav.Bar(),help1->Momentum()));
    help1->SetStatus(part_status::active);
    help2->SetStatus(part_status::active);
    AddPartonToString(help1,nhep);
    delete help1;
      unsigned int lastc(part->GetFlow(2));
      for (++i;i<blob->NInP();++i) {
      part = blob->InParticle(i);
      AddPartonToString(part,nhep);
	if (part->GetFlow(1)==lastc) {
	  lastc=0;
	  break;
	}
    }      
      if (lastc!=0)
	msg_Error()<<METHOD<<"(): Error. Open color string."<<std::endl;
    AddPartonToString(help2,nhep);
    delete help2;
      lastc=0;
  }
  else {
      for (;i<blob->NInP();i++) {
      part = blob->InParticle(i);
      AddPartonToString(part,nhep);
	if (part->GetFlow(1)==0) break;
      }  
    }  
  }
  return nhep;
}


Return_Value::code Lund_Interface::StringFragmentation(Blob *blob,Blob_List *bloblist,int nhep) 
{
  hepevt.nevhep=0;
  hepevt.nhep=nhep;
  sphepc(2);
  spdat1.mstu[70-1]=1;
  spdat1.mstu[71-1]=hepevt.nhep;
  int ip(1);
  spprep(ip);
  spstrf(ip);
  sphepc(1);
  bool goon(true), flag(false);
  if(hepevt.nhep<=nhep) { goon=false; flag=true;}
  while(goon) {
    goon=false;
    for(int i=0; i<hepevt.nhep; ++i) {
      ATOOLS::Flavour flav;
      flav.FromHepEvt(hepevt.idhep[i]);
      if(hepevt.isthep[i]==1 && (flav.Strong()||flav.IsDiQuark())) {
	goon=true;
	int ipp(i+1);
	int save(hepevt.nhep);
	spstrf(ipp);
	sphepc(1);
	if(hepevt.nhep<=save) { goon=false; flag=true;}
	break;
      }
    }
  }
  if(spdat1.mstu[24-1]!=0) {
    Vec4D cms(0.,0.,0.,0.);
    for (int i=0;i<blob->NInP();i++) cms+=blob->InParticle(i)->Momentum();
    msg_Error()<<"ERROR in "<<METHOD<<" : "<<std::endl
	       <<"   SPSTRF call results in error code : "<<spdat1.mstu[24-1]
	       <<" for "<<std::endl<<(*blob)<<"  "<<cms<<", "<<cms.Abs2()
	       <<std::endl;
    spdat1.mstu[24-1]=0;
    if(spdat1.mstu[23-1]<int(rpa.gen.NumberOfDicedEvents()/100) ||
       rpa.gen.NumberOfDicedEvents()<200) {
      msg_Error()<<"   Up to now: "<<spdat1.mstu[23-1]<<" errors, retrying..."
		 <<std::endl;
      return Return_Value::Retry_Phase;
    }
    msg_Error()<<"   Up to now: "<<spdat1.mstu[23-1]<<" errors, abort the run."
	       <<std::endl;
    THROW(critical_error,"Too many errors in lund fragmentation.");
  }
  if(flag) {
    msg_Error()<<"ERROR in "<<METHOD<<" : "<<std::endl
	       <<"   Incomplete fragmentation."<<std::endl;
    return Return_Value::Retry_Phase;
  }
  spdat1.mstu[70-1]=2;
  spdat1.mstu[72-1]=hepevt.nhep;
  return Return_Value::Success;
}

void Lund_Interface::AddPartonToString(Particle *parton,int &nhep)
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

void Lund_Interface::FillFragmentationBlob(Blob *blob)
{
  Particle *particle;
  Flavour flav;
  Vec4D momentum, position;
  for (int i=0;i<hepevt.nhep;++i) {
    if ((hepevt.isthep[i]!=2)&&(hepevt.isthep[i]!=1)&&(hepevt.isthep[i]!=149)) continue;
    if (hepevt.idhep[i]==93) flav=Flavour(kf_cluster);
                        else flav.FromHepEvt(hepevt.idhep[i]);
    if (flav==Flavour(kf_string) || 
	flav==Flavour(kf_cluster)) {
      for (int j=hepevt.jdahep[i][0]-1;j<hepevt.jdahep[i][1];j++) {
        flav.FromHepEvt(hepevt.idhep[j]);
	momentum=Vec4D(hepevt.phep[j][3],hepevt.phep[j][0],
		       hepevt.phep[j][1],hepevt.phep[j][2]);
	position=Vec4D(hepevt.vhep[j][3],hepevt.vhep[j][0],
		       hepevt.vhep[j][1],hepevt.vhep[j][2]);
	particle = new Particle(-1,flav,momentum);
	particle->SetNumber(0);
	particle->SetStatus(part_status::active);
	particle->SetInfo('P');
	particle->SetFinalMass(hepevt.phep[j][4]);
	blob->SetPosition(position);
	blob->AddToOutParticles(particle);
      }
    }
  }
  blob->SetStatus(blob_status::needs_hadrondecays);
}

void Lund_Interface::FillOutgoingParticlesInBlob(Blob *blob)
{
  Flavour    flav;
  Vec4D      momentum;
  Particle * particle;

  int n_q(0), n_g(0); Particle_Vector partons;
  for (int j=hepevt.jdahep[0][0]-1;j<hepevt.jdahep[0][1];j++) {
    int idhep = hepevt.idhep[j];
    switch (idhep) {
    case  9000111: idhep =  10111;   break;
    case  9000211: idhep =  10211;   break;
    case -9000211: idhep = -10211;   break;
    case  9010221: idhep =  10331;   break;
    case    10111: idhep =  9000111; break;
    case    10211: idhep =  9000211; break;
    case   -10211: idhep = -9000211; break;
    case    10221: idhep =  9010221; break;
    case    10331: idhep =  10221; break;
    }
    flav.FromHepEvt(idhep);
    momentum=Vec4D(hepevt.phep[j][3],hepevt.phep[j][0],
		   hepevt.phep[j][1],hepevt.phep[j][2]);
    // don't fill blob position vector here, because decay is in CMS until boosting back
    particle = new Particle(-1,flav,momentum);
    particle->SetNumber(0);
    particle->SetStatus(part_status::active);
    particle->SetInfo('P');
    particle->SetFinalMass(hepevt.phep[j][4]);
    blob->AddToOutParticles(particle);
    
    if(abs(idhep)>0 && abs(idhep)<7) {
      n_q++;
      partons.push_back(particle);
    }
    else if(abs(idhep)==21) {
      n_g++;
      partons.push_back(particle);
    }
  }
  
  size_t n=partons.size();
  if(n>0) blob->SetStatus(blob_status::needs_showers);
  if(n_q==2 && n_g==0 && n==2) {
    if(partons[0]->Flav().IsAnti()) {
      partons[0]->SetFlow(2,-1);
      partons[1]->SetFlow(1,partons[0]->GetFlow(2));
    }
    else {
      partons[0]->SetFlow(1,-1);
      partons[1]->SetFlow(2,partons[0]->GetFlow(1));
    }
  }
  else if(n_q==0 && n_g==2 && n==2) {
    partons[0]->SetFlow(2,-1);
    partons[0]->SetFlow(1,-1);
    partons[1]->SetFlow(2,partons[0]->GetFlow(1));
    partons[1]->SetFlow(1,partons[0]->GetFlow(2));
  }
  else if(n_q==0 && n_g==3 && n==3) {
    partons[0]->SetFlow(2,-1);
    partons[0]->SetFlow(1,-1);
    partons[1]->SetFlow(2,partons[0]->GetFlow(1));
    partons[1]->SetFlow(1,-1);
    partons[2]->SetFlow(2,partons[1]->GetFlow(1));
    partons[2]->SetFlow(1,partons[0]->GetFlow(2));
  }
  else if(n>0) {
    msg_Error()<<METHOD<<" wasn't able to set the color flow for"<<endl<<*blob<<endl;
  }
}

void Lund_Interface::RestoreStatus() {
  for(int i=0;i<5;i++) {
    spdatr.mrsp[i] = s_saved_mrpy[i];
  }
  for(int i=0;i<100;i++) {
    spdatr.rrsp[i] = s_saved_rrpy[i];
  }
}

void Lund_Interface::SaveStatus() {
  for(int i=0;i<5;i++) {
    s_saved_mrpy[i] = spdatr.mrsp[i];
  }
  for(int i=0;i<100;i++) {
    s_saved_rrpy[i] = spdatr.rrsp[i];
  }
}

bool Lund_Interface::ReadInStatus(const std::string &path)
{
  ReadInStatus((path+"Lund_random.dat").c_str(),0);
  return true;
}

void Lund_Interface::ReadInStatus(const std::string &filename, int mode) {
  ifstream myinstream(filename.c_str());
  if (myinstream.good()) {
    for(int i=0;i<5;i++) {
      myinstream>>spdatr.mrsp[i];
    }
    for(int i=0;i<100;i++) {
      myinstream>>spdatr.rrsp[i];
    }
    myinstream.close();
  }
  else msg_Error()<<"ERROR in "<<METHOD<<": "<<filename<<" not found!!"<<endl;
}

void Lund_Interface::WriteOutStatus(const std::string &filename)
{
  ofstream myoutstream(filename.c_str());
  if (myoutstream.good()) {
    myoutstream.precision(32);
    for(int i=0;i<5;i++) {
      myoutstream<<spdatr.mrsp[i]<<"\t";
    }
    for(int i=0;i<100;i++) {
      myoutstream<<spdatr.rrsp[i]<<"\t";
    }
    myoutstream<<endl;
    myoutstream.close();
  }
  else msg_Error()<<"ERROR in "<<METHOD<<": "<<filename<<" not found!!"<<endl;
}

void Lund_Interface::PrepareTerminate()
{
  std::string path(rpa.gen.Variable("SHERPA_STATUS_PATH"));
  if (path=="") return;
  RestoreStatus();
  WriteOutStatus((path+"/Lund_random.dat").c_str());
}

void Lund_Interface::Error(const int error)
{
  ++s_errors;
  if (s_errors>s_maxerrors) {
    THROW(critical_error,"Pythia calls PYERRM("+
	  ToString(error)+")");
  }
  else {
    msg_Error()<<"Lund_Interface::Error("<<error<<") "<<om::red
	       <<"Pythia calls PYERRM("<<error<<") in event "
	       <<rpa.gen.NumberOfDicedEvents()<<"."
	       <<om::reset<<endl;
//     if (msg_LevelIsDebugging()) {
//       msg_Tracking()<<*s_bloblist<<endl;
      splist(2);
//     }
  }
}

void Lund_Interface::NextFile(const bool newfile) 
{
  if (!m_writeout) return; 
  string oldfile;
  bool oldfileexists=false;
  ofstream *outfile=p_hepevt->GetOutStream();
  if (outfile!=NULL) {
    oldfileexists=true;
    oldfile=m_outfile+ToString(m_curfile)+string(".evts");
    if (newfile) 
      (*outfile)<<(m_outfile+ToString(++m_curfile)+string(".evts"))<<endl;
    if (m_compress) {
      system((string("gzip ")+oldfile+string(".gz ")+oldfile).c_str());
      system((string("rm ")+oldfile).c_str());
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
  string file = string(m_outfile+ToString(m_curfile)+string(".evts"));
  p_hepevt->ChangeOutStream(file,m_evtsperfile);
}

Return_Value::code Lund_Interface::OneEvent(Blob_List * const blobs,double &weight)
{
  bool okay = false;
  for (int i=0;i<200;i++) {
    spevnt();
    sphepc(1);
    //pylist(2);
    weight=1.;  //*=sppars.pari[10];
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
    if (msg_LevelIsDebugging()) splist(3);
    if (p_hepevt->HepEvt2Sherpa(blobs)) { 
      okay = true; 
      break; 
    }
  }
  if (okay) return Return_Value::Success;
  return Return_Value::Nothing;
} 
