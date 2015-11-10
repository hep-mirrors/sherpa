#include "AddOns/Root/Output_RootNtuple.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Phys/NLO_Subevt.H"
#include "PHASIC++/Process/Process_Base.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Message.H"
#include "MODEL/Main/Model_Base.H"

#include <limits>
#include <string.h>

#ifdef USING__ROOT
#include "TPluginManager.h"
#include "TROOT.h"
#endif

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

#ifdef USING__MPI
MPI_Datatype MPI_rntuple_evt2;
MPI_Datatype MPI_Vec4D;
#endif

Output_RootNtuple::Output_RootNtuple
(const Output_Arguments &args,int exact,int ftype):
  Output_Base("Root"), m_exact(exact), m_ftype(ftype)
{
  args.p_reader->SetAllowUnits(true);
  m_filesize = args.p_reader->GetValue<long int>("NTUPLE_SIZE",2147483647);
  args.p_reader->SetAllowUnits(false);
  m_mode=args.p_reader->GetValue<int>("ROOTNTUPLE_MODE",0);
  m_treename=args.p_reader->GetValue<std::string>("ROOTNTUPLE_TREENAME","t3");
  m_basename =args.m_outpath+"/"+args.m_outfile;
  m_ext = ".root";
  m_cnt2=m_cnt3=m_fcnt=m_evt=0;
  m_idcnt=0;
  m_avsize=args.p_reader->GetValue<int>("ROOTNTUPLE_AVSIZE",10000);
#ifdef USING__MPI
  int size=MPI::COMM_WORLD.Get_size();
  if (m_exact) m_avsize=Max((size_t)1,m_avsize/size);
#endif
  m_total=0;
  m_csumsqr=m_csum=m_cn=0.;
  m_sum=m_s2=m_s3=m_c1=m_c2=0.;
  m_sq=m_fsq=m_sq2=m_sq3=0.;
#ifdef USING__ROOT
  p_t3=NULL;
  p_f=NULL;
#ifdef USING__MPI
  rntuple_evt2 dummye[1];
  MPI_Datatype typee[22]={MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,
			  MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,
			  MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,
			  MPI_DOUBLE,MPI_DOUBLE,
			  MPI_LONG,MPI_INT,MPI_INT,
			  MPI_INT,MPI_INT,MPI_INT,MPI_INT,
			  MPI_INT,MPI_DOUBLE,MPI_INT,MPI_CHAR};
  int blocklene[22]={1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,18,1,2};
  MPI_Aint basee, dispe[22];
  MPI_Get_address(&dummye[0],&basee);
  MPI_Get_address(&dummye[0].weight,dispe+0);
  MPI_Get_address(&dummye[0].wgt0,dispe+1);
  MPI_Get_address(&dummye[0].x1,dispe+2);
  MPI_Get_address(&dummye[0].x2,dispe+3);
  MPI_Get_address(&dummye[0].y1,dispe+4);
  MPI_Get_address(&dummye[0].y2,dispe+5);
  MPI_Get_address(&dummye[0].fscale,dispe+6);
  MPI_Get_address(&dummye[0].rscale,dispe+7);
  MPI_Get_address(&dummye[0].alphas,dispe+8);
  MPI_Get_address(&dummye[0].kfac,dispe+9);
  MPI_Get_address(&dummye[0].psw,dispe+10);
  MPI_Get_address(&dummye[0].id,dispe+11);
  MPI_Get_address(&dummye[0].ncount,dispe+12);
  MPI_Get_address(&dummye[0].nparticle,dispe+13);
  MPI_Get_address(&dummye[0].kf1,dispe+14);
  MPI_Get_address(&dummye[0].kf2,dispe+15);
  MPI_Get_address(&dummye[0].f1,dispe+16);
  MPI_Get_address(&dummye[0].f2,dispe+17);
  MPI_Get_address(&dummye[0].nuwgt,dispe+18);
  MPI_Get_address(dummye[0].uwgt,dispe+19);
  MPI_Get_address(&dummye[0].oqcd,dispe+20);
  MPI_Get_address(dummye[0].type,dispe+21);
  for (size_t i=0;i<22;++i) dispe[i]-=basee;
  MPI_Type_create_struct(22,blocklene,dispe,typee,&MPI_rntuple_evt2);
  MPI_Type_commit(&MPI_rntuple_evt2);
  Vec4D dummyv[1];
  MPI_Datatype typev[1]={MPI_DOUBLE};
  int blocklenv[1]={4};
  MPI_Aint basev, dispv[1];
  MPI_Get_address(&dummyv[0],&basev);
  MPI_Get_address(&dummyv[0][0],dispv+0);
  for (size_t i=0;i<1;++i) dispv[i]-=basev;
  MPI_Type_create_struct(1,blocklenv,dispv,typev,&MPI_Vec4D);
  MPI_Type_commit(&MPI_Vec4D);
#endif
#endif
}

void Output_RootNtuple::Header()
{
#ifdef USING__ROOT
#ifdef USING__MPI
  int rank=MPI::COMM_WORLD.Get_rank();
  if (rank) return;
#endif
  p_f=new TFile((m_basename+m_ext).c_str(),"recreate");
  p_t3 = new TTree(m_treename.c_str(),"Reconst ntuple");
  p_t3->SetMaxTreeSize(m_filesize);
  p_t3->Branch("id",&m_id,"id/I");
  if (m_exact) p_t3->Branch("ncount",&m_ncount,"ncount/I");
  p_t3->Branch("nparticle",&m_nparticle,"nparticle/I");
  if (m_ftype) {
    p_t3->Branch("px",p_pxd,"px[nparticle]/D");
    p_t3->Branch("py",p_pyd,"py[nparticle]/D");
    p_t3->Branch("pz",p_pzd,"pz[nparticle]/D");
    p_t3->Branch("E",p_Ed,"E[nparticle]/D");
  }
  else {
    p_t3->Branch("px",p_px,"px[nparticle]/F");
    p_t3->Branch("py",p_py,"py[nparticle]/F");
    p_t3->Branch("pz",p_pz,"pz[nparticle]/F");
    p_t3->Branch("E",p_E,"E[nparticle]/F");
  }
  p_t3->Branch("alphas",&m_alphas,"alphas/D");
  if (m_exact) p_t3->Branch("kfactor",&m_kfac,"kfactor/D");
  if (m_exact) p_t3->Branch("ps_wgt",&m_psw,"ps_wgt/D");
  p_t3->Branch("kf",p_kf,"kf[nparticle]/I");
  p_t3->Branch("weight",&m_wgt,"weight/D");
  p_t3->Branch("weight2",&m_wgt2,"weight2/D");
  p_t3->Branch("me_wgt",&m_mewgt,"me_wtg/D");
  p_t3->Branch("me_wgt2",&m_mewgt2,"me_wtg2/D");
  p_t3->Branch("x1",&m_x1,"x1/D");
  p_t3->Branch("x2",&m_x2,"x2/D");
  p_t3->Branch("x1p",&m_y1,"x1p/D");
  p_t3->Branch("x2p",&m_y2,"x2p/D");
  p_t3->Branch("id1",&m_id1,"id1/I");
  p_t3->Branch("id2",&m_id2,"id2/I");
  if (m_exact) p_t3->Branch("id1p",&m_id1p,"id1p/I");
  if (m_exact) p_t3->Branch("id2p",&m_id2p,"id2p/I");
  p_t3->Branch("fac_scale",&m_fscale,"fac_scale/D");
  p_t3->Branch("ren_scale",&m_rscale,"ren_scale/D");
  p_t3->Branch("nuwgt",&m_nuwgt,"nuwgt/I");
  p_t3->Branch("usr_wgts",p_uwgt,"usr_wgts[nuwgt]/D");
  p_t3->Branch("alphasPower",&m_oqcd,"alphasPower/B");
  p_t3->Branch("part",m_type,"part[2]/C");
  ATOOLS::exh->AddTerminatorObject(this);
  gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo","*",
					"TStreamerInfo","RIO","TStreamerInfo()"); 
#endif
}

Output_RootNtuple::~Output_RootNtuple()
{
  PrepareTerminate();
#ifdef USING__MPI
  MPI_Type_free(&MPI_rntuple_evt2);
  MPI_Type_free(&MPI_Vec4D);
#endif
}

void Output_RootNtuple::Footer()
{
  PrepareTerminate();
}

void Output_RootNtuple::PrepareTerminate()
{
  StoreEvt();
#ifdef USING__ROOT
  if (p_t3==NULL) return;
  p_t3->AutoSave();
  delete p_t3;
  p_t3=NULL;
  // delete p_f;
  ATOOLS::exh->RemoveTerminatorObject(this);
#endif
  if (m_exact) {
    double xs=m_csum/m_cn;
    double err=sqrt((m_csumsqr/m_cn-sqr(m_csum/m_cn))/(m_cn-1.));
    msg_Info()<<METHOD<<"(): '"<<p_f->GetName()
	      <<"' stores "<<xs<<" pb +- ( "<<err<<" pb = "
	      <<dabs((int(10000.0*err/xs))/100.0)<<" % )\n";
    return;
  }
  msg_Info()<<"ROOTNTUPLE_OUTPUT stored: "<<m_s2/m_c2<<" +/- "<<sqrt((m_sq2/m_c2-sqr(m_s2/m_c2))/(m_c2-1.))<<" pb  (reweighted 1) \n"; 
  double c3(m_idcnt);
  msg_Info()<<"                          "<<m_s3/c3<<" +/- "<<sqrt((m_sq3/c3-sqr(m_s3/c3))/(c3-1.))<<" pb  (reweighted 2) \n"; 
  msg_Info()<<"                          "<<m_sum/m_c1<<" +/- "<<sqrt((m_sq/m_c1-sqr(m_sum/m_c1))/(m_c1-1.))<<" pb  (before reweighting) \n"<<endl; 
}

void Output_RootNtuple::AddDecayProducts(Particle *part,int &np) 
{
  DEBUG_FUNC(*part);
  if (part->DecayBlob()) {
    for (size_t i(0);i<part->DecayBlob()->NOutP();++i)
      AddDecayProducts(part->DecayBlob()->OutParticle(i),np);
    return;
  }
  DEBUG_VAR("Adding "<<*part);
  ++np;
  int kfc=part->Flav().Kfcode(); 
  if (part->Flav().IsAnti()) kfc=-kfc;
  if (m_fcnt>=m_flavlist.size()) {
    m_flavlist.resize(m_flavlist.size()+3*m_avsize);
    m_momlist.resize(m_momlist.size()+3*m_avsize);
  }
  m_flavlist[m_fcnt]=kfc;
  m_momlist[m_fcnt]=part->Momentum();
  ++m_fcnt;
}

void Output_RootNtuple::Output(Blob_List* blobs, const double weight) 
{
  Blob* signal=0, *shower=0;
  for (Blob_List::const_iterator blit=blobs->begin();blit!=blobs->end();++blit) 
    if ((*blit)->Type()==btp::Signal_Process) {
      signal=(*blit);
    }
    else if ((*blit)->Type()==btp::Shower) {
      shower=(*blit);
    }
  if (!signal || (m_mode==1 && !shower)) return;
  int ncount=(*signal)["Trials"]->Get<double>();
  m_evt+=ncount;
  m_c1+=ncount;
  m_cn+=ncount;
  m_cnt3++;
  m_idcnt++;

  Blob *blob=signal;
  if (m_mode==1) blob=shower;
  Blob_Data_Base* seinfo=(*signal)["MEWeightInfo"];
  ME_Weight_Info* wgtinfo(NULL);
  if (seinfo) wgtinfo=seinfo->Get<ME_Weight_Info*>();
  seinfo=(*signal)["NLO_subeventlist"];
  std::string type((*signal)["NLOQCDType"]->Get<std::string>());
  
  if (!seinfo) { // BVI type events
    if (m_evtlist.size()<=m_cnt2)
      m_evtlist.resize(m_evtlist.size()+m_avsize);
    m_evtlist[m_cnt2].weight=(*signal)["Weight"]->Get<double>();
    m_evtlist[m_cnt2].ncount=ncount;
    m_sum+=m_evtlist[m_cnt2].weight;
    m_csum+=m_evtlist[m_cnt2].weight;
    m_csumsqr+=sqr(m_evtlist[m_cnt2].weight);
    m_fsq+=sqr(m_evtlist[m_cnt2].weight);
    m_evtlist[m_cnt2].id=m_idcnt;
    m_evtlist[m_cnt2].fscale=sqrt((*signal)["Factorisation_Scale"]->Get<double>());
    m_evtlist[m_cnt2].rscale=sqrt((*signal)["Renormalization_Scale"]->Get<double>());
    m_evtlist[m_cnt2].alphas=MODEL::s_model->ScalarFunction("alpha_S",m_evtlist[m_cnt2].rscale*m_evtlist[m_cnt2].rscale);
    m_evtlist[m_cnt2].kfac=wgtinfo->m_K;
    m_evtlist[m_cnt2].oqcd=(*signal)["Orders"]->Get<std::vector<double> >()[0];
    if (type=="B" || type=="") strcpy(m_evtlist[m_cnt2].type,"B");
    else if (type=="V") strcpy(m_evtlist[m_cnt2].type,"V");
    else if (type=="I") strcpy(m_evtlist[m_cnt2].type,"I");
    else THROW(fatal_error,"Error in NLO type '"+type+"'");
    if (wgtinfo) {
      if      (type=="B" || type=="") {
	m_evtlist[m_cnt2].wgt0=wgtinfo->m_B;
	m_evtlist[m_cnt2].psw=wgtinfo->m_B/(*signal)["MEWeight"]->Get<double>();
      }
      else if (type=="V" || type=="I") m_evtlist[m_cnt2].wgt0=wgtinfo->m_VI;
      m_evtlist[m_cnt2].x1=wgtinfo->m_x1;
      m_evtlist[m_cnt2].x2=wgtinfo->m_x2;
      m_evtlist[m_cnt2].y1=wgtinfo->m_y1;
      m_evtlist[m_cnt2].y2=wgtinfo->m_y2;
      size_t nren(wgtinfo->m_type&mewgttype::VI?wgtinfo->m_wren.size():0),
             nfac(wgtinfo->m_type&mewgttype::KP?wgtinfo->m_wfac.size():0);
      m_evtlist[m_cnt2].nuwgt=nren+nfac;
      if (wgtinfo->m_type&mewgttype::VI) {
        for (size_t i(0);i<nren;++i)
          m_evtlist[m_cnt2].uwgt[i]=wgtinfo->m_wren[i];
      }
      if (wgtinfo->m_type&mewgttype::KP) {
        for (size_t i(0);i<nfac;++i)
          m_evtlist[m_cnt2].uwgt[nren+i]=wgtinfo->m_wfac[i];
      }
    }
    for (int inp=0, i=0;i<blob->NInP();i++) {
    Particle *part=blob->InParticle(i);
    if (part->ProductionBlob() &&
	part->ProductionBlob()->Type()==btp::Signal_Process) continue;
    int kfc=part->Flav().Kfcode(); if (part->Flav().IsAnti()) kfc=-kfc;
    if (++inp==1) m_evtlist[m_cnt2].f1=m_evtlist[m_cnt2].kf1=kfc;
    else m_evtlist[m_cnt2].f2=m_evtlist[m_cnt2].kf2=kfc;
    }
    int np=0;
    for (int i=0;i<blob->NOutP();i++) {
      Particle *part=blob->OutParticle(i);
      if (part->DecayBlob()) {
	if (m_mode==1 && part->DecayBlob()->Type()==btp::Signal_Process) continue;
	AddDecayProducts(part,np);
	continue;
      }
      ++np;
      int kfc=part->Flav().Kfcode(); 
      if (part->Flav().IsAnti()) kfc=-kfc;
      if (m_fcnt>=m_flavlist.size()) {
	m_flavlist.resize(m_flavlist.size()+3*m_avsize);
	m_momlist.resize(m_momlist.size()+3*m_avsize);
      }
      m_flavlist[m_fcnt]=kfc;
      m_momlist[m_fcnt]=part->Momentum();
      ++m_fcnt;
    }
    m_evtlist[m_cnt2].nparticle=np;
    ++m_cnt2;
  }
  else { // RS-type events
    NLO_subevtlist* nlos = seinfo->Get<NLO_subevtlist*>();
    double tweight=0.;
    for (size_t j=0;j<nlos->size();j++) {
      if ((*nlos)[j]->m_result==0.0) continue;
      if (m_evtlist.size()<=m_cnt2)
	m_evtlist.resize(m_evtlist.size()+m_avsize);
      ATOOLS::Particle_List * pl=(*nlos)[j]->CreateParticleList();
      m_evtlist[m_cnt2].weight=(*nlos)[j]->m_result;
      m_evtlist[m_cnt2].ncount=ncount;
      tweight+=m_evtlist[m_cnt2].weight;
      m_evtlist[m_cnt2].nparticle=pl->size();
      m_evtlist[m_cnt2].id=m_idcnt;
      m_evtlist[m_cnt2].wgt0=(*nlos)[j]->m_mewgt;
      m_evtlist[m_cnt2].fscale=sqrt((*nlos)[j]->m_mu2[stp::fac]);
      m_evtlist[m_cnt2].rscale=sqrt((*nlos)[j]->m_mu2[stp::ren]);
      m_evtlist[m_cnt2].alphas=MODEL::s_model->ScalarFunction("alpha_S", m_evtlist[m_cnt2].rscale*m_evtlist[m_cnt2].rscale);
      m_evtlist[m_cnt2].kfac=(*nlos)[j]->m_K;
      m_evtlist[m_cnt2].psw=m_evtlist[m_cnt2].weight/(*signal)["MEWeight"]->Get<double>();
      m_evtlist[m_cnt2].oqcd=(*signal)["Orders"]->Get<std::vector<double> >()[0];
      if (type!="RS") THROW(fatal_error,"Error in NLO type");
      if ((*nlos)[j]->p_real==(*nlos)[j]) strcpy(m_evtlist[m_cnt2].type,"R");
      else                                strcpy(m_evtlist[m_cnt2].type,"S");

      if (wgtinfo) {
	m_evtlist[m_cnt2].x1=wgtinfo->m_x1;
	m_evtlist[m_cnt2].x2=wgtinfo->m_x2;
	m_evtlist[m_cnt2].y1=wgtinfo->m_y1;
	m_evtlist[m_cnt2].y2=wgtinfo->m_y2;
      }
      m_evtlist[m_cnt2].nuwgt=0;

      Particle* part=signal->InParticle(0);
      int kfc=part->Flav().Kfcode(); if (part->Flav().IsAnti()) kfc=-kfc;
      m_evtlist[m_cnt2].kf1=kfc;
      m_evtlist[m_cnt2].f1=(long int)(*nlos)[j]->p_fl[0];
      part=signal->InParticle(1);
      kfc=part->Flav().Kfcode(); if (part->Flav().IsAnti()) kfc=-kfc;
      m_evtlist[m_cnt2].kf2=kfc;
      m_evtlist[m_cnt2].f2=(long int)(*nlos)[j]->p_fl[1];

      ++m_cnt2;
      for (ATOOLS::Particle_List::const_iterator pit=pl->begin();
	   pit!=pl->end();++pit) {
	kfc=(*pit)->Flav().Kfcode(); 
	if ((*pit)->Flav().IsAnti()) kfc=-kfc;	  
	if (m_fcnt>=m_flavlist.size()) {
	  m_flavlist.resize(m_flavlist.size()+3*m_avsize);
	  m_momlist.resize(m_momlist.size()+3*m_avsize);
	}	
	m_flavlist[m_fcnt]=kfc;
	m_momlist[m_fcnt]=(*pit)->Momentum();
	++m_fcnt;
	delete *pit;
      }      
      delete pl;
    }
    m_sum+=tweight;
    m_csum+=tweight;
    m_csumsqr+=sqr(tweight);
    m_fsq+=sqr(tweight);
  }
  
  if ((rpa->gen.NumberOfGeneratedEvents()%m_avsize)==0) StoreEvt();
}

void Output_RootNtuple::ChangeFile()
{
  StoreEvt();
#ifdef USING__ROOT
  if (p_t3==NULL) return;
  double xs=m_csum/m_cn;
  double err=sqrt((m_csumsqr/m_cn-sqr(m_csum/m_cn))/(m_cn-1.));
  msg_Info()<<METHOD<<"(): '"<<p_f->GetName()
	    <<"' stores "<<xs<<" pb +- ( "<<err<<" pb = "
	    <<dabs((int(10000.0*err/xs))/100.0)<<" % )\n";
  m_csumsqr=m_csum=m_cn=0.0;
  p_f=p_t3->ChangeFile(p_f);
#endif
}

void Output_RootNtuple::MPISync()
{
#ifdef USING__MPI
  static int s_offset=11;
  int size=MPI::COMM_WORLD.Get_size();
  if (size>1) {
    int rank=MPI::COMM_WORLD.Get_rank();
    double vals[6];
    if (rank==0) {
      m_evtlist.resize(m_cnt2);
      m_flavlist.resize(m_fcnt);
      m_momlist.resize(m_fcnt);
      for (int tag=1;tag<size;++tag) {
	MPI::COMM_WORLD.Recv(&vals,6,MPI::DOUBLE,MPI::ANY_SOURCE,s_offset*size+tag);
	std::vector<rntuple_evt2> evts((int)vals[0]);
	std::vector<int> flavs((int)vals[3]);
	std::vector<Vec4D> moms((int)vals[3]);
	MPI::COMM_WORLD.Recv(&evts.front(),(int)vals[0],MPI_rntuple_evt2,MPI::ANY_SOURCE,(s_offset+1)*size+tag);
	MPI::COMM_WORLD.Recv(&flavs.front(),(int)vals[3],MPI::INT,MPI::ANY_SOURCE,(s_offset+2)*size+tag);
	MPI::COMM_WORLD.Recv(&moms.front(),(int)vals[3],MPI_Vec4D,MPI::ANY_SOURCE,(s_offset+3)*size+tag);
	m_evtlist.insert(m_evtlist.end(),evts.begin(),evts.end());
	m_flavlist.insert(m_flavlist.end(),flavs.begin(),flavs.end());
	m_momlist.insert(m_momlist.end(),moms.begin(),moms.end());
	int oid=-1;
	for (size_t i(m_cnt2);i<m_cnt2+evts.size();++i) {
	  if (m_evtlist[i].id!=oid) {
	    oid=m_evtlist[i].id;
	    ++m_idcnt;
	  }
	  m_evtlist[i].id=m_idcnt;
	}
	m_cnt2+=vals[0];
	m_cnt3+=vals[1];
	m_evt+=vals[2];
	m_fcnt+=vals[3];
	m_fsq+=vals[4];
	m_sum+=vals[5];
	m_c1+=vals[2];
	m_csumsqr+=vals[4];
	m_csum+=vals[5];
	m_cn+=vals[2];
      }
    }
    else {
      vals[0]=m_cnt2;
      vals[1]=m_cnt3;
      vals[2]=m_evt;
      vals[3]=m_fcnt;
      vals[4]=m_fsq;
      vals[5]=m_sum;
      MPI::COMM_WORLD.Send(&vals,6,MPI::DOUBLE,0,s_offset*size+rank);
      MPI::COMM_WORLD.Send(&m_evtlist.front(),(int)vals[0],MPI_rntuple_evt2,0,(s_offset+1)*size+rank);
      MPI::COMM_WORLD.Send(&m_flavlist.front(),(int)vals[3],MPI::INT,0,(s_offset+2)*size+rank);
      MPI::COMM_WORLD.Send(&m_momlist.front(),(int)vals[3],MPI_Vec4D,0,(s_offset+3)*size+rank);
      m_cnt2=m_cnt3=m_fcnt=m_evt=0;
      m_sum=m_fsq=0.0;
    }
  }
#endif
}

void Output_RootNtuple::StoreEvt()
{
  if (m_cnt2==0) return;
  MPISync();
#ifdef USING__ROOT
  if (p_t3==NULL) return;
#endif
  size_t fc=0;
  double scale2=double(m_cnt2)/double(m_evt);
  double scale3=double(m_cnt3)/double(m_evt);
  if (m_exact) scale2=scale3=1.0;
  for (size_t i=0;i<m_cnt2;i++) {
#ifdef USING__ROOT
    m_id  = m_evtlist[i].id;
    m_wgt = m_evtlist[i].weight*scale2;
    m_wgt2= m_evtlist[i].weight*scale3;
    m_mewgt = m_evtlist[i].wgt0*scale2;
    m_mewgt2= m_evtlist[i].wgt0*scale3;
    m_ncount= m_evtlist[i].ncount;
    m_x1 = m_evtlist[i].x1;
    m_x2 = m_evtlist[i].x2;
    m_y1 = m_evtlist[i].y1;
    m_y2 = m_evtlist[i].y2;
    m_id1 = m_evtlist[i].kf1;
    m_id2 = m_evtlist[i].kf2;
    m_id1p = m_evtlist[i].f1;
    m_id2p = m_evtlist[i].f2;
    m_nuwgt = m_evtlist[i].nuwgt;
    for (int j=0;j<m_nuwgt;j++)
      p_uwgt[j]=m_evtlist[i].uwgt[j]*scale2;

    m_fscale = m_evtlist[i].fscale;
    m_rscale = m_evtlist[i].rscale;
  
    m_nparticle=m_evtlist[i].nparticle;
    m_alphas=m_evtlist[i].alphas;
    m_kfac=m_evtlist[i].kfac;
    m_psw=m_evtlist[i].psw;
    m_oqcd=m_evtlist[i].oqcd;
    strcpy(m_type,m_evtlist[i].type);
    for (size_t j=0;j<m_evtlist[i].nparticle;j++) {
      p_kf[j] = m_flavlist[fc];
      p_E[j]  = p_Ed[j]  = m_momlist[fc][0];
      p_px[j] = p_pxd[j] = m_momlist[fc][1];
      p_py[j] = p_pyd[j] = m_momlist[fc][2];
      p_pz[j] = p_pzd[j] = m_momlist[fc][3];
      fc++;
    }
    p_t3->Fill();
#endif
    m_s2+=m_evtlist[i].weight*scale2;
    m_sq2+=sqr(m_evtlist[i].weight*scale2);
    m_s3+=m_evtlist[i].weight*scale3;
    m_c2+=1.;
  }
  m_sq+=m_fsq;
  m_sq3+=m_fsq*sqr(scale3);
  m_cnt2=m_cnt3=m_fcnt=m_evt=0;
  m_fsq=0.;
}

DECLARE_GETTER(Output_RootNtuple,"Root",
	       Output_Base,Output_Arguments);

Output_Base *ATOOLS::Getter<Output_Base,Output_Arguments,Output_RootNtuple>::
operator()(const Output_Arguments &args) const
{
  return new Output_RootNtuple(args);
}

void ATOOLS::Getter<Output_Base,Output_Arguments,Output_RootNtuple>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"Root NTuple output";
}

DECLARE_GETTER(Output_ERootNtuple,"ERoot",
	       Output_Base,Output_Arguments);

Output_Base *ATOOLS::Getter<Output_Base,Output_Arguments,Output_ERootNtuple>::
operator()(const Output_Arguments &args) const
{
  return new Output_RootNtuple(args,1);
}

void ATOOLS::Getter<Output_Base,Output_Arguments,Output_ERootNtuple>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"Root NTuple output";
}

DECLARE_GETTER(Output_EDRootNtuple,"EDRoot",
	       Output_Base,Output_Arguments);

Output_Base *ATOOLS::Getter<Output_Base,Output_Arguments,Output_EDRootNtuple>::
operator()(const Output_Arguments &args) const
{
  return new Output_RootNtuple(args,1,1);
}

void ATOOLS::Getter<Output_Base,Output_Arguments,Output_EDRootNtuple>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"Root NTuple output";
}

