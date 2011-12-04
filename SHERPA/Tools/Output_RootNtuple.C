#include "SHERPA/Tools/Output_RootNtuple.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Phys/NLO_Subevt.H"
#include "PHASIC++/Process/Process_Base.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Message.H"
#include "MODEL/Main/Model_Base.H"

#ifdef USING__MPI
#include "mpi.h"
#endif

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

#ifdef USING__MPI
MPI_Datatype MPI_rntuple_evt2;
MPI_Datatype MPI_Vec4D;
#endif

Output_RootNtuple::Output_RootNtuple(std::string basename,std::string ext,
				     int precision,long int filesize)
{
  Data_Reader dr(" ",";","!","=");
  dr.AddComment("#");
  dr.AddWordSeparator("\t");
  m_mode=dr.GetValue<int>("ROOTNTUPLE_MODE",0);
  m_basename =basename;
  m_ext = ext;
  m_cnt2=m_cnt3=m_fcnt=m_evt=0;
  m_idcnt=0;
  m_avsize=10000;
#ifdef USING__MPI
  int size=MPI::COMM_WORLD.Get_size();
  int rank=MPI::COMM_WORLD.Get_rank();
#endif
  m_total=0;
  m_evtlist.resize(m_avsize+100);
  m_flavlist.resize(3*m_avsize);
  m_momlist.resize(3*m_avsize);
  m_sum=m_s2=m_s3=m_c1=m_c2=0.;
  m_sq=m_fsq=m_sq2=m_sq3=0.;
#ifdef USING__ROOT
  p_t3=NULL;
  p_f=NULL;
#ifdef USING__MPI
  rntuple_evt2 dummye[1];
  MPI_Datatype typee[15]={MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,
			  MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,
			  MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,
			  MPI_LONG,MPI_INT,MPI_INT,MPI_INT,
			  MPI_INT,MPI_DOUBLE};
  int blocklene[15]={1,1,1,1,1,1,1,1,1,1,1,1,1,1,18};
  MPI_Aint basee, dispe[15];
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
  MPI_Get_address(&dummye[0].id,dispe+9);
  MPI_Get_address(&dummye[0].nparticle,dispe+10);
  MPI_Get_address(&dummye[0].kf1,dispe+11);
  MPI_Get_address(&dummye[0].kf2,dispe+12);
  MPI_Get_address(&dummye[0].nuwgt,dispe+13);
  MPI_Get_address(dummye[0].uwgt,dispe+14);
  for (size_t i=0;i<15;++i) dispe[i]-=basee;
  MPI_Type_create_struct(15,blocklene,dispe,typee,&MPI_rntuple_evt2);
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
  if (rank) return;
#endif
  p_f=new TFile((m_basename+m_ext).c_str(),"recreate");
  p_t3 = new TTree("t3","Reconst ntuple");
  p_t3->SetMaxTreeSize(Min(filesize,Number<long int>("100GB")));
  p_t3->Branch("id",&m_id,"id/I");
  p_t3->Branch("nparticle",&m_nparticle,"nparticle/I");
  p_t3->Branch("px",p_px,"px[nparticle]/F");
  p_t3->Branch("py",p_py,"py[nparticle]/F");
  p_t3->Branch("pz",p_pz,"pz[nparticle]/F");
  p_t3->Branch("E",p_E,"E[nparticle]/F");
  p_t3->Branch("alphas",&m_alphas,"alphas/D");
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
  p_t3->Branch("fac_scale",&m_fscale,"fac_scale/D");
  p_t3->Branch("ren_scale",&m_rscale,"ren_scale/D");

  p_t3->Branch("nuwgt",&m_nuwgt,"nuwgt/I");
  p_t3->Branch("usr_wgts",p_uwgt,"usr_wgts[nuwgt]/D");
  ATOOLS::exh->AddTerminatorObject(this);
#else
  msg_Error()<<"Sherpa must be linked with root to enable ROOTNTUPLE output!"<<endl;
#endif
}

Output_RootNtuple::~Output_RootNtuple()
{
  PrepareTerminate();
#ifdef USING__MPI
  exh->MPISync();
  MPI_Type_free(&MPI_rntuple_evt2);
  MPI_Type_free(&MPI_Vec4D);
#endif
}

void Output_RootNtuple::PrepareTerminate()
{
  StoreEvt();
#ifdef USING__ROOT
  if (p_t3==NULL) return;
  p_t3->AutoSave();
  delete p_t3;
  // delete p_f;
  ATOOLS::exh->RemoveTerminatorObject(this);
#endif
  cout<<"ROOTNTUPLE_OUTPUT stored: "<<m_s2/m_c2<<" +/- "<<sqrt((m_sq2/m_c2-sqr(m_s2/m_c2))/(m_c2-1.))<<" pb  (reweighted 1) \n"; 
  double c3(m_idcnt);
  cout<<"                          "<<m_s3/c3<<" +/- "<<sqrt((m_sq3/c3-sqr(m_s3/c3))/(c3-1.))<<" pb  (reweighted 2) \n"; 
  cout<<"                          "<<m_sum/m_c1<<" +/- "<<sqrt((m_sq/m_c1-sqr(m_sum/m_c1))/(m_c1-1.))<<" pb  (before reweighting) \n"<<endl; 
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
  m_cnt3++;
  m_idcnt++;

  Blob *blob=signal;
  if (m_mode==1) blob=shower;
  Blob_Data_Base* seinfo=(*signal)["ME_wgtinfo"];
  ME_wgtinfo* wgtinfo(0);
  if (seinfo) wgtinfo=seinfo->Get<ME_wgtinfo*>();
  seinfo=(*signal)["NLO_subeventlist"];
  
  if (!seinfo) {
    m_evtlist[m_cnt2].weight=(*signal)["Weight"]->Get<double>();
    m_sum+=m_evtlist[m_cnt2].weight;
    m_fsq+=sqr(m_evtlist[m_cnt2].weight);
    m_evtlist[m_cnt2].id=m_idcnt;
    m_evtlist[m_cnt2].fscale=sqrt((*signal)["Factorisation_Scale"]->Get<double>());
    m_evtlist[m_cnt2].rscale=sqrt((*signal)["Renormalization_Scale"]->Get<double>());
    m_evtlist[m_cnt2].alphas=MODEL::s_model->ScalarFunction("alpha_S",m_evtlist[m_cnt2].rscale*m_evtlist[m_cnt2].rscale);
    if (wgtinfo) {
      m_evtlist[m_cnt2].wgt0=wgtinfo->m_w0;
      m_evtlist[m_cnt2].x1=wgtinfo->m_x1;
      m_evtlist[m_cnt2].x2=wgtinfo->m_x2;
      m_evtlist[m_cnt2].y1=wgtinfo->m_y1;
      m_evtlist[m_cnt2].y2=wgtinfo->m_y2;
      m_evtlist[m_cnt2].nuwgt=wgtinfo->m_nx;
      for (int i=0;i<m_evtlist[m_cnt2].nuwgt;i++) 
	m_evtlist[m_cnt2].uwgt[i]=wgtinfo->p_wx[i];
    }
    for (int inp=0, i=0;i<blob->NInP();i++) {
    Particle *part=blob->InParticle(i);
    if (part->ProductionBlob() &&
	part->ProductionBlob()->Type()==btp::Signal_Process) continue;
    int kfc=part->Flav().Kfcode(); if (part->Flav().IsAnti()) kfc=-kfc;
    if (++inp==1) m_evtlist[m_cnt2].kf1=kfc;
    else m_evtlist[m_cnt2].kf2=kfc;
    }
    int np=0;
    for (int i=0;i<blob->NOutP();i++) {
      Particle *part=blob->OutParticle(i);
      if (part->DecayBlob() &&
	  part->DecayBlob()->Type()==btp::Signal_Process) continue;
      ++np;
      int kfc=part->Flav().Kfcode(); 
      if (part->Flav().IsAnti()) kfc=-kfc;
      m_flavlist[m_fcnt]=kfc;
      m_momlist[m_fcnt]=part->Momentum();
      ++m_fcnt;
      if (m_fcnt>=m_flavlist.size()) {
	m_flavlist.resize(m_flavlist.size()+m_avsize);
	m_momlist.resize(m_momlist.size()+m_avsize);
      }
    }
    m_evtlist[m_cnt2].nparticle=np;
    ++m_cnt2;
  }
  else {
    NLO_subevtlist* nlos = seinfo->Get<NLO_subevtlist*>();
    double tweight=0.;
    for (size_t j=0;j<nlos->size();j++) {
      if ((*nlos)[j]->m_result==0.0) continue;
      ATOOLS::Particle_List * pl=(*nlos)[j]->CreateParticleList();
      m_evtlist[m_cnt2].weight=(*nlos)[j]->m_result;
      tweight+=m_evtlist[m_cnt2].weight;
      m_evtlist[m_cnt2].nparticle=pl->size();
      m_evtlist[m_cnt2].id=m_idcnt;
      m_evtlist[m_cnt2].wgt0=(*nlos)[j]->m_mewgt;
      m_evtlist[m_cnt2].fscale=sqrt((*nlos)[j]->m_mu2[stp::fac]);
      m_evtlist[m_cnt2].rscale=sqrt((*nlos)[j]->m_mu2[stp::ren]);
      m_evtlist[m_cnt2].alphas=MODEL::s_model->ScalarFunction("alpha_S", m_evtlist[m_cnt2].rscale*m_evtlist[m_cnt2].rscale);

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
      part=signal->InParticle(1);
      kfc=part->Flav().Kfcode(); if (part->Flav().IsAnti()) kfc=-kfc;
      m_evtlist[m_cnt2].kf2=kfc;

      ++m_cnt2;
      for (ATOOLS::Particle_List::const_iterator pit=pl->begin();
	   pit!=pl->end();++pit) {
	kfc=(*pit)->Flav().Kfcode(); 
	if ((*pit)->Flav().IsAnti()) kfc=-kfc;	  
	m_flavlist[m_fcnt]=kfc;
	m_momlist[m_fcnt]=(*pit)->Momentum();
	++m_fcnt;
	if (m_fcnt>=m_flavlist.size()) {
	  m_flavlist.resize(m_flavlist.size()+m_avsize);
	  m_momlist.resize(m_momlist.size()+m_avsize);
	}	
	delete *pit;
      }      
      delete pl;
    }
    m_sum+=tweight;
    m_fsq+=sqr(tweight);
  }
  
  if (m_cnt2>=m_avsize) StoreEvt();
}

void Output_RootNtuple::ChangeFile(std::string number)
{
  StoreEvt();
#ifdef USING__ROOT
  p_t3->AutoSave();
#endif
}

void Output_RootNtuple::MPISync()
{
#ifdef USING__MPI
  int size=MPI::COMM_WORLD.Get_size();
  if (size>1) {
    int rank=MPI::COMM_WORLD.Get_rank();
    double vals[6];
    if (rank==0) {
      m_evtlist.resize(m_cnt2);
      m_flavlist.resize(m_fcnt);
      m_momlist.resize(m_fcnt);
      for (int tag=1;tag<size;++tag) {
	if (!exh->MPIStat(tag)) continue;
	MPI::COMM_WORLD.Recv(&vals,6,MPI::DOUBLE,MPI::ANY_SOURCE,tag);
	MPI::COMM_WORLD.Send(&vals,1,MPI::DOUBLE,tag,size+tag);
	std::vector<rntuple_evt2> evts((int)vals[0]);
	std::vector<int> flavs((int)vals[3]);
	std::vector<Vec4D> moms((int)vals[3]);
	MPI::COMM_WORLD.Recv(&evts.front(),(int)vals[0],MPI_rntuple_evt2,MPI::ANY_SOURCE,tag);
	MPI::COMM_WORLD.Recv(&flavs.front(),(int)vals[3],MPI::INT,MPI::ANY_SOURCE,size+tag);
	MPI::COMM_WORLD.Recv(&moms.front(),(int)vals[3],MPI_Vec4D,MPI::ANY_SOURCE,2*size+tag);
	MPI::COMM_WORLD.Send(&vals,1,MPI::DOUBLE,tag,3*size+tag);
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
      }
    }
    else {
      vals[0]=m_cnt2;
      vals[1]=m_cnt3;
      vals[2]=m_evt;
      vals[3]=m_fcnt;
      vals[4]=m_fsq;
      vals[5]=m_sum;
      MPI::COMM_WORLD.Send(&vals,6,MPI::DOUBLE,0,rank);
      MPI::COMM_WORLD.Recv(&vals,1,MPI::DOUBLE,0,size+rank);
      MPI::COMM_WORLD.Send(&m_evtlist.front(),(int)vals[0],MPI_rntuple_evt2,0,rank);
      MPI::COMM_WORLD.Send(&m_flavlist.front(),(int)vals[3],MPI::INT,0,size+rank);
      MPI::COMM_WORLD.Send(&m_momlist.front(),(int)vals[3],MPI_Vec4D,0,2*size+rank);
      MPI::COMM_WORLD.Recv(&vals,1,MPI::DOUBLE,0,3*size+rank);
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
  for (size_t i=0;i<m_cnt2;i++) {
#ifdef USING__ROOT
    m_id  = m_evtlist[i].id;
    m_wgt = m_evtlist[i].weight*scale2;
    m_wgt2= m_evtlist[i].weight*scale3;
    m_mewgt = m_evtlist[i].wgt0*scale2;
    m_mewgt2= m_evtlist[i].wgt0*scale3;
    m_x1 = m_evtlist[i].x1;
    m_x2 = m_evtlist[i].x2;
    m_y1 = m_evtlist[i].y1;
    m_y2 = m_evtlist[i].y2;
    m_id1 = m_evtlist[i].kf1;
    m_id2 = m_evtlist[i].kf2;
    m_nuwgt = m_evtlist[i].nuwgt;
    for (int j=0;j<m_nuwgt;j++)
      p_uwgt[j]=m_evtlist[i].uwgt[j]*scale2;

    m_fscale = m_evtlist[i].fscale;
    m_rscale = m_evtlist[i].rscale;
  
    m_nparticle=m_evtlist[i].nparticle;
    m_alphas=m_evtlist[i].alphas;
    for (size_t j=0;j<m_evtlist[i].nparticle;j++) {
      p_kf[j] = m_flavlist[fc];
      p_E[j]  = m_momlist[fc][0];
      p_px[j] = m_momlist[fc][1];
      p_py[j] = m_momlist[fc][2];
      p_pz[j] = m_momlist[fc][3];
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
