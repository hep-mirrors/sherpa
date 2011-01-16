#include "SHERPA/Tools/Output_RootNtuple.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Phys/NLO_Subevt.H"
#include "PHASIC++/Process/Process_Base.H"
#include "ATOOLS/Org/Message.H"


using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Output_RootNtuple::Output_RootNtuple(std::string basename,std::string ext,int precision)
{
  m_basename =basename;
  m_ext = ext;
  m_cnt2=m_cnt3=m_fcnt=m_evt=0;
  m_idcnt=0;
  m_avsize=10000;
  m_total=0;
  m_evtlist.resize(m_avsize+100);
  m_flavlist.resize(3*m_avsize);
  m_momlist.resize(3*m_avsize);
  m_sum=m_s2=m_s3=m_c1=m_c2=0.;
  m_sq=m_fsq=m_sq2=m_sq3=0.;
#ifdef USING__ROOT
  p_f=new TFile((m_basename+m_ext).c_str(),"recreate");
  p_t3 = new TTree("t3","Reconst ntuple");
  p_t3->Branch("id",&m_id,"id/I");
  p_t3->Branch("nparticle",&m_nparticle,"nparticle/I");
  p_t3->Branch("px",p_px,"px[nparticle]/F");
  p_t3->Branch("py",p_py,"py[nparticle]/F");
  p_t3->Branch("pz",p_pz,"pz[nparticle]/F");
  p_t3->Branch("E",p_E,"E[nparticle]/F");

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
#else
  msg_Error()<<"Sherpa must be linked with root to enable ROOTNTUPLE output!"<<endl;
#endif
}

Output_RootNtuple::~Output_RootNtuple()
{
  StoreEvt();
#ifdef USING__ROOT
  p_t3->AutoSave();
  delete p_t3;
  delete p_f;
#endif
  cout<<"ROOTNTUPLE_OUTPUT stored: "<<m_s2/m_c2<<" +/- "<<sqrt((m_sq2/m_c2-sqr(m_s2/m_c2))/(m_c2-1.))<<" pb  (reweighted 1) \n"; 
  double c3(m_idcnt);
  cout<<"                          "<<m_s3/c3<<" +/- "<<sqrt((m_sq3/c3-sqr(m_s3/c3))/(c3-1.))<<" pb  (reweighted 2) \n"; 
  cout<<"                          "<<m_sum/m_c1<<" +/- "<<sqrt((m_sq/m_c1-sqr(m_sum/m_c1))/(m_c1-1.))<<" pb  (before reweighting) \n"<<endl; 
}

void Output_RootNtuple::Output(Blob_List* blobs, const double weight) 
{
  Blob* signal=0;
  for (Blob_List::const_iterator blit=blobs->begin();blit!=blobs->end();++blit) 
    if ((*blit)->Type()==btp::Signal_Process) {
      signal=(*blit);
      break;
    }
  if (!signal) return;
  int ncount=(*signal)["Trials"]->Get<double>();
  m_evt+=ncount;
  m_c1+=ncount;
  m_cnt3++;
  m_idcnt++;


  Blob_Data_Base* seinfo=(*signal)["ME_wgtinfo"];
  ME_wgtinfo* wgtinfo(0);
  if (seinfo) wgtinfo=seinfo->Get<ME_wgtinfo*>();
  seinfo=(*signal)["NLO_subeventlist"];
  
  if (!seinfo) {
    m_evtlist[m_cnt2].weight=(*signal)["Weight"]->Get<double>();
    m_sum+=m_evtlist[m_cnt2].weight;
    m_fsq+=sqr(m_evtlist[m_cnt2].weight);
    m_evtlist[m_cnt2].nparticle=signal->NOutP();
    m_evtlist[m_cnt2].id=m_idcnt;
    m_evtlist[m_cnt2].fscale=sqrt((*signal)["Factorisation_Scale"]->Get<double>());
    m_evtlist[m_cnt2].rscale=sqrt((*signal)["Renormalization_Scale"]->Get<double>());

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

    Particle* part=signal->InParticle(0);
    int kfc=part->Flav().Kfcode(); if (part->Flav().IsAnti()) kfc=-kfc;
    m_evtlist[m_cnt2].kf1=kfc;
    part=signal->InParticle(1);
    kfc=part->Flav().Kfcode(); if (part->Flav().IsAnti()) kfc=-kfc;
    m_evtlist[m_cnt2].kf2=kfc;
 
    ++m_cnt2;
    for (int i=0;i<signal->NOutP();i++) {
      part=signal->OutParticle(i);
      kfc=part->Flav().Kfcode(); 
      if (part->Flav().IsAnti()) kfc=-kfc;
      m_flavlist[m_fcnt]=kfc;
      m_momlist[m_fcnt]=part->Momentum();
      ++m_fcnt;
      if (m_fcnt>=m_flavlist.size()) {
	m_flavlist.resize(m_flavlist.size()+m_avsize);
	m_momlist.resize(m_momlist.size()+m_avsize);
      }
    }
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
      m_evtlist[m_cnt2].fscale=sqrt((*nlos)[j]->m_muf2);
      m_evtlist[m_cnt2].rscale=sqrt((*nlos)[j]->m_mur2);

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

void Output_RootNtuple::StoreEvt()
{
  if (m_cnt2==0) return;
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
