#include "Jet_Veto.H"

#ifdef PROFILE__all
#include "prof.hh"
#else 
#define PROFILE_HERE 
#define PROFILE_LOCAL(LOCALNAME) 
#endif

using namespace APACIC;
using namespace ATOOLS;

double PT_Measure::operator()(const ATOOLS::Vec4D &p1) 
{ 
  return p_jf->MTij2(Vec4D(1.,0.,0.,1.),p1); 
}

double PT_Measure::operator()
  (const ATOOLS::Vec4D &p1,const ATOOLS::Vec4D &p2) 
{ 
  return p_jf->MTij2(p1,p2); 
}

Vec4D E_Scheme::operator()(const ATOOLS::Vec4D &p1,const ATOOLS::Vec4D &p2) 
{ 
  return p1+p2; 
}

#include "Cluster_Algorithm.C"

template class Cluster_Algorithm<Vec4D,PT_Measure,E_Scheme>;

Jet_Veto::Jet_Veto(ATOOLS::Jet_Finder *const jf,
		   Timelike_Kinematics *const kin):
  p_jf(jf), p_cluster(new Cluster_Type()), 
  p_kin(kin), p_istrees(NULL), p_fstree(NULL), 
  m_mode(jv::none), m_maxjets(2), m_cmode(0), m_jmode(1), m_ljmode(1) 
{
  p_cluster->Measure().SetJetFinder(p_jf);
}

int Jet_Veto::TestKinematics(const int mode,Knot *const mo)
{
  PROFILE_HERE;
  if (mode==1 && !(m_mode&jv::mlm)) {
    p_cur=NULL;
    m_cmode=0;
  }
  if (mode==0 && !(m_mode&jv::final) && !(m_mode&jv::initial)) return 1;
  msg_Debugging()<<METHOD<<"("<<mode<<"): p_{t jet} = "
		 <<sqrt(p_jf->ShowerPt2())<<" / "
		 <<sqrt(m_ycut)*rpa.gen.Ecms()<<" {"<<std::endl;
  msg_Indent();
  size_t hard(0);
  std::vector<Vec4D> jets;
  if (mo==NULL) {
    if (m_mode&jv::initial && p_istrees!=NULL) {
      if (!CollectISMomenta(p_istrees[0]->GetRoot(),jets,hard)) return 0;
      if (!CollectISMomenta(p_istrees[1]->GetRoot(),jets,hard)) return 0;
    }
    if (m_mode&jv::final) {
      if (!CollectFSMomenta(p_fstree->GetRoot(),jets,hard)) return 0;
    }
  }
  else {
    if (!CollectFSMomenta(mo,jets,hard)) return 0;
  }
  m_cmode=0;
  std::vector<ATOOLS::Vec4D> savejets(jets);
  m_rates.resize(jets.size());
  size_t nmin(p_jf->Type()==1?1:0);
  p_cluster->SetPoints(jets);
  p_cluster->Cluster(nmin,cs::num);
  for (size_t i(0);i<jets.size();++i) {
    m_rates[jets.size()-i-1]=p_cluster->DMins()[i];
    msg_Debugging()<<"jetrate Q_{"<<(jets.size()-i)<<"->"
		   <<(jets.size()-i-1)<<"} = "
		   <<sqrt(p_cluster->DMins()[i])<<"\n";
  }
  size_t njets(nmin), nljets(nmin);
  double jcrit(p_jf->ShowerPt2()), ljcrit(m_ycut*sqr(rpa.gen.Ecms()));
  if (m_jmode==0) jcrit=std::numeric_limits<double>::max();
  if (m_ljmode==0) ljcrit=0.0;
  msg_Debugging()<<"jet veto ("<<m_jmode<<") "<<jcrit
		 <<", lose jet veto ("<<m_ljmode<<") "<<ljcrit<<"\n";
  for (;njets<m_rates.size();++njets) if (m_rates[njets]<jcrit) break;
  for (;nljets<m_rates.size();++nljets) if (m_rates[nljets]<ljcrit) break;
  msg_Debugging()<<"produced "<<njets<<" / "<<nljets
		 <<" jets out of "<<hard<<", nmax = "<<m_maxjets<<"\n";
  if (njets>hard) {
    msg_Debugging()<<"produced "<<(njets-hard)
		   <<" additional jets"<<std::endl;
    msg_Debugging()<<"}\n";
    if (mode==0) return 0;
    if (m_mode&jv::mlm) return -1;
  }
  else if (nljets<hard) {
    msg_Debugging()<<"lost "<<(hard-nljets)
		   <<" jets"<<std::endl;
    msg_Debugging()<<"}\n";
    if (mode==1) return -1;
  }
  msg_Debugging()<<"}\n";
  return 1;
}

int Jet_Veto::CollectISMomenta(Knot *knot,std::vector<Vec4D> &vecs,
			       size_t &hard)
{
  int dtest(1);
  if (knot->left && knot->left->t<dabs(knot->right->t) &&
      knot->left->E2>0.0)
    dtest=CollectFSMomenta(knot->left,vecs,hard);
  if (knot->prev==NULL || knot->stat==3) return dtest;
  dtest=CollectISMomenta(knot->prev,vecs,hard);
  return dtest;
}

int Jet_Veto::CollectFSMomenta(Knot *knot,std::vector<Vec4D> &vecs,
			       size_t &hard)
{
  if (knot->part->Flav().Strong()) {
    switch (m_cmode) {
    case 0:
      if (knot->left==NULL || knot->left->stat==3) {
	msg_Debugging()<<"take knot "<<knot->kn_no
		       <<", mom "<<knot->part->Momentum()<<"\n";
	vecs.push_back(knot->part->Momentum());
      }
      break;
    case 1:
      if (knot!=p_cur && knot->part->Info()=='H' &&
	  (knot->left==NULL || knot->left->part->Info()!='H')) 
	vecs.push_back(knot->part->Momentum());
      break;
    }
  }
  int dtest(1);
  size_t rhard(hard);
  if (knot->left!=NULL && knot->left->stat!=3) {
    dtest=CollectFSMomenta(knot->left,vecs,hard);
    if (dtest==1) 
      dtest=CollectFSMomenta(knot->right,vecs,hard);
  }
  if (rhard==hard && knot->part->Flav().Strong() && 
      knot->part->Info()=='H') ++hard;
  return dtest;
}

int Jet_Veto::TestISKinematics(Knot *const knot)
{
  PROFILE_HERE;
  if (m_mode&jv::mlm || !(m_mode&jv::initial)) return 1;
  msg_Debugging()<<METHOD<<"("<<knot->kn_no<<","<<knot->part->Info()
		 <<"): p_{t jet} = "<<sqrt(p_jf->ShowerPt2())<<"\n";
  double E2(sqr(sqrt(knot->left->E2)+sqrt(knot->right->E2)));
  double z(p_kin->LightConeZ(knot->right->z,E2,knot->t,
			     knot->right->t,knot->left->t));
  double pt2(z*(1.0-z)*knot->t-(1.0-z)*knot->right->t-z*knot->left->t);
  msg_Debugging()<<"pt_old = "<<sqrt(knot->right->pt2lcm)
 		 <<" pt = "<<sqrt(pt2)<<", z/\\tilde z-1 = "
		 <<z/knot->right->z-1.0<<"\n";
  if (knot->part->Info()!='H') {
    if (pt2<0.0 || pt2>p_jf->ShowerPt2()) return 0;
    knot->left->pt2lcm=knot->pt2lcm=pt2;
  }
  return 1;
}

int Jet_Veto::TestFSKinematics(Knot *const knot)
{
  PROFILE_HERE;
  if (m_mode&jv::global) return TestKinematics(0,NULL);
  if (knot->left==NULL ||
      m_mode&jv::mlm || !(m_mode&jv::final)) return 1;
  msg_Debugging()<<METHOD<<"("<<knot->kn_no<<","<<knot->part->Info()
		 <<"): p_{t jet} = "<<sqrt(p_jf->ShowerPt2())<<"\n";
  msg_Indent();
  if (!(knot->left->part->Info()=='H' && 
	knot->right->part->Info()=='H')) {
    double pt2(p_jf->MTij2(knot->left->part->Momentum(),
 			   knot->right->part->Momentum()));
    msg_Debugging()<<" jv  pt = "<<sqrt(pt2)<<", pt_old = "
		   <<sqrt(knot->pt2lcm)<<" <- "
		   <<knot->left->part->Momentum()<<" "
		   <<knot->right->part->Momentum()<<"\n";
    if (pt2>p_jf->ShowerPt2()) return 0;
    knot->left->pt2lcm=knot->right->pt2lcm=pt2;
  }
  return 1;
}

