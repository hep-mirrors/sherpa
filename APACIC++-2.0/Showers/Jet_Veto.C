#include "Jet_Veto.H"

using namespace APACIC;
using namespace ATOOLS;

Jet_Veto::Jet_Veto(ATOOLS::Jet_Finder *const jf,
		   Timelike_Kinematics *const kin):
  p_jf(jf), p_kin(kin), p_istrees(NULL), p_fstree(NULL), 
  m_mode(jv::none), m_maxjets(2) {}

int Jet_Veto::TestKinematics(const int mode,Knot *const mo)
{
  if (mode==1 && !(m_mode&jv::mlm)) return 1;
  msg_Debugging()<<METHOD<<"("<<mode<<"): p_{t jet} = "
		 <<sqrt(p_jf->ShowerPt2())<<" {"<<std::endl;
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
  std::vector<ATOOLS::Vec4D> savejets(jets);
  m_rates.resize(jets.size());
  size_t nmin(p_jf->Type()==1?1:0);
  while (jets.size()>nmin) {
    double pt2min(std::numeric_limits<double>::max());
    std::vector<Vec4D>::iterator jit(jets.begin());
    std::vector<Vec4D>::iterator jwin(jets.end()), kwin(jets.end());
    for (size_t j(0);j<jets.size();++j) {
      if (p_jf->Type()!=1) {
 	double pt2j(p_jf->MTij2(Vec4D(1.,0.,0.,1.),*jit));
 	if (pt2j<pt2min) {
 	  pt2min=pt2j;
 	  jwin=jit;
 	  kwin=jets.end();
 	}
      }
      std::vector<Vec4D>::iterator kit(jit);
      for (size_t k(j+1);k<jets.size();++k) {
	++kit;
	double pt2jk(p_jf->MTij2(*jit,*kit));
	if (pt2jk<pt2min) {
	  pt2min=pt2jk;
	  jwin=jit;
	  kwin=kit;
	}
      }
      ++jit;
    }
    if (jwin!=jets.end()) {
      msg_Debugging()<<"jetrate Q_{"<<jets.size()<<"->"
		     <<jets.size()-1<<"} = "<<sqrt(pt2min)<<"\n";
      m_rates[jets.size()-1]=pt2min;
      if (kwin!=jets.end()) *kwin+=*jwin;
      jets.erase(jwin);
    }
    else {
      msg.Error()<<METHOD<<"("<<mode<<"): No min p_T. Abort."<<std::endl;
      return 0;
    }
  }
  size_t njets(nmin), maxjets(m_maxjets);
  double crit(rpa.gen.Ycut()*sqr(rpa.gen.Ecms()));
  for (;njets<m_rates.size();++njets) if (m_rates[njets]<crit) break;
  msg_Debugging()<<"produced "<<njets
		 <<" jets out of "<<hard<<", nmax = "<<maxjets<<"\n";
  if (njets>hard) {
    msg_Debugging()<<"produced "<<(njets-hard)
		   <<" additional jets"<<std::endl;
    msg_Debugging()<<"}\n";
    if (hard<maxjets) {
      if (mode==0) return 0;
      else return -1;
    }
    else {
      crit=p_jf->ShowerPt2();
      for (njets=nmin;njets<m_rates.size();++njets) 
	if (m_rates[njets]<crit) break;
      if (njets>hard)
	if (mode==0) return 0;
	else return -1;
    }
  }
  else if (njets<hard) {
    msg_Debugging()<<"lost "<<(hard-njets)
		   <<" jets"<<std::endl;
    msg_Debugging()<<"}\n";
    if (mode==0) return 0;
    return -1;
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
  if ((knot->left==NULL || knot->left->stat==3) &&
      knot->part->Flav().Strong()) vecs.push_back(knot->part->Momentum());
  int dtest(1);
  size_t rhard(hard);
  if (knot->left && knot->left->stat!=3) {
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
  if (m_mode&jv::mlm || !(m_mode&jv::initial)) return 1;
  msg_Debugging()<<METHOD<<"("<<knot->kn_no<<","<<knot->part->Info()
		 <<"): p_{t jet} = "<<sqrt(p_jf->ShowerPt2())<<"\n";
  double E2(sqr(sqrt(knot->left->E2)+sqrt(knot->right->E2)));
  double z(p_kin->LightConeZ(knot->right->z,E2,knot->tout,
			       knot->right->t,knot->left->t));
  double pt2(z*(1.0-z)*knot->tout-(1.0-z)*knot->right->t-z*knot->left->t);
  msg_Debugging()<<" pt = "<<sqrt(pt2)<<"\n";
  if (knot->part->Info()!='H') {
    if (pt2<0.0 || pt2>p_jf->ShowerPt2()) return 0;
  }
  else {
    if (pt2<=p_jf->ShowerPt2()) return 0;
  }
  return 1;
}

int Jet_Veto::TestFSKinematics(Knot *const knot)
{
  if (knot->stat==0 ||
      m_mode&jv::mlm || !(m_mode&jv::final)) return 1;
  msg_Debugging()<<METHOD<<"("<<knot->kn_no<<","<<knot->part->Info()
		 <<"): p_{t jet} = "<<sqrt(p_jf->ShowerPt2())<<"\n";
  msg_Indent();
  if (!(knot->left->part->Info()=='H' && 
	knot->right->part->Info()=='H')) {
    double pt2(p_jf->MTij2(knot->left->part->Momentum(),
 			   knot->right->part->Momentum()));
    msg_Debugging()<<" pt = "<<sqrt(pt2)<<", pt_old = "
		   <<sqrt(knot->pt2lcm)<<"\n";
    knot->left->pt2lcm=knot->right->pt2lcm=pt2;
    if (pt2>p_jf->ShowerPt2()) return 0;
  }
  if (knot->part->Info()=='H' && knot->prev!=NULL) {
    Knot *si(knot->prev->left);
    if (si==knot) si=knot->prev->right;
    double pt2(p_jf->MTij2(knot->part->Momentum(),
			   si->part->Momentum()));
    msg_Debugging()<<" pt = "<<sqrt(pt2)<<", pt_old = "
		   <<sqrt(knot->pt2lcm)<<"\n";
    if (pt2<=p_jf->ShowerPt2()) return 0;
  }
  return 1;
}
