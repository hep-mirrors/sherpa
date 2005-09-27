#include "Jet_Veto.H"

using namespace APACIC;
using namespace ATOOLS;

Jet_Veto::Jet_Veto(ATOOLS::Jet_Finder *const jf,
		   Timelike_Kinematics *const kin):
  p_jf(jf), p_kin(kin), m_mode(jv::none), m_maxjets(2) {}

#ifdef DURHAM_IS_MEASURE
void BoostInCMS(const ATOOLS::Vec4D &cms,std::vector<ATOOLS::Vec4D> &moms)
{
  Poincare cmsboost(cms);
  for (std::vector<Vec4D>::iterator 
 	 mit(moms.begin());mit!=moms.end();++mit) cmsboost.Boost(*mit);
}
#endif

int Jet_Veto::TestKinematics(const int mode,Knot *const mo)
{
  msg_Debugging()<<METHOD<<"(..): p_{t jet} = "<<sqrt(p_jf->ShowerPt2())
		 <<" {"<<std::endl;
  msg_Indent();
  size_t hard(0);
  std::vector<Vec4D> jets;
  if (mo==NULL) {
    if (m_mode&jv::initial) {
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
#ifdef DURHAM_IS_MEASURE
  int oldmode(p_jf->Type());
  if (oldmode!=1) {
    BoostInCMS(p_fstree->GetRoot()->part->Momentum(),jets);
    p_jf->SetType(1);
  }
#endif
  std::vector<ATOOLS::Vec4D> savejets(jets);
  m_rates.resize(jets.size());
  size_t nmin(p_jf->Type()==1?1:0);
  while (jets.size()>nmin) {
    double pt2min(4.0*sqr(rpa.gen.Ecms()));
    std::vector<Vec4D>::iterator jwin(jets.end()), kwin(jets.end());
    for (std::vector<Vec4D>::iterator 
	   jit(jets.begin());jit!=jets.end();++jit) {
      if (p_jf->Type()!=1) {
 	double pt2j(p_jf->MTij2(Vec4D(1.,0.,0.,1.),*jit));
 	if (pt2j<pt2min) {
 	  pt2min=pt2j;
 	  jwin=jit;
 	  kwin=jets.end();
 	}
      }
      for (std::vector<Vec4D>::iterator kit(jit);kit!=jets.end();++kit) {
	if (kit==jit) 
	  if ((++kit)==jets.end()) break;
	double pt2jk(p_jf->MTij2(*jit,*kit));
	if (pt2jk<pt2min) {
	  pt2min=pt2jk;
	  jwin=jit;
	  kwin=kit;
	}
      }
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
      for (size_t i(0);i<savejets.size();++i) PRINT_INFO(i<<" "<<savejets[i]);
      return 0;
    }
  }
#ifdef DURHAM_IS_MEASURE
  p_jf->SetType(oldmode);
#endif
  size_t njets(nmin);
  double crit(rpa.gen.Ycut()*sqr(rpa.gen.Ecms()));
  for (;njets<m_rates.size();++njets) if (m_rates[njets]<crit) break;
  msg_Debugging()<<"produced "<<njets
		 <<" jets out of "<<hard<<std::endl;
  if (njets>hard) {
    msg_Debugging()<<"produced "<<(njets-hard)
		   <<" additional jets"<<std::endl;
    msg_Debugging()<<"}\n";
    if (hard<m_maxjets) {
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
    if (mode>0) {
      msg_Debugging()<<"lost "<<(hard-njets)
		     <<" jets"<<std::endl;
      msg_Debugging()<<"}\n";
      return -1;
    }
  }
  msg_Debugging()<<"}\n";
  return 1;
}

int Jet_Veto::CollectISMomenta(Knot *knot,std::vector<Vec4D> &vecs,
			       size_t &hard)
{
  msg_Debugging()<<METHOD<<"(["<<(knot->prev?knot->prev->kn_no:-1)<<","
		 <<(knot->prev?knot->prev->stat:-1)<<"]->["
		 <<knot->kn_no<<","<<knot->stat
		 <<","<<knot->part->Info()<<","<<knot->part->Flav()<<"]->["
		 <<(knot->left?knot->left->kn_no:-1)<<","
		 <<(knot->left?knot->left->stat:-1)<<"],["
		 <<(knot->right?knot->right->kn_no:-1)<<","
		 <<(knot->right?knot->right->stat:-1)<<"]): {\n";
  msg_Indent();
  int dtest(1);
  if (knot->left && knot->left->t<dabs(knot->right->t) &&
      knot->left->E2>0.0)
    dtest=CollectFSMomenta(knot->left,vecs,hard);
  if (knot->prev==NULL || knot->stat==3) {
    Vec4D mom(-1.0*knot->part->Momentum());
    if (knot->left) {
      mom=-1.0*(knot->left->part->Momentum()+
		knot->right->part->Momentum());
    }
#ifdef DURHAM_IS_MEASURE
    msg_Debugging()<<"take in mom "<<knot->kn_no<<" -> "
 		   <<knot->part->Momentum()<<"\n";
    vecs.push_back(knot->part->Momentum());
    ++hard;
#endif
    msg_Debugging()<<"}\n";
    return dtest;
  }
  dtest=CollectISMomenta(knot->prev,vecs,hard);
  msg_Debugging()<<"}\n";
  return dtest;
}

int Jet_Veto::CollectFSMomenta(Knot *knot,std::vector<Vec4D> &vecs,
			       size_t &hard)
{
  msg_Debugging()<<METHOD<<"(["<<(knot->prev?knot->prev->kn_no:-1)<<","
		 <<(knot->prev?knot->prev->stat:-1)<<"]->["
		 <<knot->kn_no<<","<<knot->stat
		 <<","<<knot->part->Info()<<","<<knot->part->Flav()<<"]->["
		 <<(knot->left?knot->left->kn_no:-1)<<","
		 <<(knot->left?knot->left->stat:-1)<<"],["
		 <<(knot->right?knot->right->kn_no:-1)<<","
		 <<(knot->right?knot->right->stat:-1)<<"]): {\n";
  msg_Indent();
  if (knot->left && knot->left->t!=knot->t) 
    p_kin->DoSingleKinematics(knot);
  if ((knot->left==NULL || knot->left->stat==3) &&
      knot->part->Flav().Strong()) {
    msg_Debugging()<<"take mom "<<knot->kn_no<<" -> "
		   <<knot->part->Momentum()<<"\n";
    vecs.push_back(knot->part->Momentum());
  }
  int dtest(1);
  size_t rhard(hard);
  if (knot->left && knot->left->stat!=3) {
    dtest=CollectFSMomenta(knot->left,vecs,hard);
    if (dtest==1) 
      dtest=CollectFSMomenta(knot->right,vecs,hard);
  }
  if (rhard==hard && knot->part->Flav().Strong() && 
      knot->part->Info()=='H') {
    msg_Debugging()<<"hard "<<knot->kn_no<<std::endl;
    ++hard;
  }
  msg_Debugging()<<"}\n";
  return dtest;
}

