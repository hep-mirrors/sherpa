#include "Jet_Veto.H"

using namespace APACIC;
using namespace ATOOLS;

Jet_Veto::Jet_Veto(ATOOLS::Jet_Finder *const jf,
		   Timelike_Kinematics *const kin):
  p_jf(jf), p_kin(kin), m_mode(jv::none), m_maxjets(2) {}

int Jet_Veto::TestKinematics(const int mode)
{
  msg_Debugging()<<METHOD<<"(..): p_{t jet} = "<<sqrt(p_jf->ShowerPt2())
		 <<" {"<<std::endl;
  msg_Indent();
  size_t hard(0);
  std::vector<Vec4D> jets;
  if (m_mode&jv::initial) {
    if (!CollectISMomenta(p_istrees[0]->GetRoot(),jets,hard)) return 0;
    if (!CollectISMomenta(p_istrees[1]->GetRoot(),jets,hard)) return 0;
  }
  if (m_mode&jv::final)
    if (!CollectFSMomenta(p_fstree->GetRoot(),jets,hard)) return 0;
  bool cluster(true);
  while (cluster) {
    cluster=false;
    double pt2min(p_jf->ShowerPt2());
    std::vector<Vec4D>::iterator jwin(jets.end()), kwin(jets.end());
    for (std::vector<Vec4D>::iterator 
	   jit(jets.begin());jit!=jets.end();++jit) {
      for (std::vector<Vec4D>::iterator kit(jit);kit!=jets.end();++kit) {
	if (kit==jit) ++kit;
	if (kit==jets.end()) break;
	double pt2jk(p_jf->MTij2(*jit,*kit));
	if (pt2jk<pt2min) {
	  pt2min=pt2jk;
	  jwin=jit;
	  kwin=kit;
	}
      }
    }
    if (jwin!=jets.end() && kwin!=jets.end()) {
      *jwin+=*kwin;
      jets.erase(kwin);
      cluster=true;
    }
  }
  if (p_jf->Type()!=1) {
    for (std::vector<Vec4D>::iterator 
	   jit(jets.begin());jit!=jets.end();++jit) {
      double pt2j(p_jf->MTij2(Vec4D(1.,0.,0.,1.),*jit));
      if (pt2j<p_jf->ShowerPt2()) jit=--jets.erase(jit);
    }
  }
  msg_Debugging()<<"produced "<<jets.size()
		 <<" jets out of "<<hard<<std::endl;
  if (jets.size()>hard) {
    if (hard<=m_maxjets) {
      msg_Debugging()<<"produced "<<(jets.size()-hard)
		     <<" additional jets"<<std::endl;
      msg_Debugging()<<"}\n";
      if (mode) return -1;
      return 0;
    }
  }
  else if (jets.size()<hard && mode) {
    msg_Debugging()<<"lost "<<(hard-jets.size())
		   <<" jets"<<std::endl;
    msg_Debugging()<<"}\n";
    return -1;
  }
  msg_Debugging()<<"}\n";
  return 1;
}

int Jet_Veto::CollectISMomenta(Knot *knot,std::vector<Vec4D> &vecs,
			       size_t &hard)
{
  msg_Debugging()<<METHOD<<"(["<<knot->part->Flav()<<","
		 <<(knot->prev?knot->prev->kn_no:-1)<<"->"
		 <<knot->kn_no<<"->("
		 <<(knot->left?knot->left->kn_no:-1)<<","
		 <<(knot->right?knot->right->kn_no:-1)<<"),"<<knot->stat
		 <<","<<knot->part->Info()<<"]): {\n";
  msg_Indent();
  int dtest(1);
  if (knot->left)
    dtest=CollectFSMomenta(knot->left,vecs,hard);
  if (knot->prev==NULL || knot->stat==3) {
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
  msg_Debugging()<<METHOD<<"(["<<knot->part->Flav()<<","<<knot->kn_no<<"->("
		 <<(knot->left?knot->left->kn_no:-1)<<","
		 <<(knot->right?knot->right->kn_no:-1)<<"),"<<knot->stat
		 <<","<<knot->part->Info()<<"]): {\n";
  msg_Indent();
  char type(toupper(knot->part->Info()));
  if (knot->left && type=='F')
    p_kin->DoSingleKinematics(knot);
  if ((knot->left==NULL || knot->left->stat==3) &&
      knot->part->Flav().Strong() && type!='G') {
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
  if (rhard==hard && knot->part->Flav().Strong() && type=='H') {
    msg_Debugging()<<"hard "<<knot->kn_no<<std::endl;
    ++hard;
  }
  msg_Debugging()<<"}\n";
  return dtest;
}

