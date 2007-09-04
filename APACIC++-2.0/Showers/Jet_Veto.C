#include "Jet_Veto.H"

#include "MyStrStream.H"
#include "Shell_Tools.H"
#include "Data_Reader.H"

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
  if (p_jf->Type()>1) {
    if (IsZero(p1.PPerp2())) return std::numeric_limits<double>::max();
    return p_jf->MTij2(Vec4D(1.,0.,0.,1.),p1); 
  }
  return std::numeric_limits<double>::max();
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

void E_Scheme::operator()(const ATOOLS::Vec4D &p1) 
{ 
}

Vec4D P_Scheme::operator()(const ATOOLS::Vec4D &p1,const ATOOLS::Vec4D &p2) 
{ 
  Vec4D p12(p1+p2);
  p12[0]=p12.PSpat();
  return p12;
}

void P_Scheme::operator()(const ATOOLS::Vec4D &p1) 
{ 
}

#include "Cluster_Algorithm.C"

template class Cluster_Algorithm<Vec4D,PT_Measure,E_Scheme>;

Jet_Veto::Jet_Veto(ATOOLS::Jet_Finder *const jf,
		   Timelike_Kinematics *const kin):
  p_jf(jf), p_cluster(new Cluster_Type()), 
  p_kin(kin), p_istrees(NULL), p_fstree(NULL), 
  m_mode(0), m_cmode(0), m_jmode(jv::none), m_ljmode(jv::none) 
{
  p_cluster->Measure().SetJetFinder(p_jf);
  std::string helps;
  Data_Reader reader;
  if (reader.ReadFromFile(helps,"PRINT_PS_RATES") && helps.length()>0) {
    m_histos.resize(5);
    for (size_t i(0);i<m_histos.size();++i)
      m_histos[i] = new Histogram(10,1.0e-5,1.0,100);
  }
  exh->AddTerminatorObject(this);
}

Jet_Veto::~Jet_Veto()
{
  PrepareTerminate();
  exh->RemoveTerminatorObject(this);
}

void Jet_Veto::PrepareTerminate()
{
  if (m_histos.empty()) return;
  std::string pname("ps_rates_"+ToString(rpa.gen.Ecms())+"/");
  if (MakeDir(pname.c_str(),448)) {
    for (size_t i(0);i<m_histos.size();++i) {
      m_histos[i]->Finalize();
      m_histos[i]->Output((pname+"r_"+ToString(i)+".dat").c_str());     
    }
  }
}

int Jet_Veto::TestKinematics(Knot *const mo)
{
  PROFILE_HERE;
  m_q2hard=0.0;
  bool jv(m_jmode&jv::global);
  bool ljv(m_ljmode&jv::global);
  if (!jv && !ljv) return 1;
  msg_Debugging()<<*p_fstree;
  msg_Debugging()<<METHOD<<"("<<mo<<"): p_{t,jv} = "
		 <<(mo?mo->qjv:p_fstree->GetRoot()->qjv)<<" ("<<jv
		 <<"), p_{t,ljv} = "<<(mo?mo->qljv:p_fstree->GetRoot()->qljv)
		 <<" ("<<ljv<<") {"<<std::endl;
  msg_Indent();
  size_t hard(0), inthard(0);
  std::vector<Knot*> next;
  std::vector<Vec4D> jets;
  if (mo==NULL) {
    if (p_istrees!=NULL) {
      if (!CollectISMomenta(p_istrees[0]->GetRoot(),jets,
			    next,hard,inthard)) return 0;
      if (!CollectISMomenta(p_istrees[1]->GetRoot(),jets,
			    next,hard,inthard)) return 0;
    }
    if (!CollectFSMomenta(p_fstree->GetRoot(),jets,
			  next,hard,inthard)) return 0;
  }
  else {
    if (!CollectFSMomenta(mo,jets,next,hard,inthard,mo)) return 0;
  }
  for (size_t i(0);i<next.size();++i) {
    int sub(TestKinematics(next[i]));
    if (sub!=1) return sub;
  }
  m_cmode=m_mode;
  std::vector<ATOOLS::Vec4D> savejets(jets);
  int nmin(p_jf->Type()==1?1:0);
  m_rates.resize(Max(int(jets.size())-nmin,0),0.0);
  p_cluster->SetPoints(jets);
  p_cluster->Cluster(nmin,cs::num);
  msg_Debugging()<<"hard scale Q_h = "<<sqrt(m_q2hard)<<"\n";
  for (size_t i(0);i<m_rates.size();++i) {
    m_rates[m_rates.size()-i-1]=p_cluster->DMins()[i];
    msg_Debugging()<<"jetrate Q_{"<<(jets.size()-i)<<"->"
		   <<(jets.size()-i-1)<<"} = "
		   <<sqrt(p_cluster->DMins()[i])<<"\n";
  }
  for (size_t i(0);i<m_histos.size() && i<m_rates.size();++i)
    m_histos[i]->Insert(m_rates[i]/m_q2hard);
  double jcrit(mo?sqr(mo->qjv):sqr(p_fstree->GetRoot()->qjv));
  double ljcrit(mo?sqr(mo->qljv):sqr(p_fstree->GetRoot()->qljv));
  if (!jv) jcrit=std::numeric_limits<double>::max();
  if (!ljv) ljcrit=0.0;
  msg_Debugging()<<"jet veto ("<<jv<<") "<<sqrt(jcrit)
		 <<", lose jet veto ("<<ljv<<") "<<sqrt(ljcrit)<<"\n";
  size_t njets(0), nljets(0);
  for (;njets<m_rates.size();++njets) if (m_rates[njets]<jcrit) break;
  for (;nljets<m_rates.size();++nljets) if (m_rates[nljets]<ljcrit) break;
  njets+=nmin;
  nljets+=nmin;
  size_t maxjets(mo!=NULL?mo->maxjets:p_fstree->GetRoot()->maxjets);
  msg_Debugging()<<(mo!=NULL)<<"?"<<(mo!=NULL?mo->maxjets:-1)<<":"
		 <<p_fstree->GetRoot()->maxjets<<"\n";
  maxjets+=inthard;
  hard+=inthard;
  msg_Debugging()<<"produced "<<njets<<" / "<<nljets
		 <<" jets out of "<<hard<<", nmax = "
		 <<maxjets<<", internal = "<<inthard
		 <<", nmin = "<<nmin<<"\n";
  if (njets>hard) {
    msg_Debugging()<<"produced "<<(njets-hard)
		   <<" additional jets"<<std::endl;
    msg_Debugging()<<"}\n";
    if (jv) return 0;
  }
  else if (nljets<hard) {
    msg_Debugging()<<"lost "<<(hard-nljets)
		   <<" jets"<<std::endl;
    msg_Debugging()<<"}\n";
    if (ljv) return -1;
  }
  msg_Debugging()<<"}\n";
  return 1;
}

int Jet_Veto::CollectISMomenta(Knot *knot,std::vector<Vec4D> &vecs,
			       std::vector<Knot*> &next,size_t &hard,
			       size_t &inthard)
{
  int dtest(1);
  if (knot->left && knot->left->t<dabs(knot->right->t) &&
      knot->left->E2>0.0)
    dtest=CollectFSMomenta(knot->left,vecs,next,hard,inthard);
  m_q2hard=Max(m_q2hard,dabs(knot->tmax));
  if (knot->prev!=NULL) 
    dtest=CollectISMomenta(knot->prev,vecs,next,hard,inthard);
  return dtest;
}

int Jet_Veto::CollectFSMomenta(Knot *knot,std::vector<Vec4D> &vecs,
			       std::vector<Knot*> &next,size_t &hard,
			       size_t &inthard,Knot *const first)
{
  if (knot->part->Flav().Strong()) {
    switch (m_cmode) {
    case 0:
      if (knot->left==NULL) {
	msg_Debugging()<<"take knot "<<knot->kn_no
		       <<", mom "<<knot->part->Momentum()<<"\n";
	vecs.push_back(knot->part->Momentum());
      }
      else if (knot->decay!=NULL && knot->left->decay!=knot->decay && 
	       knot!=first) {
	msg_Debugging()<<"split case, take knot "<<knot->kn_no
		       <<", mom "<<knot->part->Momentum()<<"\n";
	vecs.push_back(knot->part->Momentum());
	next.push_back(knot);
	if (knot->part->Flav().Strong() && knot->part->Info()=='H') {
	  msg_Debugging()<<"hard knot "<<knot->kn_no<<"\n";
	  ++inthard;
	}
	return 2;
      }
      break;
    case 1:
      if (knot!=p_cur && knot->part->Info()=='H' &&
	  (knot->left==NULL || knot->left->part->Info()!='H' || 
 	  (knot->shower==2 && knot->left->decay!=knot->decay))) 
	vecs.push_back(knot->part->Momentum());
      break;
    }
  }
  int dtest1(1), dtest2(1);
  size_t rhard(hard), rinthard(inthard);
  if (knot->left!=NULL) {
    if (knot->left->part->Flav().Strong()) m_q2hard=Max(m_q2hard,knot->tmax);
    dtest1=CollectFSMomenta(knot->left,vecs,next,hard,inthard);
    if (dtest1>0) 
      dtest2=CollectFSMomenta(knot->right,vecs,next,hard,inthard);
  }
  if (rhard==hard && rinthard==inthard && 
      knot->part->Flav().Strong() && knot->part->Info()=='H') {
    if (dtest1==2 || dtest2==2) ++inthard;
    else ++hard;
    msg_Debugging()<<"hard knot "<<knot->kn_no<<"\n";
  }
  if (dtest1>0 && dtest2>0) {
    if (dtest1==2 || dtest2==2) return 2;
    return 1;
  }
  return 0;
}

int Jet_Veto::TestISKinematics(Knot *const knot,Knot *const partner)
{
  PROFILE_HERE;
  bool jv(m_jmode&jv::initial);
  bool ljv(m_ljmode&jv::initial);
  msg_Debugging()<<METHOD<<"("<<knot->kn_no<<","<<knot->part->Info()
		 <<"): p_{t,jv} = "<<knot->qjv<<" ("<<jv
		 <<"), p_{t,ljv} = "<<knot->qljv<<" ("<<ljv<<")\n";
  Vec4D p3(knot->part->Momentum()), p2(partner->part->Momentum());
  Vec4D p4(knot->left->part->Momentum()), cms(p3+p2);
  double m3(knot->part->Flav().Mass()), m4(knot->left->part->Flav().Mass());
  Poincare boost(cms);
  boost.Boost(p2);
  boost.Boost(p3);
  boost.Boost(p4);
  Poincare rot(p3,Vec4D::ZVEC);
  rot.Rotate(p4);
  double pt2(p_jf->MTij2(p4,Vec4D(1.0,0.0,0.0,1.0),m4,m3));
  msg_Debugging()<<"kt = "<<sqrt(knot->pt2lcm)
 		 <<" pt = "<<sqrt(pt2)<<" <- knot = "
		 <<knot->part->Momentum()<<" left = "
		 <<knot->left->part->Momentum()<<"\n";
  if (knot->part->Info()!='H') {
    if (jv && pt2>sqr(knot->qjv)) return 0;
  }
  return 1;
}

int Jet_Veto::TestFSKinematics(Knot *const knot)
{
  PROFILE_HERE;
  bool jv(m_jmode&jv::final);
  bool ljv(m_ljmode&jv::final);
  if (knot->left==NULL) return 1;
  msg_Debugging()<<METHOD<<"("<<knot->kn_no<<","<<knot->part->Info()
		 <<"): p_{t,jv} = "<<knot->qjv<<" ("<<jv
		 <<"), p_{t,ljv} = "<<knot->qljv<<" ("<<ljv<<")\n";
  int type(p_jf->Type());
  Knot *d[3]={knot,knot->left,knot->right};
  for (short unsigned int i(1);i<3;++i) {
    if (d[i]->left==NULL) continue;
    bool f[3]={d[3-i]->part->Flav().Strong()&&
	       d[i]->left->part->Flav().Strong(),
	       d[3-i]->part->Flav().Strong()&&
	       d[i]->right->part->Flav().Strong(),
	       d[i]->left->part->Flav().Strong()&&
	       d[i]->right->part->Flav().Strong()};
    Vec4D p[3]={d[3-i]->part->Momentum(),
		d[i]->left->part->Momentum(),
		d[i]->right->part->Momentum()};
    double m[3]={d[3-i]->part->Flav().Mass(),
		 d[i]->left->part->Flav().Mass(),
		 d[i]->right->part->Flav().Mass()};
#ifdef BOOST_Decay_JV
    if (knot->cms!=Vec4D()) {
      Poincare cms(knot->cms);
      msg_Debugging()<<"boost in cms "<<knot->cms<<"\n";
      for (short unsigned int j(0);j<3;++j) cms.Boost(p[i]);
    }
#endif
    double sf(sqrt(sqr(p[0][0])-d[3-i]->tout)/p[0].PSpat());
    for (short unsigned int j(1);j<3;++j) p[0][j]*=sf;
    p_jf->SetType(knot->cms!=Vec4D()?1:type);
    double jpt2[3]={f[0]?p_jf->MTij2(p[0],p[1],m[0],m[1]):0.0,
		    f[1]?p_jf->MTij2(p[0],p[2],m[0],m[2]):0.0,
		    f[2]?p_jf->MTij2(p[1],p[2],m[1],m[2]):0.0};
    p_jf->SetType(type);
    int jets(f[0]+f[1]+f[2]), ljets(0);
    if (d[i]->left->part->Info()!='H' || 
	d[i]->right->part->Info()!='H') { 
      for (short unsigned int j(0);j<3;++j) {
	if (f[j]&&jpt2[j]<=sqr(knot->qjv)) --jets;
	if (f[j]&&jpt2[j]>sqr(knot->qljv)) ++ljets;
      }
      msg_Debugging()<<"  jv ("<<d[3-i]->kn_no<<","
		     <<d[3-i]->part->Flav()<<","<<d[3-i]->part->Info()
		     <<"),("<<d[i]->kn_no<<","<<d[i]->part->Flav()
		     <<","<<d[i]->part->Info()
		     <<")->("<<d[i]->left->kn_no
		     <<","<<d[i]->left->part->Flav()
		     <<","<<d[i]->left->part->Info()
		     <<"),("<<d[i]->right->kn_no
		     <<","<<d[i]->right->part->Flav()
		     <<","<<d[i]->right->part->Info()
		     <<"), jpt = {"<<sqrt(jpt2[0])<<","
		     <<sqrt(jpt2[1])<<","<<sqrt(jpt2[2])<<"} -> "
		     <<jets<<"("<<(f[0]+f[1]+f[2])<<") jets of type "
		     <<(knot->cms!=Vec4D()?1:type)<<"\n";
      if (jv && jets>Max(f[0]+f[1]+f[2]-1,0)) {
	msg_Debugging()<<"  jet veto\n";
	return 0;
      }
      if (ljv && ljets<1 && d[i]->part->Info()=='H' && 
	  d[3-i]->part->Info()=='H') {
	msg_Debugging()<<"  lose jet veto\n";
	return 0;
      }
    }
  }
  return 1;
}

