#include "AMEGIC++/Cluster/Combine_Table.H"

#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/My_Limits.H"
#include "ATOOLS/Phys/Cluster_Definitions_Base.H"
#include "ATOOLS/Phys/Blob.H"
#include <iomanip>

#define MAXD std::numeric_limits<double>::max()

using namespace AMEGIC;
using namespace AMEGIC;
using namespace MODEL;
using namespace ATOOLS;

int Combine_Table::s_all(0);

Leg::Leg(AMEGIC::Point *const point,const int anti):  
  p_point(point), m_anti(anti), 
  m_nqcd(0), m_nqed(0), m_pqcd(0), m_pqed(0), m_ext(0), 
  m_qcdjets(point!=NULL?point->fl.Strong():0), m_id(0),
  m_kt2(std::numeric_limits<double>::max()), 
  m_kt2qcd(m_kt2), m_kt2qed(m_kt2), 
  m_minkt2(m_kt2), m_minkt2qcd(m_kt2), m_minkt2qed(m_kt2), 
  p_qmin(NULL) {}

Leg::Leg(const Leg &leg): 
  p_point(leg.p_point), m_anti(leg.m_anti), 
  m_nqcd(0), m_nqed(0), m_pqcd(0), m_pqed(0), m_ext(leg.m_ext), 
  m_qcdjets(leg.m_qcdjets), m_id(leg.m_id),
  m_kt2(std::numeric_limits<double>::max()), 
  m_kt2qcd(m_kt2), m_kt2qed(m_kt2), m_minkt2(leg.m_minkt2),  
  m_minkt2qcd(leg.m_minkt2qcd), m_minkt2qed(leg.m_minkt2qed), 
  p_qmin(leg.p_qmin), m_mapfl(leg.m_mapfl) {}

std::ostream &AMEGIC::operator<<
  (std::ostream &str,const std::vector<int> &info)
{
  str<<"(";
  if (info.size()>0) str<<info[0];
  else str<<"<no entry>";
  for (size_t i=1;i<info.size();++i) str<<","<<info[i];
  return str<<")";
}

std::ostream &AMEGIC::operator<<(std::ostream &ostr,const Leg &leg)
{
  return ostr<<leg.p_point<<" "<<leg.m_anti;
}

// ============================================================
//    class Combine_Key
// ============================================================

std::ostream& AMEGIC::operator<<(std::ostream& s,const Combine_Key &ck)
{
  s<<" "<<ck.m_i<<"&"<<ck.m_j<<"%"<<ck.m_k;
  if (ck.m_flav.Kfcode()!=kf_none) s<<"["<<std::setw(6)<<ck.m_flav<<"]";
  else s<<std::string(8,' ');
  return s;
}

Combine_Key::Combine_Key(): 
  m_i(0), m_j(0), m_k(0) {}

Combine_Key::Combine_Key(const int i,const int j,const int k,const Flavour &flav) :
  m_i(i), m_j(j), m_k(k), m_flav(flav) {}

bool AMEGIC::operator<(const Combine_Key & a,const Combine_Key & b) 
{
  if (a.m_i < b.m_i) return true;
  if (a.m_i > b.m_i) return false;
  if (a.m_j < b.m_j) return true;
  if (a.m_j > b.m_j) return false;
  if (a.m_k < b.m_k) return true;
  if (a.m_k > b.m_k) return false;
  if (a.m_flav.Kfcode() > b.m_flav.Kfcode()) return true;
  return false;
}

// ============================================================
//    class Combine_Data
// ============================================================

std::ostream& AMEGIC::operator<<(std::ostream &s,const Combine_Data &cd)
{
  s<<" "<<std::setw(20)<<cd.m_pt2ij<<" "<<std::setw(10)<<sqrt(dabs(cd.m_sij))<<" "
   <<std::setw(3)<<cd.m_strong;
//   s<<" ("<<std::setw(11)<<cd.m_prop<<","<<std::setw(11)<<cd.m_coupling<<") ";
  std::string graphs;
  for (size_t k=0;k<cd.m_graphs.size();++k) graphs+=","+ToString(cd.m_graphs[k]);
  s<<std::setw(50)<<graphs.substr(1);
  if (cd.p_down) s<<" -> "<<cd.p_down->m_no;
  return s;
}

Combine_Data::Combine_Data():
  m_pt2ij(0.0), m_strong(0), p_down(NULL) {}

Combine_Data::Combine_Data(const Combine_Data &cd):
  m_pt2ij(cd.m_pt2ij), m_sij(cd.m_sij), m_prop(cd.m_prop), m_coupling(cd.m_coupling),
  m_weight(cd.m_weight), m_strong(cd.m_strong), 
  p_down(NULL), 
  m_graphs(cd.m_graphs)
{
  if (cd.p_down) THROW(fatal_error,"You are not supposed to use this constructor.");
}

Combine_Data::Combine_Data(const double pt2ij,const int ngraph):
  m_pt2ij(pt2ij), m_sij(0.), m_prop(0.), m_coupling(0.), m_weight(0.), 
  m_strong(0),
  p_down(NULL) 
{
  if (ngraph>=0) m_graphs.push_back(ngraph);
}

Combine_Data::~Combine_Data() 
{
  if (p_down!=NULL) delete p_down;
}

// ============================================================
//    class Combine_Table
// ============================================================

std::ostream& AMEGIC::operator<<(std::ostream& s ,const Combine_Table & ct) 
{
  if (&ct) {
    s<<std::endl<<" Combine_Table ("<<&ct<<") "<<ct.m_no<<" (up=";
    if (ct.p_up) s<<ct.p_up->m_no<<")"<<std::endl; else s<<"#)"<<std::endl;
    s<<" x1="<<ct.m_x1<<" x2="<<ct.m_x2<<" kt_min="<<sqrt(ct.Kt2())
     <<" kt_min(QCD)="<<sqrt(ct.Kt2QCD())<<" kt_min(QED)="
     <<sqrt(ct.Kt2QED())<<std::endl;
    s<<" min{kt_min}="<<sqrt(ct.MinKt2())<<" min{kt_min(QCD)}="
     <<sqrt(ct.MinKt2QCD())<<" min{kt_min(QED)}="
     <<sqrt(ct.MinKt2QED())<<std::endl;
    s<<" ==============="<<std::endl;
    s<<" id"<<std::setw(12)<<"content"<<std::setw(8)
     <<"flav"<<std::setw(5)<<" cut qcd qed"<<std::setw(12)
     <<"q_min"<<std::setw(12)<<"q_{min qcd}"<<std::setw(12)
     <<"q_{min qed}"<<std::setw(12)<<"Q_{cut,me}"<<std::setw(12)
     <<"Q_{cut,jv}"<<std::setw(12)<<"Q_{cut,ps}"
     <<std::setw(12)<<"\\sqrt{|t|}"<<" mom"<<std::endl;
    for (int l=0; l<ct.m_nlegs; ++l) {
      double mc(ct.GetLeg(l).MinKT2QCD()), me(ct.GetLeg(l).MinKT2QED());
      s<<std::setw(3)<<l<<std::setw(12)<<ToString(ID(ct.GetLeg(l).ID()))
       <<std::setw(8)<<ct.p_legs[0][l].Flav()<<std::setw(4)
       <<ct.p_legs[0][l].Point()->t<<" "<<ct.GetLeg(l).OrderQCD()
       <<"/"<<ct.GetLeg(l).NQCD()<<" "<<ct.GetLeg(l).OrderQED()
       <<"/"<<ct.GetLeg(l).NQED()<<std::setw(12)
       <<ct.GetLeg(l).QMin()<<std::setw(12)<<(mc!=MAXD?sqrt(mc):-1)
       <<std::setw(12)<<(me!=MAXD?sqrt(me):-1)<<std::setw(12)
       <<sqrt(dabs(ct.p_moms[l].Abs2()))<<" "<<ct.p_moms[l]<<std::endl;
    }
    s<<" ---------------"<<std::endl;
    const CD_List & cl=ct.m_combinations;
    if (cl.size()>0) {
      s<<" cmb       cos\\theta        k_t \\sqrt{|t|} qcd"
       <<std::setw(50)<<"graphs"<<std::endl;
      for (CD_List::const_iterator cit=cl.begin(); cit!=cl.end(); ++cit) {
 	s<<cit->first<<std::setw(8)
	 <<ct.p_moms[cit->first.m_i].
	  Theta(ct.p_moms[cit->first.m_j])<<cit->second<<std::endl; 
      }
      for (CD_List::const_iterator cit=cl.begin(); cit!=cl.end(); ++cit) {
	if (cit->second.p_down) {
	  s<<*cit->second.p_down<<std::endl;
	}
      }
    }
    else {
      s<<" graph"<<std::setw(8)
       <<"flav"<<std::setw(5)<<" cut qcd qed"<<std::setw(12)
       <<"q_{min qcd}"<<std::setw(12)<<"q_{min qed}"<<std::setw(12)<<std::endl;
      for (int k=0;k<ct.m_nampl;++k) {
	for (int l=0;l<2;++l) {
	  double mc(ct.p_hard[k][l].KT2QCD()), me(ct.p_hard[k][l].KT2QED());
	  s<<std::setw(3)<<k<<"("<<l<<")"<<std::setw(8)
	   <<ct.p_hard[k][l].Flav()<<std::setw(4)
	   <<ct.p_hard[k][l].Point()->t<<" "<<ct.p_hard[k][l].OrderQCD()
	   <<"/"<<ct.p_hard[k][l].NQCD()<<" "<<ct.p_hard[k][l].OrderQED()
	   <<"/"<<ct.p_hard[k][l].NQED()<<std::setw(12)<<(mc!=MAXD?sqrt(mc):-1)
	   <<std::setw(12)<<(me!=MAXD?sqrt(me):-1)<<std::endl;
	}
      }
    }
  } 
  else
    s<<"***empty Combine_Table***"<<std::endl;
  return s;
}

Combine_Table::Combine_Table(AMEGIC::Process_Base *const proc,
			     ATOOLS::Mass_Selector *const ms,
			     Cluster_Definitions_Base *clus,
			     Vec4D *moms, Combine_Table *up):
  p_ms(ms), m_nstrong(0), m_nlegs(0), m_nampl(0),
  m_graph_winner(0), 
  p_up(up), p_legs(0), p_clus(clus), p_moms(moms), p_hard(NULL), p_hardc(NULL)
{
  p_proc=proc;
  m_no=++s_all;
}

Combine_Table::~Combine_Table()
{
  delete [] p_moms;
  for (int k=0;k<m_nampl;++k) {
    delete [] p_legs[k];
    if (p_hard) delete [] p_hard[k];
    if (p_hardc) delete [] p_hardc[k];
  }
  delete [] p_legs;
  if (p_hard) delete [] p_hard;
  if (p_hardc) delete [] p_hardc;
  --s_all;
}

void Leg::DetermineCouplings(const int type) 
{
  m_nqed=m_nqcd=m_pqcd=m_pqed=0;
  AMEGIC::Point *p(p_point);
  if (type==1) p=p->prev;
  if (p->fl.Strong()) ++m_pqcd;
  else ++m_pqed;
  if (p->left->fl.Strong()) ++m_pqcd;
  else ++m_pqed;
  if (p->right->fl.Strong()) ++m_pqcd;
  else ++m_pqed;
  if (m_pqcd==3) ++m_nqcd;
  else ++m_nqed;
  if (p->Lorentz->Type()=="Triangle" ||
      p->Lorentz->Type()=="Box" ||
      p->Lorentz->Type()=="C4GS") m_nqcd+=2;
  m_type=p->Lorentz->Type();
  /*
  msg_Debugging()<<METHOD<<"("<<type<<"): "<<p->fl<<"->"
		 <<p->left->fl<<","<<p->right->fl
		 <<" => n_qcd = "<<m_nqcd<<", m_nqed = "<<m_nqed<<"\n";
  */
}

ATOOLS::Flavour Combine_Table::IsoFlip(const ATOOLS::Flavour &fl) const
{
  switch (fl.Kfcode()) {
  case kf_u: return fl.IsAnti()?Flavour(kf_d).Bar():Flavour(kf_d);
  case kf_d: return fl.IsAnti()?Flavour(kf_u).Bar():Flavour(kf_u);
  case kf_c: return fl.IsAnti()?Flavour(kf_s).Bar():Flavour(kf_s);
  case kf_s: return fl.IsAnti()?Flavour(kf_c).Bar():Flavour(kf_c);
  case kf_t: return fl.IsAnti()?Flavour(kf_b).Bar():Flavour(kf_b);
  case kf_b: return fl.IsAnti()?Flavour(kf_t).Bar():Flavour(kf_t);
  default: break;
  }
  return fl;
}

Flavour Combine_Table::MatchFlavour(const Leg &a,const Leg &b,const Leg &c,int mode) const
{
  return p_proc->ReMap(a.Point()->fl,a.Point()->GetPropID());
}

Leg Combine_Table::CombinedLeg(Leg *legs,const int i,const int j)
{
  Leg & a=legs[i], & b=legs[j], mo;
  /*
  msg_Debugging()<<METHOD<<" "<<(a.Point()->prev?
				 a.Point()->prev->fl:Flavour(kf_p_plus))
		 <<" ("<<(b.Point()->prev?
			 b.Point()->prev->fl:Flavour(kf_p_plus))
		 <<" "<<(b.Point()->prev->left?
			 b.Point()->prev->left->fl:Flavour(kf_p_plus))
		 <<" "<<(b.Point()->prev->right?
			 b.Point()->prev->right->fl:Flavour(kf_p_plus))
		 <<"), a="<<a.Point()->fl<<", b="<<b.Point()->fl<<"\n";
  */
  if ( (a.Point()->prev == b.Point()->prev) && (a.Point()->prev != 0) ) {
    // combinable-type: common mother
    mo.SetPoint(a.Point()->prev);
    mo.DetermineCouplings(0);
    mo.SetMapFlavour(MatchFlavour(mo,a,b,0));
  } 
  else if (a.Point() == b.Point()->left) {
    // combinable-type: a daughter of b
    mo.SetPoint(b.Point()->right);
    mo.DetermineCouplings(1);
    mo.SetMapFlavour(MatchFlavour(mo,b,a,1));
  } 
  else if (a.Point() == b.Point()->right) {
    // combinable-type: a daughter of b
    mo.SetPoint(b.Point()->left);
    mo.DetermineCouplings(1);
    mo.SetMapFlavour(MatchFlavour(mo,b,a,1));
  } 
  else  if (b.Point() == a.Point()->left) {
    // combinable-type: b daughter of a
    mo.SetPoint(a.Point()->right);
    mo.DetermineCouplings(1);
    mo.SetMapFlavour(MatchFlavour(mo,a,b,1));
  } 
  else  if (b.Point() == a.Point()->right) {
    // combinable-type: b daughter of a
    mo.SetPoint(a.Point()->left);
    mo.DetermineCouplings(1);
    mo.SetMapFlavour(MatchFlavour(mo,a,b,1));
  } 
  else THROW(fatal_error,"   Cannot combine legs.");
  mo.SetQCDJets((a.Point()->t<10?a.QCDJets():0)+
		(b.Point()->t<10?b.QCDJets():0));
  /*
  msg_Debugging()<<"mapped flavours: a="
		 <<a.Point()->fl<<"("<<a.MapFlavour()
		 <<")[t="<<a.Point()->t<<",j="<<a.QCDJets()<<"], b="
		 <<b.Point()->fl<<"("<<b.MapFlavour()
		 <<")[t="<<b.Point()->t<<",j="<<b.QCDJets()<<"], c="
		 <<mo.Point()->fl<<"("<<mo.MapFlavour()
		 <<")[t="<<mo.Point()->t<<",j="<<mo.QCDJets()<<"]\n";
  */
  // fix charge incase initial state has wrong
  int icharge;
  if (i<2)  icharge = a.Anti()*a.Point()->fl.IntCharge() - 
	              b.Anti()*b.Point()->fl.IntCharge();
  else      icharge = a.Point()->fl.IntCharge() + b.Point()->fl.IntCharge();
  if (icharge!=mo.Point()->fl.IntCharge()) mo.SetAnti(-1);    
  return mo;
}

void Combine_Table::SetLegScales
(Leg &leg,Leg &legi,Leg &legj,const Vec4D &pi,const Vec4D &pj,
 const double &pt2ij)
{
  /*
  msg_Debugging()<<"set leg scales "<<legi.Point()->fl<<" & "<<legj.Point()->fl
		 <<" -> "<<leg.Point()->fl<<" w/ "<<sqrt(pt2ij)<<" "
		 <<(pi+pj)<<" s = "<<sqrt((pi+pj).Abs2())<<" mt2 = "
		 <<pi.MPerp()<<" / "<<pj.MPerp()<<" / "<<(pi+pj).Mass()<<"\n"; 
  */
  leg.SetKT2(pt2ij);
  // special case: ggh vertex -> order qcd = 2
  // must therefore check for order qed as well
  if (leg.OrderQCD()>0 && leg.OrderQED()==0) leg.SetKT2QCD(pt2ij);
  else if (leg.NQCD()>0) {
    leg.SetKT2QCD(pt2ij);
    if (!leg.Point()->fl.Strong()) 
      leg.SetKT2QCD(Max(dabs((pi+pj).Abs2()),leg.KT2QCD()));
  }
  if (leg.OrderQED()>0) leg.SetKT2QED(pt2ij);
  else if (leg.NQED()>0) {
    leg.SetKT2QCD(pt2ij);
  }
  if (legi.Point()->t<10) {
    if (legj.Point()->t<10) {
      leg.SetMinKT2(Min(legi.MinKT2(),legj.MinKT2()));
      leg.SetMinKT2QCD(Min(legi.MinKT2QCD(),legj.MinKT2QCD()));
      leg.SetMinKT2QED(Min(legi.MinKT2QED(),legj.MinKT2QED()));
    }
    else {
      leg.SetMinKT2(legi.MinKT2());
      leg.SetMinKT2QCD(legi.MinKT2QCD());
      leg.SetMinKT2QED(legi.MinKT2QED());
    }
  }
  else if (legj.Point()->t<10) {
    leg.SetMinKT2(legj.MinKT2());
    leg.SetMinKT2QCD(legj.MinKT2QCD());
    leg.SetMinKT2QED(legj.MinKT2QED());
  }
}
  
Leg * Combine_Table::CombineLegs
(Leg *legs,const int i,const int j,const int nlegs,const double pt2ij) 
{
  Leg * alegs = new Leg[nlegs];
  // assume i < j 
  for (int l=0; l<j; ++l) {
    if (l==i) {
      alegs[i] = CombinedLeg(legs,i,j);
      SetLegScales(alegs[i],legs[i],legs[j],p_moms[i],p_moms[j],pt2ij);
      size_t idi(GetLeg(i).ID()), idj(GetLeg(j).ID()), id(idi|idj);
      alegs[i].SetID(id);
    }
    else {
      alegs[l] = Leg(legs[l]);
    }
  }
  for (int l=j+1; l<=nlegs; ++l) alegs[l-1] = Leg(legs[l]);
  return alegs;
}


bool Combine_Table::CombineMoms(Vec4D *moms,const int _i,const int _j,const int maxl) 
{
  Cluster_Amplitude *ampl(Cluster_Amplitude::New());
  for (int i=0;i<=maxl;++i)
    ampl->CreateLeg(i<2?-moms[i]:moms[i],p_legs[0][i].Flav(),
		    ColorID(),p_legs[0][i].ID());
  Vec4D_Vector after=p_clus->Combine
    (*ampl,m_cdata_winner->first.m_i,m_cdata_winner->first.m_j,
     m_cdata_winner->first.m_k,m_cdata_winner->second.m_mo,p_ms);
  ampl->Delete();
  if (after.empty()) return false;
  for (size_t l=0; l<after.size(); ++l) p_moms[l] = l<2?-after[l]:after[l];
  return true;
}

bool Combine_Table::CombineMoms(Vec4D *moms,const int _i,const int _j,
				     const int maxl,Vec4D *&omoms) 
{
  Cluster_Amplitude *ampl(Cluster_Amplitude::New());
  for (int i=0;i<=maxl;++i)
    ampl->CreateLeg(i<2?-moms[i]:moms[i],p_legs[0][i].Flav(),
		   ColorID(),p_legs[0][i].ID());
  Vec4D_Vector after=p_clus->Combine
    (*ampl,m_cdata_winner->first.m_i,m_cdata_winner->first.m_j,
     m_cdata_winner->first.m_k,m_cdata_winner->second.m_mo,p_ms);
  ampl->Delete();
  if (after.empty()) return false;
  omoms = new Vec4D[maxl];
  for (size_t l=0; l<after.size(); ++l) omoms[l] = l<2?-after[l]:after[l];
  return true;
}

double Combine_Table::MinKt2() const
{
  // return lowest cluster scale, which is not necessarily 
  // scale of first clustering
  double min(std::numeric_limits<double>::max());
  for (int i(0);i<m_nlegs;++i) 
    if (GetLeg(i).Point()->t<10 && 
	GetLeg(i).MinKT2()<min) min=GetLeg(i).MinKT2();
  return min;
}

double Combine_Table::MinKt2QCD(const int cpl) const
{
  // return lowest cluster scale, which is not necessarily 
  // scale of first clustering
  double min(std::numeric_limits<double>::max());
  for (int i(0);i<m_nlegs;++i) 
    if (GetLeg(i).Point()->t<10 && 
	GetLeg(i).MinKT2QCD()<min) {
      if (cpl==0 ||
	  (cpl==1 && GetLeg(i).OrderQCD()>0) ||
	  (cpl==2 && GetLeg(i).NQCD()>1)) min=GetLeg(i).MinKT2QCD();
    }
  return min;
}

double Combine_Table::MinKt2QED(const int cpl) const
{
  // return lowest cluster scale, which is not necessarily 
  // scale of first clustering
  double min(std::numeric_limits<double>::max());
  for (int i(0);i<m_nlegs;++i) 
    if (GetLeg(i).Point()->t<10 && 
	GetLeg(i).MinKT2QED()<min) {
      if (cpl==0 ||
	  (cpl==1 && GetLeg(i).OrderQED()>0) ||
	  (cpl==2 && GetLeg(i).NQED()>1)) min=GetLeg(i).MinKT2QED();
    }
  return min;
}

double Combine_Table::Kt2() const
{
  double min(std::numeric_limits<double>::max());
  for (int i(0);i<m_nlegs;++i) 
    if (GetLeg(i).Point()->t<10 && 
	GetLeg(i).KT2()<min) min=GetLeg(i).KT2();
  return min;
}

double Combine_Table::Kt2QCD(const int cpl) const
{
  double min(std::numeric_limits<double>::max());
  for (int i(0);i<m_nlegs;++i) 
    if (GetLeg(i).Point()->t<10 && 
	GetLeg(i).KT2QCD()<min) {
      if (cpl==0 ||
	  (cpl==1 && GetLeg(i).OrderQCD()>0) ||
	  (cpl==2 && GetLeg(i).NQCD()>1)) min=GetLeg(i).KT2QCD();
    }
  return min;
}

double Combine_Table::Kt2QED(const int cpl) const
{
  double min(std::numeric_limits<double>::max());
  for (int i(0);i<m_nlegs;++i) 
    if (GetLeg(i).Point()->t<10 && 
	GetLeg(i).KT2QED()<min) {
      if (cpl==0 ||
	  (cpl==1 && GetLeg(i).OrderQED()>0) ||
	  (cpl==2 && GetLeg(i).NQED()>1)) min=GetLeg(i).KT2QED();
    }
  return min;
}

double Combine_Table::Sprime() const
{
  if (!p_moms) {
    return 0;
  }
  return (p_moms[0]+p_moms[1]).Abs2();
}

int Combine_Table::AddCouplings(int &nqed,int &nqcd) const
{
  int nqedt(-1), nqcdt(-1);
  for (int i(0);i<m_nampl;++i) {
    int nqedtt(p_hard[i][0].OrderQED()+p_hard[i][1].OrderQED());
    int nqcdtt(p_hard[i][0].OrderQCD()+p_hard[i][1].OrderQCD());
    if (nqedt<0 && nqcdt<0) {
      nqedt=nqedtt;
      nqcdt=nqcdtt;
    }
    else {
      if (nqedt!=nqedtt || nqcdt!=nqcdtt) {
	msg_Tracking()<<METHOD<<"(): Warning. Ambiguous couplings."<<std::endl;
	if (nqcdtt>nqcdt) {
	  msg_Debugging()<<"n_{QCD} = "<<nqcdtt<<" in diagram "
			 <<i<<" -> reset\n";
	  nqedt=nqedtt;
	  nqcdt=nqcdtt;
	}
      }
    }
  }
  nqed=nqedt;
  nqcd=nqcdt;
  return NLegs();
}

bool Combine_Table::Combinable(const Leg &a,const Leg &b,Flavour &mo) const
{
  Leg lmo;
  if (a.Point()->prev!=NULL && a.Point()->prev==b.Point()->prev) {
    lmo.SetPoint((AMEGIC::Point*)a.Point()->prev);
    mo=MatchFlavour(lmo,a,b,0);
    return true;
  }
  if (a.Point()->prev==b.Point()) {
    lmo.SetPoint((AMEGIC::Point*)b.Point());
    mo=MatchFlavour(lmo,b,a,1);
    return true;
  }
  if (b.Point()->prev==a.Point()) {
    lmo.SetPoint((AMEGIC::Point*)a.Point());
    mo=MatchFlavour(lmo,a,b,1);
    return true;
  }
  return false;
}

double Combine_Table::Kt2Hard() const
{
  if (p_hard==NULL) THROW(fatal_error,"No hard legs");
  double min(std::numeric_limits<double>::max());
  for (int i(0);i<m_nampl;++i) {
    for (int j(0);j<2;++j) 
      if (p_hard[i][j].KT2()<min) min=p_hard[i][j].KT2();
  }
  return min;
}

double Combine_Table::Kt2QCDHard(const int cpl) const
{
  if (p_hard==NULL) THROW(fatal_error,"No hard legs");
  double min(std::numeric_limits<double>::max());
  for (int i(0);i<m_nampl;++i) {
    for (int j(0);j<2;++j) 
      if (p_hard[i][j].KT2QCD()<min) {
	if (cpl==0 ||
	    (cpl==1 && p_hard[i][j].OrderQCD()>0) ||
	    (cpl==2 && p_hard[i][j].NQCD()>1)) min=p_hard[i][j].KT2QCD();
      }
  }
  return min;
}

double Combine_Table::Kt2QEDHard(const int cpl) const
{
  if (p_hard==NULL) THROW(fatal_error,"No hard legs");
  double min(std::numeric_limits<double>::max());
  for (int i(0);i<m_nampl;++i) {
    for (int j(0);j<2;++j) 
      if (p_hard[i][j].KT2QED()<min) {
	if (cpl==0 ||
	    (cpl==1 && p_hard[i][j].OrderQCD()>0) ||
	    (cpl==2 && p_hard[i][j].NQCD()>1)) min=p_hard[i][j].KT2QED();
      }
  }
  return min;
}

double Combine_Table::GetWinner(int &i,int &j,int &k)
{ 
  i=m_cdata_winner->first.m_i; 
  j=m_cdata_winner->first.m_j;
  k=m_cdata_winner->first.m_k;
  if (m_cdata_winner->second.p_down!=NULL) {
    double kt2qcd(m_cdata_winner->second.p_down->GetLeg(i).KT2QCD());
    if (kt2qcd<std::numeric_limits<double>::max()) return sqrt(kt2qcd);
    return sqrt(m_cdata_winner->second.p_down->GetLeg(i).KT2());
  }
  THROW(fatal_error,"Legs not combined. No Scale information");
  return 0.0;
}

void Combine_Table::AddPossibility(const int i,const int j,const int k,
				   const int ngraph,const Flavour &mo) 
{
  CD_List::iterator cit=m_combinations.find(Combine_Key(i,j,k));
  if (cit!=m_combinations.end()) {
    cit->second.m_graphs.push_back(ngraph);
    cit->second.m_strong=Max(cit->second.m_strong,
			     CombinedLeg(p_legs[ngraph],i,j).OrderQCD());
  }
  else {
    Combine_Data cd(0.,ngraph);
    cd.m_mo=mo;
    cd.m_strong=CombinedLeg(p_legs[ngraph],i,j).OrderQCD();
    m_combinations[Combine_Key(i,j,k)]=cd;
  }
}

void Combine_Table::FillTable(Leg **legs,const int nlegs,const int nampl)
{
  // store information
  p_legs=legs;
  m_nlegs=nlegs;
  m_nampl=nampl;
  Flavour mo;
  // determine possible combinations and corresponding y_ij  if nlegs>4
  if (m_nlegs>4) {
    int start=0;
    // cluster initial state only if isrshower and isr_x is on. 
    if (!legs[0][0].Flav().Strong() && !legs[0][1].Flav().Strong()) start=2;
    for (int i=start; i<m_nlegs; ++i) {  
//       if (!m_isr1on && i==0) i=1;
//       if (!m_isr2on && i==1) i=2;
      for (int j=i+1; j<m_nlegs; ++j) {
	// never combine "0&1" !
	if (j==1) j=2;
	// check if leg i is combinable with leg j in any graph
	for (int k=0;k<m_nampl;++k) {
// 	  msg_Debugging()<<"start w/ "<<k<<", "
// 			 <<i<<": "<<p_legs[k][i].MapFlavour()<<"\n";
	  if (Combinable(p_legs[k][i],p_legs[k][j],mo)) {
	    int sci(p_legs[k][i].Flav().StrongCharge());
	    int scj(p_legs[k][j].Flav().StrongCharge());
	    for (int l=0;l<m_nlegs;++l)
	      if (l!=i && l!=j) {
		int sc(p_legs[k][l].Flav().StrongCharge());
		if (((sci==8 || scj==8 || sc==8) && 
		     (sci!=0 && scj!=0 && sc!=0)) ||
		    (sci!=8 && scj!=8 && sc!=8)) AddPossibility(i,j,l,k,mo);
	      }
	  }
	}
      }
    }
  }
}

CD_List::iterator Combine_Table::CalcPropagator(CD_List::iterator &cit)
{
  if (cit->first.m_flav.Kfcode()==kf_none) {
    cit->second.m_sij=(p_moms[cit->first.m_i]+p_moms[cit->first.m_j]).Abs2();
    Cluster_Amplitude *ampl(Cluster_Amplitude::New());
    for (int i=0;i<m_nlegs;++i)
      ampl->CreateLeg(i<2?-p_moms[i]:p_moms[i],p_legs[0][i].Flav(),
		     ColorID(),p_legs[0][i].ID());
    cit->second.m_pt2ij=p_clus->KPerp2
      (*ampl,cit->first.m_i,cit->first.m_j,cit->first.m_k,cit->second.m_mo,p_ms);
    msg_Debugging()<<"Calculate m_perp("<<cit->first.m_i<<"["
		   <<p_legs[0][cit->first.m_i].Flav()<<"],"
		   <<cit->first.m_j<<"["<<p_legs[0][cit->first.m_j].Flav()<<"],"
		   <<cit->first.m_k<<"["<<p_legs[0][cit->first.m_k].Flav()
		   <<"],"<<cit->second.m_mo<<") -> "<<cit->second.m_pt2ij<<std::endl;
    ampl->Delete();
    return cit;
  }
  else {
    CD_List::iterator father= 
      m_combinations.find(Combine_Key(cit->first.m_i,cit->first.m_j,
				      cit->first.m_k));
    if (father!=m_combinations.end()) {
      cit->second.m_pt2ij = father->second.m_pt2ij;
      return father;
    }
    else THROW(fatal_error,"No father.");
  }
  return m_combinations.end();
}

Combine_Table *Combine_Table::
CalcJet(int nl,const double x1,const double x2,
	ATOOLS::Vec4D * moms,const size_t mode) 
{
  if (p_up==NULL) {
    m_x1 = x1;
    m_x2 = x2;
  }
  m_rejected.clear();
  while (true) {
    bool did_boost(InitStep(moms,nl));
    if (!SelectWinner(did_boost)) {
      if (mode==1 ||
	  (nl==4 && (IdentifyHardProcess() || p_up==NULL))) {
	if (did_boost) for (size_t i=0;i<m_nl;++i) p_moms[i]=p_save_moms[i]; 
	delete [] p_save_moms;
	return this;
      }
      delete [] p_save_moms;
      delete this;
      return NULL;
    }
    // if number of legs is still greater 4 Cluster once more
    // if number of legs equals 4, determine end situation
    if (nl<4) THROW(fatal_error,"nlegs < min. Abort.");
    Combine_Table *tab(CreateNext(did_boost));
    if (tab!=NULL) {
      Combine_Table *next(NextTable(tab,x1,x2));
      if (next!=NULL) return next;
    }
    m_rejected.insert(m_cdata_winner->first);
    msg_Debugging()<<METHOD<<"(): Table "<<m_no<<": reject winner "
		   <<m_cdata_winner->first<<"\n";
  }
  return NULL;
}

bool Combine_Table::InitStep(ATOOLS::Vec4D *moms,const int nl)
{
  m_nl=nl;
  // change momenta to actual values    
  if (moms!=0) for (size_t l=0;l<m_nl;++l) p_moms[l]=moms[l];
  // boost in CMS frame and rotate to z-axis (store old moms)
  p_save_moms = new Vec4D[m_nl];
  for (size_t i=0;i<m_nl;++i) p_save_moms[i] = p_moms[i];
  bool did_boost(false);
  if (!(Vec3D(p_moms[0])==Vec3D(-1.*p_moms[1]))) {
    Poincare cms, zaxis;
    bool dir(p_moms[0][3]>0.0);
    if (p_moms[0].PSpat2()<p_moms[1].PSpat2()) dir=p_moms[1][3]<0.0;
    cms   = Poincare(p_moms[0]+p_moms[1]);
    for (size_t i=0;i<m_nl;++i) cms.Boost(p_moms[i]);
    
    if (dir) zaxis = Poincare(p_moms[0],Vec4D::ZVEC);
    else zaxis = Poincare(p_moms[1],Vec4D::ZVEC);
    
    for (size_t i=0;i<m_nl;++i) zaxis.Rotate(p_moms[i]);
    did_boost = true;
  }
  return did_boost;
}

bool Combine_Table::SelectWinner(const bool did_boost)
{
  CD_List & cl(m_combinations);
  if (cl.size()==0) return false;
  // calculate pt2ij and determine "best" combination
  m_cdata_winner = cl.end();
  CD_List::iterator qcd_winner(cl.end()), nqcd_winner(cl.end());
  double kt2(std::numeric_limits<double>::max()), kt2qcd(kt2), kt2nqcd(kt2);
  for (CD_List::iterator cit(cl.begin()); cit!=cl.end(); ++cit) {
    CD_List::iterator tit(CalcPropagator(cit));
    double pt2ij(cit->second.m_pt2ij.m_op2);
    if (cit->second.m_graphs.size()==0) continue;
    if (m_rejected.find(cit->first)==m_rejected.end()) {
      if (pt2ij<kt2) {
	m_cdata_winner=cit;
	kt2=pt2ij;
      }
      if (cit->second.m_strong && pt2ij<kt2qcd) {
	qcd_winner=cit;
	kt2qcd=pt2ij;
      }
      if (pt2ij<kt2nqcd &&
	  ((p_legs[0][cit->first.m_i].Flav().Strong() &&
	    cit->first.m_i>1)||
 	   (p_legs[0][cit->first.m_j].Flav().Strong() &&
	    cit->first.m_j>1))) {
 	nqcd_winner=cit;
 	kt2nqcd=pt2ij;
      }
    }
  }
  if (nqcd_winner!=cl.end()) m_cdata_winner=nqcd_winner;
  if (qcd_winner!=cl.end()) m_cdata_winner=qcd_winner;
  return m_cdata_winner!=cl.end();
}

bool Combine_Table::TestMomenta(const int i,const int j)
{
  // check if combined momenta are physical in their cms
  Vec4D s1 = p_save_moms[i] - p_save_moms[j];
  Vec4D s2 = p_save_moms[1-i];
  Poincare test(s1+s2);
  test.Boost(s1);
  test.Boost(s2);
  // do not check energies individually, but cms energy 
  return (s1[0]+s2[0])>0.0;
}

Combine_Table *Combine_Table::CreateNext(bool did_boost)
{
  // check if boosted  (restore saved moms)
  if (did_boost) { 
    for (size_t i=0;i<m_nl;++i) p_moms[i]=p_save_moms[i]; 
  }
  delete [] p_save_moms;
  --m_nl;
  if (!m_cdata_winner->second.p_down) {
    Vec4D * amoms;
    // generate new momenta
    if (!CombineMoms(p_moms,m_cdata_winner->first.m_i,
		     m_cdata_winner->first.m_j,m_nl,amoms)) return NULL;
    Leg ** alegs = new Leg*[m_cdata_winner->second.m_graphs.size()];
    for (size_t k=0;k<m_cdata_winner->second.m_graphs.size();++k) {
      alegs[k] = CombineLegs
	(p_legs[m_cdata_winner->second.m_graphs[k]],
	 m_cdata_winner->first.m_i,m_cdata_winner->first.m_j,m_nl,
	 m_cdata_winner->second.m_pt2ij.m_kt2);
    }
    m_cdata_winner->second.p_down = 
      new Combine_Table(p_proc,p_ms,p_clus,amoms,this);
    // initialise Combine_Table
    m_cdata_winner->second.p_down->FillTable(alegs,m_nl,m_cdata_winner->second.m_graphs.size());
  } 
  else {
    if (!((Combine_Table*)m_cdata_winner->second.p_down)->
	CombineMoms(p_moms,m_cdata_winner->first.m_i,m_cdata_winner->first.m_j,m_nl)) return NULL;
  }

  // update x1,2 (before calling CalcJet again
  Combine_Table *tab((Combine_Table*)m_cdata_winner->second.p_down);
  if (m_cdata_winner->first.m_i<2) {
    double z=tab->Sprime()/Sprime();
    if (m_cdata_winner->first.m_i==0) {
      tab->m_x1=m_x1*z;
      tab->m_x2=m_x2;
    }
    else {
      tab->m_x1=m_x1;
      tab->m_x2=m_x2*z;
    }
  }
  else {
    tab->m_x1=m_x1;
    tab->m_x2=m_x2;
  }
  return tab;
}

Combine_Table *Combine_Table::NextTable(Combine_Table *tab,
					const double x1,const double x2)
{
  Combine_Table* ct = tab->CalcJet(m_nl,x1,x2);
  if (ct!=NULL) m_graph_winner=tab->m_graph_winner;
  else m_cdata_winner->second.p_down=NULL;
  // translate back
  m_graph_winner = m_cdata_winner->second.m_graphs.front();
  return (Combine_Table*)ct;
}

bool Combine_Table::IdentifyHardProcess()
{
  msg_Debugging()<<METHOD<<"():\n";
  msg_Indent();
  m_nstrong=0;
  Flavour mo;
  if (p_hard==NULL) {
    p_hard = new Leg*[m_nampl];
    for (int i(0);i<m_nampl;++i) p_hard[i] = new Leg[2];
    p_hardc = new int*[m_nampl];
    for (int i(0);i<m_nampl;++i) p_hardc[i] = new int[4];
  }
  for (int i(0);i<m_nampl;++i) {
    if (Combinable(p_legs[i][0],p_legs[i][1],mo) &&
	Combinable(p_legs[i][2],p_legs[i][3],mo)) {
      double pt2ij1((p_moms[0]+p_moms[1]).Abs2());
      double pt2ij2((p_moms[0]+p_moms[1]).Abs2());
      msg_Debugging()<<"s-channel pt = "<<sqrt(pt2ij1)
		     <<" / "<<sqrt(pt2ij2)<<", m = "
		     <<sqrt(dabs((p_moms[0]+p_moms[1]).Abs2()))<<", "
		     <<p_legs[i][0].Flav()<<" "<<p_legs[i][1].Flav()
		     <<" -> "<<p_legs[i][2].Flav()<<" "
		     <<p_legs[i][3].Flav()<<"\n";
      p_hard[i][0]=CombinedLeg(p_legs[i],0,1);
      SetLegScales(p_hard[i][0],p_legs[i][0],p_legs[i][1],
		   p_moms[0],p_moms[1],pt2ij1);
      p_hard[i][1]=CombinedLeg(p_legs[i],2,3);
      SetLegScales(p_hard[i][1],p_legs[i][2],p_legs[i][3],
		   p_moms[2],p_moms[3],pt2ij2);
      p_hardc[i][0]=0;
      p_hardc[i][1]=0;
      p_hardc[i][2]=1;
      p_hardc[i][3]=1;
    }
    else if (Combinable(p_legs[i][0],p_legs[i][2],mo) &&
	     Combinable(p_legs[i][1],p_legs[i][3],mo)) {
      double pt2ij1((p_moms[0]-p_moms[2]).Abs2());
      double pt2ij2((p_moms[0]-p_moms[2]).Abs2());
      msg_Debugging()<<"t-channel pt = "<<sqrt(pt2ij1)
		     <<" / "<<sqrt(pt2ij2)<<", m = "
		     <<sqrt(dabs((p_moms[0]+p_moms[2]).Abs2()))<<", "
		     <<p_legs[i][0].Flav()<<" "<<p_legs[i][1].Flav()
		     <<" -> "<<p_legs[i][2].Flav()<<" "
		     <<p_legs[i][3].Flav()<<"\n";
      p_hard[i][0]=CombinedLeg(p_legs[i],0,2);
      SetLegScales(p_hard[i][0],p_legs[i][0],p_legs[i][2],
		   p_moms[0],p_moms[2],pt2ij1);
      p_hard[i][1]=CombinedLeg(p_legs[i],1,3);
      SetLegScales(p_hard[i][1],p_legs[i][1],p_legs[i][3],
		   p_moms[1],p_moms[3],pt2ij2);
      p_hardc[i][0]=0;
      p_hardc[i][1]=1;
      p_hardc[i][2]=0;
      p_hardc[i][3]=1;
    }
    else if (Combinable(p_legs[i][0],p_legs[i][3],mo) &&
	     Combinable(p_legs[i][1],p_legs[i][2],mo)) {
      double pt2ij1((p_moms[0]-p_moms[3]).Abs2());
      double pt2ij2((p_moms[0]-p_moms[3]).Abs2());
      msg_Debugging()<<"u-channel pt = "<<sqrt(pt2ij1)
		     <<" / "<<sqrt(pt2ij2)<<", m = "
		     <<sqrt(dabs((p_moms[0]+p_moms[3]).Abs2()))<<", "
		     <<p_legs[i][0].Flav()<<" "<<p_legs[i][1].Flav()
		     <<" -> "<<p_legs[i][2].Flav()<<" "
		     <<p_legs[i][3].Flav()<<"\n";
      p_hard[i][0]=CombinedLeg(p_legs[i],0,3);
      SetLegScales(p_hard[i][0],p_legs[i][0],p_legs[i][3],
		   p_moms[0],p_moms[3],pt2ij1);
      p_hard[i][1]=CombinedLeg(p_legs[i],1,2);
      SetLegScales(p_hard[i][1],p_legs[i][1],p_legs[i][2],
		   p_moms[1],p_moms[2],pt2ij2);
      p_hardc[i][0]=0;
      p_hardc[i][1]=1;
      p_hardc[i][2]=1;
      p_hardc[i][3]=0;
    }
    else THROW(fatal_error,"No match for hard process.");
    if (p_hard[i][0].Point()->t>=10) {
      msg_Debugging()<<"cut propagator "<<p_hard[i][0].Flav()<<"\n";
      return false;
    }
    m_nstrong=Max(m_nstrong,p_hard[i][0].OrderQCD()+p_hard[i][1].OrderQCD());
  }
  return true;
}

int Combine_Table::IdentifyHardPropagator() const
{
  msg_Debugging()<<METHOD<<"():\n";
  msg_Indent();
  Flavour mo;
  int channel(-1);
  for (int i(0);i<m_nampl;++i) {
    if (Combinable(p_legs[i][0],p_legs[i][1],mo) &&
	Combinable(p_legs[i][2],p_legs[i][3],mo)) {
      msg_Debugging()<<"s-channel\n";
      if (channel<0) channel=1;
      else if (channel!=1) return -1;
    }
    else if (Combinable(p_legs[i][0],p_legs[i][2],mo) &&
	     Combinable(p_legs[i][1],p_legs[i][3],mo)) {
      msg_Debugging()<<"t-channel\n";
      if (channel<0) channel=2;
      else if (channel!=2) return -1;
    }
    else if (Combinable(p_legs[i][0],p_legs[i][3],mo) &&
	     Combinable(p_legs[i][1],p_legs[i][2],mo)) {
      msg_Debugging()<<"u-channel\n";
      if (channel<0) channel=3;
      else if (channel!=3) return -1;
    }
    else THROW(fatal_error,"No match for hard process.");
  }
  return channel;
}

