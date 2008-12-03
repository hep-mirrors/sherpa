#include "Combine_Table_Base.H"

#include "Exception.H"
#include "My_Limits.H"
#include "Jet_Finder.H"
#include "Blob.H"
#include <iomanip>

#define MAXD std::numeric_limits<double>::max()

using namespace SHERPA;
using namespace AMEGIC;
using namespace MODEL;
using namespace ATOOLS;

int Combine_Table_Base::s_all(0);

Leg::Leg(AMEGIC::Point *const point,const int anti):  
  p_point(point), m_anti(anti), 
  m_nqcd(0), m_nqed(0), m_pqcd(0), m_pqed(0), m_ext(0), 
  m_qcdjets(point!=NULL?point->fl.Strong():0), m_id(0),
  m_kt2(std::numeric_limits<double>::max()), 
  m_kt2qcd(m_kt2), m_kt2qed(m_kt2), 
  m_minkt2(m_kt2), m_minkt2qcd(m_kt2), m_minkt2qed(m_kt2), 
  p_qmin(NULL) 
{
  m_q2cut[2]=m_q2cut[1]=m_q2cut[0]=ATOOLS::sqr(ATOOLS::rpa.gen.Ecms());
}

Leg::Leg(const Leg &leg): 
  p_point(leg.p_point), m_anti(leg.m_anti), 
  m_nqcd(0), m_nqed(0), m_pqcd(0), m_pqed(0), m_ext(leg.m_ext), 
  m_qcdjets(leg.m_qcdjets), m_id(leg.m_id),
  m_kt2(std::numeric_limits<double>::max()), 
  m_kt2qcd(m_kt2), m_kt2qed(m_kt2), m_minkt2(leg.m_minkt2),  
  m_minkt2qcd(leg.m_minkt2qcd), m_minkt2qed(leg.m_minkt2qed), 
  p_qmin(leg.p_qmin), m_mapfl(leg.m_mapfl) 
{
  m_q2cut[0]=leg.m_q2cut[0];
  m_q2cut[1]=leg.m_q2cut[1];
  m_q2cut[2]=leg.m_q2cut[2];
}

std::ostream &SHERPA::operator<<
  (std::ostream &str,const std::vector<int> &info)
{
  str<<"(";
  if (info.size()>0) str<<info[0];
  else str<<"<no entry>";
  for (size_t i=1;i<info.size();++i) str<<","<<info[i];
  return str<<")";
}

std::ostream &SHERPA::operator<<(std::ostream &ostr,const Leg &leg)
{
  return ostr<<leg.p_point<<" "<<leg.m_anti;
}

// ============================================================
//    class Combine_Key
// ============================================================

std::ostream& SHERPA::operator<<(std::ostream& s,const Combine_Key &ck)
{
  s<<" "<<ck.m_i<<"&"<<ck.m_j;
  if (ck.m_flav.Kfcode()!=kf_none) s<<"["<<std::setw(6)<<ck.m_flav<<"]";
  else s<<std::string(8,' ');
  return s;
}

Combine_Key::Combine_Key(): 
  m_i(0), m_j(0) {}

Combine_Key::Combine_Key(const int i,const int j,const Flavour &flav) :
  m_i(i), m_j(j), m_flav(flav) {}

bool SHERPA::operator<(const Combine_Key & a,const Combine_Key & b) 
{
  if (a.m_i < b.m_i) return true;
  if (a.m_i > b.m_i) return false;
  if (a.m_j < b.m_j) return true;
  if (a.m_j > b.m_j) return false;
  if (a.m_flav.Kfcode() > b.m_flav.Kfcode()) return true;
  return false;
}

// ============================================================
//    class Combine_Data
// ============================================================

std::ostream& SHERPA::operator<<(std::ostream &s,const Combine_Data &cd)
{
  s<<" "<<std::setw(10)<<sqrt(cd.m_pt2ij)<<" "<<std::setw(10)<<sqrt(dabs(cd.m_sij))<<" "
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
//    class Combine_Table_Base
// ============================================================

std::ostream& SHERPA::operator<<(std::ostream& s ,const Combine_Table_Base & ct) 
{
  if (&ct) {
    s<<std::endl<<" Combine_Table_Base ("<<&ct<<") "<<ct.m_no<<" (up=";
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
       <<sqrt(ct.GetLeg(l).Q2Cut())<<std::setw(12)
       <<sqrt(ct.GetLeg(l).Q2Cut(2))<<std::setw(12)
       <<sqrt(ct.GetLeg(l).Q2Cut(1))<<std::setw(12)
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
    s<<"***empty Combine_Table_Base***"<<std::endl;
  return s;
}

Combine_Table_Base::Combine_Table_Base(Jet_Finder *jf,Vec4D *moms, Combine_Table_Base *up,
				       int isrmode, int isrshoweron):
  m_nstrong(0),
  m_isr1on(isrmode&1), m_isr2on((isrmode&2)/2), m_isrshoweron(isrshoweron), 
  m_graph_winner(0), 
  p_up(up), p_legs(0), p_jf(jf), p_moms(moms), p_hard(NULL), p_hardc(NULL)
{
  m_no=++s_all;
}

Combine_Table_Base::~Combine_Table_Base()
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

ATOOLS::Flavour Combine_Table_Base::IsoFlip(const ATOOLS::Flavour &fl)
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

Flavour Combine_Table_Base::MatchFlavour(Leg &a,Leg &b,Leg &c,int mode)
{
  Flavour fla(a.Point()->fl), flb(b.Point()->fl), flc(c.Point()->fl);
  /*
  msg_Debugging()<<"have bm="<<b.MapFlavour()<<", cm="<<c.MapFlavour()
		 <<" mode "<<mode<<" -> ";
  */
  Flavour fl(fla);
  if (fla.IsQuark() && !(flb.Strong() && flc.Strong())) {
    if (flb.Kfcode()==kf_Wplus || flc.Kfcode()==kf_Wplus ||
        flb.Kfcode()==kf_Hplus || flc.Kfcode()==kf_Hplus)
      fl=IsoFlip((flb.IsQuark()?b:c).MapFlavour());
    else fl=Flavour((flb.IsQuark()?b:c).MapFlavour());
  }
  else if (fla.Kfcode()==flb.Kfcode()) {
    fl=b.MapFlavour();
    if (fla.IsAnti()^flb.IsAnti()) fl=fl.Bar();
  }
  else if (fla.Kfcode()==flc.Kfcode()) {
    fl=c.MapFlavour();
    if (fla.IsAnti()^flc.IsAnti()) fl=fl.Bar();
  }
  /*
  msg_Debugging()<<"match "<<fla<<" "<<flb<<" "
		 <<flc<<" -> "<<fl<<"\n";
  */
  return fl;
}

Leg Combine_Table_Base::CombinedLeg(Leg *legs,const int i,const int j)
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

void Combine_Table_Base::SetLegScales
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
  
void Combine_Table_Base::SetQ2Cut(const size_t &i,const size_t &j,
				  double &q2in,double &q2out,const size_t &id)
{
  if (!p_jf->YCutsSet()) return;
  size_t idi(GetLeg(i).ID()), idj(GetLeg(j).ID());
  double yci(id==1?p_jf->GetScaledYcut(idi,idi):
	     p_jf->GetScaledGlobalYcut(idi,idi));
  double ycj(id==1?p_jf->GetScaledYcut(idj,idj):
	     p_jf->GetScaledGlobalYcut(idj,idj));
  double sq2i(GetLeg(i).Q2Cut(id)/(i<2?1.0:sqr(p_jf->DeltaR())));
  double sq2j(GetLeg(j).Q2Cut(id)/(j<2?1.0:sqr(p_jf->DeltaR())));
  q2in=q2out=-1.0;
  if (yci>0.0) {
    if (IsEqual(yci,ycj)) q2out=yci*sqr(rpa.gen.Ecms());
    else if (IsEqual(sq2j,yci*sqr(rpa.gen.Ecms())) && ID(idj).size()>1) q2out=sq2j;
  }
  else if (ycj>0.0) {
    if (IsEqual(sq2i,ycj*sqr(rpa.gen.Ecms())) && ID(idi).size()>1) q2out=sq2i;
  }
  else {
    if (ID(idi).size()>1 && ID(idj).size()>1 &&	IsEqual(sq2i,sq2j)) q2out=sq2i; 
  }
  if (i>1) q2out*=sqr(p_jf->DeltaR());
  if (q2out<0.0) {
    Combine_Table_Base *ct(this);
    while (ct->Up()) ct=ct->Up();
    msg_Error()<<METHOD<<"(): Current table {\n"<<*ct<<"} -> combine "
	       <<i<<" "<<ID(idi)<<"["<<yci<<"] & "<<j<<" "<<ID(idj)
	       <<"["<<ycj<<"]"<<std::endl;
    THROW(critical_error,"Inconsistent Q_{cut} values");
  }
  q2in=q2out;
  double gyc(p_jf->GetGlobalYcut(idi|idj,idi|idj));
  if (gyc>0.0) q2in=gyc*sqr(rpa.gen.Ecms());
}

Leg * Combine_Table_Base::CombineLegs
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
      double q2in, q2out;
      SetQ2Cut(i,j,q2in,q2out,0);
      alegs[i].SetQ2Cut(q2in);
      alegs[i].SetQ2Cut(q2out,2);
      SetQ2Cut(i,j,q2in,q2out,1);
      alegs[i].SetQ2Cut(q2out,1);
    }
    else {
      alegs[l] = Leg(legs[l]);
    }
  }
  for (int l=j+1; l<=nlegs; ++l) alegs[l-1] = Leg(legs[l]);
  return alegs;
}


void Combine_Table_Base::CombineMoms(Vec4D *moms,const int i,const int j,const int maxl) 
{
  // assume i < j
  for (int l=0; l<j; ++l) {
    if (l==i) { 
      if (i<2) p_moms[i] = moms[i] - moms[j];      
      else     p_moms[i] = moms[i] + moms[j];
    }
    else       p_moms[l] = moms[l];
  }
  for (int l=j+1; l<=maxl; ++l) p_moms[l-1] = moms[l];
}

void Combine_Table_Base::CombineMoms(Vec4D *moms,const int i,const int j,
				     const int maxl,Vec4D *&omoms) 
{
  omoms = new Vec4D[maxl];
  // assume i < j
  for (int l=0; l<j; ++l) {
    if (l==i) {
      if (i<2) omoms[i] = moms[i] - moms[j];      
      else     omoms[i] = moms[i] + moms[j];
    }
    else       omoms[l] = moms[l];
  }
  for (int l=j+1; l<=maxl; ++l) omoms[l-1]=moms[l];
}

double Combine_Table_Base::MinKt2() const
{
  // return lowest cluster scale, which is not necessarily 
  // scale of first clustering
  double min(std::numeric_limits<double>::max());
  for (int i(0);i<m_nlegs;++i) 
    if (GetLeg(i).Point()->t<10 && 
	GetLeg(i).MinKT2()<min) min=GetLeg(i).MinKT2();
  return min;
}

double Combine_Table_Base::MinKt2QCD(const int cpl) const
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

double Combine_Table_Base::MinKt2QED(const int cpl) const
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

double Combine_Table_Base::Kt2() const
{
  double min(std::numeric_limits<double>::max());
  for (int i(0);i<m_nlegs;++i) 
    if (GetLeg(i).Point()->t<10 && 
	GetLeg(i).KT2()<min) min=GetLeg(i).KT2();
  return min;
}

double Combine_Table_Base::Kt2QCD(const int cpl) const
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

double Combine_Table_Base::Kt2QED(const int cpl) const
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

double Combine_Table_Base::Sprime() const
{
  if (!p_moms) {
    return 0;
  }
  return (p_moms[0]+p_moms[1]).Abs2();
}

int Combine_Table_Base::AddCouplings(int &nqed,int &nqcd) const
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

bool Combine_Table_Base::Combinable(const Leg &a,const Leg &b) const
{
  if ((a.Point()->prev!=NULL && a.Point()->prev==b.Point()->prev) ||
      a.Point()->prev==b.Point() ||
      b.Point()->prev==a.Point()) {
    return true;
  }
  return false;
}

double Combine_Table_Base::Kt2Hard() const
{
  if (p_hard==NULL) THROW(fatal_error,"No hard legs");
  double min(std::numeric_limits<double>::max());
  for (int i(0);i<m_nampl;++i) {
    for (int j(0);j<2;++j) 
      if (p_hard[i][j].KT2()<min) min=p_hard[i][j].KT2();
  }
  return min;
}

double Combine_Table_Base::Kt2QCDHard(const int cpl) const
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

double Combine_Table_Base::Kt2QEDHard(const int cpl) const
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

void Combine_Table_Base::ShuffleMomenta()
{
}
