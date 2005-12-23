#include "Combine_Table_Base.H"

#include "Exception.H"
#include "My_Limits.H"
#include <iomanip>

using namespace SHERPA;
using namespace AMEGIC;
using namespace ATOOLS;

int Combine_Table_Base::s_all(0);

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
  if (ck.m_flav.Kfcode()!=kf::none) s<<"["<<std::setw(6)<<ck.m_flav<<"]";
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
  s<<" "<<std::setw(8)<<sqrt(cd.m_pt2ij)<<" "<<std::setw(8)<<sqrt(dabs(cd.m_sij))<<" "
   <<std::setw(3)<<cd.m_strong;
  s<<" ("<<std::setw(11)<<cd.m_prop<<","<<std::setw(11)<<cd.m_coupling<<") ";
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
    s<<std::endl<<" Combine_Table_Base "<<ct.m_no<<" (up=";
    if (ct.p_up) s<<ct.p_up->m_no<<")"<<std::endl; else s<<"#)"<<std::endl;
    s<<" x1="<<ct.m_x1<<" x2="<<ct.m_x2<<" kt2_min="<<ct.m_kt2min
     <<" kt2_min(QCD)="<<ct.m_kt2QCD<<" kt2_min(QED)="<<ct.m_kt2QED<<std::endl;
    s<<" ==============="<<std::endl;
    s<<"moms="<<ct.p_moms<<std::endl;
    for (int l=0; l<ct.m_nlegs; ++l) 
      s<<" "<<l<<":"<<ct.p_moms[l]<<" "<<ct.p_moms[l].Abs2()<<" "<<ct.p_legs[0][l].Point()->fl<<std::endl;
    s<<" ---------------"<<std::endl;
    const CD_List & cl=ct.m_combinations;
    if (cl.size()>0) {
      s<<"     with "<<cl.size()<<" combinations"<<std::endl;
      for (CD_List::const_iterator cit=cl.begin(); cit!=cl.end(); ++cit) {
 	s<<cit->first<<" "<<std::setw(8)
	 <<ct.p_moms[cit->first.m_i].
	  Theta(ct.p_moms[cit->first.m_j])<<cit->second<<std::endl; 
      }
      for (CD_List::const_iterator cit=cl.begin(); cit!=cl.end(); ++cit) {
	if (cit->second.p_down) {
	  s<<*cit->second.p_down<<std::endl;
	}
      }
      // test output
      for (int k=0; k<ct.m_nampl; ++k) {
	for (int l=0; l<ct.m_nlegs; ++l) 
	  s<<" "<<ct.p_legs[k][l].Point()->fl<<"("<<ct.p_legs[k][l].Anti()<<")  ";
	s<<std::endl;
      }

    }
    else {
      for (int k=0; k<ct.m_nampl; ++k) {
	for (int l=0; l<ct.m_nlegs; ++l) 
	  s<<" "<<ct.p_legs[k][l].Point()->fl<<"("<<ct.p_legs[k][l].Anti()<<")  ";
	s<<std::endl;
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
  m_graph_winner(0), m_kt2min(std::numeric_limits<double>::max()), 
  m_kt2QCD(m_kt2min), m_kt2QED(m_kt2min),
  p_up(up), p_legs(0), p_jf(jf), p_moms(moms), p_hard(NULL)
{
  m_no=++s_all;
}

Combine_Table_Base::~Combine_Table_Base()
{
  delete [] p_moms;
  for (int k=0;k<m_nampl;++k) {
    delete [] p_legs[k];
    if (p_hard) delete [] p_hard[k];
  }
  delete [] p_legs;
  if (p_hard) delete [] p_hard;
  --s_all;
}

void Leg::DetermineCouplings(const int type) 
{
  m_nqed=m_nqcd=0;
  AMEGIC::Point *p(p_point);
  if (type==1) p=p->prev;
  if (p->fl.Strong() &&
      p->left->fl.Strong() && 
      p->right->fl.Strong()) ++m_nqcd;
  else ++m_nqed;
  switch (p->Lorentz->Type()) {
  case lf::Triangle:
  case lf::Box:
  case lf::C4GS:
    m_nqcd+=2;
    break;
  default:
    break;
  }
  m_type=p->Lorentz->Type();
  msg_Debugging()<<METHOD<<"("<<type<<"): "<<p->fl<<"->"
		 <<p->left->fl<<","<<p->right->fl
		 <<" => n_qcd = "<<m_nqcd<<", m_nqed = "<<m_nqed<<"\n";
}

Leg Combine_Table_Base::CombinedLeg(Leg *legs,const int i,const int j)
{
  Leg & a=legs[i], & b=legs[j], mo;
  msg_Debugging()<<METHOD<<" "<<(a.Point()->prev?
				 a.Point()->prev->fl:Flavour(kf::p_plus))
		 <<" ("<<(b.Point()->prev?
			 b.Point()->prev->fl:Flavour(kf::p_plus))
		 <<" "<<(b.Point()->prev->left?
			 b.Point()->prev->left->fl:Flavour(kf::p_plus))
		 <<" "<<(b.Point()->prev->right?
			 b.Point()->prev->right->fl:Flavour(kf::p_plus))
		 <<"), a="<<a.Point()->fl<<", b="<<b.Point()->fl<<"\n";
  if ( (a.Point()->prev == b.Point()->prev) && (a.Point()->prev != 0) ) {
    // combinable-type: common mother
    mo.SetPoint(a.Point()->prev);
    mo.DetermineCouplings(0);
  } 
  else if (a.Point() == b.Point()->left) {
    // combinable-type: a daughter of b
    mo.SetPoint(b.Point()->right);
    mo.DetermineCouplings(1);
  } 
  else if (a.Point() == b.Point()->right) {
    // combinable-type: a daughter of b
    mo.SetPoint(b.Point()->left);
    mo.DetermineCouplings(1);
  } 
  else  if (b.Point() == a.Point()->left) {
    // combinable-type: b daughter of a
    mo.SetPoint(a.Point()->right);
    mo.DetermineCouplings(1);
  } 
  else  if (b.Point() == a.Point()->right) {
    // combinable-type: b daughter of a
    mo.SetPoint(a.Point()->left);
    mo.DetermineCouplings(1);
  } 
  else THROW(fatal_error,"   Cannot combine legs.");

  // fix charge incase initial state has wrong
  int icharge;
  if (i<2)  icharge = a.Anti()*a.Point()->fl.IntCharge() - 
	              b.Anti()*b.Point()->fl.IntCharge();
  else      icharge = a.Point()->fl.IntCharge() + b.Point()->fl.IntCharge();
  if (icharge!=mo.Point()->fl.IntCharge()) mo.SetAnti(-1);    
  return mo;
}
  
Leg * Combine_Table_Base::CombineLegs
(Leg *legs,const int i,const int j,const int nlegs,const double pt2ij) 
{
  Leg * alegs = new Leg[nlegs];
  // assume i < j 
  for (int l=0; l<j; ++l) {
    if (l==i) {
      alegs[i] = CombinedLeg(legs,i,j);
      if (alegs[i].OrderQCD()>0) m_kt2QCD=Min(m_kt2QCD,pt2ij);
      if (alegs[i].OrderQED()>0) m_kt2QED=Min(m_kt2QED,pt2ij);
    }
    else      alegs[l] = legs[l];
  }
  for (int l=j+1; l<=nlegs; ++l) alegs[l-1] = legs[l];
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
  if (p_up) {
    double lastkt2(p_up->MinKt2());
    return Min(lastkt2,m_kt2min);
  }
  return m_kt2min;
}

double Combine_Table_Base::MinKt2QCD() const
{
  // return lowest cluster scale, which is not necessarily 
  // scale of first clustering
  if (p_up) {
    double lastkt2(p_up->MinKt2QCD());
    return Min(lastkt2,m_kt2QCD);
  }
  return m_kt2QCD;
}

double Combine_Table_Base::MinKt2QED() const
{
  // return lowest cluster scale, which is not necessarily 
  // scale of first clustering
  if (p_up) {
    double lastkt2(p_up->MinKt2QED());
    return Min(lastkt2,m_kt2QED);
  }
  return m_kt2QED;
}

double Combine_Table_Base::Kt2() const
{
  return m_kt2min;
}

double Combine_Table_Base::Kt2QCD() const
{
  return m_kt2QCD;
}

double Combine_Table_Base::Kt2QED() const
{
  return m_kt2QED;
}

double Combine_Table_Base::Sprime() const
{
  if (!p_moms) {
    return 0;
  }
  return (p_moms[0]+p_moms[1]).Abs2();
}

void Combine_Table_Base::IdentifyHardProcess()
{
  msg_Debugging()<<METHOD<<"():\n";
  msg_Indent();
  m_nstrong=0;
  if (p_hard==NULL) {
    p_hard = new Leg*[m_nampl];
    for (int i(0);i<m_nampl;++i) p_hard[i] = new Leg[2];
  }
  for (int i(0);i<m_nampl;++i) {
    if (Combinable(p_legs[i][0],p_legs[i][1]) &&
	Combinable(p_legs[i][2],p_legs[i][3])) {
      msg_Debugging()<<"s-channel\n";
      p_hard[i][0]=CombinedLeg(p_legs[i],0,1);
      p_hard[i][1]=CombinedLeg(p_legs[i],2,3);
    }
    else if (Combinable(p_legs[i][0],p_legs[i][2]) &&
	     Combinable(p_legs[i][1],p_legs[i][3])) {
      msg_Debugging()<<"t-channel\n";
      p_hard[i][0]=CombinedLeg(p_legs[i],0,2);
      p_hard[i][1]=CombinedLeg(p_legs[i],1,3);
    }
    else if (Combinable(p_legs[i][0],p_legs[i][3]) &&
	     Combinable(p_legs[i][1],p_legs[i][2])) {
      msg_Debugging()<<"u-channel\n";
      p_hard[i][0]=CombinedLeg(p_legs[i],0,3);
      p_hard[i][1]=CombinedLeg(p_legs[i],1,2);
    }
    else THROW(fatal_error,"No match for hard process.");
    m_nstrong=Max(m_nstrong,p_hard[i][0].OrderQCD()+p_hard[i][1].OrderQCD());
  }
}

int Combine_Table_Base::AddCouplings(int &nqed,int &nqcd) const
{
  if (p_up) {
    int nstrong = p_up->m_cdata_winner->second.m_strong;
    nqed+=1-nstrong;
    nqcd+=nstrong;
    return p_up->AddCouplings(nqed,nqcd);
  }
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

