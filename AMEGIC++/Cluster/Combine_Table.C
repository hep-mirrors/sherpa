#include "AMEGIC++/Cluster/Combine_Table.H"

#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/My_Limits.H"
#include "PDF/Main/Cluster_Definitions_Base.H"
#include "ATOOLS/Math/Random.H"
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
  p_qmin(NULL) {}

Leg::Leg(const Leg &leg): 
  p_point(leg.p_point), m_anti(leg.m_anti), 
  m_nqcd(0), m_nqed(0), m_pqcd(0), m_pqed(0), m_ext(leg.m_ext), 
  m_qcdjets(leg.m_qcdjets), m_id(leg.m_id),
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
  s<<" "<<cd.m_pt2ij<<" "<<cd.m_dec<<" ";
  std::string graphs;
  for (size_t k=0;k<cd.m_graphs.size();++k) graphs+=","+ToString(cd.m_graphs[k]);
  s<<graphs.substr(1);
  if (cd.p_down) s<<" -> "<<cd.p_down->m_no;
  return s;
}

Combine_Data::Combine_Data():
  m_pt2ij(0.0), m_strong(0), m_calc(0), p_down(NULL) {}

Combine_Data::Combine_Data(const double pt2ij,const int ngraph):
  m_pt2ij(pt2ij), m_strong(0), m_calc(0), p_down(NULL) 
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
    if (ct.m_decids.size()) {
      std::string ds;
      for (DecayInfo_Vector::const_iterator cit(ct.m_decids.begin());
	   cit!=ct.m_decids.end();++cit) {
        ds+=ToString(**cit)+" ";
      }
      s<<"  decs = { "<<ds<<"}\n";
    }
    s<<" id"<<std::setw(12)<<"content"<<std::setw(8)
     <<"flav"<<std::setw(5)<<"cut  qcd qed"<<std::setw(12)
     <<" mom"<<std::endl;
    for (int l=0; l<ct.m_nlegs; ++l) {
      s<<std::setw(3)<<l<<std::setw(12)<<ToString(ID(ct.GetLeg(l).ID()))
       <<std::setw(8)<<ct.p_legs[0][l].Flav()<<std::setw(4)
       <<ct.p_legs[0][l].Point()->t<<" "<<ct.GetLeg(l).OrderQCD()
       <<"/"<<ct.GetLeg(l).NQCD()<<" "<<ct.GetLeg(l).OrderQED()
       <<"/"<<ct.GetLeg(l).NQED()<<" "<<ct.p_moms[l]<<std::endl;
    }
    s<<" ---------------"<<std::endl;
    const CD_List & cl=ct.m_combinations;
    if (cl.size()>0) {
      for (CD_List::const_iterator cit=cl.begin(); cit!=cl.end(); ++cit) {
 	s<<cit->first<<std::setw(8)<<cit->second
	 <<(cit==ct.m_cdata_winner?" <-":"")<<std::endl; 
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
	  s<<std::setw(3)<<k<<"("<<l<<")"<<std::setw(8)
	   <<ct.p_hard[k][l].Flav()<<std::setw(4)
	   <<ct.p_hard[k][l].Point()->t<<" "<<ct.p_hard[k][l].OrderQCD()
	   <<"/"<<ct.p_hard[k][l].NQCD()<<" "<<ct.p_hard[k][l].OrderQED()
	   <<"/"<<ct.p_hard[k][l].NQED()<<std::setw(12)<<std::endl;
	}
      }
    }
    s<<" k_{T,min}\n";
    for (size_t i(0);i<ct.m_kt2ord.size();++i)
      s<<ID(ct.m_kt2ord[i].first)<<" -> "<<sqrt(ct.m_kt2ord[i].second)<<"\n";
  } 
  else
    s<<"***empty Combine_Table***"<<std::endl;
  return s;
}

Combine_Table::Combine_Table(AMEGIC::Process_Base *const proc,
			     ATOOLS::Mass_Selector *const ms,
			     PDF::Cluster_Definitions_Base *clus,
			     Vec4D *moms, Combine_Table *up,
			     ATOOLS::DecayInfo_Vector *const decids):
  p_ms(ms), m_nstrong(0), m_nlegs(0), m_nampl(0),
  m_graph_winner(0), 
  p_up(up), p_legs(0), p_clus(clus), p_moms(moms),
  p_hard(NULL), p_hardc(NULL), p_channel(NULL), p_scale(NULL), m_rscale(-1.0),
  p_decids(decids)
{
  p_proc=proc;
  m_no=++s_all;
  m_kt2ord=KT2Info_Vector(1,KT2_Info((1<<proc->NIn()+proc->NOut())-1,0.0));
  for (size_t i(0);i<m_decids.size();++i)
    m_kt2ord.push_back(std::make_pair(m_decids[i]->m_id,0.0));
}

Combine_Table::~Combine_Table()
{
  if (p_scale) delete [] p_scale;
  if (p_channel) delete [] p_channel;
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

Flavour Combine_Table::MatchFlavour(const Leg &a,const Leg &b,const Leg &c,int mode) const
{
  if (p_proc->Partner()==p_proc) return a.Point()->fl;
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

Leg * Combine_Table::CombineLegs
(Leg *legs,const int i,const int j,const int nlegs,const int kin) 
{
  Leg * alegs = new Leg[nlegs];
  // assume i < j 
  for (int l=0; l<j; ++l) {
    if (l==i) {
      alegs[i] = CombinedLeg(legs,i,j);
      alegs[i].SetKin(kin);
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
    ampl->CreateLeg(i<2?-moms[i]:moms[i],
		    i<2?p_legs[0][i].Flav().Bar():p_legs[0][i].Flav(),
		    ColorID(),p_legs[0][i].ID());
  Vec4D_Vector after=p_clus->Combine
    (*ampl,m_cdata_winner->first.m_i,m_cdata_winner->first.m_j,
     m_cdata_winner->first.m_k,m_cdata_winner->first.m_i<2?
     m_cdata_winner->second.m_mo.Bar():m_cdata_winner->second.m_mo,p_ms,
     m_cdata_winner->second.m_pt2ij.m_kin,m_cdata_winner->second.m_pt2ij.m_mode);
  ampl->Delete();
  if (after.empty()) {
    msg_Debugging()<<"combine moms failed\n";
    return false;
  }
  bool swap(p_legs[0][0].ID()&2);
  if ((-after[swap][0]>rpa->gen.PBeam(0)[0] &&
       !IsEqual(-after[swap][0],rpa->gen.PBeam(0)[0],1.0e-6)) ||
      (-after[1-swap][0]>rpa->gen.PBeam(1)[0] &&
       !IsEqual(-after[1-swap][0]>rpa->gen.PBeam(1)[0],1.0e-6))) {
    msg_Debugging()<<"kinematics failed\n";
    return false;
  }
  for (size_t l=0; l<after.size(); ++l) p_moms[l] = l<2?-after[l]:after[l];
  return true;
}

bool Combine_Table::CombineMoms(Vec4D *moms,const int _i,const int _j,
				     const int maxl,Vec4D *&omoms) 
{
  Cluster_Amplitude *ampl(Cluster_Amplitude::New());
  for (int i=0;i<=maxl;++i)
    ampl->CreateLeg(i<2?-moms[i]:moms[i],
		    i<2?p_legs[0][i].Flav().Bar():p_legs[0][i].Flav(),
		   ColorID(),p_legs[0][i].ID());
  Vec4D_Vector after=p_clus->Combine
    (*ampl,m_cdata_winner->first.m_i,m_cdata_winner->first.m_j,
     m_cdata_winner->first.m_k,m_cdata_winner->first.m_i<2?
     m_cdata_winner->second.m_mo.Bar():m_cdata_winner->second.m_mo,p_ms,
     m_cdata_winner->second.m_pt2ij.m_kin,m_cdata_winner->second.m_pt2ij.m_mode);
  ampl->Delete();
  if (after.empty()) {
    msg_Debugging()<<"combine moms failed\n";
    return false;
  }
  bool swap(p_legs[0][0].ID()&2);
  if ((-after[swap][0]>rpa->gen.PBeam(0)[0] &&
       !IsEqual(-after[swap][0],rpa->gen.PBeam(0)[0],1.0e-6)) ||
      (-after[1-swap][0]>rpa->gen.PBeam(1)[0] &&
       !IsEqual(-after[1-swap][0]>rpa->gen.PBeam(1)[0],1.0e-6))) {
    msg_Debugging()<<"kinematics failed\n";
    return false;
  }
  omoms = new Vec4D[maxl];
  for (size_t l=0; l<after.size(); ++l) omoms[l] = l<2?-after[l]:after[l];
  return true;
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

bool Combine_Table::Combinable(const Leg &a,const Leg &b,const int i,const int j) const
{
  Leg lmo;
  Leg tmp1 = a;
  Leg tmp2 = b;

  if ((i<2 || j<2) && (tmp1.Flav().IsSusy() || tmp2.Flav().IsSusy())) {
    return 0;
  }

  if (a.Point()->prev!=NULL && a.Point()->prev==b.Point()->prev) {
    lmo.SetPoint((AMEGIC::Point*)a.Point()->prev);
    return true;
  }
  if (a.Point()->prev==b.Point()) {
    lmo.SetPoint((AMEGIC::Point*)b.Point());
    return true;
  }
  if (b.Point()->prev==a.Point()) {
    lmo.SetPoint((AMEGIC::Point*)a.Point());
    return true;
  }
  return false;
}

double Combine_Table::GetWinner(int &i,int &j,int &k,double &mu2)
{ 
  i=m_cdata_winner->first.m_i; 
  j=m_cdata_winner->first.m_j;
  k=m_cdata_winner->first.m_k;
  mu2=m_cdata_winner->second.m_pt2ij.m_mu2;
  return m_cdata_winner->second.m_pt2ij.m_kt2;
}

void Combine_Table::AddPossibility(const int i,const int j,const int k,
				   const int ngraph) 
{
  Leg cl(CombinedLeg(p_legs[ngraph],i,j));
//   if (cl.Flav().Strong()) {
//     if (!p_legs[ngraph][k].Flav().Strong()) return;
//   }
//   else if (cl.Flav().IntCharge()==0 && cl.Flav()!=Flavour(kf_h0)) {
//     if (p_legs[ngraph][k].Flav().IntCharge()==0) return;
//   }
  CD_List::iterator cit=m_combinations.find(Combine_Key(i,j,k,cl.Flav()));
  if (cit!=m_combinations.end()) {
    cit->second.m_graphs.push_back(ngraph);
    cit->second.m_strong=Max(cit->second.m_strong,cl.OrderQCD());
  }
  else {
    Combine_Data cd(0.,ngraph);
    cd.m_strong=cl.OrderQCD();
    cd.m_mo=cl.Flav();
    cd.m_dec=cl.Point()->t;
    m_combinations[Combine_Key(i,j,k,cl.Flav())]=cd;
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
  if (m_nlegs>2+p_proc->Info().m_fi.NMinExternal()) {
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
	  if (Combinable(p_legs[k][i],p_legs[k][j],i,j)) {
	    int sci(p_legs[k][i].Flav().StrongCharge());
	    int scj(p_legs[k][j].Flav().StrongCharge());
	    Leg lmo(CombinedLeg(p_legs[k],i,j));
	    for (int l=0;l<m_nlegs;++l)
	      if (l!=i && l!=j) {
		int sc(p_legs[k][l].Flav().StrongCharge());
		if (((sci==8 || scj==8 || sc==8) && 
		     (sci!=0 && scj!=0 && sc!=0)) ||
		    (sci!=8 && scj!=8 && sc!=8) ||
		    Combinable(lmo,p_legs[k][l],i,l)) {
		  AddPossibility(i,j,l,k);
		  if (sci==8 && scj==8) AddPossibility(j,i,l,k);
		}
	      }
	  }
	}
      }
    }
  }
}

CD_List::iterator Combine_Table::CalcPropagator(CD_List::iterator &cit,int mode)
{
    if (cit->second.m_calc) return cit;
    Cluster_Amplitude *ampl(Cluster_Amplitude::New());
    for (int i=0;i<m_nlegs;++i)
      ampl->CreateLeg(i<2?-p_moms[i]:p_moms[i],
		      i<2?p_legs[0][i].Flav().Bar():p_legs[0][i].Flav(),
		     ColorID(),p_legs[0][i].ID());
    cit->second.m_calc=1;
    cit->second.m_pt2ij=p_clus->KPerp2
      (*ampl,cit->first.m_i,cit->first.m_j,cit->first.m_k,
       cit->first.m_i<2?cit->second.m_mo.Bar():cit->second.m_mo,p_ms,
       (mode&1024)||(mode&4096)?1:-1,(cit->second.m_dec>10||!cit->second.m_mo.Strong()?1:0)|
       (p_proc->Parent()->Info().m_fi.m_nloqcdtype!=PHASIC::nlo_type::lo?2:0));
    msg_Debugging()<<"Calculate m_perp("<<cit->first.m_i<<"["
		   <<p_legs[0][cit->first.m_i].Flav()<<"],"
		   <<cit->first.m_j<<"["<<p_legs[0][cit->first.m_j].Flav()<<"],"
		   <<cit->first.m_k<<"["<<p_legs[0][cit->first.m_k].Flav()
		   <<"],"<<cit->second.m_mo<<") -> "<<cit->second.m_pt2ij<<std::endl;
    ampl->Delete();
    return cit;
}

KT2Info_Vector Combine_Table::UpdateKT2(const CD_List::iterator &cdit) const
{
  KT2Info_Vector nkt2ord(m_kt2ord);
  size_t sid(p_legs[0][cdit->first.m_i].ID()|
	     p_legs[0][cdit->first.m_j].ID()), lmin(100), li(0);
  for (size_t i(0);i<nkt2ord.size();++i) {
    if ((nkt2ord[i].first&sid)==sid &&
	IdCount(nkt2ord[i].first)<lmin) {
      lmin=IdCount(nkt2ord[i].first);
      li=i;
    }
  }
  if ((cdit->second.m_dec<10 &&
       cdit->first.m_flav.Strong()) ||
      m_nlegs==2+p_proc->Info().m_fi.NMinExternal()) {
    nkt2ord[li].second=cdit->second.m_pt2ij.m_kt2;
    msg_Debugging()<<"set last k_T = "<<sqrt(nkt2ord[li].second)
		   <<" for "<<ID(nkt2ord[li].first)
		   <<" from "<<ID(sid)<<"\n";
  }
  return nkt2ord;
}

Combine_Table *Combine_Table::
CalcJet(int nl,ATOOLS::Vec4D * moms,const size_t mode,const double &kt2) 
{
  DEBUG_FUNC(mode<<" "<<nl<<" "<<sqrt(kt2));
  if (nl==3) return this;
  m_rejected.clear();
  bool invonly(true), valid(mode&512);
  while (true) {
    m_nl=nl;
    if (moms) for (size_t l=0;l<m_nl;++l) p_moms[l]=moms[l];
    if (!SelectWinner(mode)) {
      if (nl==4 && p_proc->Info().m_fi.NMinExternal()>1 &&
	  (IdentifyHardProcess() || p_up==NULL)) {
	DecayInfo_Vector decids(m_decids);
	for (int j(0);j<m_nampl;++j) {
	  size_t pid(p_hard[j][1].ID());
	  if (pid==0) continue;
	  for (size_t i(0);i<p_decids->size();++i)
	    if ((*p_decids)[i]->m_id==pid) {
	      decids.push_back((*p_decids)[i]);
	      break;
	    }
	}
	  m_graph_winner=-1;
	  double kt2min(std::numeric_limits<double>::max());
	  for (int i(0);i<m_nampl;++i)
	    if (p_scale[i]<kt2min) {
	      kt2min=p_scale[i];
	      m_graph_winner=i;
	    }
	  Cluster_Amplitude *ampl(Cluster_Amplitude::New());
	  ampl->SetProc(p_proc);
	  ampl->SetMS(p_ms);
	  for (int i=0;i<m_nlegs;++i)
	    ampl->CreateLeg(i<2?-p_moms[i]:p_moms[i],
			    i<2?p_legs[0][i].Flav().Bar():p_legs[0][i].Flav(),
			    ColorID(),p_legs[0][i].ID());
	  PDF::CParam scale((p_proc->IsMapped()?p_proc->MapProc():p_proc)
			    ->ScaleSetter()->CoreScale(ampl));
	  ampl->Delete();
	  bool ord(true);
	  for (size_t i(0);i<m_kt2ord.size();++i) {
	    if (scale.m_kt2<m_kt2ord[i].second &&
		p_hard[m_graph_winner][0].Flav().Strong()) {
	      msg_Debugging()<<"unordered configuration: "
			     <<sqrt(scale.m_kt2)<<" vs. "
			     <<sqrt(m_kt2ord[i].second)<<" "
			     <<ID(m_kt2ord[i].first)<<" "
			     <<p_hard[m_graph_winner][0].Flav()
			     <<" (mode="<<mode<<")\n";
	      if ((p_up && p_up->p_up) || !(mode&4096)) ord=false;
	      break;
	    }
	  }
	  if (!ord) {
	    if (!(mode&16)) {
	      delete this;
	      return NULL;
	    }
	  }
	  return this;
      }
      break;
    }
    m_rejected[m_cdata_winner->first]=m_cdata_winner->second;
    if (m_cdata_winner->second.m_pt2ij.m_op2<0.0) continue;
    invonly=false;
    bool ord(true);
    KT2Info_Vector nkt2ord(UpdateKT2(m_cdata_winner));
    for (size_t i(0);i<nkt2ord.size();++i) {
      if (nkt2ord[i].second<m_kt2ord[i].second) {
	msg_Debugging()<<"unordered configuration: "
		       <<sqrt(nkt2ord[i].second)<<" vs. "
		       <<sqrt(m_kt2ord[i].second)<<" "
		       <<ID(m_kt2ord[i].first)<<" (mode="<<mode<<")\n";
	if ((p_up && p_up->p_up) || !(mode&4096)) ord=false;
	break;
      }
    }
    if (!ord) {
      if (!(mode&16)) continue;
    }
    if (nl<2+p_proc->Info().m_fi.NMinExternal()) THROW(fatal_error,"nlegs < min. Abort.");
    double scale(-1.0);
    Combine_Table *tab(CreateNext());
    if (tab!=NULL) {
      scale=(tab->Momentum(0)+tab->Momentum(1)).Abs2();
      Combine_Table *next(NextTable(tab,mode,m_cdata_winner->
				    second.m_pt2ij.m_kt2));
      if (next!=NULL) return next;
    }
    if (p_up==NULL && (mode&4096)) {
      if (scale<0.0) THROW(fatal_error,"No rescue scale");
      msg_Debugging()<<"rescue scale is "<<sqrt(scale)<<"\n";
      m_rscale=scale;
      return this;
    }
    msg_Debugging()<<METHOD<<"(): Table "<<m_no<<": reject winner "
		   <<m_cdata_winner->first<<"\n";
  }
  if (valid && invonly && kt2>0.0 &&
      ((p_up && p_up->p_up) || !(mode&4096))) {
    msg_Debugging()<<"no valid configuration, retry previous\n";
    delete this;
    return NULL;    
  }
  msg_Debugging()<<"trying unordered configuration\n";
  bool zero(false);
  CD_List norejected;
  while (true) {
    double nmin(std::numeric_limits<double>::max()), pmin(nmin);
    CD_List::iterator nwin(m_rejected.end()), pwin(nwin);
    for (CD_List::iterator vit(m_rejected.begin());
	 vit!=m_rejected.end();++vit) {
      if (norejected.find(vit->first)!=norejected.end()) continue;
      if (vit->second.m_pt2ij.m_op2<0.0 && (mode&512)) continue;
      double ktdiff(vit->second.m_pt2ij.m_kt2-kt2);
      if (ktdiff<0.0) {
	if (-ktdiff<nmin) {
	  nmin=-ktdiff;
	  nwin=vit;
	}
      }
      else {
	if (ktdiff<pmin) {
	  pmin=ktdiff;
	  pwin=vit;
	}
      }
    }
    if (nwin==m_rejected.end()) nwin=pwin;
    if (nwin==m_rejected.end()) {
      if (!zero) {
	zero=true;
	norejected.clear();
	continue;
      }
      else {
	// if (valid && p_up==NULL) return CalcJet(nl,moms,mode,kt2);
	delete this;
	return NULL;
      }
    }
    norejected[nwin->first]=nwin->second;
    if (valid && nwin->second.m_pt2ij.m_op2<0.0) continue;
    m_rejected.erase(nwin);
    m_nl=nl;
    if (moms) for (size_t l=0;l<m_nl;++l) p_moms[l]=moms[l];
    if (!SelectWinner(mode)) continue;
    m_rejected[m_cdata_winner->first]=m_cdata_winner->second;
    if (nl<4) THROW(fatal_error,"nlegs < min. Abort.");
    Combine_Table *tab(CreateNext());
    if (tab!=NULL) {
      Combine_Table *next(NextTable(tab,mode,zero?0.0:m_cdata_winner->
				    second.m_pt2ij.m_kt2));
      if (next!=NULL) return next;
    }
    msg_Debugging()<<METHOD<<"(): Table "<<m_no<<": reject winner "
		   <<m_cdata_winner->first<<"\n";
  }
  return NULL;
}

bool Combine_Table::SelectWinner(const size_t &mode)
{
  CD_List & cl(m_combinations);
  if (cl.size()==0) return false;
  // calculate pt2ij and determine "best" combination
  m_cdata_winner = cl.end();
  CD_List::iterator rdata_winner(cl.end()), ldata_winner(cl.end());
  double kt2(std::numeric_limits<double>::max()), lkt2(kt2);
  double rkt2(std::numeric_limits<double>::max()), sum(0.0);
  for (CD_List::iterator cit(cl.begin()); cit!=cl.end(); ++cit) {
    CD_List::iterator tit(CalcPropagator(cit,mode));
    double pt2ij(cit->second.m_pt2ij.m_op2);
    if (cit->second.m_graphs.size()==0) continue;
    if (m_rejected.find(cit->first)==m_rejected.end()) {
      if (pt2ij>0.0) {
	// if (mode&1) {
	  if (pt2ij<kt2) {
	    m_cdata_winner=cit;
	    kt2=pt2ij;
	  }
	// }
	else {
	sum+=1.0/pt2ij;
	}
	if (cit->second.m_pt2ij.m_kt2<lkt2) {
	  ldata_winner=cit;
	  lkt2=cit->second.m_pt2ij.m_kt2;
	}
      }
      else if (cit->second.m_pt2ij.m_kt2>0.0 &&
	       cit->second.m_pt2ij.m_kt2<rkt2) {
	rdata_winner=cit;
	rkt2=cit->second.m_pt2ij.m_kt2;
      }
    }
  }
  if (!(mode&1) && !(mode&4096)) {
    double disc(sum*ran->Get()), psum(0.0);
    for (CD_List::iterator cit(cl.begin()); cit!=cl.end(); ++cit) {
      double pt2ij(cit->second.m_pt2ij.m_op2);
      if (cit->second.m_graphs.size()==0) continue;
      if (m_rejected.find(cit->first)==m_rejected.end() &&
	  pt2ij>0.0 && (psum+=1.0/pt2ij)>=disc) {
	m_cdata_winner=cit;
	break;
      }
    }
    if (sum>0.0 && m_cdata_winner==cl.end())
      THROW(fatal_error,"Internal error");
  }
  if (m_cdata_winner==cl.end() && !(mode&512)) m_cdata_winner=rdata_winner;
  // if (mode&1) m_cdata_winner=ldata_winner;
  msg_Debugging()<<*this<<"\n";
  return m_cdata_winner!=cl.end();
}

Combine_Table *Combine_Table::CreateNext()
{
  --m_nl;
  int i(m_cdata_winner->first.m_i), j(m_cdata_winner->first.m_j);
  if (i>j) std::swap<int>(i,j);
  if (!m_cdata_winner->second.p_down) {
    Vec4D * amoms;
    // generate new momenta
    if (!CombineMoms(p_moms,i,j,m_nl,amoms)) return NULL;
    Leg ** alegs = new Leg*[m_cdata_winner->second.m_graphs.size()];
    for (size_t k=0;k<m_cdata_winner->second.m_graphs.size();++k) {
      alegs[k] = CombineLegs
	(p_legs[m_cdata_winner->second.m_graphs[k]],i,j,m_nl,
	 m_cdata_winner->second.m_pt2ij.m_kin);
    }
    m_cdata_winner->second.p_down = 
      new Combine_Table(p_proc,p_ms,p_clus,amoms,this,p_decids);
    {
      Combine_Table *tab((Combine_Table*)m_cdata_winner->second.p_down);
      tab->m_kt2ord=UpdateKT2(m_cdata_winner);
      tab->m_decids=m_decids;
      size_t pid(p_legs[0][i].ID()+p_legs[0][j].ID());
      for (size_t i(0);i<p_decids->size();++i)
        if ((*p_decids)[i]->m_id==pid) {
	  tab->m_decids.push_back((*p_decids)[i]);
	  break;
	}
    }
    // initialise Combine_Table
    m_cdata_winner->second.p_down->FillTable(alegs,m_nl,m_cdata_winner->second.m_graphs.size());
  } 
  else {
    if (!((Combine_Table*)m_cdata_winner->second.p_down)->
	CombineMoms(p_moms,i,j,m_nl)) return NULL;
  }

  Combine_Table *tab((Combine_Table*)m_cdata_winner->second.p_down);
  return tab;
}

Combine_Table *Combine_Table::NextTable(Combine_Table *tab,
					const int mode,const double &kt2)
{
  Combine_Table* ct = tab->CalcJet(m_nl,NULL,mode,kt2);
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
  if (p_hard==NULL) {
    p_hard = new Leg*[m_nampl];
    for (int i(0);i<m_nampl;++i) p_hard[i] = new Leg[2];
    p_hardc = new int*[m_nampl];
    for (int i(0);i<m_nampl;++i) p_hardc[i] = new int[4];
  }
  if (p_scale==NULL) p_scale = new double[m_nampl];
  if (p_channel==NULL) p_channel = new int[m_nampl];
  for (int i(0);i<m_nampl;++i) {
    if (Combinable(p_legs[i][0],p_legs[i][1],0+2,1+2) &&
	Combinable(p_legs[i][2],p_legs[i][3],2+2,3+2)) {
      double pt2ij1((p_moms[0]+p_moms[1]).Abs2());
      double pt2ij2((p_moms[2]+p_moms[3]).Abs2());
      msg_Debugging()<<"s-channel pt = "<<sqrt(pt2ij1)
		     <<" / "<<sqrt(pt2ij2)<<", m = "
		     <<sqrt(dabs((p_moms[0]+p_moms[1]).Abs2()))<<", "
		     <<p_legs[i][0].Flav()<<" "<<p_legs[i][1].Flav()
		     <<" -> "<<p_legs[i][2].Flav()<<" "
		     <<p_legs[i][3].Flav()<<"\n";
      p_hard[i][0]=CombinedLeg(p_legs[i],0,1);
      p_hard[i][1]=CombinedLeg(p_legs[i],2,3);
      size_t idi(GetLeg(2).ID()), idj(GetLeg(3).ID());
      p_hard[i][1].SetID(idi|idj);
      p_hardc[i][0]=0;
      p_hardc[i][1]=0;
      p_hardc[i][2]=1;
      p_hardc[i][3]=1;
      p_scale[i]=sqrt(dabs(pt2ij1*pt2ij2));
      p_channel[i]=1;
    }
    else if (Combinable(p_legs[i][0],p_legs[i][2],0+2,2+2) &&
	     Combinable(p_legs[i][1],p_legs[i][3],1+2,3+2)) {
      double pt2ij1((p_moms[0]-p_moms[2]).Abs2());
      double pt2ij2((p_moms[1]-p_moms[3]).Abs2());
      msg_Debugging()<<"t-channel pt = "<<sqrt(pt2ij1)
		     <<" / "<<sqrt(pt2ij2)<<", m = "
		     <<sqrt(dabs((p_moms[0]+p_moms[2]).Abs2()))<<", "
		     <<p_legs[i][0].Flav()<<" "<<p_legs[i][1].Flav()
		     <<" -> "<<p_legs[i][2].Flav()<<" "
		     <<p_legs[i][3].Flav()<<"\n";
      p_hard[i][0]=CombinedLeg(p_legs[i],0,2);
      p_hard[i][1]=CombinedLeg(p_legs[i],1,3);
      p_hardc[i][0]=0;
      p_hardc[i][1]=1;
      p_hardc[i][2]=0;
      p_hardc[i][3]=1;
      p_scale[i]=sqrt(dabs(pt2ij1*pt2ij2));
      p_channel[i]=2;
    }
    else if (Combinable(p_legs[i][0],p_legs[i][3],0+2,3+2) &&
	     Combinable(p_legs[i][1],p_legs[i][2],1+2,2+2)) {
      double pt2ij1((p_moms[0]-p_moms[3]).Abs2());
      double pt2ij2((p_moms[1]-p_moms[2]).Abs2());
      msg_Debugging()<<"u-channel pt = "<<sqrt(pt2ij1)
		     <<" / "<<sqrt(pt2ij2)<<", m = "
		     <<sqrt(dabs((p_moms[0]+p_moms[3]).Abs2()))<<", "
		     <<p_legs[i][0].Flav()<<" "<<p_legs[i][1].Flav()
		     <<" -> "<<p_legs[i][2].Flav()<<" "
		     <<p_legs[i][3].Flav()<<"\n";
      p_hard[i][0]=CombinedLeg(p_legs[i],0,3);
      p_hard[i][1]=CombinedLeg(p_legs[i],1,2);
      p_hardc[i][0]=0;
      p_hardc[i][1]=1;
      p_hardc[i][2]=1;
      p_hardc[i][3]=0;
      p_scale[i]=sqrt(dabs(pt2ij1*pt2ij2));
      p_channel[i]=3;
    }
    else THROW(fatal_error,"No match for hard process.");
    m_nstrong=Max(m_nstrong,p_hard[i][0].OrderQCD()+p_hard[i][1].OrderQCD());
  }
  return true;
}

int Combine_Table::IdentifyHardPropagator(double &mmin) const
{
  mmin=p_scale[m_graph_winner];
  return p_channel[m_graph_winner];
}

