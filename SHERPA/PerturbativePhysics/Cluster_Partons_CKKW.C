#include "Cluster_Partons_CKKW.H"

#include "Combine_Table_CKKW.H"
#include "NLL_Branching_Probabilities.H"
#include "Data_Reader.H"
#include "Shell_Tools.H"
#include <fstream> 

using namespace SHERPA;
using namespace EXTRAXS;
using namespace AMEGIC;
using namespace ATOOLS;
using namespace PHASIC;
using namespace MODEL;

class Three_Jet_Calc: public ATOOLS::Function_Base {
private:
  NLL_Sudakov *p_sud;
  NLL_Branching_Probability_Base *p_bp;
  double m_q, m_q0;
public:
  Three_Jet_Calc(NLL_Sudakov *const sud,const Flavour &f,const double &asfac): 
    p_sud(sud),
    p_bp(new GammaQ_QG_Lambda
	 ((bpm::code)p_sud->Mode(),0.0,MODEL::as,f.Mass(),asfac)) {}
  ~Three_Jet_Calc() { delete p_bp; }
  virtual double operator()(double q);
  virtual double operator()();
  inline void SetQ(const double q)   { m_q=q;   }
  inline void SetQ0(const double q0) { m_q0=q0; }
};// end of class Three_Jet_Rate
double Three_Jet_Calc::operator()(double q) 
{
  return p_bp->Gamma(q,m_q)*p_sud->Delta(kf_gluon)(q,m_q0);
}
double Three_Jet_Calc::operator()()
{
  return m_defval; 
} 

Cluster_Partons_CKKW::
Cluster_Partons_CKKW(Matrix_Element_Handler * me,ATOOLS::Jet_Finder * jf,
		     PDF::ISR_Handler *isr,
		     const int maxjet,const int isron,const int fsron) :
  Cluster_Partons_Base(me,jf,maxjet,isr->On(),isron,fsron)
{
  p_pdf[0]=isr->PDF(0);
  p_pdf[1]=isr->PDF(1);
  std::string helps;
  Data_Reader reader(" ",";","!","=");
  if (reader.ReadFromFile(helps,"PRINT_SUDAKOV") && helps.length()>0) 
    GenerateTables(helps);
  if (!reader.ReadFromFile(m_smode,"CKKW_SHUFFLE_MODE")) m_smode=0;
  else msg_Info()<<METHOD<<"(): Set shuffle mode "<<m_smode<<".\n";
}

Cluster_Partons_CKKW::~Cluster_Partons_CKKW()
{
}

void Cluster_Partons_CKKW::GenerateTables(const std::string &path)
{
  msg_Info()<<METHOD<<"("<<path<<"): Generating sudakov tables {"<<std::endl;
    MakeDir(path,true);
  for (short unsigned int l(0);l<2;++l) {
    NLL_Sudakov *sud(l==0?p_fssud:p_issud);
    if (sud==NULL) break;
    msg_Indent();
    double ecms(rpa.gen.Ecms()), qmin(1.0), step(pow(ecms/qmin,1.0/25.0));
    for (size_t i(1);i<=7;++i) {
      if (i==7) i=21;
      Flavour f((kf_code)i);
      if (ecms>2.0*f.Mass()) {
	msg_Info()<<f<<(i<7?" quark":"")<<" sudakov ..."<<std::flush;
	std::ofstream dqout((path+"/delta_"+std::string(l==0?"":"is_")+
			     ToString(f)+"_"+ToString(ecms)+".dat").c_str());
	for (double Q(ecms);Q>qmin;Q/=step) 
	  dqout<<Q<<" "<<(sud->Delta(f))(Q,qmin)<<"\n";
	msg_Info()<<"done"<<std::endl;
      }
    }
    for (size_t i(1);i<=7;++i) {
      if (i==7) i=21;
      Flavour f((kf_code)i);
      if (ecms>2.0*f.Mass()) {
	msg_Info()<<"2-"<<f<<" rate ..."<<std::flush;
	std::ofstream r2out((path+"/r2_"+std::string(l==0?"":"is_")+
			     ToString(f)+"_"+ToString(ecms)+".dat").c_str());
	for (double Q(ecms);Q>qmin;Q/=step) 
	  r2out<<2.0*log10(Q/ecms)<<" "
	       <<sqr(sud->Delta(f)(ecms,Q))<<"\n";
	msg_Info()<<"done"<<std::endl;
      }
    }
    for (size_t i(6);i>0;--i) {
      Flavour f((kf_code)i);
      if (ecms>2.0*f.Mass()) {
	msg_Info()<<"2-"<<f<<"+1-gluon rate ..."<<std::flush;
	std::ofstream r3out((path+"/r3_"+std::string(l==0?"":"is_")+
			     ToString(f)+"_"+ToString(ecms)+".dat").c_str());
	Three_Jet_Calc *r3test(new Three_Jet_Calc(sud,f,m_fs_as_factor));
	Gauss_Integrator gauss(r3test);
	r3test->SetQ(ecms);
	for (double Q(ecms);Q>qmin;Q/=step) {
	  r3test->SetQ0(Q);
	  r3out<<2.0*log10(Q/ecms)<<" "
	       <<2.0*sqr(sud->Delta(f)(ecms,Q))*
	    gauss.Integrate(Q,ecms,1.0e-3)<<"\n";
	}
	delete r3test;
	msg_Info()<<"done"<<std::endl;
      }
    }
  }
  msg_Info()<<"}"<<std::endl;
}

Combine_Table_Base *Cluster_Partons_CKKW::
CreateTable(Jet_Finder *jf,ATOOLS::Vec4D *amoms,Combine_Table_Base *ct,
	    const int isrmode,const int isrshoweron)
{
  return new Combine_Table_CKKW(p_ajf,amoms,0,m_isrmode,m_isrshoweron);
}

EXTRAXS::XS_Base *Cluster_Partons_CKKW::GetXS(EXTRAXS::XS_Group * group, 
					      ATOOLS::Flavour * fl)
{
  XS_Base * xs(NULL);  
  const size_t nin(2), nout(2);
  if (group->XSSelector()->
      FindInGroup(group,xs,nin,nout,fl)==std::string::npos) {
    int nqed(0),nqcd(0);
    p_ct->AddCouplings(nqed,nqcd);
    xs = group->XSSelector()->GetXS
      (nin,nout,fl,false,nqed,nqcd,false);
    if (xs) group->Add(xs);
  }
  p_xs = xs;
  return xs;
}

int Cluster_Partons_CKKW::SetColours(EXTRAXS::XS_Base * xs, 
				     ATOOLS::Vec4D * p, ATOOLS::Flavour * fl)
{
  p_xs=xs;
  if (!p_xs) return Cluster_Partons_Base::SetColours(p,fl);
  
  m_hard_nqed = p_xs->OrderEWeak();
  m_hard_nqcd = p_xs->OrderStrong();
  bool test(p_xs->SetColours(p)), check(true);
  for (int i=0; i<4; ++i) {
    if (fl[i].IsQuark() || fl[i].IsSquark()) {
      if ( fl[i].IsAnti() && 
	   (p_xs->Colours()[i][0]!=0 || p_xs->Colours()[i][1]==0)) check=false;
      if (!fl[i].IsAnti() && 
	  (p_xs->Colours()[i][0]==0 || p_xs->Colours()[i][1]!=0)) check=false;
    }
    if ((fl[i].IsGluon() || fl[i].IsGluino()) && 
	(p_xs->Colours()[i][0]==0 || p_xs->Colours()[i][1]==0))   check=false;
    if (!check) {
      msg_Error()<<"Cluster_Partons_CKKW::SetColours(..): \n"
		 <<"Colour check failed for the following combination:"
		 <<std::endl;
      for (int i=0; i<4; ++i) 
	msg_Error()<<"   "<<i<<" : "<<fl[i]<<" ("
		   <<p_xs->Colours()[i][0]<<","
		   <<p_xs->Colours()[i][1]<<")"<<std::endl;
      msg_Error()<<"Abort."<<std::endl;
      abort();
    }
  }
  for (int i=0;i<4;i++) {
    m_colors[i][0] = p_xs->Colours()[i][0];
    m_colors[i][1] = p_xs->Colours()[i][1];
  }
  m_q2_iss[0]=m_q2_iss[1]=m_q2_fss=dabs(p_xs->Scale(stp::ren));
  return test;
}


void Cluster_Partons_CKKW::InitWeightCalculation()
{
  m_weight = 1.;
  m_asweight = 1.;
  m_scale  = 1.;
  m_nlegs   = p_ct->NLegs();
  m_njet    = m_nlegs-2;
  m_nstrong = 0;
  for (int l=0; l<m_nlegs; ++l) {
    if (p_ct->GetLeg(l).Flav().Strong()) ++m_nstrong;
  }
  Combine_Table_Base *ct_test(p_ct);
  while (ct_test->Up()) {
    ct_test=ct_test->Up();
    ++m_njet;
  }
  m_last_q.clear();
  m_last_i.clear();
  m_last_q.resize(m_nlegs,std::numeric_limits<double>::max());
  m_last_i.resize(m_nlegs,0);
  // scales for alpha_s and correction weight for hard interaction
  // me is reweighted w/ Min(m_is_as_factor,m_fs_as_factor), including deltar
  // and ren. scale factor see Integrable_Base::KFactor and CalculateScale
  PHASIC::Integrable_Base *proc(p_me->GetAmegic()->GetProcess());
  const Info_Key &kc(proc->KFKey());
  m_as_jet[1]=m_as_jet[0]=kc.Double(2);
  msg_Debugging()<<"scales {\n";
  msg_Debugging()<<"  is_as_fac = "<<m_is_as_factor
		 <<", fs_as_fac = "<<m_fs_as_factor
		 <<", me_as_fac = "<<m_me_as_factor
		 <<", m_kfac = "<<m_kfac<<"\n";
  msg_Debugging()<<"  amegic : "<<std::setw(12)<<sqrt(m_q2_amegic)<<"\n";
  msg_Debugging()<<"  cut is : "<<std::setw(12)<<sqrt(m_q2_isjet)
		 <<"  ->  as = "<<m_as_jet[1]<<"\n";
  msg_Debugging()<<"  cut fs : "<<std::setw(12)<<sqrt(m_q2_fsjet)
		 <<"  ->  as = "<<m_as_jet[0]<<"\n";
  msg_Debugging()<<"  is ps l: "<<std::setw(12)<<sqrt(m_q2_iss[0])<<"\n";
  msg_Debugging()<<"  is ps r: "<<std::setw(12)<<sqrt(m_q2_iss[1])<<"\n";
  msg_Debugging()<<"  fs ps  : "<<std::setw(12)<<sqrt(m_q2_fss)<<"\n";
  msg_Debugging()<<"}\n"; 
}

double Cluster_Partons_CKKW::
InternalWeight(const bool is,Leg &leg,const double iupper,
	       const double iactual,double qmin)
{
  double upper(iupper), actual(iactual);
  if (upper<actual) {
    msg_Debugging()<<METHOD<<"(): q = "<<actual<<", Q = "
		   <<upper<<". Swap."<<std::endl;
    std::swap<double>(upper,actual);
  }
  if (actual<qmin) {
    msg_Debugging()<<METHOD<<"(): q_{min} = "<<qmin<<", q = "
		   <<actual<<". Swap."<<std::endl;
    std::swap<double>(actual,qmin);
  }
  double DeltaNum(1.), DeltaDenom(1.), DeltaRatio(1.);
  if (is && (m_sud_mode&2)) {
    DeltaNum   = p_issud->Delta(leg.Flav())(upper,qmin);
    DeltaDenom = p_issud->Delta(leg.Flav())(actual,qmin);
  }
  else if (m_sud_mode&1) {
    DeltaNum   = p_fssud->Delta(leg.Flav())(upper,qmin);
    DeltaDenom = p_fssud->Delta(leg.Flav())(actual,qmin);
  }
  DeltaRatio = DeltaNum/DeltaDenom;
  msg_Debugging()<<(is?"is":"fs")<<" weight \\Delta_{"
		 <<leg.Flav()<<"}("<<upper<<","<<qmin
		 <<")/\\Delta_{"<<leg.Flav()<<"}("<<actual<<","<<qmin
		 <<") = "<<DeltaNum<<" / "<<DeltaDenom<<" = "
		 <<DeltaRatio<<"\n";
  return DeltaRatio;
}

double Cluster_Partons_CKKW::
ExternalWeight(const bool is,Leg &leg,
	       const double actual,double qmin)
{
  double DeltaNum(1.);
  if (qmin==0.0) qmin=sqrt(leg.Q2Cut());
  if (actual<qmin) {
    ++m_fails;
    DeltaNum=1.0;
  }
  else {
    if (is && (m_sud_mode&2)) DeltaNum = p_issud->Delta(leg.Flav())(actual,qmin);
    else if (m_sud_mode&1) DeltaNum = p_fssud->Delta(leg.Flav())(actual,qmin);
  }
  msg_Debugging()<<(is?"is":"fs")<<" weight \\Delta_{"
		 <<leg.Flav()<<"}("<<actual<<","<<qmin
		 <<") = "<<DeltaNum<<"\n";
  return DeltaNum;
}    

void Cluster_Partons_CKKW::StoreOldValues(const int i,const int j,
					  const int si,const double ptij)
{
  m_last_i[i] = si;
  m_last_q[i] = ptij;
  for (int l=m_nlegs-1;l>j;--l) {
    m_last_i[l] = m_last_i[l-1];
    m_last_q[l] = m_last_q[l-1];
  } 
  m_last_i[j] = si;
  m_last_q[j] = m_last_q[i];
}  

double Cluster_Partons_CKKW::CouplingWeight(const bool is,Leg &leg,
					    const double &kt)
{
  if (!(m_as_mode&1)) return 1.0;
  double asref(m_as_jet[is]);
  if (leg.Type()=="Triangle" ||
      leg.Type()=="Box" ||
      leg.Type()=="C4GS") {
    asref=(*p_runas)(sqr(Flavour(kf_h0).Mass()));
    if (m_kfac!=0.) asref*=1.+asref/(2.*M_PI)*m_kfac;
    msg_Debugging()<<"higgs vertex -> q_ref = "<<Flavour(kf_h0).Mass()
		   <<" -> asref = "<<asref<<"\n";
  }
  double as_ptij(0.);
  if (is) as_ptij=(*p_runas)(m_me_as_factor*sqr(kt)*m_is_as_factor);
  else as_ptij=(*p_runas)(m_me_as_factor*sqr(kt)*m_fs_as_factor);
  if (m_kfac!=0.) as_ptij*=1.+as_ptij/(2.*M_PI)*m_kfac;
  msg_Debugging()<<"as weight (\\alpha_s("<<m_me_as_factor<<"*sqr("<<kt<<")*"
		 <<(is?m_is_as_factor:m_fs_as_factor)
		 <<")/\\alpha_{s,ref})^O_{as} = ( "<<as_ptij<<" / "
		 <<asref<<" ) ^ "<<leg.OrderQCD()<<" = "
		 <<pow(as_ptij/asref,leg.OrderQCD())<<"\n";
  return pow(as_ptij/asref,leg.OrderQCD());
}

void Cluster_Partons_CKKW::WeightHardProcess()
{
  msg_Debugging()<<METHOD<<"(): {\n";
  msg_Indent();
  int wminqcd(-1), wminqed(-1);
  double kt2minqcd(std::numeric_limits<double>::max()), kt2minqed(kt2minqcd);
  for (int i(0);i<p_ct->NAmplitudes();++i) {
    for (int j(0);j<2;++j) {
      double kt2qcd(p_ct->GetHardLegs()[i][j].KT2QCD());
      double kt2qed(p_ct->GetHardLegs()[i][j].KT2QED());
      if (kt2qcd<kt2minqcd || 
	  (IsEqual(kt2qcd,kt2minqcd) && kt2qed<kt2minqed)) {
	wminqcd=i;
	kt2minqcd=kt2qcd;
      }
      if (kt2qed<kt2minqed || 
	  (IsEqual(kt2qed,kt2minqed) && kt2qcd<kt2minqcd)) {
	wminqed=i;
	kt2minqed=kt2qed;
      }
    }
  }
  msg_Debugging()<<"QCD: wmin_qcd = "<<wminqcd<<", ktmin_qcd = "
		 <<sqrt(kt2minqcd)<<", wmin_qed = "<<wminqed
		 <<", ktmin_qed = "<<sqrt(kt2minqed)<<"\n";
  if (wminqcd>=0) {
    // check for scale definition from Single_XS
    double mu2r[2]={p_xs?dabs(p_xs->Scale(stp::sis)):0.0,
		    p_xs?dabs(p_xs->Scale(stp::sfs)):0.0};
    if (mu2r[0]==std::numeric_limits<double>::max()) mu2r[0]=0.0; 
    if (mu2r[1]==std::numeric_limits<double>::max()) mu2r[1]=0.0; 
    double qu(sqrt(mu2r[0]!=0.0?mu2r[0]:
		   p_ct->GetHardLegs()[wminqcd][0].KT2QCD()));
    double ql(sqrt(mu2r[1]!=0.0?mu2r[1]:
		   p_ct->GetHardLegs()[wminqcd][1].KT2QCD()));
    if (qu<ql) std::swap<double>(qu,ql);
    // if intermediate particle strong, apply sudakov weight
    if (p_ct->GetHardLegs()[wminqcd][0].Flav().Strong())
      m_weight*=InternalWeight(0,p_ct->GetHardLegs()[wminqcd][0],
			       qu,ql,m_qmin[0]);
    // set possibly separate scales for the two vertices
    for (short unsigned int i(0);i<4;++i) {
      double rs(mu2r[i/2]!=0.0?mu2r[i/2]:p_ct->GetHardLegs()[wminqcd]
		[p_ct->GetHardInfo()[wminqcd][i]].KT2QCD());
      if (rs==std::numeric_limits<double>::max())
	rs=p_ct->GetHardLegs()[wminqcd]
	  [p_ct->GetHardInfo()[wminqcd][i]].KT2QED();
      m_last_q[i]=sqrt(rs);
    }
    for (short unsigned int i(0);i<2;++i) {
      Leg &leg(p_ct->GetHardLegs()[wminqcd][i]);
      if (leg.OrderQCD()>0) {
	double rs(mu2r[i]!=0.0?mu2r[i]:leg.KT2QCD());
	double asw(CouplingWeight(0,leg,sqrt(rs)));
	m_weight*=asw;
	m_asweight*=asw;
	m_scale*=sqrt(rs);
      }
    }
    if (m_q2_fss==std::numeric_limits<double>::max()) {
      m_q2_iss[0]=sqr(m_last_q[0]);
      m_q2_iss[1]=sqr(m_last_q[1]);
      m_q2_fss=m_last_q[2]*m_last_q[3];
    }
  }
  else if (wminqed>=0) {
    for (short unsigned int i(0);i<4;++i) 
      m_last_q[i]=sqrt(p_ct->GetHardLegs()[wminqed]
 		       [p_ct->GetHardInfo()[wminqed][i]].KT2QED());
    if (m_q2_fss==std::numeric_limits<double>::max()) {
      m_q2_iss[0]=sqr(m_last_q[0]);
      m_q2_iss[1]=sqr(m_last_q[1]);
      m_q2_fss=m_last_q[2]*m_last_q[3];
    }
  }
  else {
    THROW(fatal_error,"No scale in hard process");
  }
  msg_Debugging()<<"hard scales = {"<<m_last_q[0]<<","<<m_last_q[1]
		 <<","<<m_last_q[2]<<","<<m_last_q[3]<<"}\n";
  PHASIC::Integrable_Base *proc(p_me->GetAmegic()->GetProcess());
  if (proc->FactorizationScale()!="") {
    m_q2_f[1]=m_q2_f[0]=proc->Scale(stp::fac);
  }
  else {
    if (m_pdf_mode) {
      m_q2_f[0]=sqr(m_qmin[0]);
      m_q2_f[1]=sqr(m_qmin[1]);
      double qmin2me(Min(m_is_as_factor,m_fs_as_factor)*
		     proc->Scale(PHASIC::stp::fac));
      if (!IsEqual(qmin2me,m_is_as_factor*m_qmin[0]*m_qmin[1])) {
	double x[2];
	Combine_Table_Base *ct(p_ct);
	while (ct->Up()) ct=ct->Up();
	ct->GetX1X2(x[0],x[1]);
	msg_Debugging()<<"pdf reweighting: q_{fac,me} = "<<sqrt(qmin2me)
		       <<" -> q_{fac} = "<<(sqrt(m_is_as_factor)*m_qmin[0])
		       <<"/"<<(sqrt(m_is_as_factor)*m_qmin[1])
		       <<", x_1 = "<<x[0]<<", x_2 = "<<x[1]<<"\n";
	static long unsigned int cnt[2]={0,0}, psur[2]={0,0};
	for (short unsigned int i(0);i<2;++i) { 
	  if (p_pdf[i]!=NULL) {
	    ++cnt[i];
	    if (sqr(m_qmin[i])<p_pdf[i]->Q2Min()) {
	      if (++psur[i]/(double)cnt[i]>0.001) {
		msg_Error()<<METHOD<<"(): Scale under-runs minimum PDF scale"
			   <<" in "<<(psur[i]*1000/cnt[i])/10.0
			   <<" % of events."<<std::endl;
	      }
	      continue;
	    }
	    p_pdf[i]->Calculate(x[i],qmin2me);
	    double w(p_pdf[i]->GetXPDF(ct->GetLeg(i).Flav()));
	    p_pdf[i]->Calculate(x[i],m_is_as_factor*sqr(m_qmin[i]));
	    w/=p_pdf[i]->GetXPDF(ct->GetLeg(i).Flav());
	    msg_Debugging()<<"w_{"<<i<<"} = "<<(1.0/w)
			   <<" ("<<ct->GetLeg(i).Flav()<<")\n";
	    if (w>0.0) m_weight/=w;
	  }
	}
      }
    }
  }
  msg_Debugging()<<"} -> w = "<<m_weight<<"\n";
  msg_Debugging()<<"set q_{fac} = "<<sqrt(m_q2_f[0])<<"/"<<sqrt(m_q2_f[1])<<"\n";
}

void Cluster_Partons_CKKW::CalculateWeight(const double &meweight)
{
  if (m_smode&1) p_ct->ShuffleMomenta();
  msg_Debugging()<<METHOD<<"(): {\n";
  msg_Indent();
  ++m_counts;
  InitWeightCalculation();
  WeightHardProcess();
  Combine_Table_Base *ct_tmp(p_ct), *ct_down;
  int si(0), i, j;
  double ptij;
  bool strong_vertex;

  // Internal legs
  while (ct_tmp->Up()) {
    ct_down = ct_tmp;
    ct_tmp  = ct_tmp->Up();
    ++si;
    ++m_nlegs;
    m_last_q.push_back(std::numeric_limits<double>::max());
    m_last_i.push_back(si);
    ptij = ct_tmp->GetWinner(i,j);
    strong_vertex = ct_down->GetLeg(i).OrderQCD()>0;
    if (ct_down->GetLeg(i).Point()->t<10) {
      if (ct_down->GetLeg(i).Flav().Strong())
	m_weight*=InternalWeight(i<2,ct_down->GetLeg(i),m_last_q[i],
				 ptij,ct_tmp->GetLeg(i).QMin());
    }
    else {
      msg_Debugging()<<"finishing production amplitude {\n";
      if (ct_down->GetLeg(i).Flav().Strong()) {
	msg_Indent();
	m_weight*=ExternalWeight(i<2,ct_down->GetLeg(i),m_last_q[i],
				 ct_down->GetLeg(i).QMin());
	if (strong_vertex) {
	  m_last_q[i]=ptij;
	  ptij=ct_tmp->GetLeg(i).QMin();
	}
	else {
	  m_weight*=ExternalWeight(i<2,ct_down->GetLeg(i),ptij,
				   ct_tmp->GetLeg(i).QMin());
	  m_last_q[i]=ptij;
	}
      }
      msg_Debugging()<<"}\n";
    }
    if (strong_vertex) {
      double asw(CouplingWeight(i<2,ct_down->GetLeg(i),ptij));
      m_weight*=asw;
      m_asweight*=asw;
      m_scale*=ptij;
      ++m_hard_nqcd;
    }
    StoreOldValues(i,j,si,ptij);
  }

  // External legs
  for (int l=0; l<m_nlegs; ++l) {
    if (ct_tmp->GetLeg(l).Flav().Strong()) {
      m_weight*=ExternalWeight(l<2,ct_tmp->GetLeg(l),
			       m_last_q[l],ct_tmp->GetLeg(l).QMin());
    }
  }
  msg_Debugging()<<"} -> w = "<<m_weight<<"\n";
  msg_Debugging()<<"\n"<<*ct_tmp<<"\n";
  m_hard_nqed = m_nstrong-m_hard_nqcd-2;
  std::string pid(p_me->GetAmegic()->GetProcess()->Name());
  std::map<std::string,size_t>::iterator iit(m_pidmap.find(pid));
  if (iit==m_pidmap.end()) {
    m_events.resize(m_events.size()+1,0);
    m_meweight_sum.resize(m_events.size(),0.0);
    m_weight_sum.resize(m_events.size(),0.0);
    m_weight_sum_sqr.resize(m_events.size(),0.0); 
    m_asweight_sum.resize(m_events.size(),0.0);
    m_asweight_sum_sqr.resize(m_events.size(),0.0); 
    m_sweight_sum.resize(m_events.size(),0.0);
    m_sweight_sum_sqr.resize(m_events.size(),0.0); 
    m_pidmap.insert(make_pair(pid,m_pidmap.size()));
    iit=m_pidmap.find(pid);
  }
  m_events[iit->second]+=1;
  m_meweight_sum[iit->second]+=meweight;
  m_weight_sum[iit->second]+=m_weight*meweight;
  m_weight_sum_sqr[iit->second]+=sqr(m_weight)*meweight;
  m_asweight_sum[iit->second]+=m_asweight*meweight;
  m_asweight_sum_sqr[iit->second]+= sqr(m_asweight)*meweight;
  m_sweight_sum[iit->second]+=m_weight/m_asweight*meweight;
  m_sweight_sum_sqr[iit->second]+=sqr(m_weight/m_asweight)*meweight;
}

