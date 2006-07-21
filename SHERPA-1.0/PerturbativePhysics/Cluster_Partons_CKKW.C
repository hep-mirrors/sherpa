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
  return p_bp->Gamma(q,m_q)*p_sud->Delta(kf::gluon)(q,m_q0);
}
double Three_Jet_Calc::operator()()
{
  return m_defval; 
} 

Cluster_Partons_CKKW::
Cluster_Partons_CKKW(Matrix_Element_Handler * me,ATOOLS::Jet_Finder * jf,
		     const int maxjet,const int isrmode,const int isron,
		     const int fsron) :
  Cluster_Partons_Base(me,jf,maxjet,isrmode,isron,fsron),
  m_AcceptMisClusters(1.), m_LowestFromME(true)
{
  std::string helps;
  Data_Reader reader;
  if (reader.ReadFromFile(helps,"PRINT_SUDAKOV","") && helps.length()>0) 
    GenerateTables(helps);
}

Cluster_Partons_CKKW::~Cluster_Partons_CKKW()
{
}

void Cluster_Partons_CKKW::GenerateTables(const std::string &path)
{
  msg_Info()<<METHOD<<"("<<path<<"): Generating sudakov tables {"<<std::endl;
    MakeDir(path,448,true);
  for (short unsigned int l(0);l<2;++l) {
    NLL_Sudakov *sud(l==0?p_fssud:p_issud);
    if (sud==NULL) break;
    msg_Indent();
    double ecms(rpa.gen.Ecms()), qmin(1.0), step(pow(ecms/qmin,1.0/25.0));
    for (size_t i(1);i<=7;++i) {
      if (i==7) i=21;
      Flavour f((kf::code)i);
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
      Flavour f((kf::code)i);
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
      Flavour f((kf::code)i);
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
  return new Combine_Table_CKKW(p_jf,amoms,0,m_isrmode,m_isrshoweron);
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
    xs = group->XSSelector()->GetXS(nin,nout,fl,false,nqed,nqcd);
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

  m_hard_nqed = p_xs->OrderEW();
  m_hard_nqcd = p_xs->OrderStrong();
  bool test(p_xs->SetColours(p)), check(true);
  for (int i=0; i<4; ++i) {
    if (fl[i].IsQuark()) {
      if ( fl[i].IsAnti() && 
	   (p_xs->Colours()[i][0]!=0 || p_xs->Colours()[i][1]==0)) check=false;
      if (!fl[i].IsAnti() && 
	  (p_xs->Colours()[i][0]==0 || p_xs->Colours()[i][1]!=0)) check=false;
    }
    if (fl[i].IsGluon() && 
	(p_xs->Colours()[i][0]==0 || p_xs->Colours()[i][1]==0))   check=false;
    if (!check) {
      msg.Error()<<"Cluster_Partons_CKKW::SetColours(..): \n"
		 <<"Colour check failed for the following combination:"
		 <<std::endl;
      for (int i=0; i<4; ++i) 
	msg.Error()<<"   "<<i<<" : "<<fl[i]<<" ("
		   <<p_xs->Colours()[i][0]<<","
		   <<p_xs->Colours()[i][1]<<")"<<std::endl;
      msg.Error()<<"Abort."<<std::endl;
      abort();
    }
  }
  for (int i=0;i<4;i++) {
    m_colors[i][0] = p_xs->Colours()[i][0];
    m_colors[i][1] = p_xs->Colours()[i][1];
  }
  m_q2_fss  = dabs(p_xs->Scale(stp::sfs));
  m_q2_iss  = dabs(p_xs->Scale(stp::sis));
  std::cout<<"cp: found xs, fss "<<m_q2_fss<<", iss "<<m_q2_iss<<"\n";
//   Alternative choice for PS-scale
//    m_q2_fss  = dabs(p_xs->Scale(stp::as));
//    m_q2_iss  = dabs(p_xs->Scale(stp::as));
  if (m_q2_fss==std::numeric_limits<double>::max() || 
      m_q2_iss==std::numeric_limits<double>::max()) {
    m_q2_hard = m_q2_fss = m_q2_iss = p_xs->Scale(stp::fac);
  }
  else {
    m_q2_hard = dabs(p_xs->Scale(stp::as));
  }
  m_q2_qcd  = dabs(p_xs->Scale(stp::as));
  return test;
}


void Cluster_Partons_CKKW::InitWeightCalculation()
{
  m_weight = 1.;
  m_nlegs   = p_ct->NLegs();
  m_njet    = m_nlegs-2;
  m_nstrong = 0;
  for (int l=0; l<m_nlegs; ++l) {
    if (p_ct->GetLeg(l).Point()->fl.Strong()) ++m_nstrong;
  }
  Combine_Table_Base *ct_test(p_ct);
  double qmaxt(0.);
  while (ct_test->Up()) {
    ct_test=ct_test->Up();
    ++m_njet;
    //    if (qmaxt<ct_test->Kt2()) qmaxt=ct_test->Kt2();
  }
  if (qmaxt>m_q2_hard) m_q2_hard = qmaxt;
  m_qmax=sqrt(m_q2_hard);
  m_last_q.clear();
  m_last_i.clear();
  m_last_q.resize(m_nlegs,m_qmax);
  m_last_i.resize(m_nlegs,0);
  // scales for alpha_s and correction weight for hard interaction
  m_as_amegic  = (*as)(sqr(rpa.gen.Ecms()));
  double asfac(sqrt(m_is_as_factor*m_fs_as_factor));
  m_as_hard    = (*p_runas)(m_q2_hard/asfac);
  m_as_qcd     = (*p_runas)(m_q2_qcd/asfac);
  m_as_jet     = (*p_runas)(m_me_as_factor*m_q2_jet/asfac);
  std::cout<<"ct: scales: amegic "<<rpa.gen.Ecms()<<" ("<<m_as_amegic<<")"
		 <<", hard "<<sqrt(m_q2_hard)<<" ("<<m_as_hard<<")"
		 <<", qcd "<<sqrt(m_q2_qcd)<<" ("<<m_as_qcd<<")"
		 <<", jet "<<sqrt(m_q2_jet)<<" ("<<m_as_jet<<"), fac "
		 <<m_is_as_factor<<"/"<<m_fs_as_factor<<", fss "<<sqrt(m_q2_fss)
		 <<", iss "<<sqrt(m_q2_iss)<<" mefac "<<m_me_as_factor<<"\n";
  if (m_kfac!=0.) {
    m_as_hard *= 1. + m_as_hard/(2.*M_PI)*m_kfac;
    m_as_qcd  *= 1. + m_as_qcd/(2.*M_PI)*m_kfac;
    m_as_jet  *= 1. + m_as_jet/(2.*M_PI)*m_kfac;
  }
  if (m_q2_qcd<1.e-6) {
    m_as_jet   = 1.;
    m_as_hard  = 1.;
  }
  std::cout<<"ct: weight me : hard qcd "<<m_hard_nqcd<<" weight "
		 <<pow(m_as_qcd/m_as_jet,m_hard_nqcd)<<std::endl;
  if (m_hard_nqcd>0) m_weight *= pow(m_as_qcd/m_as_jet,m_hard_nqcd);
  // special treatment for effective higgs vertex
  bool found(false);
  for (int i(0);i<p_ct->NAmplitudes();++i) {
    for (int j(0);j<2;++j) {
      switch (p_ct->GetHardLegs()[i][j].Type()) {
      case lf::Triangle:
      case lf::Box:
      case lf::C4GS: {
	double mth2(0.0);
	for (int k(0);k<4;++k) {
	  if (p_ct->GetLegs()[i][k].Point()->fl==kf::h) {
	    mth2=p_ct->Momenta()[k].MPerp2();
	    found=true;
	  }
	  if (p_ct->GetLegs()[i][k].Point()->fl==kf::t) found=false;
	}
	if (!found) break;
	static double asmh2((*MODEL::as)(sqr(Flavour(kf::h).Mass())));
 	m_weight*=pow((*MODEL::as)(mth2)/asmh2,
		      p_ct->GetHardLegs()[i][j].OrderQCD());
	std::cout<<METHOD<<"(): found higgs vertex: as = "<<(*MODEL::as)(mth2)
		       <<" vs. asf = "<<asmh2<<", nqcd = "<<p_ct->GetHardLegs()[i][j].OrderQCD()<<"\n";
	break;
      }
      default:
	break;
      }
      if (found) break;
    }
    if (found) break;
  }

  m_count_startscale=-100;
  if (m_q2_hard!=m_q2_qcd && m_nstrong==3) m_count_startscale=0;

  // determine lowest scale for highest multi treatment
  m_qmin_ci = m_qmin_cf = 0.;
  if ((m_njet==m_maxjetnumber && m_njet>2) ||
      (m_njet==2 && p_ct->OrderStrong()>0)) {
    double qmin2i(0.0), qmin2f(0.0);
    JetvetoPt2(qmin2i,qmin2f);
    if (m_LowestFromME) {
      if (m_q2_qcd<qmin2i) qmin2i = m_q2_qcd;
      if (m_q2_qcd<qmin2f) qmin2f = m_q2_qcd;
    }
    if (qmin2i==0.0 || qmin2f==0.0)
      msg.Error()<<"Cluster_Partons_CKKW::InitWeightCalculation(..): "
		 <<"No minimum scale found."<<std::endl;
    m_qmin_ci=sqrt(qmin2i);
    m_qmin_cf=sqrt(qmin2f);
  }
  std::cout<<"ct: qmin "<<m_qmin_ci<<", "<<m_qmin_cf<<"\n";
}


bool Cluster_Partons_CKKW::
ApplyCombinedInternalWeight(const bool is,const Flavour & fl,
			    const double iupper,const double iactual,
			    const double asref,const int order)
{
  double upper(iupper), actual(iactual);
  if (upper<actual) {
    msg_Tracking()<<METHOD<<"(): q = "<<actual<<", Q = "
		  <<upper<<". Swap."<<std::endl;
    std::swap<double>(upper,actual);
  }
  double qmin(0.), DeltaNum(0.), DeltaDenom(1.), DeltaRatio(0.);
  double as_ptij(0.), asRatio(0.);
  if (is) {
    as_ptij = (*p_runas)(m_me_as_factor*sqr(actual)/m_is_as_factor);
    qmin = m_qmin_ci!=0.0?m_qmin_ci:m_qmin_i;
  }
  else {
    as_ptij = (*p_runas)(m_me_as_factor*sqr(actual)/m_fs_as_factor);
    qmin = m_qmin_cf!=0.0?m_qmin_cf:m_qmin_f;
  }
  if (m_kfac!=0.) as_ptij *= 1. + as_ptij/(2.*M_PI)*m_kfac;
  asRatio = as_ptij/asref;
  if (upper<actual || actual<qmin || asRatio>1.0) {
    ++m_fails;
    if (asRatio>1.0)                 asRatio    = m_AcceptMisClusters;
    if (upper<actual || actual<qmin) DeltaRatio = m_AcceptMisClusters;
  }
  else {
    if (is) {
      DeltaNum   = p_issud->Delta(fl)(upper,qmin);
      DeltaDenom = p_issud->Delta(fl)(actual,qmin);
    }
    else {
      DeltaNum   = p_fssud->Delta(fl)(upper,qmin);
      DeltaDenom = p_fssud->Delta(fl)(actual,qmin);
    }
    DeltaRatio = DeltaNum/DeltaDenom;
  }
  std::cout<<"ct: internal weight: "<<is<<" "<<fl<<" "
		 <<upper<<" -> "<<actual<<" / "<<qmin
		 <<" => delta-> "<<DeltaRatio
		 <<" ("<<DeltaNum<<"/"<<DeltaDenom<<") as-> "
		 <<asRatio<<" from "<<asref<<" "<<order<<std::endl;
  m_weight *= DeltaRatio * asRatio;
  ++m_hard_nqcd;
  return true;
}

bool Cluster_Partons_CKKW::
ApplyExternalWeight(const bool is,const Flavour & fl,
		    const double actual)
{
  double qmin(0.), DeltaNum(0.);
  if (is) qmin = m_qmin_ci!=0.0?m_qmin_ci:m_qmin_i;
     else qmin = m_qmin_cf!=0.0?m_qmin_cf:m_qmin_f;
  if (actual<qmin) {
    ++m_fails;
    DeltaNum = m_AcceptMisClusters;
  }
  else {
    if (is) DeltaNum = p_issud->Delta(fl)(actual,qmin);
       else DeltaNum = p_fssud->Delta(fl)(actual,qmin);
  }
  std::cout<<"ct: external weight: "<<is<<" "<<fl<<" "
		 <<actual<<" -> "<<qmin<<" => delta-> "<<DeltaNum<<std::endl;
  m_weight *= DeltaNum;
  return true;
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

void Cluster_Partons_CKKW::OverwriteScales(const int j)
{
  for (int l=j+1; l<m_nlegs; ++l ) {
    if (m_last_q[l]==m_qmax) m_last_q[l]=sqrt(m_q2_qcd);
  }      
}

void Cluster_Partons_CKKW::CalculateWeight()
{
  ++m_counts;
  InitWeightCalculation();
  Combine_Table_Base *ct_tmp(p_ct), *ct_down;
  int si(0), i, j;
  double ptij;
  bool strong_vertex,singlet_clustered;

  // Internal legs
  while (ct_tmp->Up()) {
    ct_down = ct_tmp;
    ct_tmp  = ct_tmp->Up();
    ++si;
    ++m_nlegs;
    m_last_q.push_back(m_qmax);
    m_last_i.push_back(si);
    ptij = ct_tmp->GetWinner(i,j);
    strong_vertex = ct_down->GetLeg(i).OrderQCD()>0;
    singlet_clustered = 
      !(ct_down->GetLeg(i).Point()->fl.Strong() && 
      ct_tmp->GetLeg(i).Point()->fl.Strong() && 
      ct_tmp->GetLeg(j).Point()->fl.Strong());
    if (strong_vertex) {
      ++m_count_startscale;
      double asref(m_as_jet);
      switch (ct_down->GetLeg(i).Type()) {
      case lf::Triangle:
      case lf::Box:
      case lf::C4GS: 
	std::cout<<"found higgs vertex\n";
	static double asmh2((*MODEL::as)(sqr(Flavour(kf::h).Mass())));
	asref=asmh2;
	break;
      default: 
	break;
      }
      ApplyCombinedInternalWeight
	(i<2,ct_down->GetLeg(i).Point()->fl,m_last_q[i],ptij,asref,
	 ct_down->GetLeg(i).OrderQCD());
  }
    if (strong_vertex || singlet_clustered) StoreOldValues(i,j,si,ptij);
    if (m_count_startscale==2) OverwriteScales(j);
  }

  // External legs
  for (int l=0; l<m_nlegs; ++l) {
    if (ct_tmp->GetLeg(l).Point()->fl.Strong()) {
      if (m_last_q[l]==m_qmax) { 
	++m_count_startscale; 
      }
      if (!ApplyExternalWeight(l<2,ct_tmp->GetLeg(l).Point()->fl,
			       m_last_q[l])) {
      }
    }
    if (m_count_startscale==2) OverwriteScales(l);
  }
  std::cout<<METHOD<<"(..): combine tables {\n"<<*ct_tmp<<"\n}\n";

  m_hard_nqed = m_nstrong-m_hard_nqcd-2;
  if (m_nstrong-m_hard_nqcd>2) 
    m_weight *= pow(m_as_amegic/m_as_jet,m_hard_nqed);
  p_events[m_njet-1]         += 1;
  p_weight_sum[m_njet-1]     += m_weight;
  p_weight_sum_sqr[m_njet-1] += sqr(m_weight);
}

