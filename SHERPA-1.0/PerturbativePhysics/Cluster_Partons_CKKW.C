#include "Cluster_Partons_CKKW.H"

#include "Combine_Table_CKKW.H"

using namespace SHERPA;
using namespace APACIC;
using namespace EXTRAXS;
using namespace AMEGIC;
using namespace ATOOLS;
using namespace PHASIC;
using namespace MODEL;


static double s_counts(0.), s_fails(0.);

Cluster_Partons_CKKW::Cluster_Partons_CKKW(Matrix_Element_Handler * me, ATOOLS::Jet_Finder * jf, 
					   int maxjetnumber, int isrmode, 
					   int isrshoweron, int fsrshoweron) :
  Cluster_Partons_Base(me,jf,maxjetnumber,isrmode,isrshoweron,fsrshoweron)
{
}

Cluster_Partons_CKKW::~Cluster_Partons_CKKW()
{
}

Combine_Table_Base *Cluster_Partons_CKKW::
CreateTable(Jet_Finder *jf,ATOOLS::Vec4D *amoms,Combine_Table_Base *ct,
	    const int isrmode,const int isrshoweron)
{
//   PRINT_INFO("ckkw "<<amoms);
  return new Combine_Table_CKKW(p_jf,amoms,0,m_isrmode,m_isrshoweron);
}

EXTRAXS::XS_Base *Cluster_Partons_CKKW::GetXS(EXTRAXS::XS_Group * group, ATOOLS::Flavour * fl)
{
  XS_Base * xs(NULL);  
  const size_t nin(2), nout(2), n(nin+nout);
  if (group->XSSelector()->FindInGroup(group,xs,nin,nout,fl)==std::string::npos) {
    int nstrong(0),nqed(0),nqcd(0);
    for (size_t i=0;i<n;++i) nstrong+=fl[i].Strong();
    p_ct->AddCouplings(nqed,nqcd);
    //      std::cout<<" Cluster_Partons::GetXS nqed="<<nqed<<" nqcd="<<nqcd<<std::endl;
    Process_Base * proc=static_cast<Process_Base*>(p_me->GetAmegic()->GetProcess());
    //      std::cout<<proc->Name()<<"("<<proc->OrderEWeak()<<","<<proc->OrderStrong()<<")"<<std::endl;
    int nrqed(proc->OrderEWeak() - nqed);
    int nrqcd(proc->OrderStrong() - nqcd);
    xs = group->XSSelector()->GetXS(nin,nout,fl,false,nrqed,nrqcd);
    if (xs) group->Add(xs);
  }
  p_xs = xs;
  return xs;
}

int Cluster_Partons_CKKW::SetColours(EXTRAXS::XS_Base * xs, ATOOLS::Vec4D * p, ATOOLS::Flavour * fl)
{
  p_xs=xs;
  if (!p_xs) return Cluster_Partons_Base::SetColours(p,fl);
  m_hard_nqed = p_xs->OrderEW();
  m_hard_nqcd = p_xs->OrderStrong();
  bool test(p_xs->SetColours(p)), check(true);
  for (int i=0; i<4; ++i) {
    if (fl[i].IsQuark()) {
      if ( fl[i].IsAnti() && (p_xs->Colours()[i][0]!=0 || p_xs->Colours()[i][1]==0)) check=false;
      if (!fl[i].IsAnti() && (p_xs->Colours()[i][0]==0 || p_xs->Colours()[i][1]!=0)) check=false;
    }
    if (fl[i].IsGluon()) {
      if (p_xs->Colours()[i][0]==0 || p_xs->Colours()[i][1]==0) check=false;
    }
    if (!check) break;
  }
  if (!check) {
    std::cout<<" colour check failed !!!\n";
    for (int i=0; i<4; ++i) {
      std::cout<<i<<" : "<<fl[i]<<" ("<<p_xs->Colours()[i][0]<<","<<p_xs->Colours()[i][1]<<") \n";
    }
    abort();
  }
//   PRINT_INFO("fac scale = "<<p_xs->Scale(stp::fac)<<", qcd scale = "<<p_xs->Scale(stp::as));
  m_q2_hard=p_xs->Scale(stp::fac);
  m_q2_qcd=p_xs->Scale(stp::as);
  return test;
}


void Cluster_Partons_CKKW::InitWeightCalculation()
{
  m_weight = 1.;
  m_nlegs=p_ct->NLegs();
//   PRINT_INFO("set nlegs "<<m_nlegs);
  m_njet=m_nlegs-2;
  Combine_Table_Base *ct_test(p_ct);
  double qmaxt(0.);
  while (ct_test->Up()) {
    int i, j;
    ct_test->Up()->GetWinner(i,j);
    bool strong_vertex(ct_test->GetLeg(i).Point()->fl.Strong() && 
		       ct_test->Up()->GetLeg(i).Point()->fl.Strong() && 
		       ct_test->Up()->GetLeg(j).Point()->fl.Strong());
    ct_test=ct_test->Up();
    ++m_njet;
    if (strong_vertex && qmaxt<ct_test->Kt2()) qmaxt=ct_test->Kt2();
  }
  if (qmaxt>m_q2_hard) {
    msg.Error()<<"Cluster_Partons_CKKW::InitWeightCalculation(): Hardest qcd scale exceeds hard scale.\n"
	       <<"   hardest = "<<qmaxt<<", hard = "<<m_q2_hard<<std::endl;
    m_q2_hard = qmaxt;
  }
//   if (qmaxt>m_q2_qcd) {
//     msg.Error()<<"Cluster_Partons_CKKW::InitWeightCalculation(): Hardest qcd scale exceeds qcd scale.\n"
// 	       <<"   hardest = "<<qmaxt<<", qcd = "<<m_q2_qcd<<std::endl;
//     m_q2_qcd = qmaxt;
//   }
  m_qmax=sqrt(m_q2_hard);
//    PRINT_INFO("init last_q with "<<m_qmax);
  m_last_q.clear();
  m_last_i.clear();
  m_last_q.resize(m_nlegs,m_qmax);
  m_last_i.resize(m_nlegs,0);
  m_as_amegic = (*as)(sqr(rpa.gen.Ecms()));
  m_as_jet    = (*p_runas)(m_q2_jet);
  m_as_hard   = (*p_runas)(m_q2_hard);
//    PRINT_INFO("as(jet) = "<<m_as_jet<<", as_hard = "<<m_as_hard<<" ^ "<<m_hard_nqcd
//  	     <<" <- "<<sqrt(m_q2_hard));
  if (m_q2_qcd<1.e-6) {
    m_as_jet  = 1.;
    m_as_hard = 1.;
  }

  // one might be tempted to remove (for e+e- -> Jets) universal 
  // 2 Quark sudakov to increase the effectivity!
  // Quark sudakov for each strong interacting particle
  m_nstrong = 0;
  for (int l=0; l<m_nlegs; ++l) {
    if (p_ct->GetLeg(l).Point()->fl.Strong()) ++m_nstrong;
  }

  // weight with correct alphas at hard scale (if needed).
  if (m_hard_nqcd>0) m_weight *= pow(m_as_hard/m_as_jet,m_hard_nqcd);
  //std::cout<<" qcd_count="<<m_hard_nqcd<<" ("<<(strong-2)<<")\n";
  //std::cout<<" hard As : "<<asscale<<"\n";

  m_count_startscale=-10;
  // make sure leg between hardscale and asscale is only taken into account once
  if (m_q2_hard!=m_q2_qcd && m_nstrong==3) {
    m_count_startscale=0;
  }

  // determine lowest scale
  m_qmin=0.;
  if (m_njet==m_maxjetnumber) {
    double qmin2(0.0);
    FixJetvetoPt2(qmin2);
    if (qmin2==0.0)
      msg.Error()<<"Cluster_Partons_CKKW::InitWeightCalculation(..): No minimum scale found."<<std::endl;
    m_qmin=sqrt(qmin2);
  }
//   PRINT_INFO("nlegs "<<m_njet<<" max "<<m_maxjetnumber<<" qmin "<<m_qmin);
}


bool Cluster_Partons_CKKW::ApplyCombinedInternalWeight(const bool is,const Flavour & fl,
						       const double upper,const double actual)
{
  double qmin, DeltaNum, DeltaDenom, DeltaRatio, as_ptij, asRatio;
  if (is) {
    qmin = Max(m_qmin,m_qmin_i);
    msg_Tracking()<<"is internal : "<<fl<<" up="<<upper<<" low="<<actual<<" jet="<<qmin<<std::endl;
    DeltaNum   = p_issud->Delta(fl)(upper,qmin);
    DeltaDenom = p_issud->Delta(fl)(actual,qmin);
    
  }
  else {
    qmin = Max(m_qmin,m_qmin_f);
    msg_Tracking()<<"fs internal : "<<fl<<" up="<<upper<<" low="<<actual<<" jet="<<qmin<<std::endl;
    DeltaNum   = p_sud->Delta(fl)(upper,qmin);
    DeltaDenom = p_sud->Delta(fl)(actual,qmin);
  }
  DeltaRatio = DeltaNum/DeltaDenom;
  
  as_ptij = (*p_runas)(sqr(actual));
  if (m_kfac!=0.) as_ptij *= 1. + as_ptij/(2.*M_PI)*m_kfac;
  asRatio = as_ptij/m_as_jet;

  msg_Tracking()<<"weights : as "<<as_ptij<<" asweight: "
	    <<asRatio<<" sud num: "<<DeltaNum<<" sud denom: "<<DeltaDenom<<std::endl;
  
  if (upper<actual || actual<qmin || asRatio>1.0) {
    if (asRatio>1.0) 
      msg.Error()<<"Cluster_Partons::ApplyCombinedInternalWeight(..): "
		 <<("\\alpha_s weight exceeds unity.\n   w = "
		    +ToString(asRatio)+", \\mu = "
		    +ToString(sqr(actual))+", \\mu_{jet} = "
		    +ToString(m_q2_jet))<<std::endl;
    else 
      msg.Error()<<"Cluster_Partons::ApplyCombinedInternalWeight(..): "
		 <<("Internal scale ordering violated.\n   \\mu_1 = "
		    +ToString(upper)+", \\mu_2 = "+ToString(actual)
		    +", \\mu_0 = "+ToString(qmin))<<std::endl;
    m_weight=0.0;
    return false;
  }
  m_weight *= DeltaRatio * asRatio;
  //std::cout<<" in  Delta("<<ct_down->GetLeg(i)->fl<<") : "
  //<<last_q[i]<<","<<ptij<<","<<qmin<<" = "<<DeltaR<<"\n";
  ++m_hard_nqcd;
  //std::cout<<" qcd_count="<<m_hard_nqcd<<"\n";
  //std::cout<<" in  As : "<<ptij<<"\n";
  return true;
}

bool Cluster_Partons_CKKW::ApplyExternalWeight(const bool is,const Flavour & fl,
					       const double actual)
{
  double qmin, DeltaNum;
  if (is) {
    qmin     = Max(m_qmin,m_qmin_i);
    msg_Tracking()<<"is external : "<<fl<<" low="<<actual<<" jet="<<qmin<<std::endl;
    DeltaNum = p_issud->Delta(fl)(actual,qmin);
  }
  else {
    qmin     = Max(m_qmin,m_qmin_f);
    msg_Tracking()<<"fs external : "<<fl<<" low="<<actual<<" jet="<<qmin<<std::endl;
    DeltaNum = p_sud->Delta(fl)(actual,qmin);
  }
  if (actual<qmin) {
    msg.Error()<<"Cluster_Partons::ApplyExternalWeight(..): "
	       <<("External scale ordering violated.\n   \\mu = "
		  +ToString(actual)+", \\mu_0 = "
		  +ToString(qmin))<<std::endl;
    m_weight = 0.0;
    return false;
  }
  m_weight *= DeltaNum;
  //std::cout<<" out Delta("<<p_ct->GetLeg(l)->fl<<") : "<<last_q[l]<<","<<qmin<<"\n";
  return true;
}    

void Cluster_Partons_CKKW::StoreOldValues(const int i,const int j,
					  const int si,const double ptij)
{
  // store old q values
  m_last_i[i] = si;
  m_last_q[i] = ptij;
  for (int l=m_nlegs-1;l>j;--l) {
    m_last_i[l] = m_last_i[l-1];
    m_last_q[l] = m_last_q[l-1];
  } 
  m_last_i[j] = si;
  m_last_q[j] = m_last_q[i];
}  

void Cluster_Partons_CKKW::CalculateWeight()
{
  InitWeightCalculation();
//   PRINT_INFO("jets "<<m_njet<<" hard: "<<sqrt(m_q2_hard)<<", qcd "<<sqrt(m_q2_qcd)<<", jet "<<sqrt(m_q2_jet));
  Combine_Table_Base *ct_tmp(p_ct), *ct_down;
  int si(0), i, j;
  double ptij;
  bool strong_vertex,singlet_clustered;
//   PRINT_INFO("--------------------");
  ++s_counts;
//   PRINT_INFO("core : "<<p_ct->GetLeg(0).Point()->fl<<" "<<p_ct->GetLeg(1).Point()->fl
// 	     <<" -> "<<p_ct->GetLeg(2).Point()->fl<<" "<<p_ct->GetLeg(3).Point()->fl);

  while (ct_tmp->Up()) {
    ct_down = ct_tmp;
    ct_tmp  = ct_tmp->Up();
    ++si;
    ++m_nlegs;
    m_last_q.push_back(m_qmax);
    m_last_i.push_back(si);
    ptij = ct_tmp->GetWinner(i,j);
    strong_vertex = 
      ct_down->GetLeg(i).Point()->fl.Strong() && 
      ct_tmp->GetLeg(i).Point()->fl.Strong() && 
      ct_tmp->GetLeg(j).Point()->fl.Strong();
    singlet_clustered = 
      !(ct_down->GetLeg(i).Point()->fl.Strong()) && 
      ct_tmp->GetLeg(i).Point()->fl.Strong() && 
      ct_tmp->GetLeg(j).Point()->fl.Strong();

//     PRINT_INFO(ct_down->GetLeg(i).Point()->fl<<" -> "<<ct_tmp->GetLeg(i).Point()->fl
// 	       <<" "<<ct_tmp->GetLeg(j).Point()->fl<<" -> "<<ptij);

//     PRINT_INFO(strong_vertex<<" "<<singlet_clustered<<" "<<m_last_q[i]<<" "<<ptij<<" "<<m_qmax);
    if (strong_vertex) {
      ++m_count_startscale;
//       PRINT_INFO("i : "<<i<<", j : "<<j<<", size = "<<m_last_q.size()<<", last q "<<m_last_q[i]);
      if (!ApplyCombinedInternalWeight(i<2,ct_down->GetLeg(i).Point()->fl,m_last_q[i],ptij)) {
	++s_fails;
	PRINT_INFO("failure rate : "<<(s_fails/s_counts));
      }
    }
    if (strong_vertex || singlet_clustered) StoreOldValues(i,j,si,ptij);

    // check that not too many hardscales are present!
    //     if (count_startscale==2) {
    //       for (int l=j+1; l<nlegs; ++l ) {
    // 	if (last_q[l]==qmax) {
    // 	  last_q[l]=sqrt(asscale);
    // 	}
    //       }      
    //     }
  }

//   PRINT_INFO(*ct_tmp);

  // go over all remaining legs
//     PRINT_INFO("external weight "<<*ct_tmp);
  for (int l=0; l<m_nlegs; ++l) {
//     PRINT_INFO("external weight "<<l<<" "<<m_nlegs);
    if (ct_tmp->GetLeg(l).Point()->fl.Strong()) {
      ++m_count_startscale;
      if (!ApplyExternalWeight(l<2,ct_tmp->GetLeg(l).Point()->fl,m_last_q[l])) {
	++s_fails;
	PRINT_INFO("failure rate : "<<(s_fails/s_counts));
      }

    }
    // check that not too many hardscales are present!
    //       if (count_startscale==2) {
    // 	for (int k=l+1; k<nlegs; ++k ) {
    // 	  if (last_q[k]==qmax) {
    // 	    last_q[k]=sqrt(asscale);
    // 	  }
    // 	}      
    //       }
  }
    

  // replace remaining strong couplings
  //std::cout<<" remaining strong couplings ( "<<m_hard_nqed<<","<<m_hard_nqcd<<") vs "
  //<<m_nstrong-2<<std::endl;

  m_hard_nqed = m_nstrong-m_hard_nqcd-2;    
  if (m_nstrong-m_hard_nqcd>2) {
    m_weight *= pow(m_as_amegic/m_as_jet,m_hard_nqed);
    //std::cout<<" cancel As : "<<as_amegic/as_jet<<"^("<<m_hard_nqed<<" ) \n";
  }
  
  //std::cout<<" weight="<<m_weight<<"\n";
  
  p_events[m_njet-1]         += 1;  // count events
  p_weight_sum[m_njet-1]     += m_weight;
  p_weight_sum_sqr[m_njet-1] += sqr(m_weight);
}

