#include "Cluster_Partons.H"
#include "Message.H"
#include "Run_Parameter.H"
#include "MyStrStream.H"
#include "Flow.H"
#include "XS_Selector.H"
#include "Matrix_Element_Handler.H"
#include "Initial_State_Shower.H"
#include "Final_State_Shower.H"
#include "Running_AlphaS.H"
#include "Amegic.H"
#include <cassert>
#include <iomanip>

using namespace SHERPA;
using namespace APACIC;
using namespace EXTRAXS;
using namespace AMEGIC;
using namespace ATOOLS;
using namespace MODEL;

Cluster_Partons::Cluster_Partons(Matrix_Element_Handler * me, Jet_Finder * jf,
				 int maxjetnumber, int isrmode,int isrshoweron, int fsrshoweron) :
  p_me(me),p_jf(jf),m_maxjetnumber(maxjetnumber),
  m_isrmode(isrmode), m_isrshoweron(isrshoweron), 
  m_fsrshoweron(fsrshoweron), m_kfac(0.), p_local_tree(NULL) 
{


  // read in some parameters
  Data_Read dr(rpa.GetPath()+"Shower.dat");     // !!!!!!!! SHOWER_DATA_FILE

  m_mode = Init(); // !!! ew merging see also Combine_Table
  m_mode = dr.GetValue<int>("EW_MERGING",m_mode);
  if (m_mode) {
    msg_Info()<<"Cluster_Partons: EW merging mode "<<m_mode<<std::endl;
  }

  m_bp_mode  = dr.GetValue<int>("SUDAKOVTYPE",0);
  m_as_order = dr.GetValue<int>("SUDAKOVASORDER",0);
  //  double as_fac = dr.GetValue<double>("SUDAKOVASFAC",1.);
  double as_fac = rpa.gen.RenormalizationScaleFactor();

  int jetratemode = dr.GetValue<int>("CALCJETRATE",-1);
  if ((m_bp_mode&(7+8+16+32+64))!=m_bp_mode) {
    msg.Error()<<"WARNING in Cluster_Partons :"<<std::endl
	       <<"   Wrong mode for NLL_Sudakovs: "<<m_bp_mode<<" vs "<<(m_bp_mode&127)<<std::endl
	       <<"   Set it to 0 = ordinary NLL_Sudakovs."<<std::endl;
    m_bp_mode=0;
  }
  //    m_kfac=  CA*(67./18.-M_PI*M_PI/6.)-10./9.*TR*Nf;  == 3.45409 (for Nf=5)
  if (m_bp_mode&16) m_kfac = 3.45409;
  msg_Tracking()<<"Cluster_Partons runs in mode : \n"
		<<"   SUDAKOVTYPE    = "<<m_bp_mode<<"\n"
		<<"   SUDAKOVASORDER = "<<m_as_order<<"\n"
		<<"   SUDAKOVASFAC   = "<<as_fac<<"\n"
		<<"   CALCJETRATE    = "<<jetratemode<<"\n"
		<<"   kfac           = "<<m_kfac<<std::endl;

  p_runas = 0;
  if (m_as_order!=-1 || as_fac!=1.) {
    double mz2 = sqr(Flavour(kf::Z).Mass());
    double as_mz = (*as)(mz2/as->ScaleFactor());
    p_runas = new Running_AlphaS(as_mz,mz2,m_as_order,as_fac);
  }

  p_combi = 0;


  /* 0 no sudakow weights, 1 alphas only, 2 full sudakov weight  (but for highest jet number) */
  /* cf. also begin of Cluster_Partons::CalculateWeight() */
  if (m_isrshoweron==0 && m_fsrshoweron==0) {
    p_sud = NULL; m_sud_mode = 0;
  }
  else {
    p_sud   = new NLL_Sudakov(m_bp_mode,p_jf->Smax(),p_jf->Smin(),p_runas,jetratemode);
    m_sud_mode = 2;  
  }

  if (p_runas==0) { p_runas = as; }
  

  p_events   = new long[maxjetnumber];
  p_weight_sum     = new double[maxjetnumber];
  p_weight_sum_sqr = new double[maxjetnumber]; 
  for (int i=0;i<maxjetnumber;++i) {
    p_events[i]=0;
    p_weight_sum[i]=p_weight_sum_sqr[i]=0.;
  }
}

int Cluster_Partons::Init() 
{
  Process_Base * procs=p_me->GetAmegic()->Processes();
  if (CheckProcess(procs)<0) return 1;
  return 0;
}

int Cluster_Partons::CheckProcess(Process_Base * procs) 
{
  int check=0;
  for (size_t i=0; i<procs->Size(); ++i) {
    int test=0;
    if ((*procs)[i]!=(*(*procs)[i])[0]) test = CheckProcess((*procs)[i]);
    else {
      int ncpl = (*procs)[i]->OrderStrong()+(*procs)[i]->OrderEWeak();
      int nexpected = (*procs)[i]->NIn()+(*procs)[i]->NOut()-2;
      test=nexpected-ncpl;

      std::cout<<" checking "<<(*procs)[i]->Name()<<" : "<<test<<std::endl;
    }
    check = Min(check,test);
  }

  return check;
}




void Cluster_Partons::WriteOutSudakovWeights() 
{
  msg_Info()<<" Statistics Sudakov Rejection "<<std::endl;

  for (int i=0;i<m_maxjetnumber;++i) {
    if (p_events[i]==0) continue;
    double w_mean  = p_weight_sum[i]/p_events[i];
    double w_delta = 1./p_events[i] * sqrt(dabs(p_weight_sum_sqr[i] - 
					   (sqr(p_weight_sum[i])-p_weight_sum_sqr[i])/(p_events[i]-1.)));
    double w_sigma =  sqrt(dabs(( sqr(p_weight_sum[i]/p_events[i])
				    - p_weight_sum_sqr[i]/p_events[i]))/(p_events[i]-1.) );
    msg_Info()<<(i+1)<<" : weight="<<w_mean<<" +- "<<w_delta<<" ("<<w_sigma<<")"<<std::endl;
    p_weight_sum[i]=w_mean;
    p_weight_sum_sqr[i]=w_sigma;
  }
  
  MyStrStream sstr;
  int ecms = int(rpa.gen.Ecms()*10.);
  double ycut = log(rpa.gen.Ycut())/log(10.);
  int mode = m_sud_mode&3 + m_bp_mode;
  sstr<<"sudweights_"<<ecms<<"_"<<mode<<".dat"<<std::endl;
  std::string filename;
  sstr>>filename;
  std::ofstream  rfile(filename.c_str(),std::ios::app);
  rfile.precision(6);
  rfile<<std::setw(10)<<ycut;
  for (int i=1;i<5;++i) {
    //    if (p_events[i]==0) continue;
    if (i<m_maxjetnumber) 
      rfile<<" "<<std::setw(10)<<p_weight_sum[i]<<" "<<std::setw(10)<<p_weight_sum_sqr[i];
    else
      rfile<<" "<<std::setw(10)<<0.<<" "<<std::setw(10)<<0.;
  }  
  rfile<<std::endl;
}


Cluster_Partons::~Cluster_Partons()
{
  if (m_as_order && p_runas) delete p_runas;


  if (p_local_tree) delete p_local_tree;
  if (p_combi) delete p_combi;
  if (p_sud)  delete p_sud;
  
  WriteOutSudakovWeights();
  
  if (p_events) delete [] p_events;
  if (p_weight_sum) delete [] p_weight_sum;
  if (p_weight_sum_sqr) delete [] p_weight_sum_sqr;
}

bool Cluster_Partons::ClusterConfiguration(Blob * blob,double x1,double x2) {
  p_blob          = blob;
  int nin         = p_blob->NInP();
  int nout        = p_blob->NOutP();

  if (nin==1) {
    if (nout<4) return 1;
    msg.Error()<<"Error in Cluster_Partons::ClusterConfiguration()"<<std::endl
	       <<"   Try to cluster decay blob, nin ="<<nin<<" with nout = "<<nout<<","<<std::endl
	       <<"   No method provided yet. Return 0."<<std::endl;
    return 0;
  }

  int nampl       = p_me->NumberOfDiagrams();
  int    nlegs    = nin + nout;
  Leg ** legs     = 0;
  bool reuse      = 0;

  // start cluster algorithm :
  if (!reuse) {
    if (p_combi) delete p_combi;
    p_combi    = 0;

    // generate a list of "legs" for each amplitude
    legs = new Leg *[nampl];
    for (int k=0;k<nampl;) {
      legs[k] = new Leg[nlegs];
      int l   = 0;
      if (FillLegs(legs[k],p_me->GetDiagram(k),l,nlegs)) {
	++k;
      } 
      else {
	delete [] legs[k];
	--nampl;
      }
    }
  }  

  p_ct = 0;
  // if no combination table exist, create it
  if (!p_combi) {
    /*
      - copy moms to insert into Combine_Table (will be delete there)
      - create new Combine_Table with given momenta and given Jet-measure
      - initialise Combine_Table
      - determine best combination sheme
    */ 
    Vec4D * amoms = new Vec4D[nlegs];
    for (int i=0;i<nin;++i)  amoms[i]     = blob->InParticle(i)->Momentum();
    for (int i=0;i<nout;++i) amoms[nin+i] = blob->OutParticle(i)->Momentum();
    if (p_me->InSwaped()) {
      // avoid flavour mismatch if using amplitudes
      Vec4D help=amoms[0];
      amoms[0]=amoms[1];
      amoms[1]=help;
    }

    p_combi = new Combine_Table(p_jf,amoms,0,m_isrmode,m_isrshoweron, m_mode);
    p_combi->FillTable(legs,nlegs,nampl);   
    p_ct = p_combi->CalcJet(nlegs,x1,x2); 
  }
  else {
    // use the existing combination table and determine best combination sheme
    Vec4D * amoms = new Vec4D[nlegs];
    for (int i=0;i<nin;++i)  amoms[i]     = blob->InParticle(i)->Momentum();
    for (int i=0;i<nout;++i) amoms[nin+i] = blob->OutParticle(i)->Momentum();
    if (p_me->InSwaped()) {
      // avoid flavour mismatch if using amplitudes
      Vec4D help=amoms[0];
      amoms[0]=amoms[1];
      amoms[1]=help;
    }
    p_ct = p_combi->CalcJet(nlegs,x1,x2,amoms);
    delete [] amoms;
  }

  CreateFlavourMap();

  //  msg_Tracking()
//      std::cout<<"========================================\n"
// 	      <<*p_combi<<std::endl;

  return 1;
}


bool Cluster_Partons::FillLegs(Leg * alegs, Point * root, int & l, int maxl) {
  if (l>= maxl) {
    msg.Error()<<" Error in Cluster_Partons::FillLegs() !!! "<<std::endl;
    return 0;
  }
  if (l==0) {
    alegs[root->number]=Leg(root);
    l++;
  }
  if (root->left) {
    if (root->middle) return 0; // four vertex 

    if (root->cpl[1]!=Complex(0.,0.)) {
      if ((root->fl.IsScalar() && root->left->fl.IsVector() && root->right->fl.IsVector()) ||
	  (root->fl.IsVector() && root->left->fl.IsScalar() && root->right->fl.IsVector()) ||
	  (root->fl.IsVector() && root->left->fl.IsVector() && root->right->fl.IsScalar())) {
// 	std::cout<<" reset VSS vertex "<<root->fl<<" "<<root->left->fl<<" "<<root->right->fl<<" \n";
	  root->cpl[1]=Complex(0.,0.);
      }
    }


    return FillLegs(alegs,root->left,l,maxl)*FillLegs(alegs,root->right,l,maxl);
  } 
  else {
    alegs[root->number]=Leg(root);
    l++;
    return 1;
  }
}

void Cluster_Partons::CalculateWeight(double hardscale,double asscale, double jetscale,double qm_i,double qm_f) 
{
//      std::cout<<" CalculateWeight("<<hardscale<<","<<asscale<<","<<jetscale<<","
//      	   <<qm_i*qm_i<<","<<qm_f*qm_f<<") called \n";

  m_qmin_i = qm_i;
  m_qmin_f = qm_f;
  double qmax = sqrt(hardscale);

  int nlegs   = p_ct->NLegs();

  // count jets
  int njet=nlegs-2;
  Combine_Table * ct_test = p_ct;
  while (ct_test->Up()) {
    ct_test=ct_test->Up();
    ++njet;
  }

  m_weight      = 1.;

  int si = 0;
  std::vector<double> last_q(nlegs,qmax);
  std::vector<int>    last_i(nlegs,si);
  double as_amegic = (*as)(sqr(rpa.gen.Ecms()));
  double as_jet    = (*as)(jetscale);
  double as_hard   = (*p_runas)(asscale);

  if (jetscale<1.e-6 || asscale<1.e-6) {
    as_jet  = 1.;
    as_hard = 1.;
  }

  // one might be tempted to remove (for e+e- -> Jets) universal 
  // 2 Quark sudakov to increase the effectivity!
  int strong=0;
  for (int l=0; l<nlegs; ++l) {
    if (p_ct->GetLeg(l)->fl.Strong()) {// Quark sudakov for each strong interacting particle
      strong++;
    }    
  }

  // weight with correct alphas at hard scale (if needed).
  if (m_hard_nqcd>0 && m_sud_mode%10>0) {
    m_weight *= pow(as_hard/as_jet,m_hard_nqcd);
//        std::cout<<" qcd_count="<<m_hard_nqcd<<" ("<<(strong-2)<<")\n";
//        std::cout<<" hard As : "<<asscale<<"\n";
  }

  int count_startscale=-10;
  // make sure leg between hardscale and asscale is only taken into account twice
  if (hardscale!=asscale && strong==3) {
    count_startscale=0;
  }

  // determine lowest scale
  if (njet==m_maxjetnumber) {
    jetscale=asscale;
    FixJetvetoPt2(jetscale);
    m_qmin_i=Max(m_qmin_i,sqrt(jetscale));
    m_qmin_f=Max(m_qmin_f,sqrt(jetscale));
  }



  Combine_Table * ct_tmp = p_ct;

  // anti-cluster and ...
  while (p_ct->Up()) {
    // ... add sudakov and alphaS factor for each clustering
    Combine_Table * ct_down = p_ct;
    p_ct=p_ct->Up();
    ++si;
    ++nlegs;

    // make space in vector:
    last_q.push_back(qmax);
    last_i.push_back(si);

    // ... determine winner: i,j, and yij
    int i,j;
    double ptij = p_ct->GetWinner(i,j);

    bool strong_vertex = ct_down->GetLeg(i)->fl.Strong()
      && p_ct->GetLeg(i)->fl.Strong() && p_ct->GetLeg(j)->fl.Strong();
    bool singlet_clustered = !(ct_down->GetLeg(i)->fl.Strong())
      && p_ct->GetLeg(i)->fl.Strong() && p_ct->GetLeg(j)->fl.Strong();
//     if (singlet_clustered) std::cout<<" singlet clustered "<<std::endl;
//     else std::cout<<" no singlet clustered "<<std::endl;

    if (m_sud_mode%10>1 && 
	( strong_vertex || 
	  (m_mode&2==0 && ct_down->GetLeg(i)->fl.Strong()))) {  // sudakov form factor
      double qmin = m_qmin_f;
      if (i<2) qmin = m_qmin_i;
      double w_in =  p_sud->Delta(ct_down->GetLeg(i)->fl)(last_q[i],qmin)/
	             p_sud->Delta(ct_down->GetLeg(i)->fl)(ptij,qmin);
      m_weight *= w_in;
//       std::cout<<" in  Delta("<<ct_down->GetLeg(i)->fl<<") : "<<last_q[i]<<","<<ptij<<","<<qmin<<" = ";
//       std::cout<<w_in<<"\n";
      ++count_startscale;
    }

    if (m_sud_mode%10>0 && strong_vertex) {    // alphaS factor
      double a = (*p_runas)(ptij*ptij);
      double w_in_as = a/as_jet;
      if (m_kfac!=0.) 	w_in_as *= 1. + a/(2. * M_PI)*m_kfac;
      m_weight *= w_in_as;
      ++m_hard_nqcd;
//          std::cout<<" qcd_count="<<m_hard_nqcd<<"\n";
//          std::cout<<" in  As : "<<ptij<<"\n";
    }

    // store old q values
    last_i[i] = si;
    if (m_mode&2==0 || strong_vertex || singlet_clustered)
      last_q[i] = ptij;
    for (int l=nlegs-1;l>j;--l) {
      last_i[l] = last_i[l-1];
      last_q[l] = last_q[l-1];
    } 
    last_i[j] = si;
    last_q[j] = last_q[i];

    // check that not too many hardscales are present!
    if (count_startscale==2) {
      for (int l=j+1; l<nlegs; ++l ) {
	if (last_q[l]==qmax) {
	  last_q[l]=sqrt(asscale);
	}
      }      
    }
  }


  // go over all remaining legs
  if (m_sud_mode%10>1) {
    for (int l=0; l<nlegs; ++l) {
      double qmin = m_qmin_f;
      if (l<2) qmin = m_qmin_i;
      if (p_ct->GetLeg(l)->fl.Strong() && last_q[l]>qmin) {//  sudakov form factor
	double w_out = p_sud->Delta(p_ct->GetLeg(l)->fl)(last_q[l],qmin);
	m_weight *= w_out;
//   	std::cout<<" out Delta("<<p_ct->GetLeg(l)->fl<<") : "<<last_q[l]<<","<<qmin<<"\n";
      }    

      // check that not too many hardscales are present!
      if (p_ct->GetLeg(l)->fl.Strong()) {
	++count_startscale;
      }
      if (count_startscale==2) {
	for (int k=l+1; k<nlegs; ++k ) {
	  if (last_q[k]==qmax) {
	    last_q[k]=sqrt(asscale);
	  }
	}      
      }
    }
  }

  // replace remaning strong couplings
//    std::cout<<" remaining strong couplings ( "<<m_hard_nqed<<","<<m_hard_nqcd<<") vs "<<m_nstrong-2<<std::endl;
  m_hard_nqed = m_nstrong-m_hard_nqcd-2;    
  if (m_nstrong-m_hard_nqcd>2) {
    m_weight *= pow(as_amegic/as_jet,m_hard_nqed);
//     std::cout<<" cancel As : "<<as_amegic/as_jet<<"^("<<m_hard_nqed<<" ) \n";
  }

//   std::cout<<" weight="<<m_weight<<"\n";

  p_ct = ct_tmp;

  p_events[njet-1]+=1;  // count events
  p_weight_sum[njet-1]+=m_weight;
  p_weight_sum_sqr[njet-1]+=sqr(m_weight);
}


int Cluster_Partons::SetDecayColours(Vec4D * p, Flavour * fl,int col1,int col2)
{
  int ncol   = 0;
  int nquark = 0;
  int ngluon = 0;
  for (int i=0; i<3; ++i) {
    if (fl[i].Strong()) {
      ++ncol;
      if (fl[i].IsQuark() || fl[i].IsSquark()) ++nquark;
      if (fl[i].IsGluon() || fl[i].IsGluino()) ++ngluon;
    }
  }  
  m_colors[0][0] = col1; m_colors[0][1] = col2; 
  for (int i=1; i<3; ++i) m_colors[i][0] = m_colors[i][1] = 0;
  m_scale = p[0].Abs2();
  m_asscale = m_scale;
  switch (ncol) {
  case 0: 
    // no colours at all.
    return 0;
  case 2:
    //3->31 8->81
    if (fl[0].Strong()) {
      for (short int i=1;i<3;i++) {
	if (fl[i].Strong()) {
	  for (short int j=0;j<2;j++) m_colors[i][j] = m_colors[0][j];
	  return 0;
	}
      }
    }
    // 1->33 1->88
    if (col1==0 && col2==0) {
      if ((fl[1].IsQuark()||fl[1].IsSquark()) && (fl[2].IsQuark()||fl[2].IsSquark())) {
	if (fl[1].IsAnti() && !(fl[2].IsAnti())) {
	  m_colors[1][1] = m_colors[2][0] = ATOOLS::Flow::Counter();
	  return 0;
	}
	if (fl[2].IsAnti() && !(fl[1].IsAnti())) {
	  m_colors[1][0] = m_colors[2][1] = ATOOLS::Flow::Counter();
	  return 0;
	}
      }
      if ((fl[1].IsGluon()||fl[1].IsGluino()) && (fl[2].IsGluon()||fl[2].IsGluino())) {
	m_colors[1][1] = m_colors[2][0] = ATOOLS::Flow::Counter();
	m_colors[1][0] = m_colors[2][1] = ATOOLS::Flow::Counter();
	return 0;
      }
    }
  case 3:
  default :
    msg.Error()<<"Error in Cluster_Partons::SetDecayColours:"<<std::endl
	       <<"   Cannot handle single color in 1 -> 2 process :"
               <<"   "<<fl[0]<<" "<<fl[1]<<" "<<fl[2]<<std::endl
	       <<"   Will abort the run."<<std::endl;
    abort();
  }
}

XS_Base * Cluster_Partons::GetXS(EXTRAXS::XS_Group * group, ATOOLS::Flavour * fl)
{
  XS_Base * xs = NULL;  
  const size_t nin=2, nout=2;
  const size_t n=nin+nout;

  Process_Base * proc=static_cast<Process_Base*>(p_me->GetAmegic()->GetProcess());
  p_lastproc=proc;

  if (group->XSSelector()->FindInGroup(group,xs,nin,nout,fl)==std::string::npos) {
    if (m_mode&1) {
      xs = group->XSSelector()->GetXS(nin,nout,fl,true);
    }
    else {
      int nstrong=0;
      for (size_t i=0;i<n;++i) {
	nstrong+=fl[i].Strong();
      }
      
      int nqed=0, nqcd=0;
      p_ct->AddCouplings(nqed,nqcd);

      //      std::cout<<" Cluster_Partons::GetXS nqed="<<nqed<<" nqcd="<<nqcd<<std::endl;

      Process_Base * proc=static_cast<Process_Base*>(p_me->GetAmegic()->GetProcess());
      //      std::cout<<proc->Name()<<"("<<proc->OrderEWeak()<<","<<proc->OrderStrong()<<")"<<std::endl;
      int nrqed = proc->OrderEWeak() - nqed;
      int nrqcd = proc->OrderStrong() - nqcd;

      xs = group->XSSelector()->GetXS(nin,nout,fl,false,nrqed,nrqcd);
    }
    if (xs) group->Add(xs);
  }
  p_xs=xs;
  return xs;
}


int Cluster_Partons::SetColours(XS_Base * xs, Vec4D * p, Flavour * fl)
{
  p_xs=xs;
  if (!xs)  {
    int test=SetColours(p,fl);
//     Process_Base * proc=static_cast<Process_Base*>(p_me->GetAmegic()->GetProcess());
//     std::cout<<" =Cluster============ \n";    
//     std::cout<<proc->Name()<<std::endl;
//     std::cout<<" ("<<m_hard_nqed<<","<<m_hard_nqcd<<")\n";
    return test;
  }

  if (!(m_mode&1) || xs->Size()==1) {
    m_hard_nqed = xs->OrderEW();
    m_hard_nqcd = xs->OrderStrong();
    // output calc remaining ew and strong order
//     Process_Base * proc=static_cast<Process_Base*>(p_me->GetAmegic()->GetProcess());
//     std::cout<<" ================== \n";    
//     std::cout<<proc->Name()<<std::endl;
//     std::cout<<xs->Name()<<std::endl;
//     std::cout<<" ("<<m_hard_nqed<<","<<m_hard_nqcd<<")\n";
    bool test=xs->SetColours(p);
    // check colour consistency 
    bool check=true;
    for (int i=0; i<4; ++i) {
      if (fl[i].IsQuark()) {
	if (fl[i].IsAnti() && (xs->Colours()[i][0]!=0 || xs->Colours()[i][1]==0)) check=false;
	if (!fl[i].IsAnti() && (xs->Colours()[i][0]==0 || xs->Colours()[i][1]!=0)) check=false;
      }
      if (fl[i].IsGluon()) {
	if (xs->Colours()[i][0]==0 || xs->Colours()[i][1]==0) check=false;
      }
      if (!check) break;
    }
    if (!check) {
      std::cout<<" colour check failed !!!\n";
      for (int i=0; i<4; ++i) {
	std::cout<<i<<" : "<<fl[i]<<" ("<<xs->Colours()[i][0]<<","<<xs->Colours()[i][1]<<") \n";
      }
    }
    return test;
  }
  else {
    // count existing ew and strong order
    int nqed=0, nqcd=0;
    int nout = p_ct->AddCouplings(nqed,nqcd);

    // calc remaining ew and strong order
    Process_Base * proc=static_cast<Process_Base*>(p_me->GetAmegic()->GetProcess());
    //    std::cout<<proc->Name()<<std::endl;
    int nrqed = proc->OrderEWeak() - nqed;
    int nrqcd = proc->OrderStrong() - nqcd;

    // select process (and colour structure)
    double xsecs[4] = {0.,0.,0.,0.};
    double s=(p[0]+p[1]).Abs2();
    double t=(p[0]-p[2]).Abs2();
    double u=(p[0]-p[3]).Abs2();

    int nxs=xs->Size();
    if (nxs>2) {
      msg.Error()<<" ERROR in Cluster_Partons::SetColours(XS_Base * xs, Vec4D * p, Flavour * fl) \n"
		 <<" ERROR more than two processes found"<<std::endl;
      for (int i=0; i<nxs; ++i) {
	msg.Error()<<(*xs)[i]->Name()<<std::endl;
      }
    }
    double sum = 0.;
    for (int i=0; i<nxs; ++i) {
      //      (*xs)[i]->Print();
      if ((int)(*xs)[i]->OrderEW()<=nrqed && (int)(*xs)[i]->OrderStrong()<=nrqcd) 
	sum+= (*xs)[i]->operator()(s,t,u);
      xsecs[i] = sum;
    }  

    double disc = sum*ran.Get();
    if (sum==0.) {
      msg.Error()<<" ERROR in Cluster_Partons::SetColours(XS_Base * xs, Vec4D * p, Flavour * fl) \n"
		 <<" ERROR no suitable process found"<<std::endl;
    }
    for (int i=0; i<nxs; ++i) {
      if (disc<=xsecs[i]) {
	m_hard_nqed = (*xs)[i]->OrderEW();
	m_hard_nqcd = (*xs)[i]->OrderStrong();
	p_xs = (*xs)[i];
	xs->SetSelected(p_xs);
        return p_xs->SetColours(s,t,u);
      }
    }
  }
  msg.Error()<<" ERROR in Cluster_Partons::SetColours(XS_Base * xs, Vec4D * p, Flavour * fl) \n";
  return 1;
}

int Cluster_Partons::SetColours(Vec4D * p, Flavour * fl)
{
  m_hard_nqed = 0;
  m_hard_nqcd = 0;

  // *** 2 -> 2 processes with unambiguous coulor structure
  // (a) no colors
  // (b) two quarks
  // (c) two quarks and one gluon
  // (d) two gluons (ADD-Model) 
  // (e) three gluons (ADD-Model)

  int ncol   = 0;
  int nquark = 0;
  int ngluon = 0;
  for (int i=0; i<4; ++i) {
    if (fl[i].Strong()) {
      ++ncol;
      if (fl[i].IsQuark() || fl[i].IsSquark()) ++nquark;
      if (fl[i].IsGluon() || fl[i].IsGluino()) ++ngluon;
    }
  }

  for (int i=0; i<4; ++i) m_colors[i][0]=m_colors[i][1]=0;
  double s = (p[0]+p[1]).Abs2();
  double t = (p[0]-p[2]).Abs2();
  double u = (p[0]-p[3]).Abs2();
  m_scale = s;
  m_asscale = m_scale;

  if (ncol==4) {
    msg.Out()<<"WARNING in Cluster_Partons::SetColours() : called for 4 coloured objects \n"
	     <<"          Don't know how to handle this ! "<<std::endl;
    for (int i=0; i<4; ++i) {
      msg.Out()<<i<<" : "<<fl[i]<<"\n";
    }
    return 1;
  }
  // (a) no colors
  if (ncol==0) {
    m_hard_nqed = 2;
    return 0;
  }

  
  int cols[3]={0,0,0};

  // (b) two quarks
  if (ncol==2 && nquark==2) {
    cols[0]=cols[1]=500;
    m_hard_nqed = 2;
  }
  // scale assuming two in or two out (cf "test" below)

  // (c) two quarks and one gluon (and one (massive) boson)
  if (ncol==3 && nquark==2 && ngluon==1) {
    m_hard_nqed = 1;
    m_hard_nqcd = 1;
    cols[0]=500;
    cols[1]=500+1;
    // in case one leg is massive we add all m^2 to the pt:
    // naive:    m_scale = (2.*s*t*u)/(s*s+t*t+u*u);
    // massive into account:
    //    m_asscale = 0.5*(sqr(p[2][1])+sqr(p[2][2])+sqr(p[3][1])+sqr(p[3][2]));
    m_asscale = 0.;
    if (fl[2].Strong()) m_asscale+=sqr(p[2][1])+sqr(p[2][2]);
    if (fl[3].Strong()) m_asscale+=sqr(p[3][1])+sqr(p[3][2]);

    // --- better determine pt in the right system  ---

    Vec4D q[4];
    for (int i=0;i<4;++i) q[i]=p[i];
    Poincare cms(q[0]+q[1]);
    for (int i=0;i<4;++i)    cms.Boost(q[i]);
    Poincare rot(q[0],Vec4D::ZVEC);
    for (int i=0;i<4;++i)      rot.Rotate(q[i]);
    //    m_asscale = 0.5*(sqr(q[2][1])+sqr(q[2][2])+sqr(q[3][1])+sqr(q[3][2]));

    // --- end ---
    m_scale = 0.5*(sqr(q[2][1])+sqr(q[2][2])+sqr(q[3][1])+sqr(q[3][2]));
    m_scale += p[2].Abs2()+p[3].Abs2();
    //    m_asscale = Min(dabs(t),dabs(u));
    
  }
  
  // (d) two gluons (ADD)
  if (ncol==2 && nquark==0 && ngluon==2) {
    cols[0]=500;
    cols[1]=500+1;
  }

  // (e) three gluons (ADD-Model)
  if (ncol==3 && nquark==0 && ngluon==3) {
    cols[0]=500;
    cols[1]=500+1;
    cols[2]=500+2;
  }



  ncol=0;
  int antis[3]={0,1,0};
  if (ngluon==3) antis[2]=2;
  int test = 3;
  for (int i=0; i<4; ++i) {
    if (fl[i].IsQuark() || fl[i].IsSquark()) {
      if (i<2) {
	antis[ncol]=fl[i].IsAnti();
	test=test && 1;
      }
      else {
	antis[ncol]=!fl[i].IsAnti();
	test=test && 2;
      }
      m_colors[i][fl[i].IsAnti()]=cols[ncol];
      ++ncol;
    }
  }
  if (!test && ncol==2) {
    // if colour flow = 1 particle in IS -> 1 in FS : m_scale = t or u.
    // naive : m_scale = (2.*s*t*u)/(s*s+t*t+u*u);
    //         mass effects ?!
    if (fl[0].Strong()) {
      if (fl[2].Strong())      m_asscale = m_scale = dabs(t);
      else if (fl[3].Strong()) m_asscale = m_scale = dabs(u);
    }
    else if (fl[1].Strong()) {
      if (fl[2].Strong())      m_asscale = m_scale = dabs(u);
      else if (fl[3].Strong()) m_asscale = m_scale = dabs(t);
    }
  }

  if (ngluon>0) {
    for (int i=0; i<4; ++i) {
      if (fl[i].IsGluon() || fl[i].IsGluino()) {
	if (i<2) {
	  m_colors[i][0] = cols[antis[1]];
	  m_colors[i][1] = cols[antis[0]];
	}
	else {
	  m_colors[i][0]=cols[antis[0]];
	  m_colors[i][1]=cols[antis[1]];
	}
	if (ngluon>=2 && nquark==0) {
	  int a=antis[0];
	  antis[0]=antis[1];
	  antis[1]=antis[2];
	  antis[2]=a;
	}
      }
    }
  }
  return 0;
}

void Cluster_Partons::FillDecayTree(Tree * fin_tree,XS_Base * xs)
{
  if (!fin_tree && m_fsrshoweron) {
    msg.Error()<<"ERROR in Cluster_Partons::FillDecayTrees: no trees!"
	       <<" No shower to be performed! "<<std::endl;
    return;
  } 
  if (p_blob->NInP()!=1 || p_blob->NOutP()!=2) {
    msg.Error()<<"ERROR in Cluster_Partons::FillDecayTrees: "<<std::endl
	       <<"   wrong number of articles in blob."<<std::endl;
    return;
  }

  if (!m_fsrshoweron) {
    // prepare dummy tree
    if (!p_local_tree)  p_local_tree=new Tree();
    p_local_tree->Reset();
    fin_tree = p_local_tree;
  }
  fin_tree->Reset();
  int flow = ATOOLS::Flow::Counter(), i=0;

  std::vector<Knot *> knots;
  knots.reserve(3);

  Knot * knot  = fin_tree->NewKnot();
  *knot->part  = *p_blob->InParticle(0);
  knot->part->SetStatus(2);
  knot->part->SetInfo('h');
  knot->stat   = 0;
  knot->z      = p_blob->OutParticle(0)->Momentum()[0]/knot->part->Momentum()[0];
  knot->E2     = sqr(knot->part->Momentum()[0]);
  knot->t      = knot->part->Momentum().Abs2();
  knot->thcrit = M_PI;
  knot->tout   = knot->t;
  knots.push_back(knot);
  if (knot->part->Flav().Strong()) {
    for (int j=0;j<2;j++) {
      if (xs && (xs->Colours()[i][j]!=0)) {	
	p_blob->InParticle(0)->SetFlow(j+1,flow+xs->Colours()[i][j]); 
      }
    }
  }
  i++;

  knot         = fin_tree->NewKnot(p_blob->OutParticle(0));
  knot->stat   = 3;
  knot->E2     = knots[0]->t;
  knot->t      = knots[0]->t;
  knot->thcrit = M_PI;
  if (knot->part->DecayBlob()) knot->tout = knot->part->Momentum().Abs2();
                          else knot->tout = Max(knot->part->Momentum().Abs2(),
						sqr(knot->part->Flav().Mass()));
  knots.push_back(knot);
  if (knot->part->Flav().Strong()) {
    for (int j=0;j<2;j++) {
      if (xs && (xs->Colours()[i][j]!=0)) {	
	knot->part->SetFlow(j+1,flow+xs->Colours()[i][j]); 
	p_blob->OutParticle(i-1)->SetFlow(j+1,flow+xs->Colours()[i][j]); 
      }
      else {
	if (m_colors[i][j]!=0) {	
	  knot->part->SetFlow(j+1,flow+m_colors[i][j]);
	  p_blob->OutParticle(i-1)->SetFlow(j+1,flow+m_colors[i][j]);
	}
      }
    }
  }
  i++;

  knot         = fin_tree->NewKnot(p_blob->OutParticle(1));
  knot->stat   = 3;
  knot->E2     = knots[0]->t;
  knot->t      = knots[0]->t;
  knot->thcrit = M_PI;
  if (knot->part->DecayBlob()) knot->tout = knot->part->Momentum().Abs2();
                          else knot->tout = Max(knot->part->Momentum().Abs2(),
						sqr(knot->part->Flav().Mass()));
  knots.push_back(knot);
  if (knot->part->Flav().Strong()) {
    for (int j=0;j<2;j++) {
      if (xs && (xs->Colours()[i][j]!=0)) {	
	knot->part->SetFlow(j+1,flow+xs->Colours()[i][j]); 
	p_blob->OutParticle(i-1)->SetFlow(j+1,flow+xs->Colours()[i][j]); 
      }
      else {
	if (m_colors[i][j]!=0) {	
	  knot->part->SetFlow(j+1,flow+m_colors[i][j]);
	  p_blob->OutParticle(i-1)->SetFlow(j+1,flow+m_colors[i][j]);
	}
      }
    }
  }
  i++;

  EstablishRelations(knots[0],knots[1],knots[2],1);
}


void Cluster_Partons::FixJetvetoPt2(double & jetvetopt2)
{
  double pt2min=p_ct->MinKt();
  if (pt2min!=0.) jetvetopt2=pt2min;

  /*
  //  std::cout<<" FixJetvetoPt2("<<jetvetopt2<<")\n"<<std::endl;
  Combine_Table * ct_tmp = p_ct;

  // anti-cluster and ...
  while (p_ct->Up()) {
    Combine_Table * ct_down = p_ct;
    p_ct=p_ct->Up();
    int i,j;
    double ptij = p_ct->GetWinner(i,j);
    if (m_sud_mode%10>0 && ct_down->GetLeg(i)->fl.Strong()  
	&& p_ct->GetLeg(i)->fl.Strong() && p_ct->GetLeg(j)->fl.Strong()) {  
      if (sqr(ptij)<jetvetopt2) jetvetopt2=sqr(ptij);
    }
  }

  p_ct = ct_tmp;
  */
  //  std::cout<<"  Cluster_Partons::FixJetvetoPt2() jetveto = "<<jetvetopt2<<")\n"<<std::endl;
}


void Cluster_Partons::FillTrees(Tree ** ini_trees,Tree * fin_tree)
{
  if ((!ini_trees && m_isrshoweron) || (!fin_tree && m_fsrshoweron)) {
    msg.Error()<<"ERROR in Cluster_Partons::FillTrees:"<<std::endl
	       <<"   No trees! no shower to be performed! "<<std::endl;
    return;
  } 
  

  if (!m_isrshoweron) {
    // prepare dummy tree
    if (!p_local_tree)  p_local_tree=new Tree();
    p_local_tree->Reset();
  }
  if (!m_fsrshoweron) {
    // prepare dummy tree
    if (!p_local_tree)  p_local_tree=new Tree();
    p_local_tree->Reset();
    fin_tree=p_local_tree;
  }

  int n[2]={0,1};
  if (p_me->InSwaped()) {
    n[0]=1;
    n[1]=0;
  }

  std::vector<Knot *> knots;
  std::vector<Knot *> ini_knots; // production points
  knots.reserve(10);
  ini_knots.reserve(10);

  // count jets
  int njet=p_ct->NLegs()-2;
  Combine_Table * ct_test = p_ct;
  while (ct_test->Up()) {
    ct_test=ct_test->Up();
    ++njet;
  }

  // generate knotlist from pointlist in Combine_Table

  // start initial state
  if (m_isrshoweron) {
    knots.push_back(Point2Knot(ini_trees[n[0]],p_ct->GetLeg(0), p_ct->Momentum(0),'G'));
    knots.push_back(Point2Knot(ini_trees[n[1]],p_ct->GetLeg(1), p_ct->Momentum(1),'G'));
  }
  else {
    knots.push_back(Point2Knot(p_local_tree,p_ct->GetLeg(0), p_ct->Momentum(0),'G'));
    knots.push_back(Point2Knot(p_local_tree,p_ct->GetLeg(1), p_ct->Momentum(1),'G'));
  }
  
  Knot * mo = 0;   

  mo   = fin_tree->NewKnot();
  knots.push_back(Point2Knot(fin_tree    ,p_ct->GetLeg(2), p_ct->Momentum(2),'H'));
  knots.push_back(Point2Knot(fin_tree    ,p_ct->GetLeg(3), p_ct->Momentum(3),'H'));

  knots[0]->part->SetDecayBlob(p_blob);  
  knots[1]->part->SetDecayBlob(p_blob);
  knots[2]->part->SetProductionBlob(p_blob);
  knots[3]->part->SetProductionBlob(p_blob);

  for (int i=0;i<4;i++) {
    for (int j=0;j<2;j++) {
      if (p_xs) {
	knots[i]->part->SetFlow(j+1,p_xs->Colours()[i][j]); 
      }
      else {
	knots[i]->part->SetFlow(j+1,m_colors[i][j]);
      }	
    }
  }

  assert(mo);
  Vec4D sum=p_ct->Momentum(2)+p_ct->Momentum(3);
  m_cms_boost=Poincare(sum);
  Vec4D p1   =m_cms_boost*sum;
  Vec4D p2   =m_cms_boost*p_ct->Momentum(2);
  Vec4D p3   =m_cms_boost*p_ct->Momentum(3);

  *(mo->part) = Particle(0,Flavour(kf::none),sum);
  mo->part->SetInfo('M');
  mo->part->SetStatus(2);
  mo->stat   = 0;
  mo->z      = p2[0]/p1[0];
  mo->E2     = sqr(p1[0]);
  mo->thcrit = M_PI;

  //we have a virtuality ordered shower, therefore:
  mo->t = mo->part->Momentum().Abs2();

  double scale = Scale();
   if (p_xs) scale = p_xs->Selected()->Scale(PHASIC::stp::as);

   scale = 4.*(m_qmin_i*m_qmin_i);
   scale = Max(scale,mo->t);

   double s=sum.Abs2();
   double t=dabs((p_ct->Momentum(0)-p_ct->Momentum(2)).Abs2());
   double u=dabs((p_ct->Momentum(1)-p_ct->Momentum(2)).Abs2());
   if (p_lastproc->OrderEWeak()==0) scale=2.0*s*t*u/(s*s+t*t+u*u);

//   std::cout<<" ew:"<<p_me->OrderEWeak()
//   	     <<" strong:"<<p_me->OrderStrong()<<"\n";
  if(p_me->OrderStrong()==0 && njet==m_maxjetnumber) scale=mo->t;
//  std::cout<<"scale="<<scale<<std::endl;
  EstablishRelations(mo,knots[0],knots[1],0,scale);
  EstablishRelations(mo,knots[2],knots[3],1);      

  // determine starting conditions for showers
  // note, that starting conditions for subsequent branches have to be 
  // evaluted during the shower evoultion (since the system, esp. for 
  // final state showers starting from the initial state shower are not
  // known.)
  DetermineColourAngles(knots);

  for (int l=0; l<4; ++l) ini_knots.push_back(mo);

  Tree * tree;
  int nlegs=4;
  int i,j;

  Combine_Table * ct_tmp = p_ct;

  p_ct=p_ct->Up();
  while (p_ct) {
    knots.push_back(0);ini_knots.push_back(0); ++nlegs;
    p_ct->GetWinner(i,j);
    for (int l=knots.size()-1;l>j;--l) knots[l] = knots[l-1];
    for (int l=ini_knots.size()-1;l>j;--l) ini_knots[l] = ini_knots[l-1];
    if (i>=2) tree = fin_tree; 
    else      tree = ini_trees[n[i]];

    Knot * d1, * d2;
    if (i>=2) {
      d1 = Point2Knot(tree, p_ct->GetLeg(i), p_ct->Momentum(i),'H');
      d2 = Point2Knot(tree, p_ct->GetLeg(j), p_ct->Momentum(j),'H');
      d1->part->SetProductionBlob(p_blob);
      d2->part->SetProductionBlob(p_blob);

      EstablishRelations(knots[i],d1,d2,1);      
    } 
    else {
      d1 = Point2Knot(tree, p_ct->GetLeg(i), p_ct->Momentum(i),'H');
      d2 = Point2Knot(tree, p_ct->GetLeg(j), p_ct->Momentum(j),'H');
      d1->part->SetDecayBlob(p_blob);  
      d2->part->SetDecayBlob(p_blob);

      EstablishRelations(d1,knots[i],d2,2+i);      
    }


    knots[i] = d1;
    knots[j] = d2;

    p_ct = p_ct->Up();
  }
          
  p_ct = ct_tmp;


  // update colours in blob
  for (int i=0; i<4 ; ++i) {
    int j=i/2;
    int k=i%2+1;
    if (p_blob->InParticle(j)->Flav()==knots[j]->part->Flav()) {
      p_blob->InParticle(j)->SetFlow(k, knots[j]->part->GetFlow(k));
    }
    else {
      p_blob->InParticle(1-j)->SetFlow(k, knots[j]->part->GetFlow(k));
    }
  }
  
  for (int i=4; i<2*nlegs ; ++i) {
    int j=i/2;
    int k=i%2+1;
    p_blob->OutParticle(j-2)->SetFlow(k, knots[j]->part->GetFlow(k));
  }


  if (msg.LevelIsDebugging()) {
    msg.Out()<<" in Cluster_Partons::FillTrees("<<m_isrshoweron<<","
	     <<m_fsrshoweron<<")"<<std::endl;
    if (ini_trees) {
      msg.Out()<<"initree[0]:"<<std::endl<<ini_trees[0]
	       <<"initree[1]:"<<std::endl<<ini_trees[1];
    }
    msg.Out()<<"fin_tree:"<<std::endl<<fin_tree
	     <<"****************************************"<<std::endl;
  }
}



Knot * Cluster_Partons::Point2Knot(Tree * tree, const Leg & po, 
				   const Vec4D & mom, char info) 
{
  Flavour flav(po->fl);
  // check in map
  Flavour_Map::const_iterator cit=m_flmap.find(flav);
  if (cit!=m_flmap.end()) flav=cit->second;

  if (po.ExtraAnti() == -1) flav = flav.Bar();

  Knot * k =0 ;

  bool found = 0;
  for (int i=0;i<p_blob->NInP();i++) {
    if ( (p_blob->InParticle(i)->Flav() == flav) &&
	 (p_blob->InParticle(i)->Momentum() == mom) ) { 
      k = tree->NewKnot(p_blob->InParticle(i));
      found = 1;
      break;
    }
  }
  if (!found) {
    for (int i=0;i<p_blob->NOutP();i++) {
      if ( (p_blob->OutParticle(i)->Flav() == flav) &&
	   (p_blob->OutParticle(i)->Momentum() == mom) ) {
	k = tree->NewKnot(p_blob->OutParticle(i));
	found = 1;
	break;
      }
    }
  }
  if (!found) {
    k = tree->NewKnot();
    *k->part = Particle(0,flav,mom);
  }


  // preliminary parton status!!!
  k->part->SetInfo(info);
  k->part->SetStatus(1);  //final
  k->tout      = sqr(flav.PSMass());
  if (flav.IsKK() || k->part->DecayBlob()) k->tout=mom.Abs2();
  k->E2        = sqr(mom[0]);
  k->costh     = 0; 
  k->stat      = 3;
  k->thcrit    = M_PI;

  return k;
}

void Cluster_Partons::EstablishRelations(Knot * mo, Knot * d1,Knot * d2,int mode,double scale)
{
  if (mode==1) {
    Vec4D p1   =m_cms_boost*mo->part->Momentum();
    Vec4D p2   =m_cms_boost*d1->part->Momentum();
    Vec4D p3   =m_cms_boost*d2->part->Momentum();

    mo->left  = d1;
    mo->right = d2;
    mo->z     = p2[0]/p1[0];
    mo->stat  = 0;
    mo->part->SetStatus(2);
    if (mo->part->Info() != 'H') mo->part->SetInfo('f');

    d1->prev  = mo;
    d2->prev  = mo;
    d1->E2    = sqr(p2[0]);
    d2->E2    = sqr(p3[0]);

    APACIC::Final_State_Shower::EstablishRelations(mo,d1,d2);

    mo->tout = mo->t;

    return;
  }
  else if (mode==0) {
    // initial state initialization
    //  status:
    //  p_blob->CMS()          - Vec4D hard event in LAB system
    //  d1->part->Momentum() - in the moment also in LAB system
    //  p_blob->InParticle(0)->Momentum() - in CMS system

    // set x1 and x2
    double x1,x2;
    p_ct->GetX1X2(x1,x2);
    d1->x=x1;
    d2->x=x2;

    // set start t
    d1->t = -scale;
    d2->t = -scale;

    // angle condition set via DetermineColourAngles called in FillTrees
  }
  else if (mode==2 || mode==3) {
    // initial state initialization
    //     mo -> d1 (IS) 
    //        -> d2 (FS)

    if (!mo || !d1 || !d2 ) {
      msg.Error()<<"ERROR: can not establish relations with less than three elements "<<std::endl;
    }
    mo->right=d1;
    mo->left =d2;
    d1->prev=mo;
    d2->prev=mo;

    double t1 = d1->part->Momentum().Abs2();
    double t0 = d1->t;
    if (1) {
      mo->t = t0;
      d1->t = t1;
      d1->tmax = t1;
      d2->t = -t1;
    }
    else {
      mo->t = t1;
      d1->t = t1;
      d1->tmax = t1;
      d2->t = -t0;
    }
    double x1,x2;
    p_ct->GetX1X2(x1,x2);
    if (mode==2) {
      mo->x = x1;
    }
    else {
      mo->x = x2;
    }
    d1->stat = 0;

    APACIC::Initial_State_Shower::SetColours(d1);
  }
}

Flavour Cluster_Partons::Flav(int i) {
  if (p_ct) {
    Flavour fl(p_ct->Flav(i));
    Flavour_Map::const_iterator cit=m_flmap.find(fl);
    if (cit!=m_flmap.end()) fl=cit->second;
    return fl;
  }
  msg.Error()<<"ERROR in Cluster_Partons::Flav. No ct."<<std::endl;
  return 0;
}

Vec4D Cluster_Partons::Momentum(int i) {
  if (p_ct) return p_ct->Momentum(i);
  msg.Error()<<"ERROR in Cluster_Partons::Momentum. No ct."<<std::endl;
  return Vec4D(0.,0.,0.,0.);
}

Vec4D Cluster_Partons::Momentum(Knot * mo, int & number) {
  if (mo->left) return Momentum(mo->left,number) + Momentum(mo->right,number);
  number++;
  if ((mo->part->Info()!='G'))  // &&(mo->part->info!='I'))
    return mo->part->Momentum();

  return Vec4D(0.,0.,0.,0.);
}


bool Cluster_Partons::IsColourConnected(Particle * a, Particle * b) {
  return (( (a->GetFlow(1)!=0) && ( (a->GetFlow(1)==b->GetFlow(1)) || 
				 (a->GetFlow(1)==b->GetFlow(2)))  ) ||
	  ( (a->GetFlow(2)!=0) && ( (a->GetFlow(2)==b->GetFlow(2)) ||
				 (a->GetFlow(2)==b->GetFlow(1)))  )    );
}

double Cluster_Partons::ColourAngle(const std::vector<Knot *> & knots, const int i) {
  if (!((knots[i]->part->Flav()).Strong())) return M_PI;

  
  double x1=1.;
  double x2=1.;
  Vec4D sum(0.,0.,0.,0.);
  for (size_t l=2;l<knots.size();++l)
    sum += knots[l]->part->Momentum();

  double sprime = sum.Abs2();

  x1=knots[0]->x;
  x2=knots[1]->x;
  sprime = (knots[0]->part->Momentum() + knots[1]->part->Momentum()).Abs2();
  Poincare lab(Vec4D(x1+x2,0.,0.,-(x1-x2)));

  double angle = 0.;
  int start = 0;
  for (int j=start;j<(int)knots.size();++j) {
    if (j!=i) {
      if (IsColourConnected(knots[i]->part,knots[j]->part)) {
	Vec3D ivec, jvec;
	Vec4D i4vec, j4vec;
	double th_crude=0.;
	if (i<2) {
	  // determine isr angles in labframe
	  i4vec=lab*knots[i]->part->Momentum();
	  j4vec=lab*knots[j]->part->Momentum();
	  ivec = Vec3D(i4vec);
	  jvec = Vec3D(j4vec);
	}
	else {

	  i4vec=knots[i]->part->Momentum();
	  j4vec=knots[j]->part->Momentum();
	  ivec = Vec3D(i4vec);
	  jvec = Vec3D(j4vec);
	  
	}
	Vec4D mvec   = i4vec+j4vec;
	double t_mo = mvec.Abs2();
	double E_mo = mvec[0];
	double z = j4vec[0]/E_mo;
	th_crude  = sqrt( t_mo/(z*(1.- z)))/E_mo;
	double test=ivec*jvec/(ivec.Abs()*jvec.Abs());
	double th_ex=M_PI;
	if (test<=-1.) {
	  test = -1.;
	}
	else if (test>=1.) {
	  th_ex=0;
	  test = 1.;
	}
	else {
	  th_ex = acos(test);
	}

	angle = Max(angle,th_ex);

	if (IsEqual(angle,M_PI)) {
	  angle=M_PI;
	}
      }
    }
  }

  return angle;
}


void Cluster_Partons::DetermineColourAngles(const std::vector<APACIC::Knot *> & knots) {
  int n=knots.size();
  // save momenta
  Vec4D * moms = new Vec4D[n];
  for (int i=0;i<n;++i) moms[i] = knots[i]->part->Momentum();
  
  Poincare cms(moms[0]+moms[1]);
  for (int i=0;i<n;++i) {
    knots[i]->part->SetMomentum(cms*knots[i]->part->Momentum());
  }
  Poincare zaxis(knots[0]->part->Momentum(),Vec4D::ZVEC);
  for (int i=0;i<n;++i) {
    knots[i]->part->SetMomentum(zaxis*knots[i]->part->Momentum());
  }
  
  for (size_t i=0;i<knots.size();++i) {
    double th= ColourAngle(knots,i);
    knots[i]->thcrit=th;
  }

  // restore momenta
  for (int i=0;i<n;++i) {
    knots[i]->part->SetMomentum(moms[i]);
  }
  delete [] moms;
}


void Cluster_Partons::CreateFlavourMap() {
  //  if (p_me->GetAmegic()->GetProcess()!=p_me->GetAmegic()->GetProcess()->Partner()) {

  // map needed if :
  //    process != partner process
  //    if HHMF is turned on (me_handler->Flavour =! proc->Flavour)
  //
  
  Process_Base * proc=static_cast<Process_Base*>(p_me->GetAmegic()->GetProcess());
  Process_Base * partner=proc->Partner();
  m_nstrong   = proc->NStrong();
  const Flavour * flavs=proc->Flavours();
  double ycut= proc->Ycut();
  if (ycut!=-1.) {
//     std::cout<<" selector ycut= "<<ycut<<" \n"<<std::endl;
    m_ycut=ycut;
  }
  else {
    m_ycut=ATOOLS::rpa.gen.Ycut();
  }

  const Flavour * partner_flavs=partner->Flavours();

  int n[2]={0,1};
  //    if (proc->InSwaped()^partner->InSwaped()) {
  if (p_me->InSwaped()^partner->InSwaped()) {
    n[0]=1;
    n[1]=0;
  }

  // create new map
  m_flmap.clear();
  for (size_t i=0;i<proc->NIn();++i) {
    if (partner_flavs[i]!=flavs[n[i]]) {
      m_flmap[partner_flavs[i]]=flavs[n[i]];
      if (partner_flavs[i]!=(Flavour(partner_flavs[i])).Bar()) {
	m_flmap[(Flavour(partner_flavs[i])).Bar()]=(Flavour(flavs[n[i]])).Bar();
      }
    }
  }
  for (size_t i=proc->NIn();i<proc->NIn()+proc->NOut();++i) {
    if (partner_flavs[i]!=flavs[i]) {
      m_flmap[partner_flavs[i]]=flavs[i];
      if (partner_flavs[i]!=(Flavour(partner_flavs[i])).Bar()) {
	m_flmap[(Flavour(partner_flavs[i])).Bar()]=(Flavour(flavs[i])).Bar();	  
      }
    }
  }

}

void   Cluster_Partons::JetvetoPt2(double & q2i, double & q2f) 
{ 
  q2i = sqr(m_qmin_i);
  q2f = sqr(m_qmin_f); 
}
