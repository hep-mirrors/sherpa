#include "Cluster_Partons_Base.H"
#include "NLL_Branching_Probabilities.H"
#include "Flow.H"
#include "Poincare.H"
#include "Message.H"
#include "Process_Base.H"
#include "Splitting_Function.H"
#include "MyStrStream.H"
#include <iomanip>

using namespace SHERPA;
using namespace EXTRAXS;
using namespace AMEGIC;
using namespace APACIC;
using namespace ATOOLS;
using namespace MODEL;

#define CA 3.0
#define TR 0.5
#define Nf 5.0

Cluster_Partons_Base::Cluster_Partons_Base(Matrix_Element_Handler * me,ATOOLS::Jet_Finder * jf,
					   const int maxjet, const int isrmode,
					   const int isron,const int fsron) :
  p_me(me), p_runas(NULL), p_jf(jf), p_fssud(NULL), p_issud(NULL), p_ct(NULL), p_combi(NULL), 
  m_njet(0), m_maxjetnumber(maxjet),m_isrmode(isrmode), m_isrshoweron(isron),m_fsrshoweron(fsron), 
  m_sud_mode(0), m_kfac(0.), m_counts(0.), m_fails(0.)
{
  // read in some parameters
  Data_Read dr(rpa.GetPath()+
	       rpa.gen.Variable("SHOWER_DATA_FILE","Shower.dat"));
  m_bp_mode  = dr.GetValue<int>("SUDAKOV_TYPE",40);
  m_sud_mode = dr.GetValue<int>("CKKW_SUDAKOV_MODE",3);
  m_as_mode  = dr.GetValue<int>("CKKW_ALPHAS_MODE",1);
  if (ToType<int>(rpa.gen.Variable("S_KFACTOR_SCHEME","1"))&2) {
    /*
      in principle we need the k factor also in the sudakovs
      however, lep data are then not reproduced, since the 
      parton shower produces different jet rates 
      -> switch it off for e+e-
    */
    if (p_jf->Type()>1) m_bp_mode=m_bp_mode|bpm::soft_kfac;
    m_kfac=CA*(67./18.-M_PI*M_PI/6.)-10./9.*TR*Nf;
  }
  else {
    if (m_bp_mode&bpm::soft_kfac) m_bp_mode-=bpm::soft_kfac; 
  }
  m_is_as_factor=ToType<double>(rpa.gen.Variable("IS_CPL_SCALE_FACTOR","1"));
  m_fs_as_factor=ToType<double>(rpa.gen.Variable("FS_CPL_SCALE_FACTOR","1"));
  m_me_as_factor=p_jf->Type()<2?0.25:1.0;
  msg_Tracking()<<"Initalize Cluster_Partons_Base with {\n"
		<<"   Sudakov type            = "<<m_bp_mode<<"\n"
		<<"   ren. scale factor       = "<<rpa.gen.RenormalizationScaleFactor()<<"\n" 
		<<"   is PS ren. scale factor = "<<m_is_as_factor<<"\n"
		<<"   fs PS ren. scale factor = "<<m_fs_as_factor<<"\n"
		<<"   ME ren. scale factor    = "<<m_me_as_factor<<"\n"
		<<"   K factor                = "<<m_kfac<<"\n}"<<std::endl;
  p_runas = MODEL::as; 
  if (m_fsrshoweron!=0 && (m_sud_mode&1)) {
    p_fssud = new NLL_Sudakov((bpm::code)(m_bp_mode+1),
			      p_jf->Smax(),p_jf->Smin(),
			      p_runas,m_fs_as_factor*m_me_as_factor);
  }
  if (m_isrshoweron!=0 && (m_sud_mode&2)) {
    p_issud = new NLL_Sudakov((bpm::code)(m_bp_mode+2),
			      p_jf->Smax(),p_jf->Smin(),
			      p_runas,m_is_as_factor*m_me_as_factor);
  }
  m_events.resize(m_maxjetnumber,0);
  m_meweight_sum.resize(m_maxjetnumber,0.0);
  m_weight_sum.resize(m_maxjetnumber,0.0);
  m_weight_sum_sqr.resize(m_maxjetnumber,0.0); 
  m_asweight_sum.resize(m_maxjetnumber,0.0);
  m_asweight_sum_sqr.resize(m_maxjetnumber,0.0); 
  m_sweight_sum.resize(m_maxjetnumber,0.0);
  m_sweight_sum_sqr.resize(m_maxjetnumber,0.0); 

  double deltar = rpa.gen.DeltaR();
  double ycut   = rpa.gen.Ycut();
  m_qmin_i = sqrt(ycut)*rpa.gen.Ecms();
  m_qmin_f = sqrt(ycut)*deltar*rpa.gen.Ecms();
  exh->AddTerminatorObject(this);
}

Cluster_Partons_Base::~Cluster_Partons_Base()
{
  if (p_combi)               { delete p_combi;      p_combi      = NULL; }
  if (p_fssud)               { delete p_fssud;    p_fssud        = NULL; }
  
  if (m_counts!=0) WriteOutWeights();
  exh->RemoveTerminatorObject(this);
}

void Cluster_Partons_Base::PrepareTerminate()
{
  if (m_counts!=0) WriteOutWeights();
}

void Cluster_Partons_Base::WriteOutWeights() 
{
  msg_Info()<<METHOD<<"(): Weight statistics {"<<std::endl;
  msg_Info()<<"  Misclusterings: "<<m_fails/m_counts<<std::endl;
  for (int i=0;i<m_maxjetnumber;++i) {
    if (m_events[i]==0) continue;
    double w_mean(m_weight_sum[i]/m_meweight_sum[i]);
    double w_sigma(sqrt((m_weight_sum_sqr[i]/m_meweight_sum[i]
			 -sqr(w_mean))/(m_events[i]-1.0)));
    double asw_mean(m_asweight_sum[i]/m_meweight_sum[i]);
    double asw_sigma(sqrt((m_asweight_sum_sqr[i]/m_meweight_sum[i]
			   -sqr(asw_mean))/(m_events[i]-1.0)));
    double sw_mean(m_sweight_sum[i]/m_meweight_sum[i]);
    double sw_sigma(sqrt((m_sweight_sum_sqr[i]/m_meweight_sum[i]
			  -sqr(sw_mean))/(m_events[i]-1.0)));
    msg_Info()<<"  <w>_{"<<(i+1)<<"-jet} = "<<std::setw(12)<<w_mean
	      <<" +- "<<std::setw(12)<<w_sigma<<" ( "<<std::setw(12)
	      <<(int(w_sigma*10000./w_mean))/100.<<" % )"<<std::endl;
//     msg_Info()<<"  <w>_{as,"<<(i+1)<<"-jet}     = "<<std::setw(12)<<asw_mean
// 	      <<" +- "<<std::setw(12)<<asw_sigma<<" ( "<<std::setw(12)
// 	      <<(int(asw_sigma*10000./asw_mean))/100.<<" % )"<<std::endl;
//     msg_Info()<<"  <w>_{\\Delta,"<<(i+1)<<"-jet} = "<<std::setw(12)<<sw_mean
// 	      <<" +- "<<std::setw(12)<<sw_sigma<<" ( "<<std::setw(12)
// 	      <<(int(sw_sigma*10000./sw_mean))/100.<<" % )"<<std::endl;
    m_weight_sum[i]=w_mean;
    m_weight_sum_sqr[i]=w_sigma;
    m_asweight_sum[i]=asw_mean;
    m_asweight_sum_sqr[i]=asw_sigma;
    m_sweight_sum[i]=sw_mean;
    m_sweight_sum_sqr[i]=sw_sigma;
  }
  msg_Info()<<"}"<<std::endl;
  MyStrStream sstr;
  sstr.precision(4);
  double ycut(log(rpa.gen.Ycut())/log(10.));
  if (rpa.gen.Beam1().IsHadron() && rpa.gen.Beam2().IsHadron())
    ycut=sqrt(rpa.gen.Ycut())*rpa.gen.Ecms();
  int mode(m_sud_mode&3+m_bp_mode);
  sstr<<"Weights_"<<rpa.gen.Ecms()<<"_"
      <<ycut<<"_"<<mode<<".dat"<<std::endl;
  std::string filename;
  sstr>>filename;
  std::ofstream  rfile(filename.c_str(),std::ios::app);
  rfile.precision(8);
  rfile<<"# n"<<std::setw(18)<<"<w>"<<std::setw(18)
       <<"\\Delta<w>"<<std::setw(18)<<"<w>_{as}"<<std::setw(18)
       <<"\\Delta<w>_{as}"<<std::setw(18)<<"<w>_{sud}"<<std::setw(18)
       <<"\\Delta<w>_{sud}"<<std::setw(18)<<"N"<<std::endl;
  for (int i=1;i<m_maxjetnumber;++i) {
    rfile<<std::setw(3)<<i+1<<std::setw(18)<<m_weight_sum[i]<<std::setw(18)
	 <<m_weight_sum_sqr[i]<<std::setw(18)<<m_asweight_sum[i]<<std::setw(18)
	 <<m_asweight_sum_sqr[i]<<std::setw(18)<<m_sweight_sum[i]<<std::setw(18)
	 <<m_sweight_sum_sqr[i]<<std::setw(18)<<m_events[i]<<std::endl;
  }  
}

Leg **Cluster_Partons_Base::CreateLegs(int &nampl,const int nlegs,const bool reuse)
{
  Leg **legs(NULL);
   // start cluster algorithm :
  if (!reuse) {
    if (p_combi) delete p_combi;
    p_combi = 0;
    // generate a list of "legs" for each amplitude
    legs = new Leg *[nampl];
    PHASIC::Integrable_Base *proc(p_me->GetAmegic()->GetProcess());
    for (int k=0;k<nampl;) {
      legs[k] = new Leg[nlegs];
      int l   = 0;
      if (FillLegs(legs[k],p_me->GetDiagram(k),l,nlegs)) ++k;
      else {
	delete [] legs[k];
	--nampl;
      }
    }
    for (int k=0;k<nampl;++k) {
      for (int i(0);i<nlegs;++i) {
	Flavour fl(proc->Flavours()[i]);
	if (i<2 && proc->InSwaped()) fl=proc->Flavours()[1-i];
	legs[k][i].SetMapFlavour(fl);
// 	msg_Debugging()<<"set mapfl: "<<k<<", "<<i<<": "<<fl<<" "
// 		       <<proc->InSwaped()<<"\n";
      }
    }
    msg_Debugging()<<"map process: "<<proc->Name()
		   <<" -> "<<static_cast<AMEGIC::Process_Base*>(proc)->
      Partner()->Name()<<"\n";
  }  
  return legs;
}

void Cluster_Partons_Base::CreateTables(Leg ** legs,const int nampl,
					const double x1,const double x2) 
{
  p_ct = 0;
  // if no combination table exist, create it
  int nin(p_blob->NInP()), nout(p_blob->NOutP()), nlegs(nin+nout);
  Vec4D * amoms = new Vec4D[nlegs];
  for (int i=0;i<nin;++i)  amoms[i]     = p_blob->InParticle(i)->Momentum();
  for (int i=0;i<nout;++i) amoms[nin+i] = p_blob->OutParticle(i)->Momentum();
  if (p_me->InSwaped()) {
    // avoid flavour mismatch if using amplitudes
    Vec4D help=amoms[0];
    amoms[0]=amoms[1];
    amoms[1]=help;
  }
  if (!p_combi) {
    /*
      - copy moms to insert into Combine_Table (will be delete there)
      - create new Combine_Table with given momenta and given Jet-measure
      - initialise Combine_Table
      - determine best combination sheme
    */ 
    p_combi = CreateTable(p_jf,amoms,0,m_isrmode,m_isrshoweron);
    p_combi->FillTable(legs,nlegs,nampl);   
    p_ct = p_combi->CalcJet(nlegs,x1,x2); 
  }
  else {
    // use the existing combination table and determine best combination sheme
    p_ct = p_combi->CalcJet(nlegs,x1,x2,amoms);
  }
  //  delete [] amoms;
}

bool Cluster_Partons_Base::ClusterConfiguration(Blob * blob,double x1,double x2) 
{
  p_blob=blob;
  int nin=p_blob->NInP();
  int nout=p_blob->NOutP();
  if (nin==1) {
    if (nout<4) return true;
    msg.Error()<<"Error in Cluster_Partons_Base::ClusterConfiguration()"<<std::endl
	       <<"   Try to cluster decay blob, nin ="<<nin<<" with nout = "<<nout<<","<<std::endl
	       <<"   No method provided yet. Return 0."<<std::endl;
    return false;
  }
  int nampl=p_me->NumberOfDiagrams();
  int nlegs=nin+nout;
  Leg **legs(CreateLegs(nampl,nlegs));
  CreateTables(legs,nampl,x1,x2);
  PHASIC::Integrable_Base *proc(p_me->GetAmegic()->GetProcess());
  m_ycut=proc->Ycut();
  msg_Debugging()<<METHOD<<"(): process = "<<proc->Name()<<", y_cut = "
		 <<m_ycut<<" ("<<sqrt(m_ycut)*rpa.gen.Ecms()<<")\n";
  m_q2_amegic=proc->Scale(PHASIC::stp::ren);
  m_qmin_i=sqrt(m_ycut)*rpa.gen.Ecms();
  m_qmin_f=rpa.gen.DeltaR()*m_qmin_i;
  m_q2_isjet=rpa.gen.Ycut()*sqr(rpa.gen.Ecms());
  m_q2_fsjet=sqr(rpa.gen.DeltaR())*m_q2_isjet;
  m_q2_iss=m_q2_fss=std::numeric_limits<double>::max();
  m_q2_hard=m_jv_pt2=std::numeric_limits<double>::max();
  return 1;
}

void Cluster_Partons_Base::SetQMin(Combine_Table_Base *const ct,
				   Combine_Table_Base *const ref,
				   size_t &decayjets,
				   size_t &hardjets,double &hardqmin)
{
  if (ct->Up()==ref || ct->Up()==NULL) {
    for (int l(2);l<ct->NLegs();++l) 
      if (ct->GetLeg(l).QMin()<0.0 && 
	  ct->GetLeg(l).Point()->fl.Strong()) ++hardjets;
    for (int l(0);l<p_ct->NLegs();++l) {
      if (p_ct->GetLeg(l).Point()->t<10) 
	hardqmin=Min(hardqmin,sqrt(p_ct->GetLeg(l).MinKT2QCD()));
    }
    msg_Debugging()<<"n_{hard jets} = "<<hardjets
		   <<", n_{max hard jets} = "<<m_maxqcdjets<<"\n";
    msg_Debugging()<<"n_{decay jets} = "<<decayjets
		   <<", n_{max decay jets} = "<<m_max_decayjets<<"\n";    
    msg_Debugging()<<"q_{min,hard} = "<<hardqmin<<"\n";    
    m_hardjets=hardjets;
    if (hardqmin==std::numeric_limits<double>::max()) {
      hardqmin=sqrt(p_ct->Kt2QCDHard());
      msg_Debugging()<<"reset q_{min,hard} -> "<<hardqmin<<"\n";    
    }
    if (sqr(m_qmin_i)<m_q2_isjet && !IsEqual(sqr(m_qmin_i),m_q2_isjet)) 
      // multi scale treatment
      for (int i(0);i<2;++i) 
	m_qmin[i]=Min(hardqmin,sqrt(m_q2_isjet));
    else if (hardjets==(size_t)m_maxqcdjets) 
      // highest mutiplicity treatment
      for (int i(0);i<2;++i) m_qmin[i]=hardqmin;
    else
      // standard ckkw
      for (int i(0);i<2;++i) m_qmin[i]=dabs(m_qmin[i]);
    if (sqr(m_qmin_f)<m_q2_fsjet && !IsEqual(sqr(m_qmin_f),m_q2_fsjet)) 
      // multi scale treatment
      for (int i(2);i<ct->NLegs();++i) {
	if (ct->GetLeg(i).QMin()<0.0)
	  ct->GetLeg(i).SetQMin(Min(hardqmin,sqrt(m_q2_fsjet)));
      }
    else if (hardjets==(size_t)m_maxqcdjets)
      // highest mutiplicity treatment
      for (int i(2);i<ct->NLegs();++i) {
	if (ct->GetLeg(i).QMin()<0.0)
	  ct->GetLeg(i).SetQMin(hardqmin);
      }
    else
      // standard ckkw
      for (int i(2);i<ct->NLegs();++i) {
	if (ct->GetLeg(i).QMin()<0.0)
	  ct->GetLeg(i).SetQMin(dabs(ct->GetLeg(i).QMin()));
      }
    return;
  }
  int i, j;
  double ptiji(ct->Up()->GetWinner(i,j)), ptijj(ptiji);
  ptijj=ptiji=sqrt(ct->GetLeg(i).KT2QCD());
  for (int l(0);l<j;++l)
    ct->Up()->GetLeg(l).SetPQMin(ct->GetLeg(l).PQMin());
  for (int l(j);l<ct->NLegs();++l)
    ct->Up()->GetLeg(l+1).SetPQMin(ct->GetLeg(l).PQMin());
  ct->Up()->GetLeg(j).SetPQMin(ct->GetLeg(i).PQMin());
  if (ct->GetLeg(i).Point()->t<10) {
    ptijj=ptiji=Min(ptiji,ct->GetLeg(i).QMin());
    if (i<2) {
      m_qmin.push_back(0.0);
      if (ptijj<0.0) ptijj=-m_qmin_f;
      ct->Up()->GetLeg(j).SetPQMin(&m_qmin.back());
    }
  }
  else {
    m_qmin.push_back(0.0);
    ct->Up()->GetLeg(i).SetPQMin(&m_qmin.back());
    if (i<2) m_qmin.push_back(0.0);
    ct->Up()->GetLeg(j).SetPQMin(&m_qmin.back());
    decayjets+=ct->GetLeg(i).QCDJets();
    m_max_decayjets+=ct->GetLeg(i).Point()->t-10;
    if (ct->GetLeg(i).QCDJets()<ct->GetLeg(i).Point()->t-10) {
      ptiji=i<2?m_qmin_i:m_qmin_f;
      ptijj=m_qmin_f;
    }
    ++m_cut;
  }
  ct->Up()->GetLeg(i).SetQMin(ptiji);
  ct->Up()->GetLeg(j).SetQMin(ptijj);
  SetQMin(ct->Up(),ref,decayjets,hardjets,hardqmin);
  for (int i(0);i<ct->NLegs();++i)
    if (ct->GetLeg(i).QMin()<0.0)
      ct->GetLeg(i).SetQMin(hardjets==(size_t)m_maxqcdjets?hardqmin:
			    Min(hardqmin,dabs(ct->GetLeg(i).QMin())));
  if (msg.LevelIsDebugging()) {
    msg_Debugging()<<"table "<<ct->Up()->Number()
		   <<" ("<<i<<"&"<<j<<") -> qmin = {";
    for (int l(0);l<ct->Up()->NLegs()-1;++l)
      msg_Debugging()<<ct->Up()->GetLeg(l).QMin()<<",";
    msg_Debugging()<<ct->Up()->GetLeg(ct->Up()->NLegs()-1).QMin()
		   <<"} <- "<<ptiji<<","<<ptijj<<" ("
		   <<sqrt(ct->GetLeg(i).KT2QCD())<<")\n";
  }
}

void Cluster_Partons_Base::SetQMin(Combine_Table_Base *const ref)
{
  msg_Debugging()<<METHOD<<"(): {\n";
  msg_Indent();
  m_cut=0;
  double hardqmin(sqrt(p_ct->Kt2QCDHard(2)));
  msg_Debugging()<<"p_t^2 core process: "<<hardqmin<<"\n";
  m_max_decayjets=0;
  Combine_Table_Base *ct(p_ct);
  while (ct->Up()) ct=ct->Up();
  m_qmin.clear();
  m_qmin.reserve(ct->NLegs());
  m_qmin.resize(p_ct->NLegs(),-m_qmin_f);
  m_qmin[1]=m_qmin[0]=-m_qmin_i;
  for (int i(0);i<p_ct->NLegs();++i) p_ct->GetLeg(i).SetPQMin(&m_qmin[i]);
  size_t hardjets(0), decayjets(0);
  SetQMin(p_ct,ref,decayjets,hardjets,hardqmin);
  if (msg.LevelIsDebugging()) {
    msg_Debugging()<<"table "<<p_ct->Number()<<"       -> qmin = {";
    for (int l(0);l<p_ct->NLegs()-1;++l)
      msg_Debugging()<<p_ct->GetLeg(l).QMin()<<",";
    msg_Debugging()<<p_ct->GetLeg(p_ct->NLegs()-1).QMin()<<"}\n";
  }
  msg_Debugging()<<"}\n";
}

bool Cluster_Partons_Base::FillLegs(Leg * alegs, Point * root, int & l, int maxl) 
{
  if (l>= maxl) {
    msg.Error()<<" Error in Cluster_Partons_Base::FillLegs() !!! "<<std::endl;
    return 0;
  }
  if (l==0) {
    alegs[root->number]=Leg(root);
    alegs[root->number].SetExternal(1);
    l++;
  }
  if (root->left) {
    if (root->middle) return 0; // four vertex 
    return FillLegs(alegs,root->left,l,maxl)*FillLegs(alegs,root->right,l,maxl);
  } 
  else {
    alegs[root->number]=Leg(root);
    alegs[root->number].SetExternal(1);
    l++;
    return 1;
  }
}

int Cluster_Partons_Base::SetDecayColours(Vec4D * p, Flavour * fl,int col1,int col2)
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
  m_q2_fss = m_q2_iss = p[0].Abs2();
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
    msg.Error()<<"Error in Cluster_Partons_Base::SetDecayColours:"<<std::endl
	       <<"   Cannot handle single color in 1 -> 2 process :"
               <<"   "<<fl[0]<<" "<<fl[1]<<" "<<fl[2]<<std::endl
	       <<"   Will abort the run."<<std::endl;
    abort();
  }
}

int Cluster_Partons_Base::SetColours(Vec4D * p,Flavour * fl)
{
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
    m_colors[i][0]=m_colors[i][1]=0;
  }

  switch (ncol) {
  case 4:
    return Set4Colours(nquark,ngluon,p,fl);
  case 3:
    return Set3Colours(nquark,ngluon,p,fl);
  case 2:
    return Set2Colours(nquark,ngluon,p,fl);
  case 1:
    msg.Out()<<"Error in Cluster_Partons_Base::SetColours() : called for 1 coloured object. \n"
	     <<"   Don't know how to handle this ! Abort the run."<<std::endl;
    for (int i=0; i<4; ++i) msg.Out()<<i<<" : "<<fl[i]<<std::endl;
    abort();
  case 0:
    m_q2_fss = m_q2_iss = (p[0]+p[1]).Abs2();
    m_hard_nqed = 2;
    m_hard_nqcd = 0;
    return 0;
  }
  return 1;
}

int Cluster_Partons_Base::Set4Colours(const int nquark,const int ngluon,Vec4D * p,Flavour * fl)
{
  int prop(p_ct->IdentifyHardPropagator());
  if (fl[0].IsGluon() || fl[1].IsGluon() || 
      fl[2].IsGluon() || fl[3].IsGluon() || prop<0) {
    msg.Out()<<METHOD<<"(): Cannot set colours for "<<std::endl;
    Combine_Table_Base *ct(p_ct);
    while (ct->Up()!=NULL) ct=ct->Up();
    msg.Error()<<*ct<<std::endl;
    return false;
  }
  switch (prop) {
  case 1:
    if (!fl[0].IsAnti()) m_colors[0][0]=m_colors[1][1]=500;
    else m_colors[0][1]=m_colors[1][0]=500;
    if (!fl[2].IsAnti()) m_colors[2][0]=m_colors[3][1]=501;
    else m_colors[2][1]=m_colors[3][0]=501;
    m_q2_hard=m_q2_fss=m_q2_iss=dabs((p[0]+p[1]).Abs2());
    break;
  case 2:
    if (!fl[0].IsAnti()) m_colors[0][0]=m_colors[2][0]=500;
    else m_colors[0][1]=m_colors[2][1]=500;
    if (!fl[1].IsAnti()) m_colors[1][0]=m_colors[3][0]=501;
    else m_colors[1][1]=m_colors[3][1]=501;
    m_q2_hard=m_q2_fss=m_q2_iss=dabs((p[0]-p[2]).Abs2());
    break;
  case 3:
    if (!fl[0].IsAnti()) m_colors[0][0]=m_colors[3][0]=500;
    else m_colors[0][1]=m_colors[3][1]=500;
    if (!fl[1].IsAnti()) m_colors[1][0]=m_colors[2][0]=501;
    else m_colors[1][1]=m_colors[2][1]=501;
    m_q2_hard=m_q2_fss=m_q2_iss=dabs((p[0]-p[3]).Abs2());
    break;
  }
  return true;
}

int Cluster_Partons_Base::Set2Colours(const int nquark,const int ngluon,Vec4D * p,Flavour * fl)
{
  m_hard_nqed = 2;
  m_hard_nqcd = 0;
  if (nquark+ngluon>2) {
    msg.Error()<<"ERROR in Cluster_Partons_Base::Set2Colours("<<nquark<<","<<ngluon<<")"<<std::endl
	       <<"   Wrong number of colours, abort."<<std::endl;
    abort();
  }
  int connected[2] = {-1,-1};
  int j(0);
  for (int i=0;i<4;i++) {
    if (!fl[i].Strong()) continue;
    if (fl[i].IsQuark() || fl[i].IsSquark()) {
      m_colors[i][0+int(fl[i].IsAnti())] = 500;
    }
    else if (fl[i].IsGluon()) {
      m_colors[i][j] = 500; m_colors[i][1-j] = 501;
    }
    connected[j++]=i;
  }    
  if (connected[0]<2 ^ connected[1]<2) {
    // t/u channel.
    m_q2_hard = m_q2_fss = m_q2_iss = 
      dabs((p[connected[0]]-p[connected[1]]).Abs2());
    msg_Debugging()<<"found t/u channel: hard = fss = iss = "
		   <<sqrt(m_q2_fss)<<"\n";
  }
  else {
    // s channel.
    m_q2_hard = m_q2_fss = m_q2_iss = 
      dabs((p[connected[0]]+p[connected[1]]).Abs2());
    msg_Debugging()<<"found s channel: hard = fss = iss = "
		   <<sqrt(m_q2_fss)<<"\n";
  }
  return 0;
}

int Cluster_Partons_Base::Set3Colours(const int nquark,const int ngluon,Vec4D * p,Flavour * fl)
{
  m_hard_nqed = 1;
  m_hard_nqcd = 1;
  
  int connected[3] = {-1,-1,-1};
  int singlet      = -1;
  int j(0);
  if (ngluon==3) {
    for (int i=0;i<4;i++) {
      if (!fl[i].Strong()) { 
	singlet = i; 
	continue; 
      }
      if (fl[i].IsGluon()) {
	connected[j] = i;
	m_colors[i][0+(i>1)] = 500+j; 
	if (j==2) j=-1;
	m_colors[i][1-(i>1)] = 501+j;
	j++;
      }    
    }
  }
  else if (ngluon==1) {
    for (int i=0;i<4;i++) {
      if (!fl[i].Strong()) { 
	singlet = i; 
	continue; 
      }
      if (fl[i].IsQuark() || fl[i].IsSquark()) {
	connected[j] = i;
	m_colors[i][0+int(fl[i].IsAnti())] = 500+j;
	j++;
      }
    }
    if (singlet<2) {
      /*
	0-------=====2
	       |
	       |  t
	       |
        1      +-----3

        shower scale is s/t
      */
      m_q2_fss = (p[0]+p[1]).Abs2();
      m_q2_iss = dabs((p[0]-p[2]).Abs2());
      msg_Debugging()<<"cp: s/t-channel, fss "<<m_q2_fss
		     <<", iss "<<m_q2_iss<<"\n";
    }
    else {
      /*
	0=====-------2
	      |
	      |  t
	      |
	1-----+      3

        shower scale is t/s
      */
      m_q2_fss = dabs((p[0]-p[2]).Abs2());
      m_q2_iss = (p[0]+p[1]).Abs2();
      msg_Debugging()<<"cp: t/s-channel, fss "<<m_q2_fss
		     <<", iss "<<m_q2_iss<<"\n";
    }
    bool tmode = (connected[0]<2 ^ connected[1]<2);
    for (int i=0;i<4;i++) {
      if (fl[i].IsGluon()) {
	if (tmode) {
	  if (i<2) {
	    m_colors[i][1-int(fl[connected[1]].IsAnti())] = 500;
	    m_colors[i][0+int(fl[connected[0]].IsAnti())] = 501;
	  }
	  else {
	    m_colors[i][1-int(fl[connected[0]].IsAnti())] = 501;
	    m_colors[i][0+int(fl[connected[1]].IsAnti())] = 500;
	  }
	}
	else {
	  for (int j=0;j<2;j++) {
	    m_colors[i][j] += m_colors[connected[0]][j] + m_colors[connected[1]][j];
	  }
	}
      }
    }    
  }
  for (int i=0;i<4;i++) {
    msg_Debugging()<<METHOD<<" "<<i<<" "<<p_ct->GetLeg(i).Point()->fl<<" "
		   <<m_colors[i][0]<<" "<<m_colors[i][1]<<std::endl;
  }
  return 0;
}

void Cluster_Partons_Base::FixJetvetoPt2(double & jetveto_pt2)
{
  msg_Debugging()<<METHOD<<"(): {\n";
  msg_Indent();
  if (m_jv_pt2!=std::numeric_limits<double>::max()) {
    jetveto_pt2=m_jv_pt2;
    msg_Debugging()<<"} -> "<<sqrt(jetveto_pt2)<<"\n";
    return;
  }
  jetveto_pt2=std::numeric_limits<double>::max();
  double pt2minqed(jetveto_pt2), pt2min(p_ct->MinKt2QCD());
  if (pt2min>0.0 && pt2min<std::numeric_limits<double>::max()) {
    jetveto_pt2=pt2min;
    msg_Debugging()<<"MinKt2QCD() : "<<sqrt(jetveto_pt2)<<"\n";
  }
  else {
    pt2min=p_ct->MinKt2QED();
    if (pt2min>0.0 && pt2min<std::numeric_limits<double>::max()) {
      pt2minqed=pt2min;
      msg_Debugging()<<"MinKt2QED() : "<<sqrt(jetveto_pt2)<<"\n";
    }
  }
  double hardq2min(p_ct->Kt2QCDHard(2));
  msg_Debugging()<<"core process: "<<sqrt(hardq2min)<<"\n";
  jetveto_pt2 = Min(jetveto_pt2,hardq2min);
  if (jetveto_pt2==std::numeric_limits<double>::max()) {
    jetveto_pt2=sqrt(p_ct->Kt2QCDHard());
    msg_Debugging()<<"reset q_{min,hard} -> "<<jetveto_pt2<<"\n";    
  }
  if (jetveto_pt2==std::numeric_limits<double>::max()) {
    msg.Error()<<METHOD<<"(): Warning. "
	       <<"Using QED scale for weight calculation."<<std::endl;
    jetveto_pt2=pt2minqed;
  }
  if (jetveto_pt2==std::numeric_limits<double>::max()) {
    msg.Error()<<METHOD<<"(..): No minimum scale found in table {\n"
	       <<*p_ct<<"}\n"<<std::endl;
  }
  msg_Debugging()<<"} -> "<<sqrt(jetveto_pt2)<<"\n";
  m_jv_pt2=jetveto_pt2;
}


Flavour Cluster_Partons_Base::Flav(int i) {
  if (p_ct) {
    return p_ct->Flav(i);
  }
  msg.Error()<<"ERROR in Cluster_Partons_Base::Flav. No ct."<<std::endl;
  return 0;
}

Vec4D Cluster_Partons_Base::Momentum(int i) {
  if (p_ct) return p_ct->Momentum(i);
  msg.Error()<<"ERROR in Cluster_Partons_Base::Momentum. No ct."<<std::endl;
  return Vec4D(0.,0.,0.,0.);
}

int Cluster_Partons_Base::Colour(const int part,const int ind) {
  if (part>-1&&part<4 && ind>-1&&ind<2) return m_colors[part][ind];
  msg.Error()<<"ERROR in Cluster_Partons_Base::Colour("<<part<<","<<ind<<"): "<<std::endl
	     <<"   Out of bounds, return -1."<<std::endl;
  return -1;
}

Combine_Table_Base * Cluster_Partons_Base::GetCombineTable() { return p_ct; }

void   Cluster_Partons_Base::JetvetoPt2(double & q2i, double & q2f,
					double &q2lji,double &q2ljf) 
{ 
  msg_Debugging()<<METHOD<<"(..): {\n";
  msg_Indent();
  double dr(rpa.gen.DeltaR());
  q2lji=sqr(m_qmin_i);
  q2ljf=sqr(m_qmin_f/dr);
  double qmin2(m_q2_fsjet/sqr(dr));
  bool maxjets(m_hardjets==m_maxqcdjets && m_cut+m_hardjets>0);
  bool twoscaleis(sqr(m_qmin_i)<m_q2_isjet && 
		  !IsEqual(sqr(m_qmin_i),m_q2_isjet));
  bool twoscalefs(sqr(m_qmin_f)<m_q2_fsjet && 
		  !IsEqual(sqr(m_qmin_f),m_q2_fsjet));
  bool qcd2jet(m_minqcdjets==m_maxqcdjets && p_ct->OrderStrong()>0);
  if (maxjets || twoscaleis || twoscalefs || qcd2jet) FixJetvetoPt2(qmin2);
  if ((twoscaleis || twoscalefs) && 
      !qcd2jet && !maxjets) qmin2=Min(m_q2_fsjet,qmin2);
  //  if (twoscale || qcd2jet) q2lj=0.0;
  q2i = Max(qmin2,sqr(m_qmin_i));
  q2f = Max(qmin2,sqr(m_qmin_f/dr));
  p_jf->SetShowerPt2(qmin2);
  msg_Debugging()<<"hard = "<<m_hardjets<<", cut = "<<m_cut<<", max = "
		 <<m_maxqcdjets<<", min = "<<m_minqcdjets<<", os = "
		 <<p_ct->OrderStrong()<<", up = "<<p_ct->Up()<<"\n";
  msg_Debugging()<<"maxjets = "<<maxjets<<", twoscaleis = "<<twoscaleis
		 <<", twoscalefs = "<<twoscalefs<<", qcd2jet = "
		 <<qcd2jet<<"\n";
  msg_Debugging()<<"} -> q_min = "<<sqrt(qmin2)<<", q_i = "<<sqrt(q2i)
		 <<", q_f = "<<sqrt(q2f)<<", q_lji = "<<sqrt(q2lji)
		 <<", q_ljf = "<<sqrt(q2ljf)<<"\n";
}

