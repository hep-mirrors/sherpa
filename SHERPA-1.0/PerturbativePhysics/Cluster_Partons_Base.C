#include "Cluster_Partons_Base.H"
#include "NLL_Branching_Probabilities.H"
#include "Flow.H"
#include "Poincare.H"
#include "Message.H"
#include "MyStrStream.H"
#include <iomanip>

using namespace SHERPA;
using namespace EXTRAXS;
using namespace AMEGIC;
using namespace ATOOLS;
using namespace MODEL;

#define CA 3.0
#define TR 0.5
#define Nf 5.0

Cluster_Partons_Base::Cluster_Partons_Base(Matrix_Element_Handler * me,ATOOLS::Jet_Finder * jf,
					   const int maxjet, const int isrmode,
					   const int isron,const int fsron) :
  p_me(me), p_runas(NULL), p_jf(jf), p_fssud(NULL), p_issud(NULL), p_ct(NULL), p_combi(NULL), 
  m_maxjetnumber(maxjet),m_isrmode(isrmode), m_isrshoweron(isron),m_fsrshoweron(fsron), 
  m_sud_mode(0), m_kfac(0.), m_perpmode(false), m_counts(0.), m_fails(0.)
{
  // read in some parameters
  Data_Read dr(rpa.GetPath()+"Shower.dat");     // !!!!!!!! SHOWER_DATA_FILE
  m_bp_mode  = dr.GetValue<int>("SUDAKOVTYPE",0);
  if ((m_bp_mode&(28))!=m_bp_mode) {
    msg.Error()<<"WARNING in Cluster_Partons_Base :"<<std::endl
	       <<"   Wrong mode for NLL_Sudakovs: "<<m_bp_mode<<" vs "<<(m_bp_mode&127)<<std::endl
	       <<"   Set it to 12 = ordinary NLL_Sudakovs."<<std::endl;
    m_bp_mode=12;
  }
  if (m_bp_mode&16) m_kfac = CA*(67./18.-M_PI*M_PI/6.)-10./9.*TR*Nf;
  m_as_order = dr.GetValue<int>("SUDAKOVASORDER",0);
  double as_fac   = dr.GetValue<double>("SUDAKOVASFAC",1.);
  int jetratemode = dr.GetValue<int>("CALCJETRATE",-1);
  msg_Info()<<"Initalize Cluster_Partons_Base, runs in mode : \n"
	    <<"   SUDAKOVTYPE      = "<<m_bp_mode<<"\n"
	    <<"   SUDAKOVASORDER   = "<<m_as_order<<"\n"
	    <<"   SUDAKOVASFAC(s)  = "<<as_fac<<"/"<<as->ScaleFactor()<<"\n"
	    <<"   REN.SCALE FACTOR = "<<rpa.gen.RenormalizationScaleFactor()<<"\n" 
	    <<"   CALCJETRATE      = "<<jetratemode<<"\n"
	    <<"   kfac             = "<<m_kfac<<"."<<std::endl;
  p_runas = MODEL::as; 
  
  /* 0 no sudakow weights, 1 alphas only, 2 full sudakov weight  (but for highest jet number) */
  /* cf. also begin of Cluster_Partons_Base::CalculateWeight() */
  if (m_fsrshoweron!=0) {
    p_fssud = new NLL_Sudakov((BPMode::code)(m_bp_mode+1),p_jf->Smax(),p_jf->Smin(),p_runas,jetratemode);
    m_sud_mode += 1;
  }
  if (m_isrshoweron!=0) {
    p_issud = new NLL_Sudakov((BPMode::code)(m_bp_mode+2),p_jf->Smax(),p_jf->Smin(),p_runas,jetratemode);
    m_sud_mode += 2;
  }
  p_events         = new long[m_maxjetnumber];
  p_weight_sum     = new double[m_maxjetnumber];
  p_weight_sum_sqr = new double[m_maxjetnumber]; 
  for (int i=0;i<m_maxjetnumber;++i) {
    p_events[i] = 0; p_weight_sum[i] = p_weight_sum_sqr[i] = 0.;
  }

  double deltar = rpa.gen.DeltaR();
  double ycut   = rpa.gen.Ycut();
  m_qmin_i = sqrt(ycut)*rpa.gen.Ecms();
  m_qmin_f = sqrt(ycut)*deltar*rpa.gen.Ecms();
  m_q2_jet = sqr(Min(m_qmin_i,m_qmin_f));
}

Cluster_Partons_Base::~Cluster_Partons_Base()
{
  if (p_combi)               { delete p_combi;      p_combi      = NULL; }
  if (p_fssud)               { delete p_fssud;    p_fssud        = NULL; }
  
  WriteOutSudakovWeights();
  
  if (p_events)         delete [] p_events;
  if (p_weight_sum)     delete [] p_weight_sum;
  if (p_weight_sum_sqr) delete [] p_weight_sum_sqr;
}

void Cluster_Partons_Base::WriteOutSudakovWeights() 
{
  msg_Info()<<" Statistics Sudakov Rejection "<<std::endl;
  msg_Info()<<" Misclusterings : "<<m_fails/m_counts<<"."<<std::endl;
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
  sstr<<"Sudweights_"<<ecms<<"_"<<mode<<".dat"<<std::endl;
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

Leg **Cluster_Partons_Base::CreateLegs(int &nampl,const int nlegs,const bool reuse)
{
  Leg **legs(NULL);
   // start cluster algorithm :
  if (!reuse) {
    if (p_combi) delete p_combi;
    p_combi = 0;
    // generate a list of "legs" for each amplitude
    legs = new Leg *[nampl];
    for (int k=0;k<nampl;) {
      legs[k] = new Leg[nlegs];
      int l   = 0;
      if (FillLegs(legs[k],p_me->GetDiagram(k),l,nlegs)) ++k;
      else {
	delete [] legs[k];
	--nampl;
      }
    }
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
  CreateFlavourMap();
  return 1;
}


bool Cluster_Partons_Base::FillLegs(Leg * alegs, Point * root, int & l, int maxl) 
{
  if (l>= maxl) {
    msg.Error()<<" Error in Cluster_Partons_Base::FillLegs() !!! "<<std::endl;
    return 0;
  }
  if (l==0) {
    alegs[root->number]=Leg(root);
    l++;
  }
  if (root->left) {
    if (root->middle) return 0; // four vertex 
    return FillLegs(alegs,root->left,l,maxl)*FillLegs(alegs,root->right,l,maxl);
  } 
  else {
    alegs[root->number]=Leg(root);
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
  m_q2_qcd = m_q2_hard = p[0].Abs2();
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
    msg.Out()<<"Error in Cluster_Partons_Base::SetColours() : called for 4 coloured objects. \n"
	     <<"   Don't know how to handle this ! Abort the run."<<std::endl;
    for (int i=0; i<4; ++i) msg.Out()<<i<<" : "<<fl[i]<<std::endl;
    abort();
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
    m_q2_qcd    = m_q2_hard = (p[0]+p[1]).Abs2();
    m_hard_nqed = 2;
    m_hard_nqcd = 0;
    return 0;
  }
  return 1;
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
    m_q2_qcd = m_q2_hard = dabs((p[connected[0]]-p[connected[1]]).Abs2());
  }
  else {
    // s channel.
    m_q2_qcd = m_q2_hard = dabs((p[connected[0]]+p[connected[1]]).Abs2());
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
    bool tmode = (connected[0]<2 ^ connected[1]<2);
    for (int i=0;i<4;i++) {
      if (fl[i].IsGluon()) {
	if (tmode) {
	  if (i<2) {
	    m_colors[i][1-int(fl[connected[1]].IsAnti())] = 500;
	    m_colors[i][0+int(fl[connected[0]].IsAnti())] = 501;
	  }
	  else {
	    m_colors[i][0+int(fl[connected[0]].IsAnti())] = 501;
	    m_colors[i][1-int(fl[connected[1]].IsAnti())] = 500;
	  }
	}
	else {
	  for (int j=0;j<2;j++) 
	    m_colors[i][j] += m_colors[connected[0]][j] + m_colors[connected[1]][j];
	}
      }
    }    
  }

  
  Vec4D q[4];
  for (int i=0;i<4;++i) q[i]=p[i];
  Poincare cms(q[0]+q[1]);
  for (int i=0;i<4;++i) cms.Boost(q[i]);
  Poincare rot(q[0],Vec4D::ZVEC);
  for (int i=0;i<4;++i) rot.Rotate(q[i]);
  
  if (singlet>=2) {
    m_q2_qcd = m_q2_hard  = 0.5*(sqr(q[2][1])+sqr(q[2][2])+sqr(q[3][1])+sqr(q[3][2]));
    m_q2_hard += q[2].Abs2()+q[3].Abs2();
    if (!m_perpmode) {
      double qz = 0;
      if (fl[2].Strong()) qz = q[2][3];
      else qz = q[3][3];
      if (q[0][3]*qz>0.) m_q2_qcd += dabs(q[0].Abs2());
      if (q[1][3]*qz>0.) m_q2_qcd += dabs(q[1].Abs2());    
    }
  }
  return 0;
}

void Cluster_Partons_Base::FixJetvetoPt2(double & jetveto_pt2)
{
  double pt2min = p_ct->MinKt2QCD();
  if (pt2min>0.0 && pt2min<std::numeric_limits<double>::max()) jetveto_pt2=pt2min;
  else {
    pt2min = p_ct->MinKt2QED();
    if (pt2min>0.0 && pt2min<std::numeric_limits<double>::max()) {
      jetveto_pt2=pt2min;
      msg.Error()<<METHOD<<"(): Warning. Using QED scale for weight calculation."<<std::endl;
    }
    else {
      msg.Error()<<"Cluster_Partons_Base::FixJetvetoPt2(..): No minimum scale found : "<<pt2min<<std::endl
		 <<*p_ct<<std::endl<<"-----------------------------------------------------------"<<std::endl;
    }
  }
}


Flavour Cluster_Partons_Base::Flav(int i) {
  if (p_ct) {
    Flavour fl(p_ct->Flav(i));
    Flavour_Map::const_iterator cit=m_flmap.find(fl);
    if (cit!=m_flmap.end()) fl=cit->second;
    return fl;
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

Flavour_Map * Cluster_Partons_Base::GetFlavourMap() { return &m_flmap; }

void Cluster_Partons_Base::CreateFlavourMap() 
{
  Process_Base * proc=static_cast<Process_Base*>(p_me->GetAmegic()->GetProcess());
  Process_Base * partner=proc->Partner();
  //m_nstrong   = proc->NStrong();
  const Flavour * flavs=proc->Flavours();
  double ycut = proc->Ycut();
  if (ycut!=-1.) {
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

void   Cluster_Partons_Base::JetvetoPt2(double & q2i, double & q2f) 
{ 
  double qmin2(m_q2_jet);
  if (m_njet==m_maxjetnumber && m_njet>2)  FixJetvetoPt2(qmin2);
  if (m_njet==m_maxjetnumber && m_njet==2) qmin2 = m_q2_hard;
  q2i = Max(qmin2,sqr(m_qmin_i));
  q2f = Max(qmin2,sqr(m_qmin_f)); 
  p_jf->SetShowerPt2(qmin2);
}

