#include "Combine_Table_CKKW.H"

#include "Run_Parameter.H"
#include "Poincare.H"
#include "Jet_Finder.H"
#include "Exception.H"
#include <iomanip>

using namespace SHERPA;
using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;

Combine_Table_CKKW::Combine_Table_CKKW(Jet_Finder * jf,Vec4D * moms, 
				       Combine_Table_CKKW * up,
				       int isrmode, int isrshoweron):
  Combine_Table_Base(jf,moms,up,isrmode,isrshoweron) {}

Combine_Table_CKKW::~Combine_Table_CKKW() 
{ 
}

double Combine_Table_CKKW::GetWinner(int &i,int &j)    
{ 
  i=m_cdata_winner->first.m_i; 
  j=m_cdata_winner->first.m_j;
  if (m_cdata_winner->second.p_down!=NULL) {
    double kt2qcd(m_cdata_winner->second.p_down->GetLeg(i).KT2QCD());
    if (kt2qcd<std::numeric_limits<double>::max()) return sqrt(kt2qcd);
    return sqrt(m_cdata_winner->second.p_down->GetLeg(i).KT2());
  }
  THROW(fatal_error,"Legs not combined. No Scale information");
  return 0.0;
}

void Combine_Table_CKKW::AddPossibility(const int i,const int j,
					const int ngraph) 
{
  CD_List::iterator cit=m_combinations.find(Combine_Key(i,j));
  if (cit!=m_combinations.end()) {
    cit->second.m_graphs.push_back(ngraph);
    if (p_legs[ngraph][i].Point()->fl.Strong() && 
	p_legs[ngraph][j].Point()->fl.Strong() &&
	p_legs[ngraph][j].Point()->prev->fl.Strong()) 
      cit->second.m_strong=true;
  }
  else {
    Combine_Data cd(0.,ngraph);
    if (p_legs[ngraph][i].Point()->fl.Strong() && 
	p_legs[ngraph][j].Point()->fl.Strong() &&
	p_legs[ngraph][j].Point()->prev->fl.Strong()) 
      cd.m_strong=true;
    m_combinations[Combine_Key(i,j)]=cd;
  }
}

void Combine_Table_CKKW::FillTable(Leg **legs,const int nlegs,const int nampl)
{
  // store information
  p_legs=legs;
  m_nlegs=nlegs;
  m_nampl=nampl;
  // determine possible combinations and corresponding y_ij  if nlegs>4
  if (m_nlegs>4) {
    int start=0;
    // cluster initial state only if isrshower and isr_x is on. 
    if (!m_isrshoweron) start=2;
    for (int i=start; i<m_nlegs; ++i) {  
      if (!m_isr1on && i==0) i=1;
      if (!m_isr2on && i==1) i=2;
      for (int j=i+1; j<m_nlegs; ++j) {
	// never combine "0&1" !
	if (j==1) j=2;
	// check if leg i is combinable with leg j in any graph
	for (int k=0;k<m_nampl;++k) {
// 	  msg_Debugging()<<"start w/ "<<k<<", "
// 			 <<i<<": "<<p_legs[k][i].MapFlavour()<<"\n";
	  if (Combinable(p_legs[k][i],p_legs[k][j])) AddPossibility(i,j,k);
	}
      }
    }
  }
}

CD_List::iterator Combine_Table_CKKW::CalcPropagator(CD_List::iterator &cit)
{
  if (cit->first.m_flav.Kfcode()==kf::none) {
    cit->second.m_sij   = (p_moms[cit->first.m_i]+p_moms[cit->first.m_j]).Abs2();
    cit->second.m_pt2ij = p_jf->MTij2(p_moms[cit->first.m_i],p_moms[cit->first.m_j]);
//     std::cout<<"Calculate m_perp("<<cit->first.m_i<<cit->first.m_j<<") : "
//     	     <<p_moms[cit->first.m_i]<<" & "<<p_moms[cit->first.m_j]
//     	     <<" -> "<<sqrt(cit->second.m_pt2ij)<<std::endl;
    return cit;
  }
  else {
    CD_List::iterator father= 
      m_combinations.find(Combine_Key(cit->first.m_i,cit->first.m_j));
    if (father!=m_combinations.end()) {
      cit->second.m_pt2ij = father->second.m_pt2ij;
      return father;
    }
    else THROW(fatal_error,"No father.");
  }
  return m_combinations.end();
}

Combine_Table_Base *Combine_Table_CKKW::
CalcJet(int nl,const double x1,const double x2,ATOOLS::Vec4D * moms) 
{
  if (p_up==NULL) {
    m_x1 = x1;
    m_x2 = x2;
  }
  bool did_boost(InitStep(moms,nl));
  if (!SelectWinner(did_boost)) return this;
  // if number of legs is still greater 4 Cluster once more
  // if number of legs equals 4, determine end situation
  if (nl<4) THROW(fatal_error,"nlegs < min. Abort.");
  return NextTable(CreateNext(did_boost),x1,x2);
}

bool Combine_Table_CKKW::InitStep(ATOOLS::Vec4D *moms,const int nl)
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

bool Combine_Table_CKKW::SelectWinner(const bool did_boost)
{
  CD_List & cl(m_combinations);
  if (cl.size()==0) {
    // check if boosted  (restore saved moms)
    if (did_boost) { 
      for (size_t i=0;i<m_nl;++i) p_moms[i] = p_save_moms[i]; 
    }
    delete [] p_save_moms;
    return false;
  }
  // calculate pt2ij and determine "best" combination
  m_cdata_winner = cl.end();
  CD_List::iterator qcd_winner(cl.end());
  double kt2(std::numeric_limits<double>::max()), kt2qcd(kt2);
  for (CD_List::iterator cit(cl.begin()); cit!=cl.end(); ++cit) {
    CD_List::iterator tit(CalcPropagator(cit));
    double pt2ij(cit->second.m_pt2ij);
    if (cit->second.m_graphs.size()==0) continue;
    if (tit==cit) {
      // Relevant initial-final clustering.
      if (cit->first.m_i<2 && pt2ij<kt2) {
	// check if this combination has right direction (clustering with correct is particle)
	double d = p_moms[cit->first.m_i][3] * p_moms[cit->first.m_j][3];
	if (d<0.) {
	  // wrong direction; make sure that this cluster will NOT be performed.
	  pt2ij*=1.001;
	  cit->second.m_pt2ij = pt2ij;
	} 
	if (!TestMomenta(cit->first.m_i,cit->first.m_j)) {
	  cit->second.m_pt2ij = pt2ij = std::numeric_limits<double>::max();
	}
      }
    }
    if (pt2ij<kt2) {
      m_cdata_winner=cit;
      kt2=pt2ij;
    }
    if (cit->second.m_strong && pt2ij<kt2qcd) {
      qcd_winner=cit;
      kt2qcd=pt2ij;
    }
  }
  if (qcd_winner!=cl.end()) m_cdata_winner=qcd_winner;
  return true;
}

bool Combine_Table_CKKW::TestMomenta(const int i,const int j)
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

Combine_Table_CKKW *Combine_Table_CKKW::CreateNext(bool did_boost)
{
  // check if boosted  (restore saved moms)
  if (did_boost) { 
    for (size_t i=0;i<m_nl;++i) p_moms[i]=p_save_moms[i]; 
  }
  delete [] p_save_moms;
  --m_nl;
  if (!m_cdata_winner->second.p_down) {
    Leg ** alegs = new Leg*[m_cdata_winner->second.m_graphs.size()];
    for (size_t k=0;k<m_cdata_winner->second.m_graphs.size();++k) {
      alegs[k] = CombineLegs
	(p_legs[m_cdata_winner->second.m_graphs[k]],
	 m_cdata_winner->first.m_i,m_cdata_winner->first.m_j,m_nl,
	 m_cdata_winner->second.m_pt2ij);
    }
    Vec4D * amoms;
    // generate new momenta
    CombineMoms(p_moms,m_cdata_winner->first.m_i,
		m_cdata_winner->first.m_j,m_nl,amoms);
    m_cdata_winner->second.p_down = 
      new Combine_Table_CKKW(p_jf,amoms,this,m_isr1on+2*m_isr2on,m_isrshoweron);
    // initialise Combine_Table
    m_cdata_winner->second.p_down->FillTable(alegs,m_nl,m_cdata_winner->second.m_graphs.size());
  } 
  else {
    ((Combine_Table_CKKW*)m_cdata_winner->second.p_down)->
      CombineMoms(p_moms,m_cdata_winner->first.m_i,m_cdata_winner->first.m_j,m_nl);
  }

  // update x1,2 (before calling CalcJet again
  Combine_Table_CKKW *tab((Combine_Table_CKKW*)m_cdata_winner->second.p_down);
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

Combine_Table_CKKW *Combine_Table_CKKW::NextTable(Combine_Table_CKKW *tab,
						  const double x1,const double x2)
{
  Combine_Table_Base* ct = tab->CalcJet(m_nl,x1,x2);
  m_graph_winner = tab->m_graph_winner;
  // translate back
  m_graph_winner = m_cdata_winner->second.m_graphs[m_graph_winner];
  return (Combine_Table_CKKW*)ct;
}

void Combine_Table_CKKW::IdentifyHardProcess()
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
  for (int i(0);i<m_nampl;++i) {
    if (Combinable(p_legs[i][0],p_legs[i][1]) &&
	Combinable(p_legs[i][2],p_legs[i][3])) {
      double pt2ij(p_jf->MTij2(p_moms[2],p_moms[3]));
      msg_Debugging()<<"s-channel pt = "<<sqrt(pt2ij)<<", m = "
		     <<sqrt(dabs((p_moms[0]+p_moms[1]).Abs2()))<<"\n";
      p_hard[i][0]=CombinedLeg(p_legs[i],0,1);
      SetLegScales(p_hard[i][0],p_legs[i][0],p_legs[i][1],
		   p_moms[0],p_moms[1],pt2ij);
      p_hard[i][1]=CombinedLeg(p_legs[i],2,3);
      SetLegScales(p_hard[i][1],p_legs[i][2],p_legs[i][3],
		   p_moms[2],p_moms[3],pt2ij);
      p_hardc[i][0]=0;
      p_hardc[i][1]=0;
      p_hardc[i][2]=1;
      p_hardc[i][3]=1;
    }
    else if (Combinable(p_legs[i][0],p_legs[i][2]) &&
	     Combinable(p_legs[i][1],p_legs[i][3])) {
      double pt2ij1(p_jf->MTij2(p_moms[0],p_moms[2]));
      double pt2ij2(p_jf->MTij2(p_moms[1],p_moms[3]));
      msg_Debugging()<<"t-channel pt = "<<sqrt(pt2ij1)
		     <<" / "<<sqrt(pt2ij2)<<", m = "
		     <<sqrt(dabs((p_moms[0]+p_moms[2]).Abs2()))<<"\n";
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
    else if (Combinable(p_legs[i][0],p_legs[i][3]) &&
	     Combinable(p_legs[i][1],p_legs[i][2])) {
      double pt2ij1(p_jf->MTij2(p_moms[0],p_moms[3]));
      double pt2ij2(p_jf->MTij2(p_moms[1],p_moms[2]));
      msg_Debugging()<<"u-channel pt = "<<sqrt(pt2ij1)
		     <<" / "<<sqrt(pt2ij2)<<", m = "
		     <<sqrt(dabs((p_moms[0]+p_moms[3]).Abs2()))<<"\n";
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
    m_nstrong=Max(m_nstrong,p_hard[i][0].OrderQCD()+p_hard[i][1].OrderQCD());
  }
}

