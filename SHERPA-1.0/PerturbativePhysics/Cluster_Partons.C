#include "Cluster_Partons.H"
#include "Message.H"
#include "Matrix_Element_Handler.H"
#include "Initial_State_Shower.H"
#include "Final_State_Shower.H"

using namespace SHERPA;
using namespace APACIC;
using namespace EXTRAXS;
using namespace AMEGIC;
using namespace ATOOLS;
using namespace MODEL;

using std::endl;

Cluster_Partons::Cluster_Partons(Matrix_Element_Handler * me, ATOOLS::Jet_Finder * jf,
				 int maxjetnumber, int isrmode,int isrshoweron, int fsrshoweron) :
  p_me(me),p_jf(jf),m_maxjetnumber(maxjetnumber),
  m_isrmode(isrmode), m_isrshoweron(isrshoweron), m_fsrshoweron(fsrshoweron), p_local_tree(0)
{
  //      p_lastproc  = 0;
  p_combi = 0;
  p_sud   = new NLL_Sudakov(p_jf->Smax(),p_jf->Smin());
  
  /* 0 no sudakow weights, 1 alphas only, 2 full sudakov weight  (but for highest jet number) */
  /* cf. also begin of Cluster_Partons::CalculateWeight() */
  m_sud_mode = 2;  

  p_events   = new long[maxjetnumber];
  p_weight_sum     = new double[maxjetnumber];
  p_weight_sum_sqr = new double[maxjetnumber]; 
  for (int i=0;i<maxjetnumber;++i) {
    p_events[i]=0;
    p_weight_sum[i]=p_weight_sum_sqr[i]=0.;
  }
}

Cluster_Partons::~Cluster_Partons()
{
  if (p_local_tree) delete p_local_tree;
  if (p_combi) delete p_combi;
  if (p_sud)  delete p_sud;
  
  msg.Debugging()<<" Statistics Sudakov Rejection "<<endl;
  for (int i=0;i<m_maxjetnumber;++i) {
    if (p_events[i]==0) continue;
    double w_mean  = p_weight_sum[i]/p_events[i];
    double w_delta = 1./p_events[i] * sqrt(p_weight_sum[i] - 
					   (sqr(p_weight_sum[i])-p_weight_sum_sqr[i])/(p_events[i]-1.));
    msg.Debugging()<<(i+1)<<" : weight="<<w_mean<<" +- "<<w_delta<<endl;
  }

  if (p_events) delete [] p_events;
  if (p_weight_sum) delete [] p_weight_sum;
  if (p_weight_sum_sqr) delete [] p_weight_sum_sqr;
}

bool Cluster_Partons::ClusterConfiguration(Blob * _blob,double _x1,double _x2) {
  int nin         = p_me->Nin();
  int nout        = p_me->Nout();
  Flavour * flavs = p_me->Flavs();
  int nampl       = p_me->NumberOfDiagrams();

  int    nlegs    = nin + nout;
  Leg ** legs     = 0;

  p_blob          = _blob;


  bool reuse=0;

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
    for (int i=0;i<nlegs;++i) amoms[i] = p_me->Momenta()[i];

    p_combi = new Combine_Table(p_jf,amoms,0,m_isrmode,m_isrshoweron);
    p_combi->FillTable(legs,nlegs,nampl);   
    p_ct = p_combi->CalcJet(nlegs,_x1,_x2); 
  }
  else {
    // use the existing combination table and determine best combination sheme
    p_ct = p_combi->CalcJet(nlegs,_x1,_x2,p_me->Momenta());
  }

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
    return FillLegs(alegs,root->left,l,maxl)*FillLegs(alegs,root->right,l,maxl);
  } 
  else {
    alegs[root->number]=Leg(root);
    l++;
    return 1;
  }
}

void Cluster_Partons::CalculateWeight(double hard,double jet) 
{
  double facycut = 1.; // FK **** rpa.test.FactorYcut();
  double facnlly = 1.; // FK **** rpa.test.FactorNLLQ();

  double qmin = facnlly*sqrt(jet);
  double qmax = facnlly*sqrt(hard);

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
  double as_jet  = (*as)(facycut*jet);
  double as_hard = (*as)(facycut*hard);
  
  if (jet<1.e-6 || hard<1.e-6) {
    as_jet  = 1.;
    as_hard = 1.;
  }


  // remove (for e+e- -> Jets) universal 2 Quark sudakov to increase the effectivity!
  int strong=0;
  for (int l=0; l<nlegs; ++l) {
    if (p_ct->GetLeg(l)->fl.Strong()) {// Quark sudakov for each strong interacting particle
      strong++;
    }    
  }
  // weight with correct alphas at hard scale (if needed).
  if (strong>2 && m_sud_mode%10>0) {
    m_weight *= pow(as_hard/as_jet,strong-2);
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

    int i,j;
    // ... determine winner: i,j, and yij
    double ptij = p_ct->GetWinner(i,j);

    if (m_sud_mode%10>1 && ct_down->GetLeg(i)->fl.IsGluon()) {  // Gluon sudakov
      double w_in_g =  p_sud->DeltaG(last_q[i],qmin)/p_sud->DeltaG(facnlly*ptij,qmin);
      m_weight *= w_in_g;
    }
    if (m_sud_mode%10>1 && ct_down->GetLeg(i)->fl.IsQuark()) {  // Quark sudakov
      double w_in_q = p_sud->DeltaQ(last_q[i],qmin)/p_sud->DeltaQ(facnlly*ptij,qmin);
      m_weight *= w_in_q;
    }
    if (m_sud_mode%10>0 && ct_down->GetLeg(i)->fl.Strong()) {   // alphaS factor
      double w_in_as = (*as)(facycut*ptij*ptij)/as_jet;
      m_weight *= w_in_as;
    }

    // store old q values
    last_i[i] = si;
    last_q[i] = facnlly*ptij;
    for (int l=j+1; l<nlegs; ++l ) {
      last_i[l] = last_i[l-1];
      last_q[l] = last_q[l-1];
    }
    last_i[j] = si;
    last_q[j] = ptij;
  }


  // go over all remaining legs
  if ((m_sud_mode%10>1 && njet!=m_maxjetnumber) || m_sud_mode>10) {

    for (int l=0; l<nlegs; ++l) {
      if (p_ct->GetLeg(l)->fl.IsGluon()) { // Gluon sudakov
	double w_out_g = p_sud->DeltaG(last_q[l],qmin);
	m_weight *= w_out_g;
      }
      if (p_ct->GetLeg(l)->fl.IsQuark()) {// Quark sudakov
	double w_out_q = p_sud->DeltaQ(last_q[l],qmin);
	m_weight *= w_out_q;
      }    
    }
  }
  msg.Debugging()<<" sudakov weight="<<m_weight<<std::endl;

  p_ct = ct_tmp;

  p_events[njet-1]+=1;  // count events
  p_weight_sum[njet-1]+=m_weight;
  p_weight_sum_sqr[njet-1]+=sqr(m_weight);
}


int  Cluster_Partons::SetColours(ATOOLS::Vec4D * p, Flavour * fl)
{
  // *** 2 -> 2 processes with unambiguous coulor structure
  // (a) no colors
  // (b) two quarks
  // (c) two quarks and one gluon
  // (d) two gluons  

  for (int i=0; i<4; ++i) m_colors[i][0]=m_colors[i][1]=0;
  double s = (p[0]+p[1]).Abs2();
  double t = (p[0]-p[2]).Abs2();
  double u = (p[0]-p[3]).Abs2();

  m_scale=s;
  
  int ncol   = 0;
  int nquark = 0;
  int ngluon = 0;
  for (int i=0; i<4; ++i) if (fl[i].Strong()) {
    ++ncol;
    if (fl[i].IsQuark() || fl[i].IsSquark()) ++nquark;
    if (fl[i].IsGluon() || fl[i].IsGluino()) ++ngluon;
  }

  // (a) no colors
  if (ncol==0) return 0;

  if (ncol==1) {
    msg.Error()<<"ERROR: Can not handle color singlet in 2 -> 2 process!"<<endl;
    abort();
    return 1;
  }  

  int cols[2]={0,0};
  // (b) two quarks
  if (ncol==2 && nquark==2)
    cols[0]=cols[1]=500;

  // (c) two quarks and one gluon (and one (massive) boson)
  if (ncol==3 && nquark==2 && ngluon==1) {
    cols[0]=500;
    cols[1]=501;
    // in case one leg is massive we add all m^2 to the pt:
    m_scale = (sqr(p[3][1])+sqr(p[3][2])+p[2].Abs2()+p[3].Abs2());
    // naive:    m_scale = (2.*s*t*u)/(s*s+t*t+u*u);
  }
  
  if (cols[0]==0) {
    msg.Error()<<"ERROR: in Cluster_Partons::SetColours "<<endl
	       <<"Can not handle color structure in 2 -> 2 process!"<<endl;
    abort();
    return 1;
  }

  ncol=0;
  int antis[2]={0,0};
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
    if (fl[0].Strong()) {
      if (fl[2].Strong()) m_scale=dabs(t);
      else if (fl[3].Strong()) m_scale=dabs(u);
    }
    else if (fl[1].Strong()) {
      if (fl[2].Strong()) m_scale=dabs(u);
      else if (fl[3].Strong()) m_scale=dabs(t);
    }
  }

  if (ngluon>0) {
    for (int i=0; i<4; ++i) {
      if (fl[i].IsGluon() || fl[i].IsGluino()) {
	if (i<2) {
	  m_colors[i][antis[0]]=cols[1];
	  m_colors[i][antis[1]]=cols[0];
	}
	else {
	  m_colors[i][antis[0]]=cols[0];
	  m_colors[i][antis[1]]=cols[1];
	}
      }
    }
  }
  return 0;

}

void Cluster_Partons::FillTrees(Tree ** ini_trees,Tree * fin_tree,XS_Base * xs)
{
  cout<< "Cluster"<<endl;
  if ((!ini_trees && m_isrshoweron) || (!fin_tree && m_fsrshoweron)) {
    msg.Error()<<"ERROR in Cluster_Partons::FillTrees: no trees! no shower to be performed! "<<endl;
    return;
  } 
  
  // **rm**  if (!m_isron && !m_fsron) return;

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

  std::vector<Knot *> knots;
  std::vector<Knot *> ini_knots; // production points
  knots.reserve(10);
  ini_knots.reserve(10);

  // generate knotlist from pointlist in Combine_Table

  // start initial state
  if (m_isrshoweron) {
    knots.push_back(Point2Knot(ini_trees[0],p_ct->GetLeg(0), p_ct->Momentum(0),'G'));
    knots.push_back(Point2Knot(ini_trees[1],p_ct->GetLeg(1), p_ct->Momentum(1),'G'));
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
      if (xs) {
	cout<<" xscols["<<i<<"]["<<j<<"]"<<xs->Colours()[i][j]<<endl;
	knots[i]->part->SetFlow(j+1,xs->Colours()[i][j]); 
      }
      else {
	cout<<" clcols["<<i<<"]["<<j<<"]"<<m_colors[i][j]<<endl;
	knots[i]->part->SetFlow(j+1,m_colors[i][j]);
      }	
    }
  }

  if (mo) {
  *(mo->part) = Parton(0,Flavour(kf::none),p_ct->Momentum(2)+p_ct->Momentum(3));
  mo->part->SetInfo('M');
  mo->part->SetStatus(2);
  mo->stat   = 0;
  mo->z      = knots[2]->part->Momentum()[0]/mo->part->Momentum()[0];
  mo->E2     = sqr(mo->part->Momentum()[0]);
  mo->t      = mo->part->Momentum().Abs2();
  mo->thcrit = M_PI;
  }

  EstablishRelations(mo,knots[0],knots[1],0);
  EstablishRelations(mo,knots[2],knots[3],1);      

  // determine starting conditions for showers
  // note, that starting conditions for subsequent branches have to be 
  // evaluted during the shower evoultion (since the system, esp. for 
  // final state showers starting from the initial state shower are not
  // known.)
  DetermineColourAngles(knots);

  for (int l=0; l<4; ++l) ini_knots.push_back(knots[l]);

  Tree * tree;
  int nlegs=4;
  int i,j;

  Combine_Table * ct_tmp = p_ct;

  p_ct=p_ct->Up();
  while (p_ct) {
    knots.push_back(0);ini_knots.push_back(0); ++nlegs;
    p_ct->GetWinner(i,j);
    for (int l=knots.size()-1;l>j;--l) knots[l] = knots[l-1];
    if (i>=2) tree = fin_tree; 
    else      tree = ini_trees[i];

    Knot * d1, * d2, * mo1, * mo2;
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
    if (knots[j]) p_blob->InParton(j)->SetFlow(k, knots[j]->part->GetFlow(k));
    else {
      cout<<" Problem in Cluster_Partons ini"<<endl;
    }
  }
  
  for (int i=4; i<2*nlegs ; ++i) {
    int j=i/2;
    int k=i%2+1;
    if (knots[j])  p_blob->OutParton(j-2)->SetFlow(k, knots[j]->part->GetFlow(k));
    else {
      cout<<" Problem in Cluster_Partons fin"<<endl;
    }
  }
}



Knot * Cluster_Partons::Point2Knot(Tree * tree, const Leg & po, 
				   const Vec4D & mom, char info='M') 
{
  Flavour flav(po->fl);
  if (po.ExtraAnti() == -1) flav = flav.Bar();

  Knot * k   = tree->NewKnot();

  bool found = 0;
  for (int i=0;i<p_blob->NInP();i++) {
    if ( (p_blob->InParton(i)->Flav() == flav) &&
	 (p_blob->InParton(i)->Momentum() == mom) ) { 
      *(k->part)   = p_blob->InParton(i);
      found = 1;
    }
  }
  for (int i=0;i<p_blob->NOutP();i++) {
    if ( (p_blob->OutParton(i)->Flav() == flav) &&
	 (p_blob->OutParton(i)->Momentum() == mom) ) {
      if (found) {
	msg.Error()<<"Blob with in- and outgoing particle identical !!!"<<std::endl
		   <<p_blob<<std::endl;
      }
      *(k->part)   = p_blob->OutParton(i);
      found = 1;
    }
  }
  if (!found) *(k->part)   = Parton(0,flav,mom);

  // preliminary parton status!!!
  double scale = sqr(mom[0]);
  k->part->SetInfo(info);
  k->part->SetStatus(1);  //final
  k->tout      = sqr(flav.PSMass());
  k->E2        = scale;
  k->costh     = 0; 
  k->stat      = 3;
  k->thcrit    = M_PI;

  return k;
}

void Cluster_Partons::EstablishRelations(Knot * mo, Knot * d1,Knot * d2,int mode)
{
  if (mode==1) {
    // fin ...
    mo->left  = d1;
    mo->right = d2;
    mo->z     = d1->part->Momentum()[0]/mo->part->Momentum()[0];
    mo->stat  = 0;
    mo->part->SetStatus(2);
    if (mo->part->Info() != 'H') mo->part->SetInfo('f');

    d1->prev  = mo;
    d2->prev  = mo;

    APACIC::Final_State_Shower::EstablishRelations(mo,d1,d2);

    return;
  }
  else if (mode==0) {
    // initial state initialization
    //  status:
    //  p_blob->CMS()          - Vec4D hard event in LAB system
    //  d1->part->Momentum() - in the moment also in LAB system
    //  p_blob->InParton(0)->Momentum() - in CMS system

    double q2      = mo->t;
    // set x1 and x2
    double x1,x2;
    p_ct->GetX1X2(x1,x2);
    d1->x=x1;
    d2->x=x2;

    // set start t
    d1->t = -q2;
    d2->t = -q2;

    // angle condition set via DetermineColourAngles called in FillTrees
  }
  else if (mode==2 || mode==3) {
    // initial state initialization
    //     mo -> d1 (IS) 
    //        -> d2 (FS)

    if (!mo || !d1 || !d2 ) {
      msg.Error()<<"ERROR: can not establish relations with less than three elements "<<endl;
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
      d2->t = -t1;
    }
    else {
      mo->t = t1;
      d1->t = t1;
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
  if (p_ct) return p_ct->Flav(i);
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


bool Cluster_Partons::IsColourConnected(Parton * a, Parton * b) {
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
  for (int l=2;l<knots.size();++l)
    sum += knots[l]->part->Momentum();

  double sprime = sum.Abs2();

  x1=knots[0]->x;
  x2=knots[1]->x;
  sprime = (knots[0]->part->Momentum() + knots[1]->part->Momentum()).Abs2();
  Poincare lab(Vec4D(x1+x2,0.,0.,-(x1-x2)));

  double angle = 0.;
  int start = 0;
  for (int j=start;j<knots.size();++j) {
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
	double th_ex=acos(ivec*jvec/(ivec.Abs()*jvec.Abs()));

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
  for (int i=0;i<n;++i) {
    moms[i]=knots[i]->part->Momentum();
  }
  
  Poincare cms(moms[0]+moms[1]);
  for (int i=0;i<n;++i) {
    knots[i]->part->SetMomentum(cms*knots[i]->part->Momentum());
  }
  Poincare zaxis(knots[0]->part->Momentum(),Vec4D::ZVEC);
  for (int i=0;i<n;++i) {
    knots[i]->part->SetMomentum(zaxis*knots[i]->part->Momentum());
  }
  
  for (int i=0;i<knots.size();++i) {
    double th= ColourAngle(knots,i);
    knots[i]->thcrit=th;
  }

  // restore momenta
  for (int i=0;i<n;++i) {
    knots[i]->part->SetMomentum(moms[i]);
  }
  delete [] moms;
}
