#include "Cluster_Partons.H"
#include "Message.H"
#include "Matrix_Element_Handler.H"
#include "Initial_State_Shower.H"
#include "Final_State_Shower.H"

using namespace SHERPA;
using namespace APACIC;
using namespace EXTRAXS;
using namespace AMEGIC;
using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace AORGTOOLS;
using namespace MODEL;

bool Cluster_Partons::ClusterConfiguration(Blob * _blob) {
  msg.Debugging()<<"In Cluster_Partons::ClusterConfiguration("<<p_me->ProcessName()<<")"<<std::endl;

  int nin         = p_me->Nin();
  int nout        = p_me->Nout();
  Flavour * flavs = p_me->Flavs();
  int nampl       = p_me->NumberOfDiagrams();

  int    nlegs    = nin + nout;
  Leg ** legs     = 0;

  blob            = _blob;


  // *AS* reusing combi does not work in the moment (mismatch of momenta!!!)
  //      would have to check if same process
  bool reuse=0;

  // start cluster algorithm :
  if (!reuse) {
    if (combi) delete combi;
    combi    = 0;

    // generate a list of "legs" for each amplitude
    legs = new Leg *[nampl];
    for (int k=0;k<nampl;) {
      legs[k] = new Leg[nlegs];
      int l   = 0;
      if (FillLegs(legs[k],p_me->GetDiagram(k),l,nlegs)) {
	// generate a "final legs"-list for all topologies
	msg.Debugging()<<(k+1)<<" graph found:"<<(p_me->GetDiagram(k))<<std::endl;
	++k;
      } 
      else {
	// generate NO "final legs"-list for topologies with four-verticies
	msg.Debugging()<<"  dropping Graph"<<(p_me->GetDiagram(k))<<std::endl;
	delete [] legs[k];
	--nampl;
      }
    }
  }  

  ct = 0;
  // if no combination table exist, create it
  if (!combi) {
    /*
      - copy moms to insert into Combine_Table (will be delete there)
      - create new Combine_Table with given momenta and given Jet-measure
      - initialise Combine_Table
      - determine best combination sheme
    */ 
    Vec4D * amoms = new Vec4D[nlegs];
    for (int i=0;i<nlegs;++i) amoms[i] = p_me->Momenta()[i];

    combi = new Combine_Table(jf,amoms,0);
    combi->FillTable(legs,nlegs,nampl);   
    ct = combi->CalcJet(nlegs); 
  }
  else {
    // use the existing combination table and determine best combination sheme
    ct = combi->CalcJet(nlegs,p_me->Momenta());
  }

  msg.Tracking()<<" graph selected ... "<<std::endl;
  msg.Tracking()<<combi<<std::endl;

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

void Cluster_Partons::CalculateWeight(double hard,double jet) {
  msg.Tracking()<<"In Cluster_Partons::CalculateWeight("<<hard<<","<<jet<<")"<<std::endl;
  double facycut=rpa.test.FactorYcut();
  double facnlly=rpa.test.FactorNLLQ();

  double qmin = facnlly*sqrt(jet);
  double qmax = facnlly*sqrt(hard);

  int nlegs   = ct->NLegs();
  // count jets

  int njet=nlegs-2;
  Combine_Table * ct_test = ct;
  while (ct_test->Up()) {
    ct_test=ct_test->Up();
    ++njet;
  }
//   msg.Out()<<" njets="<<njet<<std::endl;
//   msg.Out()<<" nlegs="<<nlegs<<std::endl;
//   if (njet==maxjetnumber) {
//     msg.Out()<<" reduced weight!!! "<<std::endl;
//   }


  weight      = 1.;

  int si = 0;
  std::vector<double> last_q(nlegs,qmax);
  std::vector<int>    last_i(nlegs,si);
  double as_jet  = (*as)(facycut*jet);
  double as_hard = (*as)(facycut*hard);  //*AS* check this!
  //  double as_hard = as->AsFixed();

  // remove (for e+e- -> Jets) universal 2 Quark sudakov to increase the effectivity!
  int strong=0;
  for (int l=0; l<nlegs; ++l) {
    if (ct->GetLeg(l)->fl.Strong()) {// Quark sudakov for each strong interacting particle
      /* *AS*
      if (njet==maxjetnumber) {
	msg.Tracking()<<"Multiply weight by (out hard)  "<<1./sud->DeltaQ(qmax,qmin)
		       <<"   "<<qmax<<" --> "<<qmin<<std::endl;
	weight /= sud->DeltaQ(qmax,qmin);
      }
      */
      strong++;
    }    
  }
  // weight with correct alphas at hard scale (if needed).
  if (strong>2) {
    msg.Tracking()<<"Multiply weight by (as hard)"<<pow(as_hard/as_jet,strong-2)<<std::endl;
    weight *= pow(as_hard/as_jet,strong-2);
  }


  Combine_Table * ct_tmp = ct;

  // anti-cluster and ...
  while (ct->Up()) {
    // ... add sudakov and alphaS factor for each clustering
    Combine_Table * ct_down = ct;
    ct=ct->Up();
    ++si;
    ++nlegs;
    // make space in vector:
    last_q.push_back(qmax);
    last_i.push_back(si);

    int i,j;
    // ... determine winner: i,j, and yij
    double ptij = ct->GetWinner(i,j);

    if (ct_down->GetLeg(i)->fl.IsGluon()) {  // Gluon sudakov
      weight *= sud->DeltaG(last_q[i],qmin)/sud->DeltaG(facnlly*ptij,qmin);
      msg.Tracking()<<"Multiply weight by (in G)  "
		     <<sud->DeltaG(last_q[i],qmin)/sud->DeltaG(facnlly*ptij,qmin)
		     <<"  "<<last_q[i]<<" -> "<<ptij<<"  ("<<qmin<<")"<<std::endl;
    }
    if (ct_down->GetLeg(i)->fl.IsQuark()) {  // Quark sudakov
      weight *= sud->DeltaQ(last_q[i],qmin)/sud->DeltaQ(facnlly*ptij,qmin);
      msg.Tracking()<<"Multiply weight by (in Q)  "
		     <<sud->DeltaQ(last_q[i],qmin)/sud->DeltaQ(facnlly*ptij,qmin)
		     <<"  "<<last_q[i]<<" -> "<<ptij<<"  ("<<qmin<<")"<<std::endl;
    }
    if (ct_down->GetLeg(i)->fl.Strong()) {   // alphaS factor
      weight *= (*as)(facycut*ptij*ptij)/as_jet;
      msg.Tracking()<<"Multiply weight by (as in)  "<<(*as)(facycut*ptij*ptij)/as_jet<<std::endl;
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
  if (njet!=maxjetnumber) {

    for (int l=0; l<nlegs; ++l) {
      if (ct->GetLeg(l)->fl.IsGluon()) { // Gluon sudakov
	weight *= sud->DeltaG(last_q[l],qmin);
	msg.Tracking()<<"Multiply weight by (out G) "<<sud->DeltaG(last_q[l],qmin)
		       <<"  "<<last_q[l]<<" -> "<<qmin<<std::endl;
      }
      if (ct->GetLeg(l)->fl.IsQuark()) {// Quark sudakov
	weight *= sud->DeltaQ(last_q[l],qmin);
	msg.Tracking()<<"Multiply weight by (out Q)"<<sud->DeltaQ(last_q[l],qmin)
		       <<"  "<<last_q[l]<<" -> "<<qmin<<std::endl;
      }    
    }
  }
  msg.Tracking()<<" sudakov weight="<<weight<<std::endl;

  ct = ct_tmp;
}


int  Cluster_Partons::SetColours(AMATOOLS::Vec4D * p, Flavour * fl)
{
  // *** 2 -> 2 processes with unambiguous coulor structure
  // (a) no colors
  // (b) two quarks
  // (c) two quarks and one gluon
  // (d) two gluons  

  for (int i=0; i<4; ++i) colors[i][0]=colors[i][1]=0;
  double s = (p[0]+p[1]).Abs2();
  double t = (p[0]-p[2]).Abs2();
  double u = (p[0]-p[3]).Abs2();

  scale=s;
  
  int ncol   = 0;
  int nquark = 0;
  int ngluon = 0;
  for (int i=0; i<4; ++i) if (fl[i].Strong()) {
    ++ncol;
    if (fl[i].IsQuark() || fl[i].IsSquark()) ++nquark;
    if (fl[i].IsGluon() || fl[i].IsGluino()) ++ngluon;
  }

  msg.Tracking()<<"SetColors : "<<ncol<<"/"<<nquark<<"/"<<ngluon<<endl;

  // (a) no colors
  if (ncol==0) return 0;

  if (ncol==1) {
    cerr<<"ERROR: Can not handle color singlet in 2 -> 2 process!"<<endl;
    abort();
    return 1;
  }  

  int cols[2]={0,0};
  // (b) two quarks
  if (ncol==2 && nquark==2)
    cols[0]=cols[1]=500;

  // (c) two quarks and one gluon
  if (ncol==3 && nquark==2 && ngluon==1) {
    cols[0]=500;
    cols[1]=501;
    scale = (2.*s*t*u)/(s*s+t*t+u*u);
  }
  
  if (cols[0]==0) {
    cout<<"ERROR: in Cluster_Partons::SetColours "<<endl;
    cout<<"Can not handle color structure in 2 -> 2 process!"<<endl;
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
      colors[i][fl[i].IsAnti()]=cols[ncol];
      ++ncol;
    }
  }
  if (!test) {
    msg.Tracking()<<" quark flow from initial to final state !!! "<<endl;
    scale = (2.*s*t*u)/(s*s+t*t+u*u);
  }
  // Note: EW Boson production schould have scale M_B !!! ?


  if (ngluon>0) {
    for (int i=0; i<4; ++i) {
      if (fl[i].IsGluon() || fl[i].IsGluino()) {
	if (i<2) {
	  colors[i][antis[0]]=cols[1];
	  colors[i][antis[1]]=cols[0];
	}
	else {
	  colors[i][antis[0]]=cols[0];
	  colors[i][antis[1]]=cols[1];
	}
      }
    }
  }
  return 0;

}

void Cluster_Partons::FillTrees(Tree ** ini_trees,Tree * fin_tree,XS_Base * xs)
{
  std::vector<Knot *> knots;
  std::vector<Knot *> ini_knots; // production points
  knots.reserve(10);
  ini_knots.reserve(10);

  // generate knotlist from pointlist in Combine_Table

  // start initial state
  knots.push_back(Point2Knot(ini_trees[0],ct->GetLeg(0), ct->Momentum(0),'G'));
  knots.push_back(Point2Knot(ini_trees[1],ct->GetLeg(1), ct->Momentum(1),'G'));

  Knot * mo   = fin_tree->NewKnot();
  knots.push_back(Point2Knot(fin_tree    ,ct->GetLeg(2), ct->Momentum(2),'H'));
  knots.push_back(Point2Knot(fin_tree    ,ct->GetLeg(3), ct->Momentum(3),'H'));

  knots[0]->part->SetDecayBlob(blob);  
  knots[1]->part->SetDecayBlob(blob);
  knots[2]->part->SetProductionBlob(blob);
  knots[3]->part->SetProductionBlob(blob);

  for (int i=0;i<4;i++) {
    for (int j=0;j<2;j++) {
      if (xs) {
	msg.Tracking()<<" [ "<<i<<" ][ "<<j<<" ] = "<<xs->Colours()[i][j]<<std::endl;
	knots[i]->part->SetFlow(j+1,xs->Colours()[i][j]);
      }
      else {
	msg.Tracking()<<" [ "<<i<<" ][ "<<j<<" ] = "<<colors[i][j]<<std::endl;
	knots[i]->part->SetFlow(j+1,colors[i][j]);
      }	
    }
  }

  *(mo->part) = Parton(0,Flavour(kf::none),ct->Momentum(2)+ct->Momentum(3));
  mo->part->SetInfo('M');
  mo->part->SetStatus(2);
  mo->stat   = 0;
  mo->z      = knots[2]->part->Momentum()[0]/mo->part->Momentum()[0];
  mo->E2     = sqr(mo->part->Momentum()[0]);
  mo->t      = mo->part->Momentum().Abs2();
  mo->thcrit = M_PI;

  EstablishRelations(mo,knots[2],knots[3],1);      
  EstablishRelations(mo,knots[0],knots[1],0);

  // determine starting conditions for showers
  // note, that starting conditions for subsequent branches have to be 
  // evaluted during the shower evoultion (since the system, esp. for 
  // final state showers starting from the initial state shower are not
  // known.
  DetermineColourAngles(knots);

  for (int l=0; l<4; ++l) ini_knots.push_back(knots[l]);

  Tree * tree;
  int nlegs=4;
  int i,j;

  Combine_Table * ct_tmp = ct;

  ct=ct->Up();
  while (ct) {
    knots.push_back(0);ini_knots.push_back(0); ++nlegs;
    ct->GetWinner(i,j);
    msg.Tracking()<<" winner i="<<i<<"  j="<<j<<std::endl;

    for (int l=knots.size()-1;l>j;--l) knots[l] = knots[l-1];
    if (i>=2) tree = fin_tree; 
    else      tree = ini_trees[i];

    Knot * d1, * d2, * mo1, * mo2;
    if (i>=2) {
      d1 = Point2Knot(tree, ct->GetLeg(i), ct->Momentum(i),'H');
      d2 = Point2Knot(tree, ct->GetLeg(j), ct->Momentum(j),'H');
      d1->part->SetProductionBlob(blob);
      d2->part->SetProductionBlob(blob);

      EstablishRelations(knots[i],d1,d2,1);      
    } 
    else {
      d1 = Point2Knot(tree, ct->GetLeg(i), ct->Momentum(i),'H');
      d2 = Point2Knot(tree, ct->GetLeg(j), ct->Momentum(j),'H');
      d1->part->SetDecayBlob(blob);  
      d2->part->SetDecayBlob(blob);

      EstablishRelations(d1,knots[i],d2,2);      
    }

    // *AS*
    knots[i] = d1;
    knots[j] = d2;

    ct = ct->Up();
  }
          
  ct = ct_tmp;

  // determine colour partners and colour angles (not here see above!)
  //   for (int i=0; i<nlegs; ++i) 
  //     knots[i]->thcrit = ColourAngle(knots,i);

  // update colours in blob
  for (int i=0; i<4 ; ++i) {
    int j=i/2;
    int k=i%2+1;
    blob->InParton(j)->SetFlow(k, knots[j]->part->GetFlow(k));
  }
  
  for (int i=4; i<2*nlegs ; ++i) {
    int j=i/2;
    int k=i%2+1;
    blob->OutParton(j-2)->SetFlow(k, knots[j]->part->GetFlow(k));
  }


  if (rpa.gen.Tracking()) {
    msg.Tracking()<<" ISR Leg 0 : "<<std::endl;
    msg.Tracking()<<ini_trees[0]<<std::endl;
    msg.Tracking()<<" ISR Leg 1 : "<<std::endl;
    msg.Tracking()<<ini_trees[1]<<std::endl;
    msg.Tracking()<<" FSR Tree : "<<std::endl;
    msg.Tracking()<<fin_tree<<std::endl;
  }

  msg.Tracking()<<" checking momentum conservation "<<std::endl;
  Vec4D sum(0.,0.,0.,0);
  int   no=0;
  sum  += Momentum(ini_trees[0]->GetInitiator(),no);
  sum  += Momentum(ini_trees[1]->GetInitiator(),no);
  sum  += Momentum(fin_tree->GetInitiator(),no);
  msg.Tracking()<<" "<<no<<" particles add to sum="<<sum<<std::endl;
  msg.Tracking()<<"out Interface_Tools::FillTrees"<<std::endl;
}



Knot * Cluster_Partons::Point2Knot(Tree * tree, const Leg & po, 
				   const Vec4D & mom, char info='M') 
{
  Flavour flav(po->fl);
  if (po.ExtraAnti() == -1) flav = flav.Bar();

  Knot * k   = tree->NewKnot();

  bool found = 0;
  for (int i=0;i<blob->NInP();i++) {
    if ( (blob->InParton(i)->Flav() == flav) &&
	 (blob->InParton(i)->Momentum() == mom) ) { 
      *(k->part)   = blob->InParton(i);
      found = 1;
    }
  }
  for (int i=0;i<blob->NOutP();i++) {
    if ( (blob->OutParton(i)->Flav() == flav) &&
	 (blob->OutParton(i)->Momentum() == mom) ) {
      if (found) {
	msg.Error()<<"Blob with in- and outgoing particle identical !!!"<<std::endl
		   <<blob<<std::endl;
      }
      *(k->part)   = blob->OutParton(i);
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

    msg.Events()<<"Established relations (FS) : "<<std::endl
		   <<*d1->prev<<std::endl<<*d1<<std::endl<<*d2<<std::endl<<std::endl;
    return;
  }
  else if (mode==0) {
    // initial state initialization
    //  status:
    //  blob->CMS()          - Vec4D hard event in LAB system
    //  d1->part->Momentum() - in the moment also in LAB system
    //  blob->InParton(0)->Momentum() - in CMS system

    double q2      = mo->t;
    Vec4D cms      = d1->part->Momentum() + d2->part->Momentum();
    Vec4D cms_blob = blob->CMS();
    msg.Tracking()<<" ======================================== "<<std::endl;
    msg.Tracking()<<" Establish relations (IS):"<<std::endl; 
    msg.Tracking()<<" ---------------------------------------- "<<std::endl;   
    msg.Tracking()<<" q2 = "<<q2;
    msg.Tracking()<<" cms= "<<cms<<" / "<<cms_blob<<std::endl;
    msg.Tracking()<<" cms_3 = "<<blob->InParton(0)->Momentum()+ blob->InParton(1)->Momentum()<<std::endl;
    // naive 
    double ebeam  = 0.5*rpa.gen.Ecms();
    double s      = sqr(2.*ebeam);
    double sprime = cms.Abs2();
    double x1,x2;
    ct->GetX1X2(x1,x2);
    d1->x=x1;
    d2->x=x2;

    // set naive start t
    d1->t = -q2;
    d2->t = -q2;


    // angle condition !!!!!! ?

    // set cms momenta?!
//     d1->part->SetMomentum(blob->InParton(0)->Momentum());
//     d2->part->SetMomentum(blob->InParton(1)->Momentum());
  }
  else if (mode==2) {
    // initial state initialization
    //     mo -> d1 (IS) 
    //        -> d2 (FS)

    if (!mo || !d1 || !d2 ) {
      cout<<" can not establish relations with less than three elements "<<endl;
    }
    mo->right=d1;
    mo->left =d2;
    d1->prev=mo;
    d2->prev=mo;

    double t1 = d1->part->Momentum().Abs2();
    d1->t = t1;
    mo->t = t1;
    d2->t = -t1;
    double ebeam  = 0.5*rpa.gen.Ecms();
    double x3=mo->part->Momentum()[0]/ebeam;
    mo->x = x3;
    d1->stat = 0;

    // set color conections
    APACIC::Initial_State_Shower::SetColours(d1);

    // angle condition !!!!!! ?


    msg.Tracking()<<" established relations for initial state partons "<<endl;
    msg.Tracking()<<" mo : "<<*mo<<endl;
    msg.Tracking()<<" d1 : "<<*d1<<endl;
    msg.Tracking()<<" d2 : "<<*d2<<endl;
  }

}

Flavour Cluster_Partons::Flav(int i) {
  if (ct) return ct->Flav(i);
  msg.Error()<<"ERROR in Cluster_Partons::Flav. No ct."<<std::endl;
  return 0;
}

Vec4D Cluster_Partons::Momentum(int i) {
  if (ct) return ct->Momentum(i);
  msg.Error()<<"ERROR in Cluster_Partons::Momentum. No ct."<<std::endl;
  return Vec4D(0.,0.,0.,0.);
}

Vec4D Cluster_Partons::Momentum(Knot * mo, int & number) {
  if (mo->left) return Momentum(mo->left,number) + Momentum(mo->right,number);
  number++;
  if ((mo->part->Info()!='G'))  // &&(mo->part->info!='I'))
    return mo->part->Momentum();

  msg.Tracking()<<" omiting "<<mo->part->Momentum()<<std::endl;
  return Vec4D(0.,0.,0.,0.);
}


bool Cluster_Partons::IsColourConnected(Parton * a, Parton * b) {
  msg.Tracking()<<"Check if "<<a->Flav()<<" and "<<b->Flav()<<" are connected."<<std::endl;
  return (( (a->GetFlow(1)!=0) && ( (a->GetFlow(1)==b->GetFlow(1)) || 
				 (a->GetFlow(1)==b->GetFlow(2)))  ) ||
	  ( (a->GetFlow(2)!=0) && ( (a->GetFlow(2)==b->GetFlow(2)) ||
				 (a->GetFlow(2)==b->GetFlow(1)))  )    );
}

double Cluster_Partons::ColourAngle(const std::vector<Knot *> & knots, const int i) {
  if (!((knots[i]->part->Flav()).Strong())) return M_PI;


  double x1=knots[0]->x;
  double x2=knots[1]->x;
  Poincare lab(Vec4D(x1+x2,0.,0.,-(x1-x2)));

  double angle = 0.;
  for (int j=0;j<knots.size();++j) {
    if (j!=i) {
      if (IsColourConnected(knots[i]->part,knots[j]->part)) {
	Vec3D ivec, jvec;
	Vec3D i4vec, j4vec;


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
	double th_ex=acos(ivec*jvec/(ivec.Abs()*jvec.Abs()));
	angle = Max(angle,th_ex);

	if (!IsEqual(angle,M_PI)) {
	  /*
	  Vec4D one, two;
	  if (i<2) one=Vec4D(knots[i]->part->Momentum()[0],Vec3D(knots[i]->part->Momentum()));
	  else one = knots[i]->part->Momentum();
	  if (j<2) two=Vec4D(knots[j]->part->Momentum()[0],Vec3D(knots[j]->part->Momentum()));
	  else two = knots[j]->part->Momentum();
	  Vec4D sum = one + two;
	  double t  = sum.Abs2();
	  double E = sum[0];
	  double z = one[0]/sum[0];
	  double th = sqrt(t/(z*(1.-z)))/E;
	  cout<<"Angle(ap) ="<<th<<"  exact ("<<th_ex<<")"<<endl;
	  */
	}
	else {
	  angle=M_PI;
	}
      }
    }
  }
  //  msg.Tracking()<<"Angle = "<<angle<<std::endl;
  return angle;
}


void Cluster_Partons::DetermineColourAngles(const std::vector<APACIC::Knot *> & knots) {
  //  cout<<" in DetermineColourAngles "<<endl;
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
    //    cout<<i<<" "<<(*knots[i])<<endl;
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
