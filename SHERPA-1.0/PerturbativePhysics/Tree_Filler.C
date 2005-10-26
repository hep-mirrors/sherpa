#include "Tree_Filler.H"
#include "Combine_Table_Base.H"


using namespace SHERPA;
using namespace APACIC;
using namespace ATOOLS;
using namespace std;

Tree_Filler::Tree_Filler(Cluster_Partons_Base * cluster,int maxjetno,int isron, int fsron) : 
  p_cluster(cluster), m_maxjetnumber(maxjetno), 
  m_isrshoweron(isron), m_fsrshoweron(fsron),
  p_local_tree(NULL)
{}

Tree_Filler::~Tree_Filler() { 
  if (p_local_tree) { delete p_local_tree; p_local_tree = NULL;}
}

void Tree_Filler::FillTrees(Blob * blob,Tree ** ini_trees,Tree * fin_tree)
{
  if ((!ini_trees && m_isrshoweron) || (!fin_tree && m_fsrshoweron)) {
    msg.Error()<<"ERROR in Tree_Filler::FillTrees:"<<std::endl
	       <<"   No trees! no shower to be performed! Aboer the run."<<std::endl;
    abort();
  } 
  
  if (!m_isrshoweron) {
    // prepare dummy tree
    if (!p_local_tree)  p_local_tree = new Tree();
    p_local_tree->Reset();
  }
  if (!m_fsrshoweron) {
    // prepare dummy tree
    if (!p_local_tree)  p_local_tree = new Tree();
    p_local_tree->Reset();
    fin_tree = p_local_tree;
  }
  
  int n[2]={0,1};
  if (p_cluster->InSwaped()) {
    n[0]=1;
    n[1]=0;
  }
  
  std::vector<Knot *> knots;
  std::vector<Knot *> ini_knots; // production points
  knots.reserve(10);
  ini_knots.reserve(10);
  
  // count jets
  Combine_Table_Base * ctb(p_cluster->GetCombineTable()), * ct_test(ctb);
  int njet(ctb->NLegs()-2);
  while (ct_test->Up()) {
    ct_test=ct_test->Up();
    ++njet;
  }
  p_flmap = p_cluster->GetFlavourMap();
  
  // generate knotlist from pointlist in Combine_Table
  
  // start initial state
  if (m_isrshoweron) {
    knots.push_back(Point2Knot(blob,ini_trees[n[0]],ctb->GetLeg(0),ctb->Momentum(0),'G'));
    knots.push_back(Point2Knot(blob,ini_trees[n[1]],ctb->GetLeg(1),ctb->Momentum(1),'G'));
  }
  else {
    knots.push_back(Point2Knot(blob,p_local_tree,ctb->GetLeg(0),ctb->Momentum(0),'G'));
    knots.push_back(Point2Knot(blob,p_local_tree,ctb->GetLeg(1),ctb->Momentum(1),'G'));
  }
  
  Knot * mo(fin_tree->NewKnot());
  knots.push_back(Point2Knot(blob,fin_tree,ctb->GetLeg(2),ctb->Momentum(2),'H'));
  knots.push_back(Point2Knot(blob,fin_tree,ctb->GetLeg(3),ctb->Momentum(3),'H'));
  
  knots[0]->part->SetDecayBlob(blob);  
  knots[1]->part->SetDecayBlob(blob);
  knots[2]->part->SetProductionBlob(blob);
  knots[3]->part->SetProductionBlob(blob);
  
  for (int i=0;i<4;i++) {
    for (int j=0;j<2;j++) knots[i]->part->SetFlow(j+1,p_cluster->Colour(i,j));
  }
  
  Vec4D sum(ctb->Momentum(2)+ctb->Momentum(3));
  m_cms_boost = Poincare(sum);
  Vec4D p1(m_cms_boost*sum);
  Vec4D p2(m_cms_boost*ctb->Momentum(2));
  Vec4D p3(m_cms_boost*ctb->Momentum(3));
  
  *(mo->part) = Particle(0,Flavour(kf::none),sum);
  mo->part->SetInfo('M');
  mo->part->SetStatus(2);
  mo->didkin = true;
  mo->stat   = 0;
  mo->zs     = mo->z = p2[0]/p1[0];
  mo->E2     = sqr(p1[0]);
  mo->thcrit = M_PI;

  //we have a virtuality ordered shower, therefore:
  mo->t = mo->part->Momentum().Abs2();
  double scale(p_cluster->HardScale());
  mo->pt2lcm = scale;
  for (int i(0);i<4;++i) knots[i]->pt2lcm=scale;
  if(!(p_cluster->OrderStrong()>0 && njet==m_maxjetnumber)) scale = Max(scale,4.*p_cluster->JetScale());

  EstablishRelations(mo,knots[0],knots[1],0,scale);
  EstablishRelations(mo,knots[2],knots[3],1);      

  // determine starting conditions for showers
  // note, that starting conditions for subsequent branches have to be 
  // evaluted during the shower evoultion (since the system, esp. for 
  // final state showers starting from the initial state shower are not
  // known.)
  DetermineColourAngles(knots);

  for (int i=0;i<4;++i) ini_knots.push_back(mo);

  Tree * tree;
  Knot * d1, * d2;
  int nlegs(4);

  ct_test = ctb;
  ct_test = ct_test->Up();
  int k,l;
  while (ct_test) {
    knots.push_back(0);ini_knots.push_back(0); ++nlegs;
    double scale(sqr(ct_test->GetWinner(k,l)));
    for (int i=knots.size()-1;i>l;--i)     knots[i]     = knots[i-1];
    for (int i=ini_knots.size()-1;i>l;--i) ini_knots[i] = ini_knots[i-1];
    if (k>=2) tree = fin_tree; 
         else tree = ini_trees[n[k]];
    if (k>=2) {
      d1 = Point2Knot(blob,tree,ct_test->GetLeg(k),ct_test->Momentum(k),'H');
      d2 = Point2Knot(blob,tree,ct_test->GetLeg(l),ct_test->Momentum(l),'H');
      d1->part->SetProductionBlob(blob);
      d2->part->SetProductionBlob(blob);
      d1->pt2lcm=scale;
      d2->pt2lcm=scale;
      EstablishRelations(knots[k],d1,d2,1);      
    } 
    else {
      d1 = Point2Knot(blob,tree,ct_test->GetLeg(k),ct_test->Momentum(k),'H');
      d2 = Point2Knot(blob,tree,ct_test->GetLeg(l),ct_test->Momentum(l),'H');
      d1->part->SetDecayBlob(blob);  
      d2->part->SetDecayBlob(blob);
      d1->pt2lcm=scale;
      d2->pt2lcm=scale;
      EstablishRelations(d1,knots[k],d2,2+k);      
    }
    knots[k] = d1;
    knots[l] = d2;
    ct_test = ct_test->Up();
  }

  // update colours in blob
  for (int i=0;i<4;++i) {
    k = i/2;
    l = i%2+1;
    if (blob->InParticle(k)->Flav()==knots[k]->part->Flav()) 
      blob->InParticle(k)->SetFlow(l,knots[k]->part->GetFlow(l));
    else blob->InParticle(1-k)->SetFlow(l,knots[k]->part->GetFlow(l));
  }
  for (int i=4;i<2*nlegs;++i) {
    k = i/2;
    l = i%2+1;
    blob->OutParticle(k-2)->SetFlow(l,knots[k]->part->GetFlow(l));
  }
  if (msg.LevelIsDebugging()) {
    msg.Out()<<" in Tree_Filler::FillTrees("<<m_isrshoweron<<","
	     <<m_fsrshoweron<<")"<<std::endl;
    if (ini_trees) {
      msg.Out()<<"initree[0]:"<<std::endl<<*ini_trees[0]
	       <<"initree[1]:"<<std::endl<<*ini_trees[1];
    }
    msg.Out()<<"fin_tree:"<<std::endl<<*fin_tree
	     <<"****************************************"<<std::endl;
  }
}

void Tree_Filler::FillDecayTree(Tree * fin_tree)
{
  /*
    if (!fin_tree && m_fsrshoweron) {
    msg.Error()<<"ERROR in Tree_Filler::FillDecayTrees: no trees!"
    <<" No shower to be performed! "<<std::endl;
    return;
    } 
    if (p_blob->NInP()!=1 || p_blob->NOutP()!=2) {
    msg.Error()<<"ERROR in Tree_Filler::FillDecayTrees: "<<std::endl
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
  */
}


Knot * Tree_Filler::Point2Knot(Blob * blob,Tree * tree,const Leg & po,const Vec4D & mom,char info) 
{
  Flavour flav(po.Point()->fl);
  // check in map
  Flavour_Map::const_iterator cit = p_flmap->find(flav);
  if (cit!=p_flmap->end()) flav = cit->second;

  if (po.Anti() == -1) flav = flav.Bar();

  Knot * k(NULL);
  bool found(false);
  for (int i=0;i<blob->NInP();i++) {
    if ( (blob->InParticle(i)->Flav() == flav) &&
	 (blob->InParticle(i)->Momentum() == mom) ) { 
      k = tree->NewKnot(blob->InParticle(i));
      found = true;
      break;
    }
  }
  if (!found) {
    for (int i=0;i<blob->NOutP();i++) {
      if ( (blob->OutParticle(i)->Flav() == flav) &&
	   (blob->OutParticle(i)->Momentum() == mom) ) {
	k = tree->NewKnot(blob->OutParticle(i));
	found = true;
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
  if (flav.IsKK() || k->part->DecayBlob()) k->tout=mom.Abs2();
  else k->tout = sqr(flav.PSMass());
  k->E2        = sqr(mom[0]);
  k->costh     = 0; 
  k->stat      = 3;
  k->didkin    = true;
  k->thcrit    = M_PI;
  return k;
}




Vec4D Tree_Filler::Momentum(Knot * mo, int & number) {
  if (mo->left) return Momentum(mo->left,number) + Momentum(mo->right,number);
  number++;
  if ((mo->part->Info()!='G'))  // &&(mo->part->info!='I'))
    return mo->part->Momentum();
  return Vec4D(0.,0.,0.,0.);
}


double Tree_Filler::ColourAngle(const std::vector<Knot *> & knots, const int i) 
{
  if (!((knots[i]->part->Flav()).Strong())) return M_PI;  

  double x1(knots[0]->x), x2(knots[1]->x);
  Vec4D sum(0.,0.,0.,0.);
  for (size_t l=2;l<knots.size();++l) sum += knots[l]->part->Momentum();
  Poincare lab(Vec4D(x1+x2,0.,0.,-(x1-x2)));

  double angle(0.),test(0.),th_ex(0.); //th_crude(0.),t_mo(0.),E_mo(0.),z(0.);
  int start(0);
  Vec3D ivec, jvec;
  Vec4D i4vec, j4vec;
  for (int j=start;j<(int)knots.size();++j) {
    if (j!=i) {
      if (IsColourConnected(knots[i]->part,knots[j]->part)) {
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
	/*
	  t_mo       = (i4vec+j4vec).Abs2();
	  E_mo       = mvec[0];
	  z          = j4vec[0]/E_mo;
	  th_crude   = sqrt( t_mo/(z*(1.- z)))/E_mo;
	*/
	test       = ivec*jvec/(ivec.Abs()*jvec.Abs());
	if (test<=-1.)     th_ex = M_PI; 
	else if (test>=1.) th_ex = 0;   
	else th_ex = acos(test);

	angle = Max(angle,th_ex);
	if (IsEqual(angle,M_PI)) angle = M_PI;
      }
    }
  }
  return angle;
}

bool Tree_Filler::IsColourConnected(Particle * a, Particle * b) 
{
  return (( (a->GetFlow(1)!=0) && ( (a->GetFlow(1)==b->GetFlow(1)) || 
				    (a->GetFlow(1)==b->GetFlow(2)))  ) ||
	  ( (a->GetFlow(2)!=0) && ( (a->GetFlow(2)==b->GetFlow(2)) ||
				    (a->GetFlow(2)==b->GetFlow(1)))  )    );
}


void Tree_Filler::DetermineColourAngles(const std::vector<APACIC::Knot *> & knots) 
{
  int n(knots.size());
  Vec4D * moms = new Vec4D[n];
  for (int i=0;i<n;++i) moms[i] = knots[i]->part->Momentum();
  
  Poincare cms(moms[0]+moms[1]);
  for (int i=0;i<n;++i) knots[i]->part->SetMomentum(cms*knots[i]->part->Momentum());
  Poincare zaxis(knots[0]->part->Momentum(),Vec4D::ZVEC);
  for (int i=0;i<n;++i) knots[i]->part->SetMomentum(zaxis*knots[i]->part->Momentum());

  for (int i=0;i<n;++i) {
    double th = ColourAngle(knots,i);
    knots[i]->thcrit=th;
  }

  for (int i=0;i<n;++i) knots[i]->part->SetMomentum(moms[i]);
  delete [] moms;
}


void Tree_Filler::EstablishRelations(APACIC::Knot * mo,APACIC::Knot * d1,APACIC::Knot * d2,
				     int mode,double scale)
{
  if (mode==1) {
    Vec4D p1(m_cms_boost*mo->part->Momentum());
    Vec4D p2(m_cms_boost*d1->part->Momentum());
    Vec4D p3(m_cms_boost*d2->part->Momentum());

    mo->left  = d1;
    mo->right = d2;
    mo->zs    = mo->z = p2[0]/p1[0];
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
    //  blob->CMS()                     - Vec4D hard event in LAB system
    //  d1->part->Momentum()              - at the moment also in LAB system
    //  blob->InParticle(0)->Momentum() - in CMS system

    // set x1 and x2
    double x1,x2;
    p_cluster->GetCombineTable()->GetX1X2(x1,x2);
    d1->x = x1;
    d2->x = x2;

    // set start t
    d1->t = -scale;
    d2->t = -scale;
  }
  else if (mode==2 || mode==3) {
    // initial state initialization
    //     mo -> d1 (IS) 
    //        -> d2 (FS)

    if (!mo || !d1 || !d2 ) {
      msg.Error()<<"ERROR in Tree_Filler::EstablishRelations: "<<std::endl
		 <<"   Can not establish relations with less than three elements, abort."<<std::endl;
      abort();
    }
    mo->right = d1;
    mo->left  = d2;
    d1->prev  = mo;
    d2->prev  = mo;

    double t1 = d1->part->Momentum().Abs2();
    double t0 = d1->t;
    mo->t     = t0;
    d1->t     = t1;
    d1->tmax  = t1;
    d2->t     = -t1;
    d1->stat  = 0;
    double x1,x2;
    p_cluster->GetCombineTable()->GetX1X2(x1,x2);
    if (mode==2) mo->x = x1;
            else mo->x = x2;

    APACIC::Initial_State_Shower::SetColours(d1);
  }
}
