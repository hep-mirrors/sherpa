#include "Tree_Filler.H"
#include "Combine_Table_Base.H"

using namespace SHERPA;
using namespace APACIC;
using namespace ATOOLS;
using namespace std;

Tree_Filler::Tree_Filler(Cluster_Partons_Base * cluster,
			 Shower_Handler *shower,int maxjetno,int showermode) : 
  p_cluster(cluster), m_maxjetnumber(maxjetno),
  m_isrshoweron(shower->ISROn()), m_fsrshoweron(shower->FSROn()),
  m_scale2factor(shower->PSScale2Factor()),  
  p_local_tree(NULL), p_fsrshower(shower->GetApacic()->FinShower()),
  m_showermode(showermode)
{
  m_iss_scale_fac=1.0;
  m_fss_scale_fac=1.0;
}

Tree_Filler::~Tree_Filler() { 
  if (p_local_tree) { delete p_local_tree; p_local_tree = NULL;}
}

void Tree_Filler::FillTrees(Blob * blob,Tree ** ini_trees,Tree * fin_tree)
{
  if ((!ini_trees && m_isrshoweron) || (!fin_tree && m_fsrshoweron)) {
    msg_Error()<<"ERROR in Tree_Filler::FillTrees:"<<std::endl
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

  // generate knotlist from pointlist in Combine_Table
  
  // start initial state
  Particle *mep(NULL);
  if (m_isrshoweron) {
    knots.push_back(Point2Knot(blob,ini_trees[n[0]],ctb->GetLeg(0),ctb->Momentum(0),'G',mep));
    if (mep!=NULL) {
      mep->SetFlow(1,p_cluster->Colour(0,0));
      mep->SetFlow(2,p_cluster->Colour(0,1));
    }
    knots.push_back(Point2Knot(blob,ini_trees[n[1]],ctb->GetLeg(1),ctb->Momentum(1),'G',mep));
    if (mep!=NULL) {
      mep->SetFlow(1,p_cluster->Colour(1,0));
      mep->SetFlow(2,p_cluster->Colour(1,1));
    }
  }
  else {
    knots.push_back(Point2Knot(blob,p_local_tree,ctb->GetLeg(0),ctb->Momentum(0),'G',mep));
    knots.push_back(Point2Knot(blob,p_local_tree,ctb->GetLeg(1),ctb->Momentum(1),'G',mep));
  }
  
  Knot * mo(fin_tree->NewKnot());
  knots.push_back(Point2Knot(blob,fin_tree,ctb->GetLeg(2),ctb->Momentum(2),'H',mep));
  if (mep!=NULL) {
    mep->SetFlow(1,p_cluster->Colour(2,0));
    mep->SetFlow(2,p_cluster->Colour(2,1));
  }
  knots.push_back(Point2Knot(blob,fin_tree,ctb->GetLeg(3),ctb->Momentum(3),'H',mep));
  if (mep!=NULL) {
    mep->SetFlow(1,p_cluster->Colour(3,0));
    mep->SetFlow(2,p_cluster->Colour(3,1));
  }
  
  knots[0]->part->SetDecayBlob(blob);  
  knots[1]->part->SetDecayBlob(blob);
  knots[2]->part->SetProductionBlob(blob);
  knots[3]->part->SetProductionBlob(blob);
  
  for (int i=0;i<4;i++) {
    knots[i]->part->SetFlow(1,p_cluster->Colour(i,0));
    knots[i]->part->SetFlow(2,p_cluster->Colour(i,1));
  }
  
  Vec4D sum(ctb->Momentum(2)+ctb->Momentum(3));
  m_cms_boost = Poincare(sum);
  Vec4D p1(m_cms_boost*sum);
  Vec4D p2(m_cms_boost*ctb->Momentum(2));
  Vec4D p3(m_cms_boost*ctb->Momentum(3));
  
  *(mo->part) = Particle(0,Flavour(kf_none),sum);
  mo->part->SetInfo('M');
  mo->part->SetStatus(part_status::decayed);
  mo->didkin = true;
  mo->stat   = 0;
  mo->shower = 0;
  mo->zs     = mo->z = p2[0]/p1[0];
  mo->E2     = sqr(p1[0]);
  mo->thcrit = M_PI;

  msg_Debugging()<<METHOD<<"(): PS scales {"
		 <<"\n  is l: "<<sqrt(p_cluster->ISShowerScale(0))
		 <<"\n  is r: "<<sqrt(p_cluster->ISShowerScale(1))
		 <<"\n  fs  : "<<sqrt(p_cluster->FSShowerScale())<<"\n}\n";
  //we have a virtuality ordered shower, therefore:
  mo->t = m_scale2factor*p_cluster->FSShowerScale();
  // set jet veto scale for each emission
  double q2j(p_cluster->FSJetScale());
  if (p_cluster->OrderStrong()>0) p_cluster->FixJetvetoPt2(q2j);
  mo->pt2lcm = p_cluster->FSShowerScale();
  mo->maxpt2 = m_scale2factor*
    (m_ckkwon?q2j/m_fss_scale_fac:
              sqr(sqrt(mo->t)-sqrt(knots[2]->tout)-sqrt(knots[3]->tout)));
  double scale[2]={p_cluster->ISShowerScale(0),p_cluster->ISShowerScale(1)};
  /*
  //  the setting below aims at generating enough shower radiation
  //  to fill below the onset of the hard me.
  //  this is unphysical, since it does not respect the factorisation scale !
  if (p_cluster->OrderStrong()==0) {
    scale[0]=Max(scale[0],4.*p_cluster->ISJetScale());
    scale[1]=Max(scale[1],4.*p_cluster->ISJetScale());
  }
  */
  double x1,x2;
  p_cluster->GetCombineTable()->GetX1X2(x1,x2);
  EstablishRelations(mo,knots[0],knots[1],0,x1,x2,scale[0],scale[1]);
  EstablishRelations(mo,knots[2],knots[3],1,x1,x2);      
  for (int i(0);i<2;++i) {
    knots[i]->pt2lcm=p_cluster->FactorizationScale(i);
    knots[i]->maxpt2=m_scale2factor*(m_ckkwon?q2j/m_iss_scale_fac:mo->t);
  }
  for (int i(2);i<4;++i) {
    knots[i]->pt2lcm=mo->pt2lcm;
    knots[i]->maxpt2=m_scale2factor*
      (m_ckkwon?q2j/m_fss_scale_fac:
       sqr(sqrt(mo->t)-sqrt(knots[2]->tout)-sqrt(knots[3]->tout)));
  }
  // determine starting conditions for showers
  // note, that starting conditions for subsequent branches have to be 
  // evaluted during the shower evoultion (since the system, esp. for 
  // final state showers starting from the initial state shower are not
  // known.)
  DetermineColourAngles(knots);

  if (dabs(ct_test->Momentum(0)[3])>dabs(ct_test->Momentum(1)[3]))
    knots[1]->dir=-(knots[0]->dir=ct_test->Momentum(0)[3]>0?1:-1);
  else knots[0]->dir=-(knots[1]->dir=ct_test->Momentum(1)[3]>0?1:-1);

  for (int i=0;i<4;++i) ini_knots.push_back(mo);

  Tree * tree;
  Knot * d1, * d2;
  int nlegs(4);

  ct_test = ctb;
  ct_test = ct_test->Up();
  int k,l;
  Particle *mepk(NULL), *mepl(NULL);
  while (ct_test) {
    knots.push_back(0);ini_knots.push_back(0); ++nlegs;
    double scale(sqr(ct_test->GetWinner(k,l)));
    for (int i=knots.size()-1;i>l;--i)     knots[i]     = knots[i-1];
    for (int i=ini_knots.size()-1;i>l;--i) ini_knots[i] = ini_knots[i-1];
    if (k>=2) tree = fin_tree; 
         else tree = ini_trees[n[k]];
    if (k>=2) {
      d1 = Point2Knot(blob,tree,ct_test->GetLeg(k),ct_test->Momentum(k),'H',mepk);
      d2 = Point2Knot(blob,tree,ct_test->GetLeg(l),ct_test->Momentum(l),'H',mepl);
      d1->part->SetProductionBlob(blob);
      d2->part->SetProductionBlob(blob);
      double x1,x2;
      p_cluster->GetCombineTable()->GetX1X2(x1,x2);
      EstablishRelations(knots[k],d1,d2,1,x1,x2);      
    } 
    else {
      d1 = Point2Knot(blob,tree,ct_test->GetLeg(k),ct_test->Momentum(k),'H',mepk);
      d2 = Point2Knot(blob,tree,ct_test->GetLeg(l),ct_test->Momentum(l),'H',mepl);
      d1->part->SetDecayBlob(blob);  
      d2->part->SetDecayBlob(blob);
      double x1,x2;
      p_cluster->GetCombineTable()->GetX1X2(x1,x2);
      EstablishRelations(d1,knots[k],d2,2+k,x1,x2);      
    }
    d1->pt2lcm=d2->pt2lcm=knots[k]->pt2lcm;
    // set max kt2 scale for each emission
    d1->maxpt2=knots[k]->maxpt2;
    d2->maxpt2=knots[k]->maxpt2;
    double asfac(k<2?m_iss_scale_fac:m_fss_scale_fac);
    if (scale*asfac<knots[k]->maxpt2 && 
	!IsEqual(scale*asfac,knots[k]->maxpt2) &&
	d1->part->Flav().Strong() && d2->part->Flav().Strong() &&
	mo->part->Flav().Strong()) {
      msg_Error()<<METHOD<<"(): scale ordering violated in knot "<<k<<".\n"
		 <<"   last scale = "<<knots[k]->maxpt2
		 <<", new scale = "<<scale<<std::endl;
      d1->maxpt2=scale*asfac;
      d2->maxpt2=scale*asfac;
    }
    msg_Debugging()<<"set st "<<d1->kn_no<<" & "<<d2->kn_no
		   <<" -> "<<sqrt(scale)<<" * "<<sqrt(asfac)<<" = "
		   <<sqrt(scale*asfac)<<"\n";
    knots[k] = d1;
    knots[l] = d2;
    if (mepk!=NULL) {
      mepk->SetFlow(1,d1->part->GetFlow(1));
      mepk->SetFlow(2,d1->part->GetFlow(2));
    }
    if (mepl!=NULL) {
      mepl->SetFlow(1,d2->part->GetFlow(1));
      mepl->SetFlow(2,d2->part->GetFlow(2));
    }
    ct_test = ct_test->Up();
  }
  for (size_t i(0);i<knots.size();++i) {
    knots[i]->tout=sqr(knots[i]->part->Flav().PSMass());
    if (knots[i]->part->Flav().IsMassive() &&
	knots[i]->part->Momentum().Abs2()>rpa.gen.Accu()) 
      knots[i]->tout=knots[i]->part->Momentum().Abs2();
    if (dabs(knots[i]->tout)>dabs(knots[i]->t)) knots[i]->t = Sign(knots[i]->t)*knots[i]->tout;
    //std::cout<<METHOD<<" sets tout = "<<knots[i]->tout
    //	     <<" for "<<knots[i]->part->Flav()
    //	     <<" from "<<knots[i]->part->Momentum().Abs2()
    //	     <<" and "<<sqr(knots[i]->part->Flav().PSMass())<<"."<<std::endl;
  }

  if (msg_LevelIsDebugging()) {
    msg_Out()<<" in Tree_Filler::FillTrees("<<m_isrshoweron<<","
	     <<m_fsrshoweron<<")"<<std::endl;
    if (ini_trees) {
      msg_Out()<<"initree[0]:"<<std::endl<<*ini_trees[0]
	       <<"initree[1]:"<<std::endl<<*ini_trees[1];
    }
    msg_Out()<<"fin_tree:"<<std::endl<<*fin_tree
	     <<"****************************************"<<std::endl;
  }
}

void Tree_Filler::FillDecayTree(Tree * fin_tree)
{
  /*
    if (!fin_tree && m_fsrshoweron) {
    msg_Error()<<"ERROR in Tree_Filler::FillDecayTrees: no trees!"
    <<" No shower to be performed! "<<std::endl;
    return;
    } 
    if (p_blob->NInP()!=1 || p_blob->NOutP()!=2) {
    msg_Error()<<"ERROR in Tree_Filler::FillDecayTrees: "<<std::endl
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
    knot->part->SetStatus(part_status::decayed);
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

    double x1,x2;
    p_cluster->GetCombineTable()->GetX1X2(x1,x2);
    EstablishRelations(knots[0],knots[1],knots[2],1,x1,x2);
  */
}


Knot * Tree_Filler::Point2Knot(Blob * blob,Tree * tree,const Leg & po,const Vec4D & mom,char info,Particle *&mep) 
{
  mep=NULL;
  Flavour flav(po.MapFlavour());

  if (po.Anti() == -1) flav = flav.Bar();

  Knot * k(NULL);
  bool found(false);
  for (int i=0;i<blob->NInP();i++) {
    if ( (blob->InParticle(i)->Flav() == flav) &&
	 (blob->InParticle(i)->Momentum() == mom) ) { 
      k = tree->NewKnot(mep=blob->InParticle(i));
      found = true;
      break;
    }
  }
  if (!found) {
    for (int i=0;i<blob->NOutP();i++) {
      if ( (blob->OutParticle(i)->Flav() == flav) &&
	   (blob->OutParticle(i)->Momentum() == mom) ) {
	k = tree->NewKnot(mep=blob->OutParticle(i));
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
  k->part->SetStatus(part_status::active);  //final
  k->tout=sqr(flav.PSMass());
  if (k->part->Flav().IsMassive() && 
      mom.Abs2()>rpa.gen.Accu()) k->tout=mom.Abs2(); 
  if (dabs(k->tout)>dabs(k->t)) k->t = Sign(k->t)*k->tout;
  //std::cout<<METHOD<<" sets tout = "<<k->tout<<" mass = "<<k->tout
  //	   <<" ("<<mom.Abs2()<<")"<<" for "<<flav<<"."<<std::endl;
  k->E2        = sqr(mom[0]);
  k->costh     = 0; 
  k->shower=po.External();
  if ((m_showermode&1) && po.Point()->t-10>=0) {
    k->shower+=2;
    if (po.QCDJets()<po.Point()->t-10) {
      if (po.Point()->t-10+po.Point()->fl.Strong()>2) 
	k->qjv=sqrt(po.Q2Cut(2));
      else k->qjv=sqrt(po.MinKT2QCD());
    }
    else k->qjv=sqrt(po.MinKT2QCD());
    k->qljv=sqrt(po.Q2Cut(1));
    k->maxjets=po.Point()->t-10;
    msg_Debugging()<<METHOD<<"(): ("<<k->kn_no<<") -> n_jets = "
		   <<po.QCDJets()<<"+"<<po.Point()->fl.Strong()
		   <<" ("<<po.Point()->t-10<<"+"<<po.Point()->fl.Strong()
		   <<") at m_kt = "<<sqrt(po.MinKT2QCD())
		   <<" ("<<sqrt(po.Q2Cut(2))<<","<<sqrt(po.Q2Cut(1))
		   <<") -> qjv = "<<k->qjv<<", qljv = "<<k->qljv
		   <<", n_max = "<<k->maxjets<<"\n";
  }
  k->stat      = 3*(k->shower!=0);
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
    if (j!=i && (i>1 || j>1)) {
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
  if (i<2 && angle==0.0) {
    int j(1-i);
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
      test       = ivec*jvec/(ivec.Abs()*jvec.Abs());
      if (test<=-1.)     th_ex = M_PI; 
      else if (test>=1.) th_ex = 0;   
      else th_ex = acos(test);
      
      angle = Max(angle,th_ex);
      if (IsEqual(angle,M_PI)) angle = M_PI;
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
  
  bool dir(moms[0][3]>0.0);
  if (moms[0].PSpat2()<moms[1].PSpat2()) dir=moms[1][3]<0.0;
  Poincare cms(moms[0]+moms[1]);
  for (int i=0;i<n;++i) knots[i]->part->SetMomentum(cms*knots[i]->part->Momentum());
  
  Poincare zaxis;
  if (dir) zaxis = Poincare(knots[0]->part->Momentum(),Vec4D::ZVEC);
  else zaxis = Poincare(knots[1]->part->Momentum(),Vec4D::ZVEC);
  
  for (int i=0;i<n;++i) knots[i]->part->SetMomentum(zaxis*knots[i]->part->Momentum());
  for (int i=0;i<n;++i) {
    double th = ColourAngle(knots,i);
    knots[i]->thcrit=th;
  }

  for (int i=0;i<n;++i) knots[i]->part->SetMomentum(moms[i]);
  delete [] moms;
}


void Tree_Filler::EstablishRelations(APACIC::Knot * mo,APACIC::Knot * d1,APACIC::Knot * d2,
				     int mode,double x1,double x2,double scale1,double scale2)
{
  if (mode==1) {
    Vec4D p1(m_cms_boost*mo->part->Momentum());
    Vec4D p2(m_cms_boost*d1->part->Momentum());
    Vec4D p3(m_cms_boost*d2->part->Momentum());

    mo->left  = d1;
    mo->right = d2;
    mo->zs    = mo->z = p2[0]/p1[0];
    //     mo->stat  = 0;
    mo->part->SetStatus(part_status::decayed);
    if (mo->part->Info() != 'H') mo->part->SetInfo('f');

    d1->prev  = mo;
    d2->prev  = mo;
    d1->E2    = sqr(p2[0]);
    d2->E2    = sqr(p3[0]);

    p_fsrshower->EstablishRelations(mo,d1,d2);

    mo->tout = mo->t;
    if (m_showermode&1 && mo->prev!=NULL && mo->shower==2) {
      // is mother
      if (mo->prev->dir!=0) mo->t=dabs(mo->prev->right->t);
      // fs mother
      else mo->t=mo->prev->t; 
    }
    return;
  }
  else if (mode==0) {
    // initial state initialization
    //  status:
    //  blob->CMS()                     - Vec4D hard event in LAB system
    //  d1->part->Momentum()              - at the moment also in LAB system
    //  blob->InParticle(0)->Momentum() - in CMS system

    // set x1 and x2
    d1->x = x1;
    d2->x = x2;
    // set start t
    d1->t = -m_scale2factor*scale1;
    d2->t = -m_scale2factor*scale2;
  }
  else if (mode==2 || mode==3) {
    // initial state initialization
    //     mo -> d1 (IS) 
    //        -> d2 (FS)

    if (!mo || !d1 || !d2 ) {
      msg_Error()<<"ERROR in Tree_Filler::EstablishRelations: "<<std::endl
		 <<"   Can not establish relations with less than three elements, abort."<<std::endl;
      abort();
    }
    mo->right = d1;
    mo->left  = d2;
    mo->dir   = d1->dir;
    d1->prev  = mo;
    d2->prev  = mo;

    double t1 = d1->part->Momentum().Abs2();
    double t0 = d1->t;
    mo->t     = t0;
    d1->t     = t1;
    d1->tmax  = t1;
    d2->t     = -t1;
    // account for correct amount of radiation 
    // from is/fs color dipole q->gq
    if (d1->part->Flav().IsGluon()) d2->t=-t0;
    if (mode==2) mo->x = x1;
            else mo->x = x2;

    APACIC::Initial_State_Shower::SetColours(d1);
    
    double th12(d1->part->Momentum().Theta(d2->part->Momentum()));
    double th1((-1.*mo->part->Momentum()).Theta(d1->part->Momentum()));
    double th2((-1.*mo->part->Momentum()).Theta(d2->part->Momentum()));
    if (mo->part->Flav().IsGluon()) {
      mo->thcrit=Max(th1,th2);
      if (d1->part->Flav().IsGluon()) d2->thcrit=Max(th12,th2);
      else d2->thcrit=th2;
    }
    else {
      if (d1->part->Flav().IsGluon()) {
	mo->thcrit=th1;
	d2->thcrit=th12;
      }
      else {
	mo->thcrit=th2;
	d2->thcrit=Max(th12,th2);
      }
    }
    msg_Debugging()<<mo->part->Flav()<<"->"<<d1->part->Flav()<<","
		   <<d2->part->Flav()<<" ("<<th1<<","<<th2<<","
		   <<th12<<") -> "<<mo->thcrit<<","<<d2->thcrit<<std::endl;
  }
}

