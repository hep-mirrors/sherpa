#include"Cluster_Partons.H"
#include"Message.H"

using namespace SHERPA;
using namespace APACIC;
using namespace EXTRAXS;
using namespace AMEGIC;
using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace AORGTOOLS;


bool Cluster_Partons::ClusterConfiguration(Process_Base * proc,Blob * _blob) {
  msg.Debugging()<<"In Cluster_Partons::ClusterConfiguration("<<proc->Name()<<")"<<std::endl;
  if (proc==0) {
    msg.Error()<<" ERROR: in Cluster_Partons::ClusterConfiguration("<<proc<<")"<<std::endl
	       <<"   no process found."<<std::endl;
    abort();
  }

  int nin         = proc->Nin();
  int nout        = proc->Nout();
  Flavour * flavs = proc->Flavs();
  int nampl       = proc->NumberOfDiagrams();

  int    nlegs    = nin + nout;
  Leg ** legs     = 0;

  blob            = _blob;

  // start cluster algorithm :
  if (proc!=lastproc) {
    lastproc = proc;
    if (combi) delete combi;
    combi    = 0;

    // generate a list of "legs" for each amplitude
    legs = new Leg *[nampl];
    for (int k=0;k<nampl;) {
      legs[k] = new Leg[nlegs];
      int l   = 0;
      if (FillLegs(legs[k],proc->Diagram(k),l,nlegs)) {
	// generate a "final legs"-list for all topologies
	msg.Debugging()<<(k+1)<<" graph found:"<<(proc->Diagram(k))<<std::endl;
	++k;
      } 
      else {
	// generate NO "final legs"-list for topologies with four-verticies
	msg.Debugging()<<"  dropping Graph"<<(proc->Diagram(k))<<std::endl;
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
    for (int i=0;i<nlegs;++i) amoms[i] = proc->Momenta()[i];

    combi = new Combine_Table(jf,amoms,0);
    combi->FillTable(legs,nlegs,nampl);   
    ct = combi->CalcJet(nlegs); 
  }
  else {
    // use the existing combination table and determine best combination sheme
    ct = combi->CalcJet(nlegs,proc->Momenta());
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
  msg.Debugging()<<"In Cluster_Partons::CalculateWeight("<<hard<<","<<jet<<")"<<std::endl;
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
//   cout<<" njets="<<njet<<endl;
//   cout<<" nlegs="<<nlegs<<std::endl;
//   if (njet==maxjetnumber) {
//     cout<<" reduced weight!!! "<<endl;
//   }


  weight      = 1.;

  int si = 0;
  std::vector<double> last_q(nlegs,qmax);
  std::vector<int>    last_i(nlegs,si);
  double as_jet  = (*as)(facycut*jet);
  double as_hard = (*as)(hard);

  // remove (for e+e- -> Jets) universal 2 Quark sudakov to increase the effectivity!
  int strong=0;
  for (int l=0; l<nlegs; ++l) {
    if (ct->GetLeg(l)->fl.Strong()) {// Quark sudakov for each strong interacting particle
      /* *AS*
      if (njet==maxjetnumber) {
	msg.Debugging()<<"Multiply weight by (out hard)  "<<1./sud->DeltaQ(qmax,qmin)
		       <<"   "<<qmax<<" --> "<<qmin<<std::endl;
	weight /= sud->DeltaQ(qmax,qmin);
      }
      */
      strong++;
    }    
  }
  // weight with correct alphas at hard scale (if needed).
  if (strong>2) {
    msg.Debugging()<<"Multiply weight by (as hard)"<<pow(as_hard/as_jet,strong-2)<<std::endl;
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
      msg.Debugging()<<"Multiply weight by (in G)  "
		     <<sud->DeltaG(last_q[i],qmin)/sud->DeltaG(facnlly*ptij,qmin)
		     <<"  "<<last_q[i]<<" -> "<<ptij<<"  ("<<qmin<<")"<<std::endl;
    }
    if (ct_down->GetLeg(i)->fl.IsQuark()) {  // Quark sudakov
      weight *= sud->DeltaQ(last_q[i],qmin)/sud->DeltaQ(facnlly*ptij,qmin);
      msg.Debugging()<<"Multiply weight by (in Q)  "
		     <<sud->DeltaQ(last_q[i],qmin)/sud->DeltaQ(facnlly*ptij,qmin)
		     <<"  "<<last_q[i]<<" -> "<<ptij<<"  ("<<qmin<<")"<<std::endl;
    }
    if (ct_down->GetLeg(i)->fl.Strong()) {   // alphaS factor
      weight *= (*as)(facycut*ptij*ptij)/as_jet;
      msg.Debugging()<<"Multiply weight by (as in)  "<<(*as)(facycut*ptij*ptij)/as_jet<<std::endl;
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
	msg.Debugging()<<"Multiply weight by (out G) "<<sud->DeltaG(last_q[l],qmin)
		       <<"  "<<last_q[l]<<" -> "<<qmin<<std::endl;
      }
      if (ct->GetLeg(l)->fl.IsQuark()) {// Quark sudakov
	weight *= sud->DeltaQ(last_q[l],qmin);
	msg.Debugging()<<"Multiply weight by (out Q)"<<sud->DeltaQ(last_q[l],qmin)
		       <<"  "<<last_q[l]<<" -> "<<qmin<<std::endl;
      }    
    }
  }
  msg.Tracking()<<" sudkov weight="<<weight<<std::endl;

  ct = ct_tmp;
}

//int   Cluster_Partons::SetColours(APHYTOOLS::Vec4D * );


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

  knots[0]->part->SetDec(blob);  
  knots[1]->part->SetDec(blob);
  knots[2]->part->SetProd(blob);
  knots[3]->part->SetProd(blob);

  for (int i=0;i<4;i++) {
    for (int j=0;j<2;j++) {
      if (xs) {
	msg.Debugging()<<" [ "<<i<<" ][ "<<j<<" ] = "<<xs->Colours()[i][j]<<std::endl;
	knots[i]->part->SetFlow(j+1,xs->Colours()[i][j]);
      }
      else {
	msg.Debugging()<<" [ "<<i<<" ][ "<<j<<" ] = "<<colors[i][j]<<std::endl;
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

  for (int l=0; l<4; ++l) ini_knots.push_back(knots[l]);

  Tree * tree;
  int nlegs=4;
  int i,j;

  Combine_Table * ct_tmp = ct;

  ct=ct->Up();
  while (ct) {
    knots.push_back(0);ini_knots.push_back(0); ++nlegs;
    ct->GetWinner(i,j);
    msg.Debugging()<<" winner i="<<i<<"  j="<<j<<std::endl;

    for (int l=knots.size()-1;l>j;--l) knots[l] = knots[l-1];
    if (i>=2) tree = fin_tree; 
    else      tree = ini_trees[i];

    Knot * d1, * d2, * mo1, * mo2;
    if (i>=2) {
      d1 = Point2Knot(tree, ct->GetLeg(i), ct->Momentum(i),'H');
      d2 = Point2Knot(tree, ct->GetLeg(j), ct->Momentum(j),'H');

      EstablishRelations(knots[i],d1,d2,1);      
    } 
    else {
      msg.Error()<<"Initial state clustering not fully implemented yet !!"<<std::endl;
      abort();
    }

    // *AS*
    knots[i] = d1;
    knots[j] = d2;


    /*
      //========================================
      // determine mothers for start conditions
      //========================================
      if (knots[i]->part->Flav().IsBoson()) {
      if (knots[i]->part->Flav().IsGluon()) {
      // g -> g/q   g/q
      // i    d1    d2
      //     knots[i] prod[i] or vice verca	  .... choose the more energetic;
      if (d1->part->Momentum()[0] > d2->part->Momentum()[0]) {
      mo1 = ini_knots[i]; mo2 = knots[i];
      }
      else {
      mo1 = knots[i]; mo2 = ini_knots[i];
      }
      } 
      else {
      // ph/Z/W/h -> ...
      mo1 = knots[i]; mo2 = knots[i];
      }
      } 
      else if (knots[i]->part->Flav().IsFermion()) {
      if (d1->part->Flav().IsFermion()) {
      // f -> f + g/ph/Z
      mo1 = ini_knots[i]; mo2 = knots[i];
      }
      else if (d2->part->Flav().IsFermion()) {
      // f -> g/ph/Z + f
      mo1 = knots[i]; mo2 = ini_knots[i];
      }
      else {
      msg.Error()<<" ERROR: unexpected branch in Interface_Tools::FillTrees() "<<std::endl
      <<"        fermion decay without a fermion?! "<<std::endl;
      abort();
      } 
      }
      else {
      // should never happen!
      msg.Error()<<" ERROR: unexpected branch in Interface_Tools::FillTrees() "<<std::endl
      <<"        mother not a boson and not a fermion ?"<<std::endl;
      abort();
      } 
      
      knots[i]       = d1;
      ini_knots[i]   = mo1;
      for (int l=j+1; l<nlegs; ++l ) {
      knots[l]     = knots[l-1];
      ini_knots[l] = ini_knots[l-1];
      }
      knots[j]       = d2;
      ini_knots[j]   = mo2;
    */

    ct = ct->Up();
  }
      
  /*
    // set start scale and crit angle of outgoing particle according to ini_knots:
    for (int l=0; l<nlegs; ++l ) {
    knots[l]->t      = ini_knots[l]->t;
    knots[l]->thcrit = ini_knots[l]->thcrit;
    }
  */
    
  ct = ct_tmp;

  msg.Tracking()<<" ISR Leg 0 : "<<std::endl;
  msg.Tracking()<<ini_trees[0]<<std::endl;
  msg.Tracking()<<" ISR Leg 1 : "<<std::endl;
  msg.Tracking()<<ini_trees[1]<<std::endl;
  msg.Tracking()<<" FSR Tree : "<<std::endl;
  msg.Tracking()<<fin_tree<<std::endl;

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

void Cluster_Partons::EstablishRelations(Knot * mo, Knot * d1,Knot * d2,bool mode)
{
  if (mode) {
    // fin ...
    mo->left  = d1;
    mo->right = d2;
    mo->z     = d1->part->Momentum()[0]/mo->part->Momentum()[0];
    mo->stat  = 0;
    mo->part->SetStatus(2);
    if (mo->part->Info() != 'H') mo->part->SetInfo('f');

    d1->prev  = mo;
    d2->prev  = mo;

    Vec3D vec1 = Vec3D(d1->part->Momentum());
    Vec3D vec2 = Vec3D(d2->part->Momentum());
    //    double th_old  = acos(vec1*vec2/(vec1.Abs()*vec2.Abs()));
    double t_mo = mo->part->Momentum().Abs2();
    double E_mo= mo->part->Momentum()[0];

    // using shower variables
    //    double th  = 
    //    cout<<" th_old "<<th_old<<std::endl;
    //    cout<<" th="<<sqrt( t_mo/(mo->z*(1.- mo->z)*sqr(E_mo) ))<<std::endl;
    double th  = sqrt( t_mo/(mo->z*(1.- mo->z)))/E_mo;

    // thcrit : approximation (for small t):
    if (mo->part->Flav().Strong()) {
      if ((d1->part->Flav().Strong()) && (d2->part->Flav().Strong())) {
	if ((d1->E2) > (d2->E2)) {
	  d1->t      = mo->t;
	  d1->thcrit = mo->thcrit;
	  d2->t      = t_mo;
	  d2->thcrit = th;
	}
	else {
	  d1->t      = t_mo;
	  d1->thcrit = th;
	  d2->t      = mo->t;
	  d2->thcrit = mo->thcrit;
	}
      }
      else if (!(d1->part->Flav().Strong()) && (d2->part->Flav().Strong())) {
	d1->t      = t_mo;
	d1->thcrit = M_PI;
	d2->t      = mo->t;
	d2->thcrit = mo->thcrit;
      }
      else if ((d1->part->Flav().Strong()) && !(d2->part->Flav().Strong())) {
	d1->t      = mo->t;
	d1->thcrit = mo->thcrit;
	d2->t      = t_mo;
	d2->thcrit = M_PI;
      }
    }
    else {
      if ((d1->part->Flav().Strong()) && (d2->part->Flav().Strong())) {
	d1->t      = t_mo;
	d1->thcrit = th;
	d2->t      = t_mo;
	d2->thcrit = th;
      }
      else {
	d1->t      = t_mo;
	d1->thcrit = M_PI;
	d2->t      = t_mo;
	d2->thcrit = M_PI;
      }
    }
    mo->t = t_mo;

    msg.Events()<<"Established relations (FS) : "<<std::endl
		   <<*d1->prev<<std::endl<<*d1<<std::endl<<*d2<<std::endl<<std::endl;
    return;
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
  msg.Debugging()<<"Check if "<<a->Flav()<<" and "<<b->Flav()<<" are connected."<<std::endl;
  return (( (a->GetFlow(1)!=0) && ( (a->GetFlow(1)==b->GetFlow(1)) || 
				 (a->GetFlow(1)==b->GetFlow(2)))  ) ||
	  ( (a->GetFlow(2)!=0) && ( (a->GetFlow(2)==b->GetFlow(2)) ||
				 (a->GetFlow(2)==b->GetFlow(1)))  )    );
}

double Cluster_Partons::ColourAngle(const std::vector<Knot *> & knots, const int i) {
  if (!((knots[i]->part->Flav()).Strong())) return M_PI;

  double angle = 0.;
  for (int j=0;j<knots.size();++j) {
    if (j!=i) {
      if (IsColourConnected(knots[i]->part,knots[j]->part)) {
	Vec3D ivec = Vec3D(knots[i]->part->Momentum());
	Vec3D jvec = Vec3D(knots[j]->part->Momentum());

	angle = Max(angle,acos(ivec*jvec/(ivec.Abs()*jvec.Abs())));
      }
    }
  }
  msg.Debugging()<<"Angle = "<<angle<<std::endl;
  return angle;
}


