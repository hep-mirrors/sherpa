#include "Final_State_Shower.H"


#include "Primitive_Analysis.H"
#include "Shower_Observables.H"

using namespace APACIC;
using namespace AMATOOLS;
using namespace APHYTOOLS;
using namespace AORGTOOLS;

//-----------------------------------------------------------------------
//--------------------------- Constructors ------------------------------
//----------------------------------------------------------------------- 

Final_State_Shower::Final_State_Shower() 
{
  kin  = new Timelike_Kinematics();
  sud  = new Timelike_Sudakov(kin);
};

Final_State_Shower::~Final_State_Shower() 
{
  if (sud) delete sud;
  if (kin) delete kin;
}

//-----------------------------------------------------------------------
//----------------------- Performing the Shower -------------------------
//----------------------------------------------------------------------- 

int Final_State_Shower::PerformShower(Tree * tree,bool jetveto) 
{
  msg.Debugging()<<"----------------------------------------------------------"<<std::endl
		 <<"Final_State_Shower::PerformShower "<<std::endl;

  kin->SetJetVeto(jetveto);

  int stat=InitializeJets(tree,tree->GetRoot());

  /*  
  int statb=ExtraJetCheck();
  if (stat==1 && statb==0)
    stat=3;
  */
  if (stat) {
    msg.Tracking()<<" Now DoKinematics "<<std::endl;
    if (!kin->DoKinematics(tree->GetRoot())) {
      msg.Error()<<"Error in Final_State_Shower : "<<std::endl
		 <<"Final_State_Shower::PerformShower : "
		 <<"Kinematics did not work out."<<std::endl;
      std::cout<<tree<<std::endl;

      return 0;
    }
    int statb=ExtraJetCheck();
    if (stat==1 && statb==0)
      stat=3;

    return stat;
  }
  else {
    msg.Error()<<"Error in Final_State_Shower : "<<std::endl
	       <<"Final_State_Shower::PerformShower : "
	       <<"Initializing the jets did not work out !!!"<<std::endl;
    return 0;
  }
}

void Final_State_Shower::FirstTimelikeFromSpacelike(Tree * tree,Knot* mo,bool jetveto)
{
  msg.Debugging()<<"-----------------------------------------"<<std::endl
		 <<"Final_State_Shower::FirstTimelikeFromSpacelike for Knot "
		 <<mo->kn_no<<":"<<std::endl
		 <<"    Knot has E2/t : "<<mo->E2<<"/"<<mo->t<<std::endl;
  kin->SetJetVeto(jetveto);

  Flavour flavs[2];
  for (;;) {
    if (sud->Dice(mo)) {
      flavs[0]  = sud->GetFlB();
      flavs[1]  = sud->GetFlC();
      InitDaughters(tree,mo,flavs,1);
      if (EvolveJet(tree,mo)) return;
    }
    msg.Debugging()<<"     no branch found ! Sud->Dice yielded a zero result."<<std::endl; 
    
    Reset(mo);
    return;
  }
}

//-----------------------------------------------------------------------
//---------------------------- After the Shower -------------------------
//----------------------------------------------------------------------- 

bool Final_State_Shower::SetColours(Knot * mo)
{
  if (!mo) {
    msg.Error()<<"ERROR in Final_State_Shower::SetColours()"<<std::endl
	       <<"    void mother knot."<<std::endl;
    return 0;
  }

  if ( (!mo->left) && (!mo->right) ) return 1;
  if ( !mo->left ) {
    mo->right->part->SetFlow(1,mo->part->GetFlow(1));
    mo->right->part->SetFlow(2,mo->part->GetFlow(2));
    return 1;
  }
  if ( !mo->right ) {
    mo->left->part->SetFlow(1,mo->part->GetFlow(1));
    mo->left->part->SetFlow(2,mo->part->GetFlow(2));
    return 1;
  }

  Knot * d1 = mo->left;
  Knot * d2 = mo->right;
  Knot * partner, * nopart;
  if (mo->part->Flav().Strong()) {
    if (mo->part->Flav().IsQuark()) {
      partner = d1; nopart = d2;
      if (d2->part->Flav().IsQuark()) {
	partner = d2;
	nopart  = d1;
      }
      if ((partner->part->Flav().IsQuark()) && (partner->part->Flav().IsAnti()) && 
	  (nopart->part->Flav().IsGluon())) {
	partner->part->SetFlow(2,-1);
	nopart->part->SetFlow(1,partner->part->GetFlow(2));
	nopart->part->SetFlow(2,mo->part->GetFlow(2));
      } 
      if ((partner->part->Flav().IsQuark()) && (!partner->part->Flav().IsAnti()) && 
	  (nopart->part->Flav().IsGluon())) {
	partner->part->SetFlow(1,-1);
	nopart->part->SetFlow(1,mo->part->GetFlow(1));
	nopart->part->SetFlow(2,partner->part->GetFlow(1));
      } 
      if ( (partner->part->Flav().IsQuark()) && (!(nopart->part->Flav().Strong())) ) {
	partner->part->SetFlow(1,mo->part->GetFlow(1));
	partner->part->SetFlow(2,mo->part->GetFlow(2));
      }
    } 
    else if (mo->part->Flav().IsGluon()) {
      if (mo->prev) {
	if ( (d1->part->Flav().IsQuark()) && (d2->part->Flav().IsQuark())) {
	  if (d1->part->Flav().IsAnti()) {
	    d1->part->SetFlow(2,mo->part->GetFlow(2));
	    d2->part->SetFlow(1,mo->part->GetFlow(1));
	  }
	  else if (d2->part->Flav().IsAnti()) {
	    d2->part->SetFlow(2,mo->part->GetFlow(2));
	    d1->part->SetFlow(1,mo->part->GetFlow(1));
	  }
	}
	else if ( (d1->part->Flav().IsGluon()) && (d2->part->Flav().IsGluon())) {
	  Knot * au = mo->prev->left;
	  if (mo->prev->left == mo) au = mo->prev->right;

	  partner = d1; nopart = d2;
	  if (kin->ArrangeColourPartners(au,d1,d2)) { partner = d2; nopart = d1; }
	  for (int i=1;i<3;i++) {
	    if (au->part->GetFlow(i) == mo->part->GetFlow(3-i)) {
	      partner->part->SetFlow(3-i,mo->part->GetFlow(3-i));
	      partner->part->SetFlow(i,-1);
	      nopart->part->SetFlow(3-i,partner->part->GetFlow(i));
	      nopart->part->SetFlow(i,mo->part->GetFlow(i));
	      break;
	    }
	  }
	}
      }
    }
    else {
      msg.Error()<<"ERROR in Final_State_Shower::SetColours()"<<std::endl
		 <<"    case not covered : colourful "<<mo->part->Flav()
		 <<" -> "<<d1->part->Flav()<<" + "<<d2->part->Flav()<<std::endl;
      return 0;
    }
  }
  else {
    // colour neutral mother
    if ((d1->part->Flav().Strong()) && (d2->part->Flav().Strong())) {      
      if ( (d1->part->Flav().IsQuark()) && (d2->part->Flav().IsQuark())) {
	partner = d1; nopart = d2;
	if ((d1->part->Flav().IsQuark()) && (d1->part->Flav().IsAnti()) && 
	    (d2->part->Flav().IsQuark()) && (!d2->part->Flav().IsAnti())) {
	  partner = d2; nopart = d1;
	}
	partner->part->SetFlow(1,-1);
	partner->part->SetFlow(2,0);
	nopart->part->SetFlow(2,partner->part->GetFlow(1));
	nopart->part->SetFlow(1,0);
      }
      else if ( (d1->part->Flav().IsGluon()) && (d2->part->Flav().IsGluon())) {
	d1->part->SetFlow(1,-1);
	d1->part->SetFlow(2,-1);
	d2->part->SetFlow(2,d1->part->GetFlow(1));
	d2->part->SetFlow(1,d1->part->GetFlow(2));
      }
      else {
	msg.Error()<<"ERROR in Final_State_Shower::SetColours()"<<std::endl
		   <<"    case not covered : colourless -> "
		   <<d1->part->Flav()<<" + "<<d2->part->Flav()<<std::endl;
	return 0;
      }
    }
  }

  msg.Events()<<"SetColours : "
	      <<mo->part->Flav()<<" ("<<mo->part->GetFlow(1)<<" "<<mo->part->GetFlow(2)<<")"
	      <<d1->part->Flav()<<" ("<<d1->part->GetFlow(1)<<" "<<d1->part->GetFlow(2)<<")"
	      <<d2->part->Flav()<<" ("<<d2->part->GetFlow(1)<<" "<<d2->part->GetFlow(2)<<")"<<std::endl;

  return ( SetColours(d1) && SetColours(d2) );
}


void Final_State_Shower::ExtractPartons(Knot * kn,Blob * jet,Blob_List * bl,Parton_List * pl) 
{
  if (!kn) return;
  msg.Debugging()<<"----------------------------------------------------------"<<std::endl
		 <<"Final_State_Shower::ExtractPartons for Knot "<<kn->kn_no<<std::endl;
  if (kn->part->Info() == 'H') {
    /* 
       New jet : kn = hard parton from ME info = 'HF'
                 and kn outgoing
		 or kn->left or kn->right not from ME
    */
    if (!(kn->left)) {
      pl->push_back(kn->part);
      jet = new Blob();
      jet->AddToInPartons(new Parton(kn->part));
      jet->AddToOutPartons(new Parton(kn->part));
      jet->SetId(bl->size());
      jet->SetType(std::string("FS Parton Shower (APACIC++2.0)"));
      jet->SetPosition(kn->part->XProd() + Vec4D(kn->part->LifeTime(),kn->part->Distance()));
      bl->push_back(jet);
      return;
    }
    else {
      if ((kn->left->part->Info() != 'H') || (kn->right->part->Info() != 'H')) {
	jet = new Blob();
      	kn->part->SetDec(jet);
      	kn->part->SetStatus(2);
	pl->push_back(kn->part);
	jet->AddToInPartons(new Parton(kn->part));
	jet->SetId(bl->size());
	jet->SetType(std::string("FS Parton Shower (APACIC++2.0)"));
	jet->SetPosition(kn->part->XProd() + Vec4D(kn->part->LifeTime(),kn->part->Distance()));
	bl->push_back(jet);
      }
    }
  }
  else {
    if (!kn->left) {
      if (!jet) {
	msg.Error()<<"ERROR in Final_State_Shower ::ExtractPartons :"<<std::endl
		   <<" No jet for Parton : "<<kn->part->Number()<<std::endl;
	abort();
      }
      kn->part->SetNumber(pl->size());
      kn->part->SetProd(jet);
      kn->part->SetStatus(1);
      pl->push_back(kn->part);
      jet->AddToOutPartons(new Parton(kn->part));
    }
  }
  ExtractPartons(kn->left,jet,bl,pl); 
  ExtractPartons(kn->right,jet,bl,pl); 
}


void Final_State_Shower::ExtractPartons(Knot * kn, Parton_List * pl) 
{
  if (!kn) 
    return;
  if (kn->left==0) {
    pl->push_back(new Parton(kn->part));
    return;
  }

  ExtractPartons(kn->left,pl); 
  ExtractPartons(kn->right,pl); 
}

//-----------------------------------------------------------------------
//---------------------------- Helpers ----------------------------------
//----------------------------------------------------------------------- 

bool Final_State_Shower::TestShower(Tree * tree) 
{
  bool do_ana=1;

  //----------------------------------------
  if (do_ana) {
    // Initialise Histos
  }
  Primitive_Analysis ana;
  //  ana.AddObservable(new Shower_Observables(11,1.e-6,1.,180,0));

  //  Jetrates * sobs =new Jetrates(11,1.e-6,1.,60,0);
   Jetrates * sobs =new Jetrates(11,1.e-6,1.,180,0);
   ana.AddObservable(sobs);
   ana.AddObservable(new Multiplicity(00,-0.5,50.5,51,0));
  //----------------------------------------

  int number;
  double E2 = sqr(rpa.gen.Ecms());
  double E  =rpa.gen.Ecms()*0.5;




  //  Ran.ReadInStatus("RandomA.dat",20);
  //  Ran.ReadInStatus("RandomB.dat",23);
  //  Ran.ReadInStatus("RandomC.dat",31); // nan Daughters
  //  Ran.ReadInStatus("RandomD.dat",35); // nan 
  //  Ran.ReadInStatus("RandomE.dat",339); // ycut!
  for (int i=1;i<=rpa.gen.NumberOfEvents();i++) {
    if (i%2500==0) {
      msg.Out()<<" "<<i<<" th event "<<std::endl;
      //      msg.Out()<<" ran"<<Ran.WriteOutStatus("RandomE.dat")<<std::endl;
    }
   
    msg.Events()<<"++++++++++++++++++++++++++++"<<i<<" th event +++++++++++++++++++++++++"<<std::endl;
    tree->Reset();
    InitTwojetTree(tree,E2);
    if (!PerformShower(tree,0)) return 0;     //*AS*   (tree,0) = no jetveto! (default)
    if (rpa.gen.Tracking()) OutputTree(tree);
    if (do_ana) {
      // Fill Histos
      Parton_List pl;
      pl.push_back(new Parton(0,Flavour(kf::e),Vec4D(E,0,0,E)));
      pl.push_back(new Parton(1,Flavour(kf::e).Bar(),Vec4D(E,0,0,-E)));
      ExtractPartons(tree->GetRoot(),&pl);
      for (int k=0; k<pl.size(); ++k) {
	msg.Events()<<pl[k];
      }
      // Cluster
      ana.DoAnalysis(pl,1.);


      for (int i=0; i<pl.size();++i) {
	delete pl[i];
      }
      pl.clear();

      /*
      if (sobs->ymax>0.2) {
	std::cout<<" ymax="<<sobs->ymax<<std::endl;
	std::cout<<" Tree:"<<tree<<std::endl;
	abort();
      }
      */
    }
    

  }
  ana.FinishAnalysis("testout_shower_FB",0);
  msg.Events()<<"Final_State_Shower::TestShower : "
	      <<"Terminated loops over events successfully."<<std::endl;
  return 1;
}

//-----------------------------------------------------------------------
//------------------- Initialisation of the Shower ----------------------
//----------------------------------------------------------------------- 

int Final_State_Shower::InitializeJets(Tree * tree,Knot * mo)
{
  if (!mo) {
    msg.Error()<<"Error in Final_State_Shower : "<<std::endl
	       <<"Final_State_Shower::InitializeJets : No mother found !"<<std::endl;
    return 0;
  }
  if ((!mo->left) || (!mo->right)) {
    msg.Error()<<"Error in Final_State_Shower : "<<std::endl
	       <<"Final_State_Shower::InitializeJets : No daughters found !"<<std::endl
	       <<" for "<<*mo<<std::endl;
    return 0;
  }


  double z_tmp = mo->z;
  Knot * d1    = mo->left;
  Knot * d2    = mo->right;
  Vec4D p1_tmp = d1->part->Momentum();
  Vec4D p2_tmp = d2->part->Momentum();

  bool decay1 = (d1->stat>0),decay2 = (d2->stat>0);

  // first==2  means "pythia kinematics"
  int first=1;
 //  
  if (mo==tree->GetRoot() && decay1 && decay2) first=2;  // two jet event


  msg.Tracking()<<"Final_State_Shower::InitializeJets : Try to evolve daughters :"<<std::endl
		<<"    "<<d1->kn_no<<" "<<d1->part->Flav()<<" "<<decay1<<" "<<d1->t<<" / "
		<<d2->kn_no<<" "<<d2->part->Flav()<<" "<<decay2<<" "<<d2->t<<std::endl;

  if ((decay1) || (decay2)) {
    bool accept = SmearDaughters(mo);
    for (;;) {
      if (FillBranch(tree,mo,first)) {
	accept = 1;
	if (decay1) accept = EvolveJet(tree,d1) && accept;
	if (decay2) accept = EvolveJet(tree,d2) && accept;
	if (accept) break; //return 1;
	else {
	  mo->z  = z_tmp;
	  d1->E2 = mo->z*mo->z*mo->E2;
	  d2->E2 = (1.-mo->z)*(1.-mo->z)*mo->E2;
	  d1->part->SetMomentum( p1_tmp );
	  d2->part->SetMomentum( p2_tmp );
	  msg.Debugging()<<"Final_State_Shower::InitializeJets : "
			 <<"      Could not evolve the first branch further."<<std::endl
			 <<"      Reset daughters on ME kinematics."<<std::endl;
	  // should be set on their t_out's !!!! ? 


	  // delete grand daughters
	  /*
	  d1->left=0;
	  d1->right=0;
	  d2->left=0;
	  d2->right=0;
	  if (d1->stat!=3 && d1->t>d1->tout) d1->stat=1;
	  if (d2->stat!=3 && d2->t>d2->tout) d2->stat=1;
	  */
	}
      }
      if (!(d1->stat) && !(d2->stat)) break;
    }
  }

  int ok=1; int ej=0;
  if (!decay1) { 
    int ok1=InitializeJets(tree,d1);
    if (ok1==3) ej=3;
    ok = ok && ok1;
  }
  else {
    ini_partons.push_back(d1);
  }
  if (!decay2) {
    int ok2=InitializeJets(tree,d2);
    if (ok2==3) ej=3;
    ok = ok && ok2;
  }
  else {
    ini_partons.push_back(d2);
  }

  // *AS*  if (!ExtraJetCheck(mo,d1,d2)) ej=3;

  if (ok==0) return 0;  // did not work out!
  if (ej==3) return 3;  // jetnumber reduced after shower! need new kinematics from (same) ME!
  return 1;             // everything ok.
}

bool  Final_State_Shower::ExtraJetCheck(Knot * mo, Knot * d1, Knot * d2) {
  return kin->ExtraJetCheck(0,d1,d2); // *AS* test E2 dependence
}

bool  Final_State_Shower::ExtraJetCheck() {
  bool test=1;
  //  cout<<" (A) "<<ini_partons.size()<<endl;
  for (int i=0;i<ini_partons.size()-1;++i) {
    for (int j=i+1;j<ini_partons.size();++j) {
      test=test & kin->ExtraJetCheck(0,ini_partons[i],ini_partons[j]);
      if (test==0) break;
    }
    if (test==0) break;
  }
  ini_partons.clear();
  return test;
}

bool Final_State_Shower::SmearDaughters(Knot * mo) {
  Knot * d1 = mo->left;
  Knot * d2 = mo->right;
  if ((d1->part->Flav().IsStable()) && (d2->part->Flav().IsStable())) return 0;
  if ((d1->stat==0) && (d2->stat==0)) return 0;

  double smax    = mo->t;

  double mass1   = (d1->part->Flav()).PSMass();
  double width1  = (d1->part->Flav()).Width();
  double mass12  = mass1*mass1;
  double upper1  = (smax-mass12)/mass1/width1;
  double lower1  = -mass1/width1;
  double ymin1   = atan(lower1);
  double yrange1 = atan(smax/mass1/width1/(1.+lower1*upper1));
  if (lower1*upper1<-1.) {
    if (upper1>0) yrange1 = yrange1 + M_PI;
    if (upper1<0) yrange1 = yrange1 - M_PI;
  }     

  double mass2   = (d2->part->Flav()).PSMass();
  double width2  = (d2->part->Flav()).Width();
  double mass22  = mass2*mass2;
  double upper2  = (smax-mass22)/mass2/width2;
  double lower2  = -mass2/width2;
  double ymin2   = atan(lower2);
  double yrange2 = atan(smax/mass2/width2/(1.+lower2*upper2));
  if (lower2*upper2<-1.) {
    if (upper2>0) yrange2 = yrange2 + M_PI;
    if (upper2<0) yrange2 = yrange2 - M_PI;
  }     

  double t1 = mass12;
  double t2 = mass22;
  for(;;) {
    if (!(d1->part->Flav().IsStable()) && (d1->stat))
      t1 = mass12 + mass1*width1 * tan(ran.Get()*yrange1 + ymin1);
    if (!(d2->part->Flav().IsStable()) && (d2->stat))
      t2 = mass22 + mass2*width2 * tan(ran.Get()*yrange2 + ymin2);
    if ( t1+t2+sqrt(2.*t1*t2) < (mo->t)) {
      msg.Debugging()<<"Final_State_Shower::SmearDaughters :"<<std::endl
		     <<"    Shifted t of "<<d1->part->Flav()
		     <<" from "<<d1->tout<<" to "<<t1<<std::endl
		     <<"    Shifted t of "<<d2->part->Flav()
		     <<" from "<<d2->tout<<" to "<<t2<<std::endl
		     <<"    Mother :"<<mo->t<<std::endl
		     <<"    Knot 1 :"<<mass1<<", "<<width1
		     <<" -> "<<yrange1<<", "<<ymin1<<std::endl
		     <<"    Knot 2 :"<<mass2<<", "<<width2
		     <<" -> "<<yrange2<<", "<<ymin2<<std::endl;
      d1->tout = t1;
      d2->tout = t2;
      return 1;
    }
  }
}

void Final_State_Shower::InitTwojetTree(Tree * tree,double scale) {
  double start_th=200;

  Knot * mo   = tree->NewKnot();
  *(mo->part) = Parton(1,Flavour(kf::photon),Vec4D(sqrt(scale),0,0,0));
  mo->part->SetStatus(2);
  mo->part->SetInfo('M');
  mo->t       = scale;
  mo->E2      = scale;
  mo->maxpt2  = 0.;
  mo->z       = 0.5;
  mo->costh   = -1.; 
  mo->thcrit  = start_th;
  mo->phi     = 2.*M_PI*ran.Get();
  mo->stat    = 1;  // *AS* ? 0

  //  The daughters are defined and distributed according to
  //  1+cos^2 theta w.r.t. the z-axis.

  Flavour mo_flavs[2];
  for(;;) {
    //    mo_flavs[0] = Flavour(kf::code(1+int(ran.Get()*4.)));   // *AS* default *5.
    //    mo_flavs[0] = Flavour(kf::b);
    mo_flavs[0] = Flavour(kf::d);
    if (4.*sqr(mo_flavs[0].PSMass()) < scale) break;
  }
  mo_flavs[1] = mo_flavs[0].Bar();

  double E    = 0.5 * sqrt(scale);
  double p    = sqrt(E*E-sqr(mo_flavs[0].PSMass()));
  double cth;
  for (;;) {
    cth = 1.-2.*ran.Get();
    if (1.+cth*cth < 2.*ran.Get()) break;
  }
  double sth  = sqrt(1.-cth*cth);
  Vec3D p1(p*sth*cos(mo->phi),p*sth*sin(mo->phi),p*cth);

  mo->left             = tree->NewKnot();
  mo->left->prev       = mo;
  mo->left->stat       = 3;    
  *(mo->left->part)    = Parton(2,mo_flavs[0],Vec4D(E,p1));
  mo->left->part->SetStatus(1);
  mo->left->part->SetInfo('H');
  mo->left->part->SetFlow(1,-1);
  mo->left->t          = mo->t;
  mo->left->tout       = sqr(mo_flavs[0].PSMass());
  //  std::cout<<" da set on="<<mo->left->tout<<std::endl;

  mo->left->E2         = E*E;
  mo->left->thcrit     = start_th;
 
  mo->right            = tree->NewKnot();
  mo->right->prev      = mo;
  mo->right->stat      = 3;     
  *(mo->right->part) = Parton(3,mo_flavs[1],Vec4D(E,(-1.)*p1)); 
  mo->right->part->SetStatus(1);
  mo->right->part->SetInfo('H');
  mo->right->part->SetFlow(2,mo->left->part->GetFlow(1));
  mo->right->t         = mo->t;
  mo->right->tout      = sqr(mo_flavs[1].PSMass());
  //  std::cout<<" da set on="<<mo->right->tout<<std::endl;
  mo->right->E2        = E*E;
  mo->right->thcrit    = start_th;
};

//-----------------------------------------------------------------------
//------------------------ Evolution of the Shower ----------------------
//----------------------------------------------------------------------- 

bool Final_State_Shower::EvolveJet(Tree * tree,Knot* mo)
{
  if (!mo->stat) {
    // it is not possible to dice this knot further down
    // delete already known daughters, and set mother on its "t_out"
    ResetDaughters(mo);
    // Q: When is it neccessary to reset daughters?
    return 1;
  }

  // save old status
  double z_tmp = mo->z;

  msg.Debugging()<<"Final_State_Shower::EvolveJet for ("<<mo->kn_no<<"), "
		 <<"      "<<mo->t<<", "<<mo->E2<<", "<<mo->z<<std::endl
		 <<"      "<<mo->left->kn_no<<" "<<mo->left->part->Flav()<<" / "
		 <<mo->right->kn_no<<" "<<mo->right->part->Flav()<<std::endl;

  // start evolution
  bool accept1 ,accept2;
  for (;;) {
    if (FillBranch(tree,mo,0)) {
      // select one/two daughter systems and try to evolve them 
      accept1 = EvolveJet(tree,mo->left);
      accept2 = EvolveJet(tree,mo->right);
      if (accept1 && accept2) return 1; // successful
      else {
	// reset to original (massless) kinematics
	mo->z         = z_tmp;
	mo->left->E2  = mo->z*mo->z*mo->E2;
	mo->right->E2 = (1.-mo->z)*(1.-mo->z)*mo->E2;
      }      
      // try again until daughters evolve successfully, or ...
    }
    else {
      // ... no new daughter system could be found

      // *AS* if you want to disable "second chance" change here
      bool enable2ndchance=1;
      if (enable2ndchance) {
	mo->left  = 0;  
	mo->right = 0;
	mo->stat  = 3;  // means, it needs to be diced again
      }
      else {
	Reset(mo);
      }
      return 0; // fail
    }
  }
}

bool Final_State_Shower::FillBranch(Tree * tree,Knot* mo,int first)
{
  Knot * d1 = mo->left;
  Knot * d2 = mo->right;

  msg.Debugging()<<"Final_State_Shower::FillBranch for ("<<mo->kn_no<<"), "
		 <<"      "<<mo->t<<", "<<mo->E2<<", "<<mo->z<<std::endl
		 <<"      ("<<d1->kn_no<<") "<<d1->part->Flav()<<" st"<<d1->stat<<" /  ("
		 <<d2->kn_no<<") "<<d2->part->Flav()<<" st"<<d2->stat<<std::endl;
  

  if (!(first) && (mo->t <= mo->tout) ) return 0; // failed


  Flavour d1_flavs[2];
  Flavour d2_flavs[2];
  
  bool do12,accept;
  bool diced1=0, diced2=0;

  double z_tmp = mo->z;

  for (;;) {
    Knot * g =0;
    // select one daugther to dice down
    do12 = ChooseDaughter(mo);
    if (first==2) g=mo;
    //    g=mo;
    if (!do12) {
      ResetDaughters(d1);  // if grand children already determined, delete them
      if (sud->Dice(d1,g)) { // determine t,z, and flavours
	d1_flavs[0]  = sud->GetFlB();
	d1_flavs[1]  = sud->GetFlC();
	d1->stat=1;   // *AS*
	diced1=1;
      }   
      else Reset(d1);
    }
    else {
      ResetDaughters(d2);  // if grand children already determined, delete them
      if (sud->Dice(d2,g)) { // determine t,z, and flavours
	d2_flavs[0]  = sud->GetFlB();
	d2_flavs[1]  = sud->GetFlC();
	d2->stat=1;   // *AS*
	diced2=1;
      }    
      else Reset(d2);
    }

    // check zrange first
    if (first==2) 
      kin->CheckZRange(mo);
    //    std::cout<<"b st"<<d1->stat<<" st"<<d2->stat<<std::endl;

      //    std::cout<<"a st"<<d1->stat<<" st"<<d2->stat<<std::endl;
    if ((d1->stat != 3) && (d2->stat != 3)) { // *AS*
      if (kin->Shuffle(mo,first)) { // check (and adjust) kinematics 
                                  //  (usually both children have to be diced atleast once)
	// if sucessfull set all grand children as determined above and exit
	if (d1->stat) {
	  msg.Debugging()<<"InitDaughters for "<<d1->kn_no<<" "<<d1->part<<std::endl;
      	  if (d1->left) { 
	    msg.Debugging()<<"InitDaughters for "<<d1->kn_no<<" "
			<<std::endl<<d1->left->part<<std::endl<<d1->right->part<<std::endl;
	  }
	  InitDaughters(tree,d1,d1_flavs,diced1); 
	  diced1=0;
	}
	if (d2->stat) {
	  msg.Debugging()<<"InitDaughters for "<<d2->kn_no<<" "<<d2->part<<std::endl;
	  if (d2->left) { 
	    msg.Debugging()<<"InitDaughters for "<<d2->kn_no<<" "
			<<std::endl<<d2->left->part<<std::endl<<d2->right->part<<std::endl;
	  }
	  InitDaughters(tree,d2,d2_flavs,diced2); 
	  diced2=0;
	}
	return 1;  // successfull
      }
      else 
	msg.Debugging()<<"Shuffle for "<<mo->kn_no<<" failed."<<std::endl;
    }
    msg.Debugging()<<"Final_State_Shower::FillBranch"<<std::endl
		   <<"      Daughters did not fit, "<<d1->t<<"("<<d1->tout<<")"
		   <<" / "<<d2->t<<"("<<d2->tout<<")"<<std::endl;
 
    // exit if both daughters can not be diced down further
    if (!(d1->stat) && !(d2->stat)) return 0; 
  }
}

//-----------------------------------------------------------------------
//------------------- Service for the Branchings ------------------------
//----------------------------------------------------------------------- 

bool Final_State_Shower::ChooseDaughter(Knot * mo)
{
  if (rpa.gen.Debugging()) {
    Knot * d1 = mo->left;
    Knot * d2 = mo->right;
  }

  if (mo->left->stat==3)  return 0;
  if (mo->right->stat==3) return 1;

  if ((mo->left->stat) && 
      ((!mo->right->stat) || 
       (mo->right->t < mo->right->tout) )) return 0; // left leg
  if ((mo->right->stat) && 
      ((!mo->left->stat) || 
       (mo->left->t < mo->left->tout) )) return 1; // right leg

  
  double tm1   = Min(mo->t,mo->left->E2);
  double tm2   = Min(mo->t,mo->right->E2);
  // *AS* test:
  // *AS*   if ((mo->left->t) > (mo->right->t)) return 1; // right leg
  // *AS*   return 0; // right leg
 
  if ((mo->left->t/tm1) > (mo->right->t/tm2)) return 0; // left leg
  return 1; // right leg
 
}

void Final_State_Shower::InitDaughters(Tree * tree,Knot * mo,Flavour * mo_flavs, bool diced) 
{ 
  if (!mo->left) {
    // Initialize new knots, provide them with flavours and link them
    mo->left        = tree->NewKnot();

    mo->right       = tree->NewKnot();
  }
  if (diced) {
    mo->left->prev  = mo;
    mo->left->part->SetFlav(mo_flavs[0]);
    mo->left->part->SetInfo('F');
    mo->left->part->SetStatus(1);
    mo->left->tout    = sqr(mo_flavs[0].PSMass());
    //    std::cout<<" res da set on="<<mo->left->tout<<std::endl;
    mo->left->stat    = 3;  // *AS*

    mo->right->prev = mo;
    mo->right->part->SetFlav(mo_flavs[1]);
    mo->right->part->SetInfo('F');
    mo->right->part->SetStatus(1);
    mo->right->tout   = sqr(mo_flavs[1].PSMass());
    //    std::cout<<" res da set on="<<mo->left->tout<<std::endl;
    mo->right->stat   = 3;  // *AS*

    if (mo->part->Info() != 'H') mo->part->SetInfo('f');
    mo->part->SetStatus(2);

    msg.Debugging()<<"Final_State_Shower::InitDaughters New "
		   <<mo->kn_no<<" --> "<<mo->left->kn_no<<", "<<mo->right->kn_no<<std::endl
		   <<"      "<<mo->left->part->Flav()<<", "<<mo->right->part->Flav()<<std::endl
		   <<"      mothers t = "<<mo->t<<", "<<mo->z<<std::endl;
  }

  // Reset kinematics
  mo->left->t       = mo->t; 

  mo->thcrit           = sqrt( mo->t/(mo->z*(1.- mo->z)*mo->E2) );
  double th            = M_PI;

  if ((mo->left->part->Flav().Strong()) && (mo->right->part->Flav().Strong())) th = mo->thcrit;
  else {
    if (mo->prev) mo->thcrit = th = mo->prev->thcrit;
    else mo->thcrit = th;
  }

  mo->left->t       = mo->t; 
  mo->left->E2      = (mo->z)*(mo->z)*(mo->E2); 
  mo->left->thcrit  = th;
  mo->left->maxpt2  = mo->maxpt2;

  mo->right->t      = mo->t; 
  mo->right->E2     = (1. - mo->z)*(1. - mo->z)*(mo->E2);
  mo->right->thcrit = th;
  mo->right->maxpt2 = mo->maxpt2; 
}




void Final_State_Shower::ResetDaughters(Knot * mo)
{
  msg.Debugging()<<"Final_State_Shower::ResetDaughters(Knot "<<mo->kn_no<<")"<<std::endl; 
  Reset(mo->left);
  Reset(mo->right);

  if (mo->part->Info() != 'H') mo->part->SetInfo('F');
  mo->part->SetStatus(1);
} 

void Final_State_Shower::Reset(Knot * mo) 
{ 
  if (!mo) return;
  msg.Debugging()<<"Final_State_Shower::Reset(Knot "<<mo->kn_no
		 <<", "<<mo->part->Info()<<")"<<std::endl; 

  mo->left  = 0;  
  mo->right = 0;

  mo->t     = mo->tout;
  mo->stat  = 0;
  mo->part->SetStatus(1);
  if (mo->part->Info() != 'H') mo->part->SetInfo('F');
  msg.Debugging()<<"   Now: (Knot "<<mo->kn_no<<", "<<mo->part->Info()<<")"<<std::endl; 
}

void Final_State_Shower::OutputTree(Tree * tree) 
{
  int number = 0;
  msg.Out()<<"final Tree:"<<std::endl<<tree<<std::endl;
  msg.Out()<<"Total 4 Mom = "<<GetMomentum(tree->GetRoot(),number);
  msg.Out()<<" for "<<number<<" FS particles."<<std::endl;
}


Vec4D  Final_State_Shower::GetMomentum(Knot * mo, int & number) {
  if (mo->left) {
    Vec4D p     = GetMomentum(mo->left,number) + GetMomentum(mo->right,number);
    Vec4D ptest = mo->left->part->Momentum() + mo->right->part->Momentum();
    if (!(ptest==mo->part->Momentum())) 
      msg.Tracking()<<" momentum conservation violation in Knot :"<<mo->kn_no<<std::endl
		    <<"          is: "<<ptest<<"  "<<p.Abs2()<<std::endl
		    <<"   should be:"<<mo->part->Momentum()
		    <<"  "<<mo->part->Momentum().Abs2()<<std::endl;
    return p;
  }
  number++;
  return mo->part->Momentum();
}
