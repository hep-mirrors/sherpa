#include "Initial_State_Shower.H"
#include "PDF_Handler.H"
#include "Data_Read.H"

using namespace APACIC;
using namespace AMATOOLS;
using namespace APHYTOOLS;
using namespace AORGTOOLS;

//-----------------------------------------------------------------------
//--------------------------- Constructors ------------------------------
//----------------------------------------------------------------------- 

Initial_State_Shower::Initial_State_Shower(ISR::ISR_Handler * _isr, Final_State_Shower * _fin) : 
  fin(_fin) {// i think isr doesn't have to be modified no more when i'm through with this
  
  if (_isr) {
    msg.Debugging()<<"-----------------------------------------------------------"<<std::endl;
    msg.Debugging()<<"Initial_State_Shower::Initial_State_Shower : Beams : "
		   <<rpa.gen.Beam1()<<", "<<rpa.gen.Beam2()<<std::endl
		   <<"ISR_Handler : "<<_isr->Type()<<std::endl;
    msg.Debugging()<<"-----------------------------------------------------------"<<std::endl;
  }
  else {
    // test shower only
    msg.Error()<<"Initialise Initial_State_Shower + ISR_Handler for testing purposes only !!"<<endl;

    msg.Debugging()<<"Open file "<<string("./ISR.dat")<<endl;
    Data_Read dr(string("./ISR.dat"));

    int    * isrtypes = new int[2];
    isrtypes[0]       = dr.GetValue<ISR_Type::code>("ISR1");
    isrtypes[1]       = dr.GetValue<ISR_Type::code>("ISR2");
    double s          = sqr(AORGTOOLS::rpa.gen.Ecms());
    double * splimits = new double[2];
    splimits[0]       = s*sqr(dr.GetValue<double>("SMIN"));
    splimits[1]       = s*sqr(dr.GetValue<double>("SMAX"));

    Flavour * beams, * partons;
    beams        = new Flavour[2];
    partons      = new Flavour[2];
    beams[0]     = rpa.gen.Beam1();
    beams[1]     = rpa.gen.Beam2();
//     beams[0]     = Flavour(kf::p_plus);
//     beams[1]     = Flavour(kf::p_plus);
    rpa.gen.SetBeam1(beams[0]);
    rpa.gen.SetBeam2(beams[1]);

    partons[0]   = Flavour(kf::gluon);
    partons[1]   = Flavour(kf::gluon);

    _isr         = new ISR::ISR_Handler(isrtypes,beams,partons,splimits);

    delete [] isrtypes;
    delete [] splimits;
    delete [] beams;
    delete [] partons;
  }

  if (_isr->On()) {
    t0        = rpa.pshower.InitialQ02();
    tools     = new Sudakov_Tools(1,t0,(rpa.gen.Ecms())*(rpa.gen.Ecms()));
    kin       = new Spacelike_Kinematics();
    suds      = new Spacelike_Sudakov*[2];
    suds[0]   = new Spacelike_Sudakov(_isr->PDF(0),tools);
    suds[1]   = new Spacelike_Sudakov(_isr->PDF(1),tools);
  
    allowed   = 20;
  }
  else {
    msg.Error()<<"Initialise Initial_State_Shower according to file ISR.dat"<<endl;

    tools = 0;
    suds  = 0;
    kin   = 0;
  }
}

Initial_State_Shower::~Initial_State_Shower()
{
  if (tools) delete tools;
  if (kin)   delete kin;
  if (suds) {
    for (int i=0;i<2;i++) delete suds[i];
    delete suds;
  }
}

//-----------------------------------------------------------------------
//----------------------- Performing the Shower -------------------------
//----------------------------------------------------------------------- 

bool Initial_State_Shower::PerformShower(Tree ** trees,bool _jetveto) {
  jetveto = _jetveto;

  msg.Debugging()<<"----------------------------------------------------------"<<std::endl
		 <<"Initial_State_Shower::PerformShower : for Trees "<<trees<<std::endl;

  if (InitializeSystem(trees,trees[0]->GetRoot(),trees[1]->GetRoot())) {
    double x1,x2;
    Vec4D  cms;
    double E2 = 4.*sqr(rpa.gen.Ecms());
    double E  = rpa.gen.Ecms();

    kin->BoostInCMS(trees,GetInitiator(trees[0]),GetInitiator(trees[1]));
    cms = trees[0]->GetInitiator()->part->Momentum() +
          trees[1]->GetInitiator()->part->Momentum();
    x1  = trees[0]->GetInitiator()->x;
    x2  = trees[1]->GetInitiator()->x;
    if ((dabs(sprime - cms.Abs2())/sprime > rpa.gen.Accu()) ||
	(dabs(sprime - x1*x2*E2)/sprime > rpa.gen.Accu())) {
      msg.Error()<<"Error in Initial_State_Shower : "<<std::endl
		 <<"Initial_State_Shower::PerformShower : "
		 <<"Mismatch of sprimes in CMS !"<<std::endl
		 <<"   "<<sprime<<" / "<<x1*x2*E2<<std::endl
		 <<"   "<<cms<<" / "<<cms.Abs2()<<std::endl;
      return 0;
    }
    msg.Tracking()<<" In CMS "<<std::endl
		  <<" s after shower : "<<cms.Abs2()<<" =?= "<<std::endl
		  <<"Internally : "<<sprime<<std::endl;
    lab = kin->BoostInLab(trees);
    cms = trees[0]->GetInitiator()->part->Momentum() +
      trees[1]->GetInitiator()->part->Momentum();
    x1  = trees[0]->GetInitiator()->x;
    x2  = trees[1]->GetInitiator()->x;
    if ((dabs(sprime - cms.Abs2())/sprime > rpa.gen.Accu()) ||
	(dabs(sprime - x1*x2*E2)/sprime > rpa.gen.Accu())) {
      msg.Error()<<"Error in Initial_State_Shower : "<<std::endl
		 <<"Initial_State_Shower::PerformShower : "
		 <<"Mismatch of sprimes in Lab !"<<std::endl
		 <<"   "<<sprime<<" / "<<x1*x2*E2<<std::endl
		 <<"   "<<cms<<" / "<<cms.Abs2()<<std::endl;
      return 0;
    }
    if (rpa.gen.Events()) {
      OutputTree(trees[0]);
      OutputTree(trees[1]);
    }
    msg.Events()<<" s after shower : "<<cms.Abs2()<<" =?= "
		<<x1*x2*E2<<", "<<"Internally : "<<sprime<<std::endl;
    return 1;
  }
    
  msg.Error()<<"Error in Initial_State_Shower : "<<std::endl
	     <<"Initial_State_Shower::PerformShower : "
	     <<"  Could not initialize any system and perform the shower."<<std::endl;

  return 0;
}

//-----------------------------------------------------------------------
//---------------------------- Before the Shower ------------------------
//----------------------------------------------------------------------- 

void   Initial_State_Shower::InitShowerPT(double pt2max) {
  pt2_1   = pt2_2 = pt2max;
  th_1    = th_2  = M_PI;
}

//-----------------------------------------------------------------------
//---------------------------- After the Shower -------------------------
//----------------------------------------------------------------------- 

void Initial_State_Shower::ExtractPartons(Knot * kn,Blob * jet,
					  Blob_List * bl,Parton_List * pl,int beam) 
{
  if (!kn) return;
  if (!kn->prev) {
    /* 
       New jet : kn = incoming parton from hadron info = 'I'
    */
    if (kn->part->Info() != 'G') {
      kn->part->SetNumber(pl->size());
      kn->part->SetStatus(2);
      pl->push_back(kn->part);
    }
    jet = new Blob();
    jet->AddToInPartons(kn->part);
    jet->SetId(bl->size());
    jet->SetType(std::string("IS Parton Shower (APACIC++2.0)"));
    jet->SetBeam(beam);
    kn->part->SetDec(jet);
    bl->push_back(jet);
    if (!(kn->left)) {
      if (kn->part->Info() != 'G') {
	msg.Error()<<"Error in Initial_State_Shower : "<<std::endl
		   <<"Initial_State_Shower::ExtractPartons : "<<std::endl
		   <<" Mislabelling of parton "<<kn->part->Number()<<std::endl;
	kn->part->SetNumber(pl->size());
      }
      jet->AddToOutPartons(kn->part);
      return;
    }
  }
  else if (kn->part->Info() == 'H') {
    /* 
       New jet : kn = hard parton from ME info = 'H'
                 and kn outgoing
		 or kn->left or kn->right not from ME
    */
    jet = new Blob();
    jet->AddToInPartons(kn->part);
    jet->SetId(bl->size());
    jet->SetType(std::string("IS Parton Shower (APACIC++2.0)"));
    bl->push_back(jet);
    if (kn->left) {
      kn->part->SetDec(jet);
      kn->part->SetStatus(2);
    }
    else {
      kn->part->SetStatus(1);
      jet->AddToOutPartons(kn->part);
      return;
    }
  }
  else {
    if (!kn->left) {
      if (!jet) {
	msg.Error()<<"Error in Initial_State_Shower : "<<std::endl
		   <<"Initial_State_Shower::ExtractPartons : "<<std::endl
		   <<"No jet for Parton "<<kn->part->Number()<<std::endl;
      }
      if (kn->part->Info() != 'G') {
	kn->part->SetNumber(pl->size());
	pl->push_back(kn->part);
      }
      kn->part->SetProd(jet);
      kn->part->SetStatus(1);
      jet->AddToOutPartons(kn->part);
    }
  }
  ExtractPartons(kn->left,jet,bl,pl,beam); 
  ExtractPartons(kn->right,jet,bl,pl,beam); 
}

//-----------------------------------------------------------------------
//---------------------------- Helpers ----------------------------------
//----------------------------------------------------------------------- 

bool Initial_State_Shower::TestShower(Tree ** trees) 
{
  int number;
  double x1,x2;
  Vec4D  cms;
  double E2 = 4.*sqr(rpa.gen.Ecms());
  double E  = rpa.gen.Ecms();
  for (long int n=1;n<=rpa.gen.NumberOfEvents();n++) {
    msg.Events()<<"++++++++++++++++++++++++++++"<<n
		<<" th event +++++++++++++++++++++++++"<<std::endl;
    for (int i=0;i<2;i++) trees[i]->Reset();
    InitTwoTrees(trees,E2);
    if (!PerformShower(trees,0)) return 0;
  }
  msg.Events()<<"Initial_State_Shower::TestShower : "
	      <<"Terminated loops over events successfully."<<std::endl;
  return 1;
}


//-----------------------------------------------------------------------
//------------------- Initialisation of the Shower ----------------------
//----------------------------------------------------------------------- 

bool Initial_State_Shower::InitializeSystem(Tree ** trees,Knot * k1,Knot * k2){
  if ( (!k1) || (!k2) ) {
    msg.Error()<<"Error in Initial_State_Shower : "<<std::endl
	       <<"Initial_State_Shower::InitializeSystem : No trees found !"<<std::endl;
    return 0;
  }

  Flavour k1_flavs[2];
  Flavour k2_flavs[2];
  bool decay1 = (k1->stat==1),decay2 = (k2->stat==1);

  Knot * k1save = new Knot(k1);
  Knot * k2save = new Knot(k2);

  int mismatch = 0;
  bool accepted; 
  sprime      = (k1->part->Momentum()+k2->part->Momentum()).Abs2();

  msg.Debugging()<<"Initial_State_Shower::InitializeSystem with s' = "<<sprime<<std::endl
		 <<"   Compare with scale : "<<4.*k1->x*k2->x*sqr(rpa.gen.Ecms())<<std::endl
		 <<"   x1,2 = "<<k1->x<<"("<<k1->stat<<"), "<<k2->x<<"("<<k2->stat<<")"<<std::endl
		 <<"Vecs :"<<std::endl<<k1->part->Momentum()<<std::endl<<k2->part->Momentum()<<std::endl;

  for (;;) {
    sprime      = (k1->part->Momentum()+k2->part->Momentum()).Abs2();
    accepted = 1;  
    // Parton 1/Tree 1 is the one to decay.
    if (!decay1) {
      if (!(InitializeSystem(trees,k1->prev,k2))) accepted = 0;
    }
    else {
      for (;;) {
	if (FillBranch(trees,k1,k2,0)) {
	  if (k1->z > 0.) sprime = sprime/k1->z;
	  if ((k1->maxpt2 < pt2_1) && (k1->thcrit < th_1)) break;
	  if (k1->stat==0) break;
	}
	else {
	  accepted = 0;
	  break;
	}
      }
    }
    // Parton 2/Tree 2 is the one to decay.    
    if (!decay2) { if (!(InitializeSystem(trees,k1,k2->prev))) accepted = 0; }
    else {
      for (;;) {
	if (FillBranch(trees,k2,k1,1)) {
	  if (k2->z > 0.) sprime = sprime/k2->z;
	  if ((k2->maxpt2 < pt2_2) && (k2->thcrit < th_2)) break;
	  if (k2->stat==0) break;
	}
	else {
	  accepted = 0;
	  break;
	}
      }
    }
    
    if (accepted) {
      kin->InitKinematics(k1,k2);
      if (EvolveSystem(trees,k1,k2)) {
	if (k1save) delete k1save;
	if (k2save) delete k2save;
	msg.Events()<<"Initial_State_Shower::InitializeSystem complete after "
		    <<mismatch<<" trials"<<std::endl;
	return 1;
      }
    }
    mismatch++;
    if (mismatch > allowed) {
      msg.Error()<<"Error in Initial_State_Shower : "<<std::endl
		 <<"---------------------------------------"<<std::endl
		 <<"Initial_State_Shower::InitializeSystem failed. "
		 <<mismatch<<" trials."<<std::endl
		 <<" Knots :"<<std::endl
		 <<*(k1save)<<std::endl<<*(k2save)<<std::endl
		 <<"---------------------------------------"<<std::endl;
      if (k1save) delete k1save;
      if (k2save) delete k2save;
      return 0;
    }
    msg.Debugging()<<"---------------------------------------"<<std::endl
		   <<"Initial_State_Shower::InitializeSystem failed. "
		   <<mismatch<<" trials."<<std::endl
		   <<"---------------------------------------"<<std::endl;

    trees[0]->Restore(k1);
    trees[1]->Restore(k2);
    k1->Copy(k1save);
    k2->Copy(k2save);
  }
}

//-----------------------------------------------------------------------
//------------------------ Evolution of the Shower ----------------------
//----------------------------------------------------------------------- 

bool Initial_State_Shower::EvolveSystem(Tree ** trees,Knot * k1,Knot * k2)
{
  msg.Debugging()<<"===================================================="<<std::endl
		 <<"Initial_State_Shower::EvolveSystem for Knots:"<<std::endl
		 <<k1->kn_no<<"("<<k1->stat<<") with "<<k1->t<<" and "
		 <<k2->kn_no<<"("<<k2->stat<<") with "<<k2->t<<std::endl;
  if (!(k1->stat) && !(k2->stat)) return 1;

  if ((k1->t) < (k2->t)) {
    k1->stat           = 0;
    k1->E2             = sqr(k1->part->Momentum()[0]);
    k1->prev->E2       = k1->E2/sqr(k1->z);
    k1->prev->left->E2 = k1->prev->E2*sqr(1.-k1->z);
    if (!FillBranch(trees,k1->prev,k2,0)) return 0;
    msg.Debugging()<<"Calling FirstTimelikeFromSpacelike1 for "
		   <<k1->prev->left->kn_no<<std::endl<<"    "
		   <<k1->prev->left->E2<<" from mother "<<k1->prev->kn_no
		   <<" "<<k1->prev->E2<<" and "<<k1->z<<" : "<<std::endl
		   <<"    sister : "<<k1->E2<<" : "<<k1->part->Momentum()<<std::endl;
    fin->FirstTimelikeFromSpacelike(trees[0],k1->prev->left,jetveto);
    if (!kin->DoKinematics(trees,k1,k2,0)) {
      msg.Debugging()<<"Initial_State_Shower::EvolveSystem for knots"
		     <<k1->kn_no<<", "<<k2->kn_no<<std::endl
		     <<"   Couldn't do the kinematics."<<std::endl;
      return 0;
    }
    if (k1->prev->z > 0.) sprime = sprime/k1->prev->z;
    return EvolveSystem(trees,k1->prev,k2);
  }
  else {
    k2->stat = 0;
    k2->E2             = sqr(k2->part->Momentum()[0]);
    k2->prev->E2       = k2->E2/sqr(k2->z);
    k2->prev->left->E2 = k2->prev->E2*sqr(1.-k2->z);
    if (!FillBranch(trees,k2->prev,k1,1)) return 0;
    msg.Debugging()<<"Calling FirstTimelikeFromSpacelike2 for "
		   <<k2->prev->left->kn_no<<std::endl<<"    "
		   <<k2->prev->left->E2<<" from mother "<<k2->prev->kn_no
		   <<" "<<k2->prev->E2<<" and "<<k2->z<<" : "<<std::endl
		   <<"    sister : "<<k2->E2<<" : "<<k2->part->Momentum()<<std::endl;
    fin->FirstTimelikeFromSpacelike(trees[1],k2->prev->left,jetveto);
    if (!kin->DoKinematics(trees,k2,k1,1)) {
      msg.Debugging()<<"Initial_State_Shower::EvolveSystem for knots"
		     <<k1->kn_no<<", "<<k2->kn_no<<std::endl
		     <<"   Couldn't do the kinematics."<<std::endl;
      return 0;
    }
    if (k2->prev->z > 0.) sprime = sprime/k2->prev->z;
    return EvolveSystem(trees,k1,k2->prev);
  };
}  

bool Initial_State_Shower::FillBranch(Tree ** trees,Knot * active,Knot * partner,int leg) {
  Flavour flavs[2];
  if (suds[leg]->Dice(active,sprime)) {
    if (kin->KinCheck(active,jetveto)) return 0;
    flavs[0] = suds[leg]->GetFlA();
    flavs[1] = suds[leg]->GetFlC();    
    FillMotherAndSister(trees[leg],active,flavs);
    double maxt = kin->CalculateMaxT(active,partner);
    if (maxt<active->prev->left->tout) {
      msg.Debugging()<<"Initial_State_Shower::FillBranch : for "<<active->kn_no<<std::endl
		     <<"    No timelike branch possible here : "
		     <<maxt<<" < "<<active->prev->left->tout<<std::endl;
      trees[leg]->Restore(active);
      return 0;
    }
    else {
      active->prev->left->t = maxt;
      return 1;
    }
  }
  active->stat   = 0;
  active->t      = active->tout;
  active->thcrit = 0.;
  active->maxpt2 = 0.;
  active->part->SetStatus(1);
  return 1;
}

//-----------------------------------------------------------------------
//------------------- Service for the Branchings ------------------------
//----------------------------------------------------------------------- 

void Initial_State_Shower::FillMotherAndSister(Tree * tree,Knot * k,Flavour * k_flavs)
{
  Knot * mother  = tree->NewKnot();
  k->prev        = mother;
  mother->right  = k;
  mother->part->SetFlav(k_flavs[0]);
  mother->part->SetInfo('I');
  mother->part->SetStatus(1);
  mother->t      = k->t;
  mother->tout   = rpa.pshower.InitialQ02();
  mother->x      = k->x/k->z;
  mother->stat   = 1;
  mother->E2     = 0.;

  Knot * sister  = tree->NewKnot();
  sister->prev   = mother;
  mother->left   = sister;
  sister->part->SetFlav(k_flavs[1]);
  sister->part->SetInfo('F');
  sister->part->SetStatus(1);
  sister->t      = 0.;
  sister->tout   = sqr(k_flavs[1].PSMass());
  sister->x      = (mother->x)*(1.-k->z);
  sister->stat   = 1;
  sister->E2     = 0.;
  
  if (k->part->Info() != 'G') k->part->SetInfo('i');
  k->part->SetStatus(2);
  SetColours(k);
}

void Initial_State_Shower::SetColours(Knot * k)
{
  Knot * mother = k->prev;
  Knot * sister = mother->left;
  if (mother->part->Flav().IsQuark()) {
    if (mother->part->Flav().IsAnti()) {
      if (k->part->Flav().IsQuark()) {
	/*
	  mother             k = antiquark
	    -----   ---------
	        |   |
	        |   |________
	        -------------
		        sister
	 */
	mother->part->SetFlow(1,0);
	mother->part->SetFlow(2,-1);
	sister->part->SetFlow(1,k->part->GetFlow(2));
	sister->part->SetFlow(2,mother->part->GetFlow(2));
      }
      else {
	/*
	  mother        sister = antiquark
	    -----   ---------
	        |   |
	        |   |________
	        -------------
		             k
	 */
	mother->part->SetFlow(1,0);
	mother->part->SetFlow(2,k->part->GetFlow(2));
	sister->part->SetFlow(1,0);
	sister->part->SetFlow(2,k->part->GetFlow(1));
      }
    }
    else {
      if (k->part->Flav().IsQuark()) {
	/*
	  mother             k = quark
	    -----   ---------
	        |   |
	        |   |________
	        -------------
		        sister
	*/
	mother->part->SetFlow(1,-1);
	mother->part->SetFlow(2,0);
	sister->part->SetFlow(1,mother->part->GetFlow(1));
	sister->part->SetFlow(2,k->part->GetFlow(1));
      }
      else {
	/*
	  mother        sister = quark
	    -----   ---------
	        |   |
	        |   |________
	        -------------
		             k
	 */
	mother->part->SetFlow(1,k->part->GetFlow(1));
	mother->part->SetFlow(2,0);
	sister->part->SetFlow(1,k->part->GetFlow(2));
	sister->part->SetFlow(2,0);
      }
    }
  }
  else {
    if (k->part->Flav().IsQuark()) {
      if (k->part->Flav().IsAnti()) {
	/*
	  mother        sister = quark
	         |------------
	    =====   
	         |------------
		             k = antiquark
	 */
	mother->part->SetFlow(1,-1);
	mother->part->SetFlow(2,k->part->GetFlow(2));
	sister->part->SetFlow(1,mother->part->GetFlow(1));
	sister->part->SetFlow(2,0);
      }
      else {
	/*
	  mother        sister = antiquark
	         |------------
	    =====   
	         |------------
		             k = quark
	 */
	mother->part->SetFlow(1,k->part->GetFlow(1));
	mother->part->SetFlow(2,-1);
	sister->part->SetFlow(1,0);
	sister->part->SetFlow(2,mother->part->GetFlow(2));
      }
    }
    else {
      /*
	mother        sister = gluon
             |------------
	     | |---------
	=====  | 
	     | |__________
	     |------------

	                 k = gluon
      */
      if (ran.Get() > 0) {
	mother->part->SetFlow(1,k->part->GetFlow(1));
	mother->part->SetFlow(2,-1);
	sister->part->SetFlow(1,k->part->GetFlow(2));
	sister->part->SetFlow(2,mother->part->GetFlow(2));
      }
      else {
	mother->part->SetFlow(1,-1);
	mother->part->SetFlow(2,k->part->GetFlow(2));
	sister->part->SetFlow(1,mother->part->GetFlow(1));
	sister->part->SetFlow(2,k->part->GetFlow(1));
      }
    }
  }
  return;
}


void Initial_State_Shower::InitTwoTrees(Tree ** trees,double E2) {
  double x1     = 0.005+ran.Get()*0.295;
  double x2     = 0.005+ran.Get()*0.295;
  double scale  = x1*x2*E2;
  double E      = 0.5 * sqrt(E2);
  pt2_1 = pt2_2 = scale;
  th_1  = th_2  = M_PI;

  Knot * d1   = trees[0]->NewKnot();
  *(d1->part) = Parton(1,Flavour(kf::gluon),x1*E*Vec4D(1.,0.,0.,1.));
  d1->part->SetStatus(1);
  d1->part->SetInfo('G');
  d1->part->SetFlow(1,500);
  d1->part->SetFlow(2,501);
  d1->t       = -scale;
  d1->tout    = rpa.pshower.InitialQ02();
  d1->x       = x1;
  d1->E2      = sqr(x1*E);
  d1->maxpt2  = scale;
  d1->costh   = -1.; 
  d1->thcrit  = M_PI;
  d1->stat    = 1;
  
  Knot * d2   = trees[1]->NewKnot();
  *(d2->part) = Parton(2,Flavour(kf::gluon),x2*E*Vec4D(1.,0.,0.,-1.));
  d2->part->SetStatus(1);
  d2->part->SetInfo('G');
  d2->part->SetFlow(1,502);
  d2->part->SetFlow(2,d1->part->GetFlow(1));
  d2->t       = -scale;
  d2->tout    = rpa.pshower.InitialQ02();
  d2->x       = x2;
  d2->E2      = sqr(x2*E);
  d2->maxpt2  = scale;
  d2->costh   = -1.; 
  d2->thcrit  = M_PI;
  d2->stat    = 1;

  msg.Debugging()<<"Initial_State_Shower::InitTwoTrees :"<<std::endl
		 <<"    Daughter1 :"<<d1->part->Momentum()<<std::endl
		 <<"    Daughter2 :"<<d2->part->Momentum()<<std::endl
		 <<"    Check :"<<scale<<", "
		 <<(d1->part->Momentum()+d2->part->Momentum()).Abs2()<<std::endl;
};

void Initial_State_Shower::OutputTree(Tree * tree) 
{
  int number=0;
  msg.Out()<<"final Tree:"<<std::endl<<tree<<std::endl
	   <<"Total 4 Mom = "<<GetMomentum(tree->GetInitiator(),number);
  msg.Out()<<" for "<<number<<" FS particles."<<std::endl;
}







