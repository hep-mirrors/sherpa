#include "Initial_State_Shower.H"
#include "PDF_Handler.H"
#include "Data_Read.H"

#include <iomanip>

using namespace APACIC;
using namespace AMATOOLS;
using namespace APHYTOOLS;
using namespace AORGTOOLS;
using namespace std;

//-----------------------------------------------------------------------
//--------------------------- Constructors ------------------------------
//----------------------------------------------------------------------- 

Initial_State_Shower::Initial_State_Shower(ISR::ISR_Handler * _isr, 
					   Final_State_Shower * _fin,
					   Data_Read * _dataread) : 
  p_fin(_fin), m_fsron(0) 
{
  if (p_fin) m_fsron      = 1;
  if (_isr->On()) {
    double pt2fin     = 0.;
    if (p_fin) pt2fin = p_fin->PT2Min();  
    m_t0              = _dataread->GetValue<double>("IS PT2MIN",4.);
    p_tools           = new Sudakov_Tools(1,m_t0,(rpa.gen.Ecms())*(rpa.gen.Ecms()));
    p_kin             = new Spacelike_Kinematics(pt2fin);
    p_suds            = new Spacelike_Sudakov*[2];
    p_suds[0]         = new Spacelike_Sudakov(_isr->PDF(0),p_tools,p_kin,m_t0,_dataread);
    p_suds[1]         = new Spacelike_Sudakov(_isr->PDF(1),p_tools,p_kin,m_t0,_dataread);
    m_allowed         = 20;
    return;
  }

  msg.Error()<<"Error in Initial_State_Shower::Initial_State_Shower."<<endl
	     <<"   No ISR_Handler/ISR is switched off. Abort."<<endl;
  abort();
}

Initial_State_Shower::~Initial_State_Shower()
{
  if (p_tools) delete p_tools;
  if (p_kin)   delete p_kin;
  if (p_suds) {
    for (int i=0;i<2;i++) delete p_suds[i];
    delete p_suds;
  }
}

//-----------------------------------------------------------------------
//----------------------- Performing the Shower -------------------------
//----------------------------------------------------------------------- 

bool Initial_State_Shower::PerformShower(Tree ** trees,bool _jetveto) {
  m_jetveto = _jetveto;

  msg.Events()<<"----------------------------------------------------------"<<std::endl
		 <<"Initial_State_Shower::PerformShower : for Trees "<<trees<<std::endl;

  if (InitializeSystem(trees,trees[0]->GetRoot(),trees[1]->GetRoot())) {
    double x1,x2;
    Vec4D  cms;
    double E2 = sqr(rpa.gen.Ecms());   // *AS*
    //    double E2 = 4.*sqr(rpa.gen.Ecms());  // *FK*

    p_kin->BoostInCMS(trees,GetInitiator(trees[0]),GetInitiator(trees[1]));
    cms = trees[0]->GetInitiator()->part->Momentum() +
          trees[1]->GetInitiator()->part->Momentum();
    x1  = trees[0]->GetInitiator()->x;
    x2  = trees[1]->GetInitiator()->x;
    if ((dabs(m_sprime - cms.Abs2())/m_sprime > 1.e-6 ) ||
	(dabs(m_sprime - x1*x2*E2)/m_sprime > 1.e-6 )) {
      msg.Out()<<"ERROR in Initial_State_Shower : "<<std::endl;
      msg.Error()<<setprecision(12)
		 <<"ERROR in Initial_State_Shower : "<<std::endl
		 <<"Initial_State_Shower::PerformShower : "
		 <<"Mismatch of sprimes in CMS !"<<std::endl
		 <<"   "<<m_sprime<<" / "<<x1*x2*E2<<std::endl
		 <<"   "<<cms<<" / "<<cms.Abs2()<<std::endl;
      OutputTree(trees[0]);
      OutputTree(trees[1]);

      return 1;  // should return 0; but not supported by wrapper !!
    }
    else if ((dabs(m_sprime - cms.Abs2())/m_sprime > rpa.gen.Accu()) ||
	(dabs(m_sprime - x1*x2*E2)/m_sprime > rpa.gen.Accu())) {
      /* *AS*
      msg.Error()<<setprecision(12)
		 <<"WARNING in Initial_State_Shower : "<<std::endl
		 <<"Initial_State_Shower::PerformShower : "
		 <<"Mismatch of sprimes in CMS !"<<std::endl
		 <<"   "<<m_sprime<<" / "<<x1*x2*E2<<std::endl
		 <<"   "<<cms<<" / "<<cms.Abs2()<<std::endl;
      //      return 0;
      */
    }
    msg.Tracking()<<" In CMS "<<std::endl
		  <<" s after shower : "<<cms.Abs2()<<" == "<<m_sprime<<std::endl;
    m_lab = p_kin->BoostInLab(trees);
    cms   = trees[0]->GetInitiator()->part->Momentum() +
      trees[1]->GetInitiator()->part->Momentum();
    x1    = trees[0]->GetInitiator()->x;
    x2    = trees[1]->GetInitiator()->x;
    if ((dabs(m_sprime - cms.Abs2())/m_sprime > 1.e-6 ) ||
	(dabs(m_sprime - x1*x2*E2)/m_sprime > 1.e-6 )) {
      msg.Error()<<setprecision(12)
		 <<"ERROR in Initial_State_Shower : "<<std::endl
		 <<"Initial_State_Shower::PerformShower : "
		 <<"Mismatch of sprimes in Lab !"<<std::endl
		 <<"   "<<m_sprime<<" / "<<x1*x2*E2<<std::endl
		 <<"   "<<cms<<" / "<<cms.Abs2()<<std::endl;
      return 1;  // should return 0; but not supported by wrapper !!
    }
    else if ((dabs(m_sprime - cms.Abs2())/m_sprime > rpa.gen.Accu()) ||
	(dabs(m_sprime - x1*x2*E2)/m_sprime > rpa.gen.Accu())) {
      /* *AS*
      msg.Error()<<setprecision(12)
		 <<"WARNING in Initial_State_Shower : "<<std::endl
		 <<"Initial_State_Shower::PerformShower : "
		 <<"Mismatch of sprimes in LAB !"<<std::endl
		 <<"   "<<m_sprime<<" / "<<x1*x2*E2<<std::endl
		 <<"   "<<cms<<" / "<<cms.Abs2()<<std::endl;
      //      return 0;
      */
    }

    /*
    if (rpa.gen.Events()) {
      OutputTree(trees[0]);
      OutputTree(trees[1]);
    }
    msg.Events()<<" s after shower : "<<cms.Abs2()<<" == "<<x1*x2*E2<<", "<<"Internally : "<<m_sprime<<std::endl;
    */
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
  m_pt2_1   = m_pt2_2 = pt2max;
  m_th_1    = m_th_2  = M_PI;
}

//-----------------------------------------------------------------------
//---------------------------- After the Shower -------------------------
//----------------------------------------------------------------------- 

void Initial_State_Shower::ExtractPartons(Knot * kn,int beam,Blob * jet,
					  Blob_List * bl,Parton_List * pl) 
{
  if (!kn) return;
  int number;
  if (!kn->prev) {
    /* 
       New jet : kn = incoming parton from hadron info = 'I'
    */
    if (kn->part->Info() != 'G') {
      if (pl) number = pl->size();
         else number = int(kn->part);
      kn->part->SetNumber(number);
      kn->part->SetStatus(2);
      if (pl) pl->push_back(kn->part);
    }
    jet = new Blob();
    jet->SetStatus(1);
    jet->AddToInPartons(new Parton(kn->part));
    jet->SetId(bl->size());
    jet->SetType(std::string("IS Shower (APACIC++2.0)"));
    jet->SetBeam(beam);
    kn->part->SetDecayBlob(jet);
    bl->insert(bl->begin(),jet);
    if (!(kn->left)) {
      if (kn->part->Info() != 'G') {
	msg.Error()<<"Error in Initial_State_Shower : "<<std::endl
		   <<"Initial_State_Shower::ExtractPartons : "<<std::endl
		   <<" Mislabelling of parton "<<kn->part->Number()<<std::endl;
	if (pl) number = pl->size();
	else    number = int(kn->part);
	kn->part->SetNumber(number);
      }
      jet->AddToOutPartons(new Parton(kn->part));
      jet->SetStatus(0);
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
    jet->SetStatus(1);
    jet->AddToInPartons(new Parton(kn->part));
    jet->SetId(bl->size());
    jet->SetType(std::string("IS Shower (APACIC++2.0)"));
    bl->insert(bl->begin(),jet);
    if (kn->left) {
      kn->part->SetDecayBlob(jet);
      kn->part->SetStatus(2);
    }
    else {
      kn->part->SetStatus(1);
      jet->AddToOutPartons(new Parton(kn->part));
      jet->SetStatus(0);
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
	if (pl) number = pl->size();
	else    number = int(kn->part);
	kn->part->SetNumber(number);
	if (pl) pl->push_back(kn->part);
      }
      kn->part->SetProductionBlob(jet);
      kn->part->SetStatus(1);
      jet->AddToOutPartons(new Parton(kn->part));
    }
  }
  ExtractPartons(kn->left,beam,jet,bl,pl); 
  ExtractPartons(kn->right,beam,jet,bl,pl); 
}

//-----------------------------------------------------------------------
//---------------------------- Helpers ----------------------------------
//----------------------------------------------------------------------- 

bool Initial_State_Shower::TestShower(Tree ** trees) 
{
  int number;
  double x1,x2;
  Vec4D  cms;
  double E2 = sqr(rpa.gen.Ecms());
  //  double E  = rpa.gen.Ecms();

  ran.ReadInStatus("RandomA.dat",2735);

  msg.Out()<<" Starting Test IS Shower :"<<endl;
  for (long int n=1;n<=rpa.gen.NumberOfEvents();n++) { 
    //    int no=ran.WriteOutStatus("RandomA.dat");
    if (n%2500==0) msg.Out()<<" "<<n<<" events"<<endl;

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
//------------------- Initialization of the Shower ----------------------
//----------------------------------------------------------------------- 

bool Initial_State_Shower::InitializeSystem(Tree ** trees,Knot * k1,Knot * k2){
  msg.Events()<<"In Initial_State_Shower::InitializeSystem "<<std::endl;
  if (k1) msg.Events()<<(*k1)<<endl; else msg.Events()<<"###"<<endl;
  if (k2) msg.Events()<<(*k2)<<endl; else msg.Events()<<"###"<<endl;


  if ( (!k1) || (!k2) ) {
    msg.Error()<<"ERROR: Initial_State_Shower::InitializeSystem : No trees found !"<<std::endl;
    exit(1);
    return 0;
  }

  Flavour k1_flavs[2];
  Flavour k2_flavs[2];
  bool decay1 = (k1->stat>=1), decay2 = (k2->stat>=1);
  int first = 0;
  if ((!decay1 && decay2)||(decay1 && !decay2)) first=1;
  if (!decay1 && !decay2) first=2;


  //  Knot * k1save = new Knot(k1);
  //  Knot * k2save = new Knot(k2);
  trees[0]->Store();
  trees[1]->Store();


  int mismatch  = 0;
  bool accepted =1; 

  int asc1 = 0;

  for (;;) {
    msg.Tracking()<<"=============================="<<endl;
    msg.Tracking()<<" main loop run "<<++asc1<<endl;

    m_sprime      = (k1->part->Momentum()+k2->part->Momentum()).Abs2();

    accepted = 1;  
    // Parton 1/Tree 1 is the one to decay.
    if (decay1) {
      //      for (;;) {  // *AS* angle and pt vetos only done in dice!!!
	msg.Debugging()<<"------------------------------"<<endl;
	msg.Debugging()<<"  * calling FillBranch  I ("<<k1->kn_no<<")"<<endl;
	if (FillBranch(trees,k1,k2,0)) {
	  if (k1->z>0.) m_sprime = m_sprime/k1->z;

	  //	  if ((k1->maxpt2 < m_pt2_1) && (k1->thcrit < m_th_1)) break;   // *AS* !!!!?
	  //	  if (k1->stat==0) break;
	}
	else {
	  accepted = 0;
	  //	  break;
	}
	//      }
    }
    // Parton 2/Tree 2 is the one to decay.    
    if (decay2) {
      //      for (;;) {
	msg.Debugging()<<"------------------------------"<<endl;
	msg.Debugging()<<"  * calling FillBranch II ("<<k2->kn_no<<")"<<endl;
	if (FillBranch(trees,k2,k1,1)) {
	  if (k2->z > 0.) m_sprime = m_sprime/k2->z;

	  //	  if ((k2->maxpt2 < m_pt2_2) && (k2->thcrit < m_th_2)) break;   // *AS* !!!!?
	  //	  if (k2->stat==0) break;
	}
	else {
	  accepted = 0;
	  //	  break;
	}
	//    }
    }
    
    if (accepted) {

      p_kin->InitKinematics(trees,k1,k2,first);
      if (!decay1) SetColours(k1);
      if (!decay2) SetColours(k2);

      msg.Debugging()<<"------------------------------"<<endl;
      msg.Debugging()<<"  * calling EvolveSystem ("<<k1->kn_no<<")  ("<<k2->kn_no<<")"<<endl;
	
      if (EvolveSystem(trees,k1,k2)) {
	//	  if (k1save) delete k1save;
	//	  if (k2save) delete k2save;
	msg.Events()<<"Initial_State_Shower::InitializeSystem complete after "
		    <<mismatch<<" trials"<<std::endl;
	return 1;
      }
    }
    else {

	  
      
    }
    ++mismatch;
    if (mismatch > m_allowed) {
      msg.Error()<<"Error in Initial_State_Shower : "<<std::endl
		 <<"---------------------------------------"<<std::endl
		 <<"Initial_State_Shower::InitializeSystem failed. "
		 <<mismatch<<" trials."<<std::endl
	//		 <<" Knots :"<<std::endl
	//		 <<*(k1save)<<std::endl<<*(k2save)<<std::endl
		 <<"---------------------------------------"<<std::endl;
      //      if (k1save) delete k1save;
      //      if (k2save) delete k2save;
      return 0;
    }
    msg.Debugging()<<"---------------------------------------"<<std::endl
		   <<"Initial_State_Shower::InitializeSystem failed "
		   <<mismatch<<" trials."<<std::endl
		   <<"---------------------------------------"<<std::endl;
    /*
    //    trees[0]->Restore(k1);
    trees[0]->Restore(trees[0]->GetInitiator(),k1);
    //    trees[1]->Restore(k2);
    trees[1]->Restore(trees[1]->GetInitiator(),k2);
    k1->Copy(k1save);
    k2->Copy(k2save);
    */
    //    msg.Out()<<" Restoring Trees "<<endl;
    trees[0]->Restore();
    trees[1]->Restore();
    if (rpa.gen.Events()) {
      OutputTree(trees[0]);
      OutputTree(trees[1]);
    }

  }
}

//-----------------------------------------------------------------------
//------------------------ Evolution of the Shower ----------------------
//----------------------------------------------------------------------- 

bool Initial_State_Shower::EvolveSystem(Tree ** trees,Knot * k1,Knot * k2)
{
  msg.Debugging()<<" Initial_State_Shower::EvolveSystem( ("<<k1->kn_no<<"), ("<<k2->kn_no<<") );"<<endl;

  msg.Debugging()<<"===================================================="<<std::endl
		 <<"Initial_State_Shower::EvolveSystem for Knots:"<<std::endl
		 <<k1->kn_no<<"("<<k1->stat<<") with "<<k1->t<<" and "
		 <<k2->kn_no<<"("<<k2->stat<<") with "<<k2->t<<std::endl;


  bool decay1 = (k1->stat>=1), decay2 = (k2->stat>=1);
  int first = 0;
  if ((!decay1 && decay2)||(decay1 && !decay2)) first=1;
  if (!decay1 && !decay2) first=2;

  //  if (!(k1->stat) && !(k2->stat)) return 1;   // ----- *AS* 4jet ?!!!
  if ((k1->t == k1->tout) && (k2->t == k2->tout)) return 1;   // evolution finished


  int ntree0=0, ntree1=1;
  if (((k1->t) > (k2->t)) && (k2->t != k2->tout)) {  // make shure winner is still on (stat>0) !!! 
    Knot * kh=k1;
    k1 =k2;
    k2 =kh;
    ntree0=1;
    ntree1=0;
    msg.Debugging()<<"    - calling FillBranch  b "<<endl;
  }
  else {
    msg.Debugging()<<"    - calling FillBranch  a "<<endl;
  }


  if (k1->stat>0) {
    k1->stat           = 0;
    k1->E2             = sqr(k1->part->Momentum()[0]);
    k1->prev->E2       = k1->E2/sqr(k1->z);                        // ? exact ?
    k1->prev->left->E2 = k1->prev->E2*sqr(1.-k1->z);
  }  
  else {
    //    determine z and update x
    double sprime_a = (k1->part->Momentum()+k2->part->Momentum()).Abs2();
    double sprime_b = (k1->prev->part->Momentum()+k2->part->Momentum()).Abs2();
    //    cout<<" z,x old = "<<k1->z<<","<<k1->prev->x<<endl;
    k1->z=sprime_a/sprime_b;
    k1->prev->x=k1->x/k1->z;
//     double s = sqr(rpa.gen.Ecms());
//     double sprime_c = k1->x  * k2 ->x * s;
//     double sprime_d = k1->prev->x  * k2 ->x * s;
//     cout<<" sprime I  mom/x1x2 : "<<sprime_a<<"/"<<sprime_c<<endl;
//     cout<<" sprime II mom/x1x2 : "<<sprime_b<<"/"<<sprime_d<<endl;
//     cout<<" z,x fac = "<<k1->z<<","<<k1->prev->x<<endl;
    m_sprime/=k1->z;


    double pt2max = sqr(rpa.gen.Ecms());
    double th     = 4.*k1->z*k1->z*k1->t/(4.*k1->z*k1->z*k1->t-(1.-k1->z)*k1->x*k1->x*pt2max);

    k1->prev->thcrit       = k1->thcrit;
    k1->prev->t            = k1->t;
    k1->prev->left->thcrit = th;  // will be updated in FirstTimelike...
    k1->prev->left->t      = k1->prev->part->Momentum().Abs2();
  }

  if (k1->prev->stat>0) {
    if (!FillBranch(trees,k1->prev,k2,ntree0)) return 0; // try from beginning
  }
  double maxt = p_kin->CalculateMaxT(k1,k2);
  if (maxt<k1->prev->left->tout) {
    msg.Debugging()<<"Initial_State_Shower::FillBranch : for "<<k1->kn_no<<std::endl
		   <<"    No timelike branch possible here : "
		   <<maxt<<" < "<<k1->prev->left->tout<<std::endl;
    return 0; // try from beginning
  }
  else {
    if (k1->prev->left->stat>0)
      k1->prev->left->t = maxt;
    else {
      msg.Debugging()<<" timelike daughter already known "<<endl;
    }
  } 
  
  double sprime_a = (k1->part->Momentum()+k2->part->Momentum()).Abs2();
  p_fin->FirstTimelikeFromSpacelike(trees[0],k1->prev->left,m_jetveto,sprime_a,k1->z);

  if (!p_kin->DoKinematics(trees,k1,k2,ntree0,first)) {
    msg.Debugging()<<"Initial_State_Shower::EvolveSystem for knots"
		   <<k1->kn_no<<", "<<k2->kn_no<<std::endl
		   <<"   Couldn't do the kinematics."<<std::endl;
    return 0; // try from beginning
  }
  // *AS* extra output
//   OutputTree(trees[0]);
//   OutputTree(trees[1]);

  /* *AS* calling Color routine should be called with (k1->prev) */
  p_fin->SetAllColours(k1->prev->left);
  
  if (k1->prev->z>0.) m_sprime = m_sprime/k1->prev->z;

  if (ntree0==0) {
    msg.Debugging()<<"    - calling EvolveSystem  a"<<endl;
    return EvolveSystem(trees,k1->prev,k2);
  }
  else {  
    msg.Debugging()<<"    - calling EvolveSystem  b"<<endl;
    return EvolveSystem(trees,k2,k1->prev);
  }
   
}  

bool Initial_State_Shower::FillBranch(Tree ** trees,Knot * active,Knot * partner,int leg) {
  if (leg==0) msg.Debugging()<<" Initial_State_Shower::FillBranch( I ";
  else msg.Debugging()<<" Initial_State_Shower::FillBranch( II ";
  msg.Debugging()<<",("<<active->kn_no<<"), <"<<partner->kn_no<<"> );"<<endl;

  Flavour flavs[2];
  if (p_suds[leg]->Dice(active,m_sprime,m_jetveto)) {

    // *AS* ? jetveto hier ? active not active->prev->left? moved to dice !!!
    if (p_kin->KinCheck(active,m_jetveto)) return 0;

    flavs[0] = p_suds[leg]->GetFlA();
    flavs[1] = p_suds[leg]->GetFlC();    

    FillMotherAndSister(trees[leg],active,flavs);

    return 1;
  }
  active->prev   = 0;
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
  msg.Debugging()<<"Initial_State_Shower::FillMotherAndSister(";
  msg.Debugging()<<" ("<<k->kn_no<<") )   t="<<k->t<<endl;

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
  mother->thcrit = k->thcrit;

  //  mother->E2     = sqr(1./k->z) * k->E2;

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
  sister->thcrit = k->thcrit; 
  //  sister->E2     = sqr(1./k->z -1.) * k->E2;
  
//   msg.Debugging()<<"      mother : ("<<mother->kn_no<<") E2="<<mother->E2<<endl;
//   msg.Debugging()<<"      sister : ("<<sister->kn_no<<") E2="<<sister->E2<<endl;

  if (k->part->Info() != 'G') k->part->SetInfo('i');
  k->part->SetStatus(2);
  SetColours(k);
}

void Initial_State_Shower::SetColours(Knot * k)
{

  if (!k) return;
  Knot * mother = k->prev;
  if (!mother) return;
  Knot * sister = mother->left;
  if (!sister) return;

    //  check if already enough colours
  Knot * test=0;
  
  int all_colors_known=1;
  for (int i=0;i<3;++i) {
    if (i==0) test = k;
    if (i==1) test = mother;
    if (i==2) test = sister;
    
    if (test->part->Flav().Strong()) {
      int nc=0;
      if (test->part->GetFlow(1)) ++nc;
      if (test->part->GetFlow(2)) ++nc;
      if (test->part->Flav().IsQuark()) {
	if (nc!=1) {
	  all_colors_known=0;
	  break;
	  // NO
	}
      }
      else if (test->part->Flav().IsGluon()) {
	if (nc!=2) {
	  all_colors_known=0;
	  break;
	  // NO
	}
      }
      else {
	cout<<" ERROR: strong particle "<<test->part->Flav()<<" not covered by SetColours "<<endl;
      }
    }      
  }

  if (all_colors_known) {
    return;
  }

  // set colors
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
//   x1=0.1;
//   x2=0.3;
  x1=0.1;
  x2=0.04;

  double scale  = x1*x2*E2;
  double E      = 0.5 * sqrt(E2);
  m_pt2_1 = m_pt2_2 = scale;
  m_th_1  = m_th_2  = M_PI;

  Knot * d1   = trees[0]->NewKnot();
  *(d1->part) = Parton(1,Flavour(kf::u),x1*E*Vec4D(1.,0.,0.,1.));
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
  *(d2->part) = Parton(2,Flavour(kf::u).Bar(),x2*E*Vec4D(1.,0.,0.,-1.));
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
  if (tree->GetInitiator()==0) {
    msg.Out()<<"empty Tree"<<endl;
  }
  else {
    msg.Out()<<"final Tree:"<<std::endl<<tree<<std::endl
	     <<"Total 4 Mom = "<<GetMomentum(tree->GetInitiator(),number);
    msg.Out()<<" for "<<number<<" FS particles."<<std::endl;
    // Note: we NEED two "msg.Out()" since otherwise "number" is print
    //       before it is calculated!!
  }
}







