#include "Initial_State_Shower.H"
#include "PDF_Handler.H"
#include "Data_Read.H"

#include <iomanip>

using namespace APACIC;
using namespace ATOOLS;
using namespace std;

//-----------------------------------------------------------------------
//--------------------------- Constructors ------------------------------
//----------------------------------------------------------------------- 

Initial_State_Shower::Initial_State_Shower(PDF::ISR_Handler * _isr, 
					   Final_State_Shower * _fin,
					   MODEL::Model_Base * _model,
					   Data_Read * const _dataread) : 
  p_fin(_fin), m_fsron(0) 
{
  if (p_fin) m_fsron      = 1;
  if (_isr->On()) {
    double pt2fin     = 0.;
    if (p_fin) pt2fin = p_fin->PT2Min();  
    m_t0              = _dataread->GetValue<double>("IS PT2MIN",4.);
    m_t0 = dabs(m_t0);
    m_jetveto_scheme  = _dataread->GetValue<int>("IS JETVETOSCHEME",1);

    p_tools           = new Sudakov_Tools(1,_model,m_t0,(rpa.gen.Ecms())*(rpa.gen.Ecms()));
    p_kin             = new Spacelike_Kinematics(pt2fin, _dataread);
    p_suds            = new Spacelike_Sudakov*[2];
    p_suds[0]         = new Spacelike_Sudakov(_isr->PDF(0),p_tools,p_kin,m_t0,_dataread);
    p_suds[1]         = new Spacelike_Sudakov(_isr->PDF(1),p_tools,p_kin,m_t0,_dataread);

    m_allowed         = 200;

    m_extra_pdf[0]    = 1;
    m_extra_pdf[1]    = 1;
    m_to_be_diced[0]  = 1;
    m_to_be_diced[1]  = 1;
    m_t0 = - dabs(m_t0);
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

bool Initial_State_Shower::PerformShower(Tree ** trees,int _jetveto) {
  m_jetveto = (_jetveto>0);
  if (_jetveto<0) {
    m_extra_pdf[0]    = 0;
    m_extra_pdf[1]    = 0;
  }
  else {
    m_extra_pdf[0]    = 1;
    m_extra_pdf[1]    = 1;
  }

  if (InitializeSystem(trees,trees[0]->GetRoot(),trees[1]->GetRoot())) {
    double x1,x2;
    Vec4D  cms;
    p_kin->BoostInCMS(trees,GetInitiator(trees[0]),GetInitiator(trees[1]));
    cms = trees[0]->GetInitiator()->part->Momentum() +
          trees[1]->GetInitiator()->part->Momentum();
    x1  = trees[0]->GetInitiator()->x;
    x2  = trees[1]->GetInitiator()->x;
    m_lab = p_kin->BoostInLab(trees);
    cms   = trees[0]->GetInitiator()->part->Momentum() +
      trees[1]->GetInitiator()->part->Momentum();
    x1    = trees[0]->GetInitiator()->x;
    x2    = trees[1]->GetInitiator()->x;
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
					  Blob_List * bl,Particle_List * pl) 
{
  if (!kn) return;

  // fetch last PSME blobs
  m_bl_meps_is=0;
  m_bl_meps_fs=0;
  for (Blob_Iterator blit=bl->begin();blit!=bl->end();++blit) {
    if ((*blit)->Type()==btp::ME_PS_Interface_IS) {
      m_bl_meps_is=(*blit);
    }
    if ((*blit)->Type()==btp::ME_PS_Interface_FS) {
      m_bl_meps_fs=(*blit);
    }
  }
  if (m_bl_meps_is==NULL) {
    ATOOLS::msg.Error()<<"Initial_State_Shower::ExtractPartons(..): "
		       <<"No ME PS Interface found!"<<std::endl
		       <<"   Cannot proceed. Abort."<<std::endl;
    exit(126);
  }
  m_bl_meps_is->SetStatus(0);

  int nr=1000;
  SingleExtract(kn,beam,jet,bl,nr);
}


void Initial_State_Shower::SingleExtract(Knot * kn,int beam,Blob * jet,
					  Blob_List * bl,int & nr) 
{
  if (!kn) return;

  Particle * p=NULL;

  bool newblob  = false;
  bool lastknot = false;
  bool is_is    = true;
  bool ignore   = false;

  if (!kn->prev) {
    newblob=true;
  }
  else {
    for (Knot * k=kn;is_is && k->prev; k=k->prev) {
      if (k==k->prev->left) is_is=false;
    }

    if (is_is && (kn->prev->part->Info()=='G' || kn->prev->part->Info()=='H'))
      ignore = true;
  }
  if (!kn->left)   lastknot=true;

  if (!is_is && !lastknot) {
    if (kn->left->part->Info()=='H') ignore=true;
  }
  if (!ignore && !is_is && kn->part->Info()=='H') newblob=true;

  // --- create new blob ---
  if (newblob) {
    jet = new Blob();
    jet->SetStatus(1);
    jet->SetId(bl->size());
    if (is_is) {
      jet->SetType(btp::IS_Shower);
      jet->SetTypeSpec("APACIC++2.0");
    }
    else {
      jet->SetType(btp::FS_Shower);
      jet->SetTypeSpec("APACIC++-2.0");
    }
    bl->insert(bl->begin(),jet);

    if (!kn->prev) jet->SetBeam(beam);
    p = new Particle(kn->part);
    p->SetStatus(2);
    jet->AddToInParticles(p);
  }

  // --- add to MEPS blob ---
  if (!ignore && (kn->part->Info()=='H' || kn->part->Info()=='G')) {
    if  (is_is && m_bl_meps_is) {
      p = new Particle(kn->part);
      jet->AddToOutParticles(p);
      p->SetStatus(2);
      m_bl_meps_is->AddToInParticles(p);
    }
    else if (!is_is && m_bl_meps_fs) {
      if (!p) {
	p = new Particle(kn->part);
	jet->AddToInParticles(p);
      }
      p->SetStatus(2);
      m_bl_meps_fs->AddToOutParticles(p);
    }
  }

  // --- add final state particle ---
  if (lastknot && !is_is) {
    p = new Particle(kn->part);
    p->SetStatus(1);
    jet->AddToOutParticles(p);
  }

  SingleExtract(kn->left,beam,jet,bl,nr); 
  SingleExtract(kn->right,beam,jet,bl,nr); 
}
//-----------------------------------------------------------------------
//---------------------------- Helpers ----------------------------------
//----------------------------------------------------------------------- 

bool Initial_State_Shower::TestShower(Tree ** trees) 
{
  double E2 = sqr(rpa.gen.Ecms());

  ran.ReadInStatus("RandomA.dat",2735);

  msg.Out()<<" Starting Test IS Shower :"<<endl;
  for (long int n=1;n<=rpa.gen.NumberOfEvents();n++) { 
    if (n%2500==0) msg.Out()<<" "<<n<<" events"<<endl;

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
  m_to_be_diced[0]  = 1;
  m_to_be_diced[1]  = 1;

  if ( (!k1) || (!k2) ) {
    msg.Error()<<"ERROR: Initial_State_Shower::InitializeSystem : No trees found !"<<std::endl;
    exit(1);
    return 0;
  }

  if (msg.LevelIsDebugging()) {
    cout<<" Initialize System : "<<endl;
    cout<<*k1<<*k2<<endl;;
  }

  bool decay1 = (k1->stat>=1), decay2 = (k2->stat>=1);
  int first = 0;
  if ((!decay1 && decay2)||(decay1 && !decay2)) first=1;
  if (!decay1 && !decay2) first=2;

  trees[0]->Store();
  trees[1]->Store();


  int mismatch  = 0;
  bool accepted =1; 
  int caught_jetveto=0;

  for (;;) {
    m_sprime      = (k1->part->Momentum()+k2->part->Momentum()).Abs2();

    accepted = 1;  
    // Parton 1/Tree 1 is the one to decay.
    if (decay1 && caught_jetveto!=3) {
      m_to_be_diced[0]=0;
      msg.Debugging()<<" call A "<<endl;
      if (FillBranch(trees,k1,k2,0)) {
	if (k1->z>0.) {
	  m_sprime = m_sprime/k1->z;
	  msg.Tracking()<<" z="<<k1->z<<"   ("<<k1->kn_no<<") "<<endl;
	}
      }
      else {
	accepted = 0;
      }
    }
    // Parton 2/Tree 2 is the one to decay.    
    if (decay2 && caught_jetveto!=2) {
      m_to_be_diced[1]=0;
      msg.Debugging()<<" call B "<<endl;
      if (FillBranch(trees,k2,k1,1)) {
	if (k2->z > 0.) {
	  m_sprime = m_sprime/k2->z;
	  msg.Tracking()<<" z="<<k2->z<<"   ("<<k2->kn_no<<") "<<endl;
	}
      }
      else {
	accepted = 0;
      }
    }
    
    if (accepted) {

      p_kin->InitKinematics(trees,k1,k2,first);
      if (!decay1) SetColours(k1);
      if (!decay2) SetColours(k2);

      int stat=EvolveSystem(trees,k1,k2);
      if (stat==1) {
	return 1;
      }
      else if (stat==2 || stat==3) {
	if (msg.LevelIsDebugging()) {
	  cout<<" EvolveSystem returned with : "<<stat<<endl;
	  cout<<trees[0];
	  cout<<trees[1];
	}
	if (stat==2) { 
	  m_sprime = m_sprime*k1->z;
	  msg.Tracking()<<" uz="<<k1->z<<" ("<<k1->kn_no<<") "<<endl;
	}
	else {
	  m_sprime = m_sprime*k2->z;
	  msg.Tracking()<<" uz="<<k2->z<<" ("<<k2->kn_no<<") "<<endl;
	}
	caught_jetveto=stat;
      }
      else {
	caught_jetveto=0;
      }
	
    }
    else {
      caught_jetveto=0;
	  
      
    }
    if (caught_jetveto==0) {
      ++mismatch;
      if (mismatch > m_allowed) {
	msg.Error()<<"Error in Initial_State_Shower : "<<std::endl
		   <<"---------------------------------------"<<std::endl
		   <<"Initial_State_Shower::InitializeSystem failed. "
		   <<mismatch<<" trials."<<std::endl
		   <<"---------------------------------------"<<std::endl;
	return 0;
      }
      trees[0]->Restore();
      trees[1]->Restore();
      if (msg.LevelIsTracking()) {
	OutputTree(trees[0]);
	OutputTree(trees[1]);
      }
    }
  }
}

//-----------------------------------------------------------------------
//------------------------ Evolution of the Shower ----------------------
//----------------------------------------------------------------------- 

int Initial_State_Shower::EvolveSystem(Tree ** trees,Knot * k1,Knot * k2)
{
  msg.Debugging()<<"EvolveSystem("<<k1->kn_no<<","<<k2->kn_no<<")  ["
		 <<m_to_be_diced[0]<<","<<m_to_be_diced[1]<<"]"<<endl;
  double sprime_aa = (k1->part->Momentum()+k2->part->Momentum()).Abs2();
  msg.Debugging()<<" sprime="<<m_sprime<<"  "<<sprime_aa<<endl;


  if ((k1->t == k1->tout) && (k2->t == k2->tout)) return 1;  

  bool decay1 = (k1->stat>=1), decay2 = (k2->stat>=1);
  int first = 0;
  if ((!decay1 && decay2 && (k1->t != k1->tout))||(decay1 && !decay2 && (k2->t == k2->tout))) first=1;
  if (!decay1 && !decay2 && ((k1->t == k1->tout)||(k2->t == k2->tout))) first=1;
  if (!decay1 && !decay2) first=2;



  int ntree0=0, ntree1=1;
  if (((k1->t) > (k2->t)) && (k2->t != k2->tout)) {  
    Knot * kh=k1;
    k1 =k2;
    k2 =kh;
    ntree0=1;
    ntree1=0;
  }

  int caught_jetveto=0;

  for (;;) {
    msg.Debugging()<<" ntree="<<ntree0<<endl;
    if (k1->stat>0 || caught_jetveto) {
      msg.Debugging()<<" prepare (I) "<<ntree0<<endl;
      k1->stat           = 0;
      k1->E2             = sqr(k1->part->Momentum()[0]);
      k1->prev->E2       = k1->E2/sqr(k1->z);                      
      k1->prev->left->E2 = k1->prev->E2*sqr(1.-k1->z);
    }  
    else {

      double sprime_a = (k1->part->Momentum()+k2->part->Momentum()).Abs2();
      double sprime_b = (k1->prev->part->Momentum()+k2->part->Momentum()).Abs2();
      k1->z=sprime_a/sprime_b;
      msg.Tracking()<<" detz="<<k1->z<<" ("<<k1->kn_no<<")"<<endl;
      k1->prev->x=k1->x/k1->z;
      m_sprime/=k1->z;
      if (msg.LevelIsDebugging()) {
	cout<<" prepare (II) "<<ntree0<<endl;
	cout<<" ac : "<<k1->part->Momentum()<<" ("<<k1->part->Momentum().Abs2()<<")"<<endl;
	cout<<" mo : "<<k1->prev->part->Momentum()<<" ("<<k1->prev->part->Momentum().Abs2()<<")"<<endl;
	cout<<" pa : "<<k2->part->Momentum()<<" ("<<k2->part->Momentum().Abs2()<<")"<<endl;
	cout<<" sa, sb, z, x "<<sprime_a<<","<<sprime_b<<"  ,"<<k1->z<<","<<k1->prev->x<<endl;
      }

      //double pt2max = sqr(rpa.gen.Ecms());
      //double th     = 4.*k1->z*k1->z*k1->t/(4.*k1->z*k1->z*k1->t-(1.-k1->z)*k1->x*k1->x*pt2max);

      /*  // *AS*
      if (m_to_be_diced[ntree0]) {
	if (msg.LevelIsDebugging()) {
	  cout<<" set scales: "<<k1->t<<" vs. "<<k1->prev->part->Momentum().Abs2()<<endl;
	  cout<<" from "<<endl<<*k1;
	  cout<<" fixing sister "<<endl<<*k1->prev->left;
	  cout<<" fixing mother "<<endl<<*k1->prev<<endl;
	}
	k1->prev->thcrit       = k1->thcrit;
	k1->prev->t            = k1->t;

	k1->prev->left->thcrit = th;  
	k1->prev->left->t      = k1->prev->part->Momentum().Abs2();

	if (msg.LevelIsDebugging()) {
	  cout<<" to "<<endl<<*k1;
	  cout<<" fixing sister "<<endl<<*k1->prev->left;
	  cout<<" fixing mother "<<endl<<*k1->prev<<endl;
	}
      }
      */
    }
    // *AS*    if (k1->prev->stat>0 || m_to_be_diced[ntree0]) {
    if (k1->prev->stat>0) {
      m_to_be_diced[ntree0]=0;
      bool fill=true;
      if (k1->prev->part->Info()=='H') {
	if (k1->prev->prev) {
	  if (k1->prev->prev->part->Info()=='H') {
	    msg.Debugging()<<" internal line "<<endl;
	    fill=false;
	  }
	}
      }
	  

      if (fill) {
	if (!FillBranch(trees,k1->prev,k2,ntree0)) {
	  msg.Debugging()<<" fillbranch failed "<<endl;
	  return 0; 
	}
      }
    }
    double maxt = p_kin->CalculateMaxT(k1,k2);
    if (maxt<k1->prev->left->tout) {
      msg.Debugging()<<" calc max t failed "<<endl;
      // *AS* what about massless ME and PSMass >0 ? 
      return 0;
      if (k1->prev->t==k1->prev->tout) return 0; 
      p_kin->ResetMomenta(k1,trees[ntree0]);
      return (ntree0+2);
    }
    else {
      if (k1->prev->left->stat>0) {
	//cout<<" setting startscale for sister ("<<k1->prev->left->kn_no<<") to : "<<dabs(k1->t)<<" vs. "<<maxt<<endl;
	k1->prev->left->t = dabs(k1->t);
	k1->prev->left->tmax = maxt;
      }
    } 
  
    double sprime_a = (k1->part->Momentum()+k2->part->Momentum()).Abs2();
    // *AS* should always be done 
    if (caught_jetveto==0) {
      if (p_fin) {

	// in case mass already known from ME:
	// determine real daughter momenta first!
	if (!p_kin->DoKinematics(trees,k1,k2,ntree0,first,true)) {
	  msg.Debugging()<<" C2 do kinematics failed "<<endl;
	  return 0;
	}
	p_fin->FirstTimelikeFromSpacelike(trees[ntree0],k1->prev->left,m_jetveto,sprime_a,k1->z);
      }
      else {
	Knot * mo = k1->prev->left;
	mo->t     = mo->tout;
	mo->stat  = 0;
	mo->part->SetStatus(1);
      }
    }


    if (!p_kin->DoKinematics(trees,k1,k2,ntree0,first,false)) {
      msg.Debugging()<<" do kinematics failed "<<endl;
      return 0;
      if (k1->prev->t==k1->prev->tout) return 0; 
      p_kin->ResetMomenta(k1,trees[ntree0]);
      return (ntree0+2);
      //      return 0; 
    }
    if (m_jetveto && m_jetveto_scheme==2 && k1->prev->left->part->Info()!='H') {
    //      if (p_kin->JetVeto(k1->prev->left->part->Momentum()) {
      if (p_kin->JetVeto(k1,k2)) {
	//	if (k2->stat==1) k2->stat=2; 
	//	msg.Debugging()<<" JetVeto : "<<ntree0+2<<endl;
	//	cout<<" JetVeto : "<<ntree0+2<<endl;

	// caught_jetveto=ntree0+2;
	m_to_be_diced[ntree0]=1;
	p_kin->ResetMomenta(k1,trees[ntree0]);
	return (ntree0+2);
      }
    }
    
    if (p_fin) p_fin->SetAllColours(k1->prev->left);
  
    if (k1->prev->z>0.) {
      m_sprime = m_sprime/k1->prev->z;
      msg.Tracking()<<" z="<<k1->prev->z<<"   ("<<k1->prev->kn_no<<") "<<endl;
    }
    
    if (ntree0==0) {
      int stat = EvolveSystem(trees,k1->prev,k2);
      if (stat==0 || stat ==1) return stat;
      caught_jetveto=stat;
    }
    else {  
      int stat=EvolveSystem(trees,k2,k1->prev);
      if (stat==0 || stat ==1) return stat;
      caught_jetveto=stat;
    }
    
    if (caught_jetveto && k1->prev->z>0.) {
      m_sprime = m_sprime*k1->prev->z;
      msg.Tracking()<<" uz="<<k1->prev->z<<"   ("<<k1->prev->kn_no<<") "<<endl;
    }

    if (caught_jetveto!=ntree0+2) {
      if (msg.LevelIsDebugging()) {
	cout<<" jetveto=="<<caught_jetveto<<"  ntree0="<<ntree0<<"  leaving Evolve"<<endl;
	cout<<" k1="<<*k1;
	cout<<" k2="<<*k2;
      }
      m_to_be_diced[ntree0]=1;
      msg.Tracking()<<" rtstart="<<k1->prev->tmax<<" <- "<<k1->prev->t<<" ("<<k1->prev->kn_no<<" )"<<endl;
      k1->prev->t=k1->prev->tmax;
      p_kin->ResetMomenta(k1,trees[ntree0]);
      return caught_jetveto;
    }
    /*
    else {
      //      if (k1->prev && k2->prev)
      if (((k1->prev->t) > (k2->t)) && (k2->t != k2->tout)) {  
	cout<<"SWAP"<<endl;
	Knot * kh=k1;
	k1 =k2;
	k2 =kh;
	int nh=ntree0;
        ntree0=ntree1;
	ntree1=nh;
	
      }
      else {
	cout<<"DONT SWAP"<<endl;
	cout<<" t1="<<k1->t<<" t2="<<k2->t<<endl;
      }
    }
    */
  }
}  

int Initial_State_Shower::FillBranch(Tree ** trees,Knot * active,Knot * partner,int leg) {
  if (msg.LevelIsDebugging()) {
    cout<<" FillBranch "<<endl;
    cout<<*active<<*partner<<endl;
  }
  Flavour flavs[2];
  if (p_suds[leg]->Dice(active,m_sprime,m_jetveto && m_jetveto_scheme==1  ,m_extra_pdf[leg])) {

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
  double pt2max = sqr(rpa.gen.Ecms());
  double th     = 4.*k->z*k->z*k->t/(4.*k->z*k->z*k->t-(1.-k->z)*k->x*k->x*pt2max);

  if (!k->part->Flav().Strong() || !k_flavs[0].Strong() || !k_flavs[1].Strong() ) th=k->thcrit;
  Knot * mother = 0;
  if (k->prev) {
    mother = k->prev;
    // delete old color informations if flav content changed
    if (mother->part->Flav()!=k_flavs[0]) {
      mother->part->SetFlow(1,0);
      mother->part->SetFlow(2,0);
      //      cout<<" mo flav changed from "<<mother->part->Flav()<<" to "<<k_flavs[0]<<endl;
    }
    if (mother->prev) {
      if (msg.LevelIsDebugging()) {
	cout<<" deleting mother tree : "<<endl;
	cout<<*mother<<endl;
      }
      mother->prev=0;
      // *AS* check sure no 'H' can be lost

    }
  }
  else {
    mother  = tree->NewKnot();
    k->prev = mother;
  }

  mother->right  = k;
  mother->part->SetFlav(k_flavs[0]);
  mother->part->SetInfo('I');
  mother->part->SetStatus(1);
  mother->t      = k->t;
  mother->tout   = sqr(k_flavs[0].PSMass()); 
  mother->x      = k->x/k->z;
  mother->stat   = 1;
  mother->E2     = 0.;
  //  mother->thcrit = k->thcrit;
  mother->thcrit = th;

  Knot * sister = 0;
  if (mother->left) {
    sister = mother->left;
    // delete old color informations if flav content changed
    if (sister->part->Flav()!=k_flavs[1]) {
      sister->part->SetFlow(1,0);
      sister->part->SetFlow(2,0);
      //      cout<<" si flav changed from "<<sister->part->Flav()<<" to "<<k_flavs[1]<<endl;
    }
    if (sister->left) {
      msg.Debugging()<<" deleting sister tree : "<<endl<<*sister<<endl;
      sister->left=0;
      sister->right=0;

    }
  }
  else {
    sister=tree->NewKnot();
    sister->prev   = mother;
    mother->left   = sister;
  }
  sister->part->SetFlav(k_flavs[1]);
  sister->part->SetInfo('F');
  sister->part->SetStatus(1);
  sister->t      = 0.;
  sister->tout   = sqr(k_flavs[1].PSMass());
  sister->x      = (mother->x)*(1.-k->z);
  sister->stat   = 1;
  sister->E2     = 0.;
  //  sister->thcrit = k->thcrit; 
  sister->thcrit = th; 

  if (k->part->Info() != 'G' && k->part->Info() != 'H') k->part->SetInfo('i');
  k->part->SetStatus(2);
  SetColours(k);

  if (msg.LevelIsDebugging()) {
    cout<<" FillMotherAndSister: "<<endl;
    cout<<" k  : "<<*k;
    cout<<" mo : "<<*mother;
    cout<<" si : "<<*sister<<endl;
  }
}

void Initial_State_Shower::SetColours(Knot * k)
{

  if (!k) return;
  Knot * mother = k->prev;
  if (!mother) return;
  Knot * sister = mother->left;
  if (!sister) return;

  Knot * test=0;
  
  int all_colors_known=1;
  int nquark=0;
  int ngluon=0;
  for (int i=0;i<3;++i) {
    if (i==0) test = k;
    if (i==1) test = mother;
    if (i==2) test = sister;
    
    if (test->part->Flav().Strong()) {
      int nc=0;
      if (test->part->GetFlow(1)) ++nc;
      if (test->part->GetFlow(2)) ++nc;
      if (test->part->Flav().IsQuark()) {
	++nquark;
	if (nc!=1) {
	  all_colors_known=0;
	}
      }
      else if (test->part->Flav().IsGluon()) {
	++ngluon;
	if (nc!=2) {
	  all_colors_known=0;
	}
      }
      else {
	msg.Error()<<" ERROR: strong particle "<<test->part->Flav()<<" not covered by SetColours "<<endl;
      }
    }      
  }

  if (all_colors_known) {
    return;
  }


  if (nquark+ngluon==2) {
    if (ngluon==2) {
      // a) g -> g + X
      // b) g -> X + g
      // c) X -> g + g
      int col[2];
      bool swap=false;
      if (k->part->Flav().Strong()) {
	swap=true;
	for (size_t i=0; i<2; ++i) col[i]=k->part->GetFlow(i+1);
      }
      else {
	for (size_t i=0; i<2; ++i) col[i]=Flow::Counter();
      }
      if (mother->part->Flav().Strong()) {
	swap=false;
	mother->part->SetFlow(1,col[0]);
	mother->part->SetFlow(2,col[1]);	
      }
      if (sister->part->Flav().Strong()) {
	if (swap) {
	  sister->part->SetFlow(1,col[1]);
	  sister->part->SetFlow(2,col[0]);	
	}
	else {
	  sister->part->SetFlow(1,col[0]);
	  sister->part->SetFlow(2,col[1]);	
	}
      }
    }
    else {
      // d) q/qbar -> q/qbar + X
      // e) q/qbar -> X + qbar
      // f) X -> q/qbar + qbar/q
      int col;
      if (k->part->Flav().Strong()) {
	int anti=k->part->Flav().IsAnti();
	col=k->part->GetFlow(anti+1);
      }
      else {
	col=Flow::Counter();
      }
      if (mother->part->Flav().Strong()) {
	int anti=mother->part->Flav().IsAnti();
	mother->part->SetFlow(anti+1,col);	
      }
      if (sister->part->Flav().Strong()) {
	int anti=sister->part->Flav().IsAnti();
	sister->part->SetFlow(anti+1,col);	
      }
    }
    return;
  }


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
  m_pt2_1 = m_pt2_2 = scale;
  m_th_1  = m_th_2  = M_PI;

  Knot * d1   = trees[0]->NewKnot();
  *(d1->part) = Particle(1,Flavour(kf::u),x1*E*Vec4D(1.,0.,0.,1.));
  d1->part->SetStatus(1);
  d1->part->SetInfo('G');
  d1->part->SetFlow(1,500);
  d1->part->SetFlow(2,501);
  d1->t       = -scale;
  d1->tout    = sqr(Flavour(kf::u).PSMass()); 
  d1->x       = x1;
  d1->E2      = sqr(x1*E);
  d1->maxpt2  = scale;
  d1->costh   = -1.; 
  d1->thcrit  = M_PI;
  d1->stat    = 1;
  
  Knot * d2   = trees[1]->NewKnot();
  *(d2->part) = Particle(2,Flavour(kf::u).Bar(),x2*E*Vec4D(1.,0.,0.,-1.));
  d2->part->SetStatus(1);
  d2->part->SetInfo('G');
  d2->part->SetFlow(1,502);
  d2->part->SetFlow(2,d1->part->GetFlow(1));
  d2->t       = -scale;
  d2->tout    = d1->tout;  
  d2->x       = x2;
  d2->E2      = sqr(x2*E);
  d2->maxpt2  = scale;
  d2->costh   = -1.; 
  d2->thcrit  = M_PI;
  d2->stat    = 1;
}

ATOOLS::Vec4D Initial_State_Shower::GetMomentum(Knot * mo,int & number) 
{
  if (mo->left) return GetMomentum(mo->left,number) + GetMomentum(mo->right,number);
  ++number;
  return mo->part->Momentum();
}


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
  }
}


/*
void Initial_State_Shower::ExtractPartons(Knot * kn,int beam,Blob * jet,
					  Blob_List * bl,Particle_List * pl) 
{
  if (!kn) return;

  // fetch last PSME blobs
  m_bl_meps_is=0;
  m_bl_meps_fs=0;
  for (Blob_Iterator blit=bl->begin();blit!=bl->end();++blit) {
    if ((*blit)->Type().find(string("ME PS Interface (Sherpa, IS)"))!=string::npos) {
      m_bl_meps_is=(*blit);
    }
    if ((*blit)->Type().find(string("ME PS Interface (Sherpa, FS)"))!=string::npos) {
      m_bl_meps_fs=(*blit);
    }
  }
  if (m_bl_meps_is==NULL) {
    ATOOLS::msg.Error()<<"Initial_State_Shower::ExtractPartons(..): "
		       <<"No ME PS Interface found!"<<std::endl
		       <<"   Cannot proceed. Abort."<<std::endl;
    exit(126);
  }
  m_bl_meps_is->SetStatus(0);


  int number=-1;;
  Particle * p=NULL;
  if (!kn->prev) {
// 
//       New jet : kn = incoming parton from hadron info = 'I'
//
    if (kn->part->Info() != 'G') {
      if (pl) number = pl->size();
         else number = int(kn->part);
      kn->part->SetNumber(number);
      kn->part->SetStatus(2);
      if (pl) pl->push_back(kn->part);
    }
    jet = new Blob();
    jet->SetStatus(1);
    p = new Particle(kn->part);
    p->SetStatus(2);
    jet->AddToInParticles(p);
    jet->SetId(bl->size());
    jet->SetType(std::string("IS Shower (APACIC++2.0)"));
    jet->SetBeam(beam);
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
      p = new Particle(kn->part);
      p->SetStatus(2);
      jet->AddToOutParticles(p);
      if (m_bl_meps_is) {
	m_bl_meps_is->AddToInParticles(p);
      }
      jet->SetStatus(1);
      return;
    }
    else if (kn->part->Info()=='H') {
      p = new Particle(kn->part);
      p->SetStatus(2);
      jet->AddToOutParticles(p);
      if (m_bl_meps_is) {
	m_bl_meps_is->AddToInParticles(p);
      }
      jet->SetStatus(1);
    }
  }
  else if (kn->part->Info() == 'H') {
// 
//        New jet : kn = hard parton from ME info = 'H'
//                  and kn outgoing
// 		 or kn->left or kn->right not from ME

//        a) fs  and left and right not from ME add to bl_meps_fs
//              create new blob (FSR like ISR blob)
//        b) fs  and left also ME ignore
//        c) is  and prev and prev->left not from ME add to bl_meps_is->in and to jet->out
//              do NOT create any new blob
//        d) is  AND prev also ME ignore kn

//
    bool is_fs=false;
    Knot * k=kn;
    do {
      if (k==k->prev->left) {
	is_fs=true;
	break;
      }
      k=k->prev;
    } while (k->prev);

    int fillin=0;
    if (is_fs) {
      // fs
      if (!kn->left) fillin=1;
      else if (kn->left->part->Info()!='H') fillin=1;
    }
    else {
      // is
      if (kn->prev->part->Info()!='G' && kn->prev->part->Info()!='H') fillin=2;
    }

    if (fillin) {
      if (fillin==1) {
	jet = new Blob();
	jet->SetStatus(1);
	jet->SetId(bl->size());
	jet->SetType(std::string("IS Shower (APACIC++2.0)"));
	bl->insert(bl->begin(),jet);
      }
      p = new Particle(kn->part);
      p->SetStatus(2);
      if (fillin==1 && m_bl_meps_fs) {
	jet->AddToInParticles(p);
	m_bl_meps_fs->AddToOutParticles(p);
      }
      else if (fillin==2 && m_bl_meps_is) {
	jet->AddToOutParticles(p);
	m_bl_meps_is->AddToInParticles(p);
      }
      if (kn->left) {
	kn->part->SetStatus(2);
      }
      else {
	kn->part->SetStatus(1);
	p = new Particle(kn->part);
	jet->AddToOutParticles(p);
	jet->SetStatus(1);
	return;
      }
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
	p = new Particle(kn->part);
      } 
      else {
	int
	p = new Particle(kn->part);
	if (m_bl_meps_is) {
	  m_bl_meps_is->AddToInParticles(p);
	}
      }
      if (p->Info() == 'G') p->SetStatus(2);
                       else p->SetStatus(1);
      jet->AddToOutParticles(p);
    }
    else {
      kn->part->SetStatus(2);
    }
  }
  ExtractPartons(kn->left,beam,jet,bl,pl); 
  ExtractPartons(kn->right,beam,jet,bl,pl); 
}
*/







