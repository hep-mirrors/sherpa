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

bool Initial_State_Shower::PerformShower(Tree ** trees,bool _jetveto) {
  m_jetveto = _jetveto;

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
					  Blob_List * bl,Parton_List * pl) 
{
  if (!kn) return;

  // fetch PSME blob
  Blob * bl_meps=0;
  for (Blob_Iterator blit=bl->begin();blit!=bl->end();++blit) {
    int pos = (*blit)->Type().find(string("ME PS Interface (Sherpa)"));
    if (pos>-1) {
      bl_meps=(*blit);
      bl_meps->SetStatus(0);
      break;
    }
  }

  int number;
  Parton * p;
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
    p = new Parton(kn->part);
    p->SetDecayBlob(jet);
    p->SetStatus(2);
    jet->AddToInPartons(p);
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
      p = new Parton(kn->part);
      p->SetStatus(2);
      jet->AddToOutPartons(p);
      p->SetProductionBlob(jet);
      if (bl_meps) {
	p->SetDecayBlob(bl_meps);
	bl_meps->AddToInPartons(p);
      }
      jet->SetStatus(1);
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
    p = new Parton(kn->part);
    p->SetStatus(2);
    jet->AddToInPartons(p);
    p -> SetDecayBlob(jet);
    if (bl_meps) {
      p -> SetProductionBlob(bl_meps);
      bl_meps->AddToOutPartons(p);
    }
    jet->SetId(bl->size());
    jet->SetType(std::string("IS Shower (APACIC++2.0)"));
    bl->insert(bl->begin(),jet);
    if (kn->left) {
      kn->part->SetDecayBlob(jet);
      kn->part->SetStatus(2);
    }
    else {
      kn->part->SetStatus(1);
      p = new Parton(kn->part);
      p->SetProductionBlob(jet);
      p->SetDecayBlob(NULL);
      jet->AddToOutPartons(p);
      jet->SetStatus(1);
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
	p = new Parton(kn->part);
      } 
      else {
	p = new Parton(kn->part);
	if (bl_meps) {
	  p -> SetDecayBlob(bl_meps);
	  bl_meps->AddToInPartons(p);
	}
	else p -> SetDecayBlob(NULL);
      }
      p->SetProductionBlob(jet);
      if (p->Info() == 'G') p->SetStatus(2);
                       else p->SetStatus(1);
      jet->AddToOutPartons(p);
    }
    else {
      kn->part->SetStatus(2);
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
  m_extra_pdf[0]    = 1;
  m_extra_pdf[1]    = 1;
  m_to_be_diced[0]  = 1;
  m_to_be_diced[1]  = 1;

  if ( (!k1) || (!k2) ) {
    msg.Error()<<"ERROR: Initial_State_Shower::InitializeSystem : No trees found !"<<std::endl;
    exit(1);
    return 0;
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
	m_to_be_diced[1]=0;
	if (FillBranch(trees,k1,k2,0)) {
	  if (k1->z>0.) m_sprime = m_sprime/k1->z;
	}
	else {
	  accepted = 0;
	}
    }
    // Parton 2/Tree 2 is the one to decay.    
    if (decay2 && caught_jetveto!=2) {
	m_to_be_diced[1]=0;
	if (FillBranch(trees,k2,k1,1)) {
	  if (k2->z > 0.) m_sprime = m_sprime/k2->z;
	}
	else {
	  accepted = 0;
	}
    }
    
    if (accepted) {

      p_kin->InitKinematics(trees,k1,k2,first||caught_jetveto);
      if (!decay1) SetColours(k1);
      if (!decay2) SetColours(k2);

      int stat=EvolveSystem(trees,k1,k2);
      if (stat==1) {
	return 1;
      }
      else if (stat==2 || stat==3) {
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
      if (rpa.gen.Events()) {
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
  bool decay1 = (k1->stat>=1), decay2 = (k2->stat>=1);
  int first = 0;
  if ((!decay1 && decay2)||(decay1 && !decay2)) first=1;
  if (!decay1 && !decay2) first=2;

  if ((k1->t == k1->tout) && (k2->t == k2->tout)) return 1;  


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
    if (k1->stat>0 || caught_jetveto) {
      k1->stat           = 0;
      k1->E2             = sqr(k1->part->Momentum()[0]);
      k1->prev->E2       = k1->E2/sqr(k1->z);                      
      k1->prev->left->E2 = k1->prev->E2*sqr(1.-k1->z);
    }  
    else {
      double sprime_a = (k1->part->Momentum()+k2->part->Momentum()).Abs2();
      double sprime_b = (k1->prev->part->Momentum()+k2->part->Momentum()).Abs2();
      k1->z=sprime_a/sprime_b;
      k1->prev->x=k1->x/k1->z;
      m_sprime/=k1->z;


      double pt2max = sqr(rpa.gen.Ecms());
      double th     = 4.*k1->z*k1->z*k1->t/(4.*k1->z*k1->z*k1->t-(1.-k1->z)*k1->x*k1->x*pt2max);

      if (m_to_be_diced[ntree0]) {
	k1->prev->thcrit       = k1->thcrit;
	k1->prev->t            = k1->t;
	k1->prev->left->thcrit = th;  
	k1->prev->left->t      = k1->prev->part->Momentum().Abs2();
      }
    }
    if (k1->prev->stat>0 || caught_jetveto) {
      m_to_be_diced[ntree0]=0;
      if (k1->prev->stat!=2) {  
	if (!FillBranch(trees,k1->prev,k2,ntree0)) return 0; 
      }
    }
    double maxt = p_kin->CalculateMaxT(k1,k2);
    if (maxt<k1->prev->left->tout) {
      return 0; 
    }
    else {
      if (k1->prev->left->stat>0)
	k1->prev->left->t = maxt;
    } 
  
    double sprime_a = (k1->part->Momentum()+k2->part->Momentum()).Abs2();
    if (caught_jetveto==0) {
      if (k1->prev->stat!=2) { 
	p_fin->FirstTimelikeFromSpacelike(trees[ntree0],k1->prev->left,m_jetveto,sprime_a,k1->z);
      }
    }

    if (k1->prev->stat==2) {
      k1->prev->stat=1;
    }

    if (!p_kin->DoKinematics(trees,k1,k2,ntree0,first)) {
      return 0; 
    }
    if (m_jetveto && m_jetveto_scheme==2 && p_kin->JetVeto(k1->prev->left->part->Momentum())) {
      if (k2->stat==1) k2->stat=2; 
      return (ntree0+2);
    }

    p_fin->SetAllColours(k1->prev->left);
  
    if (k1->prev->z>0.) m_sprime = m_sprime/k1->prev->z;

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

    if (caught_jetveto!=ntree0+2) return caught_jetveto;
  }
}  

int Initial_State_Shower::FillBranch(Tree ** trees,Knot * active,Knot * partner,int leg) {
  Flavour flavs[2];
  if (p_suds[leg]->Dice(active,m_sprime,m_jetveto && m_jetveto_scheme==1  ,m_extra_pdf[leg])) {

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
  Knot * mother = 0;
  if (k->prev) {
    mother = k->prev;
    if (mother->prev) {
      mother->prev=0;
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
  mother->thcrit = k->thcrit;

  Knot * sister = 0;
  if (mother->left) {
    sister = mother->left;
    if (sister->left) {
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
  sister->thcrit = k->thcrit; 

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
	}
      }
      else if (test->part->Flav().IsGluon()) {
	if (nc!=2) {
	  all_colors_known=0;
	  break;
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
  *(d1->part) = Parton(1,Flavour(kf::u),x1*E*Vec4D(1.,0.,0.,1.));
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
  *(d2->part) = Parton(2,Flavour(kf::u).Bar(),x2*E*Vec4D(1.,0.,0.,-1.));
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







