#include "Final_State_Shower.H"
#include "Blob_List.H"
#include "Blob.H"
#include "Particle_List.H"
#include "Tree.H"
#include "Timelike_Sudakov.H"
#include "Timelike_Kinematics.H"
#include "Run_Parameter.H"

#include "Random.H"

using namespace APACIC;
using namespace ATOOLS;
using namespace std;

//-----------------------------------------------------------------------
//--------------------------- Constructors ------------------------------
//----------------------------------------------------------------------- 

Final_State_Shower::Final_State_Shower(MODEL::Model_Base * _model,Data_Read * const _dataread) 
{
  m_pt2min = _dataread->GetValue<double>("FS PT2MIN",1.);
  p_kin    = new Timelike_Kinematics(m_pt2min,_dataread);
  p_sud    = new Timelike_Sudakov(p_kin,m_pt2min,_model,_dataread);
}

Final_State_Shower::~Final_State_Shower() 
{
  if (p_sud) delete p_sud;
  if (p_kin) delete p_kin;
}

//-----------------------------------------------------------------------
//----------------------- Performing the Shower -------------------------
//----------------------------------------------------------------------- 
void Final_State_Shower::SetJetvetoPt2(const double pt2) { 
  p_kin->SetJetvetoPt2(pt2); 
}

int Final_State_Shower::PerformShower(Tree * tree,int _jetveto) 
{
  bool jetveto = (_jetveto>0);

  p_kin->SetJetVeto(jetveto);

  m_ini_partons.clear();
  int stat=InitializeJets(tree,tree->GetRoot());

  if (stat) {
    if (!p_kin->DoKinematics(tree->GetRoot())) {
      msg.Error()<<"Error in Final_State_Shower : "<<std::endl
		 <<"Final_State_Shower::PerformShower : "
		 <<"Kinematics did not work out."<<std::endl;

      return 0;
    }
    return 1;
  }
  else {
    msg.Error()<<"Error in Final_State_Shower : "<<std::endl
	       <<"Final_State_Shower::PerformShower : "
	       <<"Initializing the jets did not work out !!!"<<std::endl;
    return 0;
  }
}

void Final_State_Shower::FirstTimelikeFromSpacelike(Tree * tree,Knot* mo,bool jetveto,double sprime,double z)
{
  p_kin->SetJetVeto(jetveto);

  Knot * gr = mo ->prev;
  Knot * si = gr ->right;
  if (mo->thcrit!=M_PI) {
    // updateing theta crit
    double th1c   = sqrt( dabs(si->t)/((1.- si->z)*si->E2));
    mo->thcrit=th1c;
  }

  if (mo->left && mo->right) {
    EstablishRelations(mo,mo->left,mo->right);

    int stat = InitializeJets(tree,mo,1);
  }
  else {
    if (mo->part->Info()=='H') {
      int found=0;
      for (int i=0;i<m_ini_partons.size();++i) {
	if (m_ini_partons[i]==mo) found=1;
      }
      if (!found) {
	m_ini_partons.push_back(mo);
      }
    }

    Flavour flavs[2];
    for (;;) {
      if (p_sud->Dice(mo)) {
	double test_e4  =((1./z-1.)*sprime - mo->t)/(2.*sqrt(sprime));
	mo->E2=sqr(test_e4);

	// init daughters
	flavs[0]  = p_sud->GetFlB();
	flavs[1]  = p_sud->GetFlC();
	InitDaughters(tree,mo,flavs,1);
	if (EvolveJet(tree,mo)) {
	  return;
	}
      }
    
      Reset(mo);

      return;
    }
    
  }
}

//-----------------------------------------------------------------------
//---------------------------- After the Shower -------------------------
//----------------------------------------------------------------------- 

bool Final_State_Shower::SetAllColours(Knot * mo) {
  return SetColours(mo,p_kin);
}


bool Final_State_Shower::SetColours(Knot * mo, Timelike_Kinematics * kin)
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

  //  check if already enough colours
  Knot * test=0;
  
  int all_colors_known=1;
  for (int i=0;i<3;++i) {
    if (i==0) test = mo;
    if (i==1) test = mo->left;
    if (i==2) test = mo->right;
    
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
	msg.Error()<<" ERROR: strong particle "<<test->part->Flav()<<" not covered by SetColours "<<endl;
      }
    }      
  }

  Knot * d1 = mo->left;
  Knot * d2 = mo->right;


  if (all_colors_known) return ( SetColours(d1,kin) && SetColours(d2,kin) );

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
	  Particle * aup = FindAuntParton(mo);

	  partner = d1; nopart = d2;
	  if (aup->Flav().Strong() && kin) {
	    if (kin->ArrangeColourPartners(aup,d1,d2)) { partner = d2; nopart = d1; }
	    int cset=0;
 	    for (int i=1;i<3;i++) {
	      if (aup->GetFlow(i) == mo->part->GetFlow(3-i)) {
 		partner->part->SetFlow(3-i,mo->part->GetFlow(3-i));
		partner->part->SetFlow(i,-1);
		nopart->part->SetFlow(3-i,partner->part->GetFlow(i));
		nopart->part->SetFlow(i,mo->part->GetFlow(i));
		break;
	      }
	    }
	    if (!cset) {
	      for (int i=1;i<3;i++) {
		if (aup->GetFlow(i) == mo->part->GetFlow(i)) {
		  partner->part->SetFlow(i,mo->part->GetFlow(i));
		  partner->part->SetFlow(3-i,-1);
		  nopart->part->SetFlow(i,partner->part->GetFlow(3-i));
		  nopart->part->SetFlow(3-i,mo->part->GetFlow(3-i));
		  break;
		}
	      }
	    }
	  }
	  else {  // ie. single gluon (from hard event - connected to initial states)
	    if (ran.Get()<0.5) {
	      partner=d1;
	      nopart=d2;
	    }
	    else {
	      partner=d2;
	      nopart=d1;
	    }

	    partner->part->SetFlow(1,mo->part->GetFlow(1));
	    partner->part->SetFlow(2,-1);
	    nopart->part->SetFlow(1,partner->part->GetFlow(2));
	    nopart->part->SetFlow(2,mo->part->GetFlow(2));
	  }
	}
      }
    }
    else {
      msg.Out()<<"ERROR in Final_State_Shower::SetColours()"<<std::endl
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
	msg.Out()<<"ERROR in Final_State_Shower::SetColours()"<<std::endl
		   <<"    case not covered : colourless "<<mo->part->Flav()<<" -> "
		   <<d1->part->Flav()<<" + "<<d2->part->Flav()<<std::endl;
	return 0;
      }
    }
  }

  return ( SetColours(d1,kin) && SetColours(d2,kin) );
}


void Final_State_Shower::EstablishRelations(Knot * mo, Knot * d1,Knot * d2) {
  if (!d1 || !d2 || !mo) {
    msg.Error()<<" WARNING:Final_State_Shower::EstablishRelations() called with"<<endl;
    if (mo) msg.Error()<<"mo :"<<*mo<<endl; else msg.Error()<<"mo : 0x0"<<endl;
    if (d1) msg.Error()<<"d1 :"<<*d1<<endl; else msg.Error()<<"d1 : 0x0"<<endl;
    if (d2) msg.Error()<<"d2 :"<<*d2<<endl; else msg.Error()<<"d2 : 0x0"<<endl;
    return;
  }
  // set color connections (if not jet known)
  APACIC::Final_State_Shower::SetColours(mo,0);
  
  double t_mo = mo->part->Momentum().Abs2();
  double E_mo= mo->part->Momentum()[0];
  double th  = sqrt( t_mo/(mo->z*(1.- mo->z)))/E_mo;
  if (mo->part->Flav().IsQuark() && d1->part->Flav().Strong() && d2->part->Flav().Strong()) {
    if (d1->part->Flav().IsQuark()) {
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
  else  if (mo->part->Flav().Strong()) {
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
}


void Final_State_Shower::ExtractPartons(Knot * kn,Blob * jet,Blob_List * bl,Particle_List * pl) 
{
  // fetch last ME PS blob
  Blob *bl_meps=NULL;
  for (Blob_Iterator blit=bl->begin();blit!=bl->end();++blit) {
    if((*blit)->Type()==btp::ME_PS_Interface_FS) {
      bl_meps=*blit;
    }
  }
  if (bl_meps==NULL) {
    ATOOLS::msg.Error()<<"Final_State_Shower::ExtractPartons(..): "
		       <<"   No ME PS Interface found. Cannot proceed. Abort."<<std::endl;
    exit(127);
  }
  // deactivate in partons!
  for (int i=0;i<bl_meps->NInP();++i) {
    bl_meps->InParticle(i)->SetStatus(2);
  }

  if (!kn) return;
  int number;
  Particle * p = 0;
  if (kn->part->Info()=='H') {
    /* 
       New jet : kn = hard parton from ME info = 'HF'
                 and kn outgoing
		 or kn->left or kn->right not from ME
    */
    if (!(kn->left)) {
      if (pl) pl->push_back(kn->part);
      jet = new Blob();
      jet->SetStatus(1);
      p=new Particle(kn->part);
      jet->AddToInParticles(p);
      if (bl_meps) {
	bl_meps->AddToOutParticles(p);
	bl_meps->SetStatus(0);
      }

      p=new Particle(kn->part);
      jet->AddToOutParticles(p);
      if (pl) number = pl->size();
      else number = (long int)(kn->part);
      p->SetNumber(number);

      kn->part->SetNumber(number);
      jet->SetId(bl->size());
      jet->SetType(btp::FS_Shower);
      jet->SetTypeSpec("APACIC++2.0");
      jet->SetPosition(p->XProd() + Vec4D(p->LifeTime(),p->Distance()));
      bl->push_back(jet);
      return;
    }
    else {
      if ((kn->left->part->Info() != 'H') || (kn->right->part->Info() != 'H')) {
	jet = new Blob();
	jet->SetStatus(1);
	p = new Particle(kn->part);
      	p->SetStatus(2);
	if (pl) pl->push_back(p);
	jet->AddToInParticles(p);
	if (bl_meps) {
	  bl_meps->AddToOutParticles(p);
	  bl_meps->SetStatus(0);
	}
	if (pl) number = pl->size();
	else number = (long int)(kn->part);
	p->SetNumber(number);
	kn->part->SetNumber(number);
	jet->SetId(bl->size());
	jet->SetType(btp::FS_Shower);
	jet->SetTypeSpec("APACIC++2.0");
	jet->SetPosition(p->XProd() + Vec4D(p->LifeTime(),p->Distance()));
	bl->push_back(jet);
      }
    }
  }
  else {
    if (!kn->left) {
      if (!jet) {
	msg.Error()<<"ERROR in Final_State_Shower ::ExtractPartons :"<<std::endl
		   <<"    No jet for Parton : "<<kn->part->Number()<<std::endl;
	abort();
      }
      if (pl) number = pl->size();
         else number = (long int)(kn->part);
      kn->part->SetNumber(number);
      kn->part->SetStatus(1);
      if (pl) pl->push_back(kn->part);
      jet->AddToOutParticles(new Particle(kn->part));
    }
  }
  ExtractPartons(kn->left,jet,bl,pl); 
  ExtractPartons(kn->right,jet,bl,pl); 
}


void Final_State_Shower::ExtractPartons(Knot * kn, Particle_List * pl) 
{
  if (!kn) return;
  if (kn->left==0) {
    if (pl) pl->push_back(new Particle(kn->part));
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
  double E2 = sqr(rpa.gen.Ecms());
  double E  =rpa.gen.Ecms()*0.5;
  for (int i=1;i<=rpa.gen.NumberOfEvents();i++) {
    if (i%2500==0) {
      msg.Out()<<" "<<i<<" th event "<<std::endl;
    }
    
    tree->Reset();
    InitTwojetTree(tree,E2);
    if (!PerformShower(tree,0)) return 0;     // *AS*   (tree,0) = no jetveto! (default)
    if (msg.LevelIsTracking()) OutputTree(tree);
  }

  return 1;
}

//-----------------------------------------------------------------------
//------------------- Initialisation of the Shower ----------------------
//----------------------------------------------------------------------- 
int Final_State_Shower::InitializeJets(Tree * tree,Knot * mo,int init_rel)
{
  if (!mo) {
    msg.Error()<<"Error in Final_State_Shower : "<<std::endl
	       <<"   Final_State_Shower::InitializeJets : No mother found !"<<std::endl;
    return 0;
  }
  if ((!mo->left) || (!mo->right)) {
    msg.Error()<<"Error in Final_State_Shower : "<<std::endl
	       <<"   Final_State_Shower::InitializeJets : No daughters found !"<<std::endl
	       <<" for "<<*mo<<std::endl;
    return 0;
  }


  double z_tmp = mo->z;
  Knot * d1    = mo->left;
  Knot * d2    = mo->right;
  Vec4D p1_tmp = d1->part->Momentum();
  Vec4D p2_tmp = d2->part->Momentum();

  bool decay1 = (d1->stat>0), decay2 = (d2->stat>0);

  int first=1;
  if (mo==tree->GetRoot() && decay1 && decay2) first=2; 
  if ((decay1) || (decay2)) {
    bool accept = SmearDaughters(mo);
    for (;;) {
      if (FillBranch(tree,mo,first)) {
	accept = 1;
	if (decay1) accept = EvolveJet(tree,d1) && accept;
	if (decay2) accept = EvolveJet(tree,d2) && accept;
	if (accept) break; 
	else {
	  mo->z  = z_tmp;
	  d1->E2 = mo->z*mo->z*mo->E2;
	  d2->E2 = (1.-mo->z)*(1.-mo->z)*mo->E2;
	  d1->part->SetMomentum( p1_tmp );
	  d2->part->SetMomentum( p2_tmp );
	}
      }
      if (!(d1->stat) && !(d2->stat)) break;
    }
  }

  int ok=1; 
  if (!decay1) {
    if (init_rel) EstablishRelations(d1,d1->left,d1->right);

    ok = ok && InitializeJets(tree,d1);
  }
  else {
    if (d1->part->Flav().Strong()) {
      int found=0;
      for (int i=0;i<m_ini_partons.size();++i) {
	if (m_ini_partons[i]==d1) found=1;
      }
      if (!found) {
	m_ini_partons.push_back(d1);
      }
    }
  }
  if (!decay2) {
    if (init_rel) EstablishRelations(d2,d2->left,d2->right);

    ok = ok && InitializeJets(tree,d2);
  }
  else {
    if (d2->part->Flav().Strong()) {
      int found=0;
      for (int i=0;i<m_ini_partons.size();++i) {
	if (m_ini_partons[i]==d2) found=1;
      }
      if (!found) {
	m_ini_partons.push_back(d2);
      }
    }
  }

  if (ok==0) return 0;  // did not work out!
  return 1;             // everything ok.
}

bool  Final_State_Shower::ExtraJetCheck(Knot * mo, Knot * d1, Knot * d2) {
  return p_kin->ExtraJetCheck(0,d1,d2); 
}

bool  Final_State_Shower::ExtraJetCheck() const {
  bool test=1;

  if (m_ini_partons.size()>=2) {
    for (int i=0;i<m_ini_partons.size()-1;++i) {
      for (int j=i+1;j<m_ini_partons.size();++j) {
	test=test & p_kin->ExtraJetCheck(0,m_ini_partons[i],m_ini_partons[j]);
	if (test==0) break;
      }
      if (test==0) break;
    }
  }
  else if (m_ini_partons.size()==1) {
    test=p_kin->ExtraJetCheck(0,0,m_ini_partons[0]);
  }

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
      d1->tout = t1;
      d2->tout = t2;
      return 1;
    }
  }
}

void Final_State_Shower::InitTwojetTree(Tree * tree,double scale) {
  double start_th=200;

  Knot * mo   = tree->NewKnot();
  *(mo->part) = Particle(1,Flavour(kf::photon),Vec4D(sqrt(scale),0,0,0));
  mo->part->SetStatus(2);
  mo->part->SetInfo('M');
  mo->t       = scale;
  mo->E2      = scale;
  mo->maxpt2  = 0.;
  mo->z       = 0.5;
  mo->costh   = -1.; 
  mo->thcrit  = start_th;
  mo->phi     = 2.*M_PI*ran.Get();
  mo->stat    = 1;  

  //  The daughters are defined and distributed according to
  //  1+cos^2 theta w.r.t. the z-axis.

  Flavour mo_flavs[2];
  for(;;) {
    mo_flavs[0] = Flavour(kf::code(1+int(ran.Get()*4.)));   
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
  *(mo->left->part)    = Particle(2,mo_flavs[0],Vec4D(E,p1));
  mo->left->part->SetStatus(1);
  mo->left->part->SetInfo('H');
  mo->left->part->SetFlow(1,-1);
  mo->left->t          = mo->t;
  mo->left->tout       = sqr(mo_flavs[0].PSMass());
  mo->left->maxpt2     = 0.;


  mo->left->E2         = E*E;
  mo->left->thcrit     = start_th;
 
  mo->right            = tree->NewKnot();
  mo->right->prev      = mo;
  mo->right->stat      = 3;     
  *(mo->right->part) = Particle(3,mo_flavs[1],Vec4D(E,(-1.)*p1)); 
  mo->right->part->SetStatus(1);
  mo->right->part->SetInfo('H');
  mo->right->part->SetFlow(2,mo->left->part->GetFlow(1));
  mo->right->t         = mo->t;
  mo->right->tout      = sqr(mo_flavs[1].PSMass());
  mo->right->maxpt2    = 0.;
  mo->right->E2        = E*E;
  mo->right->thcrit    = start_th;
};

//-----------------------------------------------------------------------
//------------------------ Evolution of the Shower ----------------------
//----------------------------------------------------------------------- 

bool Final_State_Shower::EvolveJet(Tree * tree,Knot* mo)
{
  if (!mo->stat) {
    ResetDaughters(mo);
    return 1;
  }

  double z_tmp = mo->z;
  bool accept1 ,accept2;
  for (;;) {
    if (FillBranch(tree,mo,0)) {
      // select one/two daughter systems and try to evolve them 
      accept1 = EvolveJet(tree,mo->left);
      accept2 = EvolveJet(tree,mo->right);
      if (accept1 && accept2) return 1; 
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
	mo->stat  = 3;  
      }
      else {
	Reset(mo);
      }
      return 0; 
    }
  }
}

bool Final_State_Shower::FillBranch(Tree * tree,Knot* mo,int first)
{
  Knot * d1 = mo->left;
  Knot * d2 = mo->right;
  if (!(first) && (mo->t <= mo->tout) ) return 0; 


  Flavour d1_flavs[2];
  Flavour d2_flavs[2];
  
  bool do12;
  bool diced1=0, diced2=0;

  for (;;) {
    Knot * g =0;
    // select one daugther to dice down
    do12 = ChooseDaughter(mo);
    if (first==2) g = mo;
    if (!do12) {
      ResetDaughters(d1);  // if grand children already determined, delete them
      if (p_sud->Dice(d1,g)) { // determine t,z, and flavours
	d1_flavs[0]  = p_sud->GetFlB();
	d1_flavs[1]  = p_sud->GetFlC();
	d1->stat=1;
	diced1=1;
      }   
      else Reset(d1);
    }
    else {
      ResetDaughters(d2);  // if grand children already determined, delete them
      if (p_sud->Dice(d2,g)) { // determine t,z, and flavours
	d2_flavs[0]  = p_sud->GetFlB();
	d2_flavs[1]  = p_sud->GetFlC();
	d2->stat=1;
	diced2=1;
      }    
      else Reset(d2);
    }

    // check zrange first
    if (first==2) 
      p_kin->CheckZRange(mo);

    if ((d1->stat != 3) && (d2->stat != 3)) {
      if (p_kin->Shuffle(mo,first)) { 
	if (d1->stat) {
	  InitDaughters(tree,d1,d1_flavs,diced1); 
	  diced1=0;
	}
	if (d2->stat) {
	  InitDaughters(tree,d2,d2_flavs,diced2); 
	  diced2=0;
	}
	return 1; 
      }
    }
    // exit if both daughters can not be diced down further
    if (!(d1->stat) && !(d2->stat)) return 0; 
  }
}

//-----------------------------------------------------------------------
//------------------- Service for the Branchings ------------------------
//----------------------------------------------------------------------- 

bool Final_State_Shower::ChooseDaughter(Knot * mo)
{
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
  if ((mo->left->t/tm1) > (mo->right->t/tm2)) return 0; // left leg
  return 1; // right leg
 
}

void Final_State_Shower::InitDaughters(Tree * tree,Knot * mo,Flavour * mo_flavs, bool diced) 
{ 
  if (!mo->left) {
    mo->left        = tree->NewKnot();
    mo->right       = tree->NewKnot();
  }
  if (diced) {
    mo->left->prev  = mo;
    mo->left->part->SetFlav(mo_flavs[0]);
    mo->left->part->SetInfo('F');
    mo->left->part->SetStatus(1);
    mo->left->tout  = sqr(mo_flavs[0].PSMass());
    mo->left->stat  = 3;  

    mo->right->prev = mo;
    mo->right->part->SetFlav(mo_flavs[1]);
    mo->right->part->SetInfo('F');
    mo->right->part->SetStatus(1);
    mo->right->tout = sqr(mo_flavs[1].PSMass());
    mo->right->stat = 3;  

    if (mo->part->Info() != 'H') mo->part->SetInfo('f');
    mo->part->SetStatus(2);
  }

  // Reset kinematics
  mo->left->t          = mo->t; 
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
  Reset(mo->left);
  Reset(mo->right);

  if (mo->part->Info() != 'H') mo->part->SetInfo('F');
  mo->part->SetStatus(1);
} 

void Final_State_Shower::Reset(Knot * mo) 
{ 
  if (!mo) return;
  mo->left  = 0;  
  mo->right = 0;

  mo->t     = mo->tout;
  mo->stat  = 0;
  mo->part->SetStatus(1);
  if (mo->part->Info() != 'H') mo->part->SetInfo('F');
}

void Final_State_Shower::OutputTree(Tree * tree) 
{
  int number = 0;
  if (tree->GetRoot()==0) {
    msg.Out()<<"empty Tree"<<endl;
  }
  else {
    msg.Out()<<"final Tree:"<<std::endl<<tree<<std::endl
	     <<"Total 4 Mom = "<<GetMomentum(tree->GetRoot(),number);
    msg.Out()<<" for "<<number<<" FS particles."<<std::endl;
  }
}


Vec4D  Final_State_Shower::GetMomentum(Knot * mo, int & number) {
  if (mo->left) {
    Vec4D p     = GetMomentum(mo->left,number) + GetMomentum(mo->right,number);
    Vec4D ptest = mo->left->part->Momentum() + mo->right->part->Momentum();
    if (!(ptest==mo->part->Momentum())) 
      msg.Out()<<"WARNING in Final_State_Shower :"<<std::endl
	       <<"   Momentum conservation violation in Knot :"<<mo->kn_no<<std::endl
	       <<"          is: "<<ptest<<"  "<<p.Abs2()<<std::endl
	       <<"   should be:"<<mo->part->Momentum()<<"  "<<mo->part->Momentum().Abs2()<<std::endl;
    return p;
  }
  number++;
  return mo->part->Momentum();
}


Particle * Final_State_Shower::FindAuntParton(Knot * mo) 
{
  Knot * au = mo->prev->left;
  if (au == mo) au = mo->prev->right;
  bool found=0;
  for (int k1=1;k1<=2;++k1)
    for (int k2=1;k2<=2;++k2)
      if ((mo->part->GetFlow(k1) > 0 ) &&
	  (au->part->GetFlow(3-k1)==mo->part->GetFlow(k1))) found=1;
	  
  if (found) return au->part;

  Blob * bl = mo->part->ProductionBlob();
  if (!bl) {
    return au->part;
  }
  
  Particle * aup=0;
  for (int i=0; i<bl->NInP();++i) {
    aup=bl->InParticle(i);
    for (int k1=1;k1<=2;++k1)
      for (int k2=1;k2<=2;++k2)
	if ((mo->part->GetFlow(k1) > 0 ) &&
	    (aup->GetFlow(k2)==mo->part->GetFlow(k1))) found=1;
    if (found) return aup;
  }

  return au->part;
}
