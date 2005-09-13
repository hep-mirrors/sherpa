#include "Final_State_Shower.H"
#include "Run_Parameter.H"
#include "Random.H"

#ifdef PROFILE__all
#include "prof.hh"
#else 
#define PROFILE_HERE
#define PROFILE_LOCAL(LOCALNAME)
#endif

using namespace APACIC;
using namespace ATOOLS;

Final_State_Shower::
Final_State_Shower(MODEL::Model_Base * model,ATOOLS::Jet_Finder * jf,
		   Data_Read *dataread) :
  p_kin(NULL), p_sud(NULL), p_jf(jf),
  m_jetmode(1), m_losejetmode(2), 
  m_maxjets(2)
{ 
  m_pt2min=dataread->GetValue<double>("FS_PT2MIN",1.0);
  m_losejetmode=dataread->GetValue<int>("FS_LOSEJETVETO",0);
  p_kin=new Timelike_Kinematics(jf,m_pt2min);
  p_sud=new Timelike_Sudakov(p_kin,model,m_pt2min);
  p_kin->SetPTCalcScheme(dataread->GetValue<int>("FS_PT_DEFINITION",2),
			 dataread->GetValue<double>("FS_PT_FACTOR",0.25)); 
  p_kin->SetZRangeScheme(dataread->GetValue<int>("FS_Z_RANGES",1));      
  p_kin->SetAngleScheme(dataread->GetValue<int>("FS_ANGLE_DEF",1));  
  p_sud->SetOrderingScheme(dataread->GetValue<int>("FS_ORDERING",1));  
  p_sud->SetCouplingScheme(dataread->GetValue<int>("FS_COUPLINGS",1));
  p_sud->SetMassScheme(dataread->GetValue<int>("FS_MASS_SCHEME",1));   
  p_sud->SetWidthScheme(dataread->GetValue<int>("FS_WIDTH_SCHEME",0));   
  p_sud->
    SetMECorrectionScheme(dataread->GetValue<int>("FS_ME_CORRECTIONS",0)); 
  p_sud->SetQEDScheme(dataread->GetValue<int>("FS_PHOTONS",0));        
  p_sud->SetCorrelationScheme(dataread->GetValue<int>("FS_ANGLE_CORR",0));
  p_sud->Init();
}

Final_State_Shower::~Final_State_Shower() 
{
  delete p_sud;
  delete p_kin;
}

int Final_State_Shower::PerformShower(Tree *tree,int jetveto) 
{
  PROFILE_HERE;
  // m_jetmode=jetveto>0
  switch(InitializeJets(tree,tree->GetRoot())) {
  case -1:
    return -1;
  case 1: {
    int jetcheck(TestKinematics(tree,1));
    if (jetcheck!=1) return jetcheck;
    if (p_kin->DoKinematics(tree->GetRoot())) return 1;
    msg.Error()<<METHOD<<"(..,"<<jetveto<<"): "
	       <<"Kinematics failed."<<std::endl;
    return 0;
  }
  case 0:
    msg.Error()<<METHOD<<"(..,"<<jetveto<<"): "
	       <<"Shower evolution failed."<<std::endl;
    return 0;
  }
  return 0;
}

int Final_State_Shower::
TimelikeFromSpacelike(Tree *tree,Knot *mo,bool jetveto,
		      double sprime,double z)
{
  msg_Debugging()<<METHOD<<"(..,"<<mo->kn_no<<","<<jetveto<<","
		 <<sprime<<","<<z<<"): {\n";
  msg_Indent();
  m_jetmode=2;
  if (mo->left && mo->right) {
    EstablishRelations(mo,mo->left,mo->right);
    return InitializeJets(tree,mo,1);
  }
  else {
    Flavour flavs[2];
    Simple_Polarisation_Info polinfos[2];
    if (p_sud->Dice(mo)) {
      // set temporary kinematics, real kinematics set in is shower
      mo->E2 = sqr(((1./z-1.)*sprime - mo->t)/(2.*sqrt(sprime)));
      mo->part->SetMomentum(Vec4D(sqrt(mo->E2),0.,0.,sqrt(mo->E2-mo->t)));
      switch(TestKinematics(tree)) {
      case 1: {
	flavs[0] = p_sud->GetFlB();
	flavs[1] = p_sud->GetFlC();
	InitDaughters(tree,mo,flavs,polinfos,1);
	mo->left->t=mo->left->tout;
	mo->right->t=mo->right->tout;
	int stat(EvolveJet(tree,mo));
	msg_Debugging()<<"}\n";
	return stat;
      }
      case -1: 
	msg_Debugging()<<"}\n";
	return -1;
      }
    } 
    msg_Debugging()<<"reset knot "<<mo->kn_no<<"\n";
    Reset(mo);
    msg_Debugging()<<"}\n";
    return 1;
  }
}

int Final_State_Shower::InitializeJets(Tree *tree,Knot *mo,int init_rel)
{
  msg_Debugging()<<METHOD<<"(["<<mo->kn_no<<"],"<<init_rel<<"): {\n";
  msg_Indent();
  m_jetmode=1;
  if (!mo || !mo->left || !mo->right) {
    msg.Error()<<METHOD<<"(..): Error. No knots."<<std::endl;
    return 0;
  }
  Knot *d1(mo->left), *d2(mo->right);
  // store old values
  double z_tmp(mo->z);
  Vec4D p1_tmp(d1->part->Momentum()), p2_tmp(d2->part->Momentum());
  // decide, whether to branch 
  bool mustbranch1(d1->stat>0), mustbranch2(d2->stat>0), cont(true);
  int first(1);
  if (mo==tree->GetRoot() && mustbranch1 && mustbranch2) first=2; 
  if (mustbranch1 || mustbranch2) {
    SmearDaughters(mo);
    int accept1(1), accept2(1);
    do {
      // dice
      switch (FillBranch(tree,mo,first)) {
      case 1:
	// evolve daughters
	if (d1->t > d2->t) {
	  if (mustbranch1) accept1=EvolveJet(tree,d1);
	  if (mustbranch2 && accept1!=-1) accept2=EvolveJet(tree,d2);
	}
	else {
	  if (mustbranch2) accept2=EvolveJet(tree,d2);
	  if (mustbranch1 && accept2!=-1) accept1=EvolveJet(tree,d1);
	}
	if (accept1==1 && accept2==1) {
	  // finished branching
	  msg_Debugging()<<"accept jets\n";
	  cont = false;
	  break;
	} 
	else if (accept1==-1 || accept2==-1) {
	  // caught losejet veto
	  return -1;
	}
	else {
	  // caught jet veto
	  mo->z  = z_tmp;
	  d1->E2 = mo->z*mo->z*mo->E2;
	  d2->E2 = (1.-mo->z)*(1.-mo->z)*mo->E2;
	  d1->part->SetMomentum( p1_tmp );
	  d2->part->SetMomentum( p2_tmp );
	  msg_Debugging()<<"restore jets "<<d1->kn_no
			 <<" and "<<d2->kn_no<<"\n";
	  msg_Debugging()<<"  p_d1 = "<<p1_tmp<<" "<<p1_tmp.Abs2()
			 <<" vs. "<<d1->t<<"\n";
	  msg_Debugging()<<"  p_d2 = "<<p2_tmp<<" "<<p2_tmp.Abs2()
			 <<" vs. "<<d2->t<<"\n";
	}
      case 0:
	if (d1->stat==0 && d2->stat==0) {
	  // reset branching
	  p_kin->Shuffle(mo,1);
	  p_kin->UpdateDaughters(mo);
	  cont = false;
	}
	break;
      case -1:
	msg_Debugging()<<"}\n";
	return -1;
      }
    } while(cont);
  }
  msg_Debugging()<<"out ("<<d1->kn_no<<","<<d1->stat<<") "
		 <<mustbranch1<<", ("<<d2->kn_no<<","<<d2->stat<<") "
		 <<mustbranch2<<std::endl;
  if (!mustbranch1) {
    if (init_rel) EstablishRelations(d1,d1->left,d1->right);
    switch(InitializeJets(tree,d1)) {
    case -1: 
      msg_Debugging()<<"}\n";
      return -1;
    case 0:  
      msg_Debugging()<<"}\n";
      return 0;
    }
  }
  if (!mustbranch2) {
    if (init_rel) EstablishRelations(d2,d2->left,d2->right);
    switch(InitializeJets(tree,d2)) {
    case -1: 
      msg_Debugging()<<"}\n";
      return -1;
    case 0:  
      msg_Debugging()<<"}\n";
      return 0;
    }
  }
  msg_Debugging()<<"}\n";
  return 1;
}

int Final_State_Shower::FillBranch(Tree *tree,Knot *mo,int first)
{
  msg_Debugging()<<METHOD<<"(["<<mo->kn_no<<","<<mo->stat<<"],"
		 <<first<<"): p_mo = "<<mo->part->Momentum()<<" {\n";
  msg_Indent();
  if (!first && mo->t<=mo->tout) {
    msg_Debugging()<<"!first & "<<mo->t<<"<"<<mo->tout<<std::endl; 
    return 0; 
  }
  Knot *d1(mo->left), *d2(mo->right);
  msg_Debugging()<<"p_d1 = "<<d1->part->Momentum()<<std::endl;
  msg_Debugging()<<"p_d2 = "<<d2->part->Momentum()<<std::endl;
  Flavour d1_flavs[2], d2_flavs[2];
  Simple_Polarisation_Info d1_polinfos[2], d2_polinfos[2];
  bool do12, diced1(false), diced2(false);
  Vec4D p1(d1->part->Momentum()), p2(d2->part->Momentum());
  Knot *g(NULL);
  while (true) {
    g=NULL;
    do12=ChooseDaughter(mo);
    msg_Debugging()<<"selected daughter "<<(do12?"2":"1")<<": "
		   <<(do12?d2->part->Flav():d1->part->Flav())<<std::endl;
    if (first==2) g=mo;
    if (!do12) {
      ResetDaughters(d1);
      if (p_sud->Dice(d1,g)) { 
	d1_flavs[0]    = p_sud->GetFlB();
	d1_flavs[1]    = p_sud->GetFlC();
	d1_polinfos[0] = p_sud->GetPolB();
	d1_polinfos[1] = p_sud->GetPolC();
	d1->stat       = 1;
	diced1         = true;
      }   
      else Reset(d1);
    }
    else {
      ResetDaughters(d2);  
      if (p_sud->Dice(d2,g)) {
	d2_flavs[0]    = p_sud->GetFlB();
	d2_flavs[1]    = p_sud->GetFlC();
	d2_polinfos[0] = p_sud->GetPolB();
	d2_polinfos[1] = p_sud->GetPolC();
	d2->stat       = 1;
	diced2         = true;
      }    
      else Reset(d2);
    }
    msg_Debugging()<<"test emission at ("<<(do12?d2->t:d1->t)
		   <<","<<(do12?d2->z:d1->z)<<")\n";
    p_kin->CheckZRange(mo,d1_flavs,d2_flavs);
    if (d1->stat!=3 && d2->stat!=3) {
      if (p_kin->Shuffle(mo,first)) {
	switch (TestKinematics(tree)) { 
	case 1:
	  msg_Debugging()<<"kinematics check passed"<<std::endl;
	  if (d1->stat>0)
	    InitDaughters(tree,d1,d1_flavs,d1_polinfos,diced1);
	  if (d2->stat>0) 
	    InitDaughters(tree,d2,d2_flavs,d2_polinfos,diced2);	  
	  mo->stat=0;
	  msg_Debugging()<<"}\n";
	  return 1;
	case -1: 
	  msg_Debugging()<<"kinematics check failed\n";
	  msg_Debugging()<<"}\n";
	  return -1;
	}
	msg_Debugging()<<"kinematics vetoed\n";
      }
      if (d1->didkin) d1->part->SetMomentum(p1);
      if (d2->didkin) d2->part->SetMomentum(p2);
    }
    if (d1->stat==0 && d2->stat==0) break;
  }
  msg_Debugging()<<"found no possible branch\n";
  msg_Debugging()<<"}\n";
  p_kin->Shuffle(mo,first);
  p_kin->UpdateDaughters(mo);
  return 0;
}

int Final_State_Shower::EvolveJet(Tree *tree,Knot *mo)
{
  msg_Debugging()<<METHOD<<"(["<<mo->kn_no<<","<<mo->part->Flav()<<"]): "
		 <<"p_mo = "<<mo->part->Momentum()<<" "
		 <<mo->part->Momentum().Abs2()<<", t_mo = "<<mo->t<<" {\n";
  msg_Indent();
  if (!mo->stat) {
    msg_Debugging()<<"stat = 0\n";
    ResetDaughters(mo);
    msg_Debugging()<<"}\n";
    return true;
  }
  double z_tmp(mo->z);
  Knot *d1(mo->left), *d2(mo->right);
  int evolve1, evolve2;
  while (true) {
    switch (FillBranch(tree,mo,0)) {
    case 1:
      if (d1->t > d2->t) { 
	evolve1 = EvolveJet(tree,d1);
	evolve2 = EvolveJet(tree,d2);
      }
      else { 
	evolve1 = EvolveJet(tree,d2);
	evolve2 = EvolveJet(tree,d1);
      }
      if (evolve1==1 && evolve2==1) {
	msg_Debugging()<<"}\n";
	return 1; 
      }
      else if (evolve1==-1 || evolve2==-1) {
	msg_Debugging()<<"}\n";
	return -1;
      } 
      // reset to original (massless) kinematics
      mo->z  = z_tmp;
      d1->E2 = mo->z*mo->z*mo->E2;
      d2->E2 = (1.-mo->z)*(1.-mo->z)*mo->E2;
      break;
    case 0:
      Reset(mo);
      msg_Debugging()<<"reset knot "<<mo->kn_no<<std::endl;
      msg_Debugging()<<"}\n";
      return 0;
    case -1:
      msg_Debugging()<<"terminate EvolveJets."<<std::endl;
      msg_Debugging()<<"}\n";
      return -1;
    }
  }
  msg_Debugging()<<"}\n";
  return 0;
}

bool Final_State_Shower::SetAllColours(Knot *mo) 
{ 
  return SetColours(mo,p_kin); 
}

bool Final_State_Shower::SetColours(Knot *mo,Timelike_Kinematics *kin)
{
  if (!mo) {
    msg.Error()<<METHOD<<"(..): Error. Void mother knot."<<std::endl;
    return false;
  }
  if ( (!mo->left) && (!mo->right) ) return true;
  if ( !mo->left ) {
    mo->right->part->SetFlow(1,mo->part->GetFlow(1));
    mo->right->part->SetFlow(2,mo->part->GetFlow(2));
    return true;
  }
  if ( !mo->right ) {
    mo->left->part->SetFlow(1,mo->part->GetFlow(1));
    mo->left->part->SetFlow(2,mo->part->GetFlow(2));
    return true;
  }
  //  check if already enough colours
  Knot *test(NULL);
  int all_colors_known(1);
  for (int i=0;i<3;++i) {
    if (i==0) test = mo;
    if (i==1) test = mo->left;
    if (i==2) test = mo->right;
    if (test->part->Flav().Strong()) {
      int nc(0);
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
	msg.Error()<<METHOD<<"(..): Error.\n  Strong particle "
		   <<test->part->Flav()<<" (n_c = "<<nc
		   <<") not covered by SetColours. "<<std::endl;
      }
    }      
  }
  Knot *d1(mo->left), *d2(mo->right);
  if (!all_colors_known) {
    Knot * partner(NULL), * nopart(NULL);
    if (mo->part->Flav().Strong()) {
      if (mo->part->Flav().IsQuark()) {
	partner = d1; nopart = d2;
	if (d2->part->Flav().IsQuark()) {
	  partner = d2;
	  nopart  = d1;
	}
	if (partner->part->Flav().IsQuark() && 
	    partner->part->Flav().IsAnti() && 
	    nopart->part->Flav().IsGluon()) {
	  partner->part->SetFlow(2,-1);
	  nopart->part->SetFlow(1,partner->part->GetFlow(2));
	  nopart->part->SetFlow(2,mo->part->GetFlow(2));
	} 
	if (partner->part->Flav().IsQuark() && 
	    !partner->part->Flav().IsAnti() && 
	    nopart->part->Flav().IsGluon()) {
	  partner->part->SetFlow(1,-1);
	  nopart->part->SetFlow(1,mo->part->GetFlow(1));
	  nopart->part->SetFlow(2,partner->part->GetFlow(1));
	} 
	if (partner->part->Flav().IsQuark() && 
	    !nopart->part->Flav().Strong()) {
	  partner->part->SetFlow(1,mo->part->GetFlow(1));
	  partner->part->SetFlow(2,mo->part->GetFlow(2));
	}
      } 
      else if (mo->part->Flav().IsGluon()) {
	if (mo->prev) {
	  if (d1->part->Flav().IsQuark() && 
	      d2->part->Flav().IsQuark()) {
	    if (d1->part->Flav().IsAnti()) {
	      d1->part->SetFlow(2,mo->part->GetFlow(2));
	      d2->part->SetFlow(1,mo->part->GetFlow(1));
	    }
	    else if (d2->part->Flav().IsAnti()) {
	      d2->part->SetFlow(2,mo->part->GetFlow(2));
	      d1->part->SetFlow(1,mo->part->GetFlow(1));
	    }
	  }
	  else if (d1->part->Flav().IsGluon() && 
		   d2->part->Flav().IsGluon()) {
	    Particle *aup(FindAuntParton(mo));
	    partner = d1; 
	    nopart = d2;
	    if (aup->Flav().Strong() && kin) {
	      if (kin->ArrangeColourPartners(aup,d1,d2)) { 
		partner = d2; 
		nopart = d1; 
	      }
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
	  else {
	    if (d1->part->Flav().IsGluon() && 
		!d2->part->Flav().Strong()) {
	      d1->part->SetFlow(1,mo->part->GetFlow(1));
	      d1->part->SetFlow(2,mo->part->GetFlow(2));
	    }
	    else if (d2->part->Flav().IsGluon() && 
		     !d1->part->Flav().Strong()) {
	      d2->part->SetFlow(1,mo->part->GetFlow(1));
	      d2->part->SetFlow(2,mo->part->GetFlow(2));
	    }
	  }
	}
      }
      else {
	msg.Out()<<METHOD<<"(..): Error.\n  Coloured "<<mo->part->Flav()
		 <<" -> "<<d1->part->Flav()<<" + "
		 <<d2->part->Flav()<<std::endl;
	return 0;
      }
    }
    else {
      // colour neutral mother
      if (d1->part->Flav().Strong() && d2->part->Flav().Strong()) {      
	if (d1->part->Flav().IsQuark() && d2->part->Flav().IsQuark()) {
	  partner = d1; 
	  nopart = d2;
	  if (d1->part->Flav().IsQuark() && d1->part->Flav().IsAnti() && 
	      d2->part->Flav().IsQuark() && !d2->part->Flav().IsAnti()) {
	    partner = d2; 
	    nopart = d1;
	  }
	  partner->part->SetFlow(1,-1);
	  partner->part->SetFlow(2,0);
	  nopart->part->SetFlow(2,partner->part->GetFlow(1));
	  nopart->part->SetFlow(1,0);
	}
	else if (d1->part->Flav().IsGluon() && d2->part->Flav().IsGluon()) {
	  d1->part->SetFlow(1,-1);
	  d1->part->SetFlow(2,-1);
	  d2->part->SetFlow(2,d1->part->GetFlow(1));
	  d2->part->SetFlow(1,d1->part->GetFlow(2));
	}
	else {
	  msg.Out()<<METHOD<<"(..): Error.\n  Colourless "<<mo->part->Flav()
		   <<" -> "<<d1->part->Flav()<<" + "
		   <<d2->part->Flav()<<std::endl;
	  return 0;
	}
      }
    }
  }
  return ( SetColours(d1,kin) && SetColours(d2,kin) );
}

void Final_State_Shower::ExtractPartons(Knot *kn,Blob *jet,
					Blob_List *bl,Particle_List *pl) 
{
  // fetch last ME PS blob
  Blob *bl_meps=NULL;
  for (Blob_List::iterator blit=bl->begin();blit!=bl->end();++blit) {
    if((*blit)->Type()==btp::ME_PS_Interface_FS) {
      bl_meps=*blit;
    }
  }
  if (bl_meps==NULL) {
    ATOOLS::msg.Error()<<METHOD<<"(..): Error.\n  No ME PS Interface. "
		       <<"Abort."<<std::endl;
    abort();
  }
  // deactivate in partons!
  for (int i=0;i<bl_meps->NInP();++i) bl_meps->InParticle(i)->SetStatus(2);
  if (!kn) return;
  Particle * p(NULL);
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
      p = new Particle(*kn->part);
      jet->AddToInParticles(p);
      if (bl_meps) {
	bl_meps->AddToOutParticles(p);
	bl_meps->SetStatus(0);
      }

      p = new Particle(*kn->part);
      jet->AddToOutParticles(p);
      if (pl) p->SetNumber(pl->size());
         else p->SetNumber(0);
      kn->part->SetNumber(p->Number());
      jet->SetId();
      jet->SetType(btp::FS_Shower);
      jet->SetTypeSpec("APACIC++2.0");
      //      jet->SetPosition(p->XProd() + Vec4D(p->LifeTime(),p->Distance()));
      bl->push_back(jet);
      return;
    }
    else {
      if ((kn->left->part->Info() != 'H') || (kn->right->part->Info() != 'H')) {
	jet = new Blob();
	jet->SetStatus(1);
	p = new Particle(*kn->part);
      	p->SetStatus(2);
	if (pl) pl->push_back(p);
	jet->AddToInParticles(p);
	if (bl_meps) {
	  bl_meps->AddToOutParticles(p);
	  bl_meps->SetStatus(0);
	}
	if (pl) p->SetNumber(pl->size());
	   else p->SetNumber(0);
        kn->part->SetNumber(p->Number());
	jet->SetId();
	jet->SetType(btp::FS_Shower);
	jet->SetTypeSpec("APACIC++2.0");
	//	jet->SetPosition(p->XProd() + Vec4D(p->LifeTime(),p->Distance()));
	bl->push_back(jet);
      }
    }
  }
  else {
    if (!kn->left) {
      if (!jet) {
	msg.Error()<<"ERROR in Final_State_Shower ::ExtractPartons :\n"
		   <<"    No jet for Parton : "<<kn->part->Number()<<std::endl;
	abort();
      }
      if (pl) kn->part->SetNumber(pl->size());
	 else kn->part->SetNumber(0);
      kn->part->SetStatus(1);
      if (pl) pl->push_back(kn->part);
      jet->AddToOutParticles(new Particle(*kn->part));
    }
  }
  ExtractPartons(kn->left,jet,bl,pl); 
  ExtractPartons(kn->right,jet,bl,pl); 
}


//-----------------------------------------------------------------------
//------------------- Helpers for initialisation of the Shower-----------
//----------------------------------------------------------------------- 

bool Final_State_Shower::SmearDaughters(Knot *mo) 
{
  Knot * d1(mo->left), * d2(mo->right);
  if (d1->part->Flav().IsStable() && 
      d2->part->Flav().IsStable()) return false;
  if (d1->stat==0 && d2->stat==0) return false;
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
  while (true) {
    if (!d1->part->Flav().IsStable() && d1->stat)
      t1 = mass12 + mass1*width1 * tan(ran.Get()*yrange1 + ymin1);
    if (!d2->part->Flav().IsStable() && d2->stat)
      t2 = mass22 + mass2*width2 * tan(ran.Get()*yrange2 + ymin2);
    if (t1+t2+sqrt(2.*t1*t2) < mo->t) {
      d1->tout = t1;
      d2->tout = t2;
      return true;
    }
  }
}

void Final_State_Shower::EstablishRelations(Knot *mo, Knot *d1,Knot *d2) 
{
  if (!d1 || !d2 || !mo) {
    msg.Error()<<METHOD<<"(..): Warning. called with\n";
    if (mo) msg.Error()<<"mo :"<<*mo<<"\n"; else msg.Error()<<"mo : 0x0\n";
    if (d1) msg.Error()<<"d1 :"<<*d1<<"\n"; else msg.Error()<<"d1 : 0x0\n";
    if (d2) msg.Error()<<"d2 :"<<*d2<<std::endl; else msg.Error()<<"d2 : 0x0"<<std::endl;
    return;
  }
  // set color connections (if not jet known)
  APACIC::Final_State_Shower::SetColours(mo,0);
  
  double t_mo(mo->part->Momentum().Abs2()), 
    E_mo(mo->part->Momentum()[0]), z_mo(d1->part->Momentum()[0]/E_mo); 
  double th(sqrt(dabs(t_mo)/(z_mo*(1.- z_mo)))/E_mo),
    maxpt2(sqr(ATOOLS::Min(z_mo,1.-z_mo))/(z_mo*(1.-z_mo))*t_mo);
  if (mo->part->Flav().IsQuark() && 
      d1->part->Flav().Strong() && d2->part->Flav().Strong()) {
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
    d1->maxpt2 = d2->maxpt2 = maxpt2;
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
      d1->maxpt2 = d2->maxpt2 = maxpt2;
    }
    else if (!(d1->part->Flav().Strong()) && (d2->part->Flav().Strong())) {
      d1->t      = t_mo;
      d1->thcrit = M_PI;
      d2->t      = mo->t;
      d2->thcrit = mo->thcrit;
      d1->maxpt2 = d2->maxpt2 = mo->maxpt2;
    }
    else if ((d1->part->Flav().Strong()) && !(d2->part->Flav().Strong())) {
      d1->t      = mo->t;
      d1->thcrit = mo->thcrit;
      d2->t      = t_mo;
      d2->thcrit = M_PI;
      d1->maxpt2 = d2->maxpt2 = mo->maxpt2;
    }
  }
  else {
    if ((d1->part->Flav().Strong()) && (d2->part->Flav().Strong())) {
      d1->t      = t_mo;
      d1->thcrit = th;
      d2->t      = t_mo;
      d2->thcrit = th;
      d1->maxpt2 = d2->maxpt2 = t_mo;
    }
    else {
      d1->t      = t_mo;
      d1->thcrit = M_PI;
      d2->t      = t_mo;
      d2->thcrit = M_PI;
      d1->maxpt2 = d2->maxpt2 = t_mo;
    }
  }
  mo->t = t_mo; 
}

bool Final_State_Shower::ChooseDaughter(Knot * mo)
{
  if (mo->left->stat==3)  return 0;
  if (mo->right->stat==3) return 1;

  if (mo->left->stat && 
      (!mo->right->stat || mo->right->t<mo->right->tout)) return 0;// left
  if (mo->right->stat && 
      (!mo->left->stat || mo->left->t<mo->left->tout)) return 1;// right
  double tm1(Min(mo->t,mo->left->E2)), tm2(Min(mo->t,mo->right->E2));
  if (mo->left->t/tm1 > mo->right->t/tm2) return 0;// left
  return 1;// right
 
}

void Final_State_Shower::
InitDaughters(Tree * tree,Knot * mo,Flavour * mo_flavs, 
	      Simple_Polarisation_Info * mo_pols, bool diced) 
{ 
  if (!mo->left) {
    mo->left          = tree->NewKnot();
    mo->right         = tree->NewKnot();
  }
  if (diced) {
    mo->left->prev     = mo;
    mo->left->polinfo  = mo_pols[0];
    mo->left->part->SetFlav(mo_flavs[0]);
    mo->left->part->SetInfo('F');
    mo->left->part->SetStatus(1);
    mo->left->tout     = sqr(mo_flavs[0].PSMass());
    mo->left->stat     = 3;  

    mo->right->prev    = mo;
    mo->right->polinfo = mo_pols[1];
    mo->right->part->SetFlav(mo_flavs[1]);
    mo->right->part->SetInfo('F');
    mo->right->part->SetStatus(1);
    mo->right->tout    = sqr(mo_flavs[1].PSMass());
    mo->right->stat    = 3;  

    if (mo->part->Info()!='H') mo->part->SetInfo('f');
    mo->part->SetStatus(2);
  }

  // Reset kinematics
  double th, maxpt2;
  if ((mo->left->part->Flav().Strong()) && 
      (mo->right->part->Flav().Strong())) {
    th     = mo->thcrit = p_kin->CalculateAngle(mo);
    maxpt2 = p_kin->CalcKt2(mo->z,mo->E2,mo->t,
			    mo->left->part->Flav().PSMass(),
			    mo->right->part->Flav().PSMass());
  }
  else {
    if (mo->prev) {
      mo->thcrit = th  = mo->prev->thcrit;
      maxpt2           = mo->prev->maxpt2;
    }
    else {
      mo->thcrit = th = M_PI;
      maxpt2     = mo->maxpt2;
    }
  }

  mo->left->t       = mo->t; 
  mo->left->E2      = (mo->z)*(mo->z)*(mo->E2); 
  mo->left->thcrit  = th;
  mo->left->maxpt2  = maxpt2;

  mo->right->t      = mo->t; 
  mo->right->E2     = (1. - mo->z)*(1. - mo->z)*(mo->E2);
  mo->right->thcrit = th;
  mo->right->maxpt2 = maxpt2; 
}

//#####################################################################
//
// Tester methods
//
//#####################################################################

bool Final_State_Shower::TestShower(Tree * tree) 
{
  double E2 = sqr(rpa.gen.Ecms());
  for (int i=1;i<=rpa.gen.NumberOfEvents();i++) {
    if (i%2500==0) {
      msg.Out()<<" "<<i<<" th event "<<std::endl;
    }
    
    tree->Reset();
    InitTwojetTree(tree,E2);
    if (!PerformShower(tree,0)) return 0; 
    if (msg.LevelIsTracking()) OutputTree(tree);
  }

  return 1;
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
}

//#####################################################################
//
// Helper methods
//
//#####################################################################


void Final_State_Shower::ResetDaughters(Knot * mo)
{
  Reset(mo->left);
  Reset(mo->right);
  if (mo->part->Info()!='H') mo->part->SetInfo('F');
  mo->part->SetStatus(1);
} 

void Final_State_Shower::Reset(Knot * mo) 
{ 
  if (!mo) return;
  mo->left  = NULL;  
  mo->right = NULL;
  mo->t     = mo->tout;
  mo->stat  = 0;
  mo->part->SetStatus(1);
  if (mo->part->Info()!='H') mo->part->SetInfo('F');
}

void Final_State_Shower::OutputTree(Tree * tree) 
{
  if (tree->GetRoot()==0) {
    msg.Out()<<"Empty Tree\n";
  }
  else {
    int number(0);
    msg.Out()<<"Final Tree:\n"<<tree<<std::endl
	     <<"Total 4 Mom = "<<GetMomentum(tree->GetRoot(),number);
    msg.Out()<<"   for "<<number<<" FS particles.\n";
  }
}


Vec4D  Final_State_Shower::GetMomentum(Knot * mo, int & number) {
  if (mo->left) {
    Vec4D p     = GetMomentum(mo->left,number) + GetMomentum(mo->right,number);
    Vec4D ptest = mo->left->part->Momentum() + mo->right->part->Momentum();
    if (!(ptest==mo->part->Momentum())) {
      number -= 10000;
      msg.Error()<<"--------------------------------------------------------"<<std::endl
		 <<"WARNING in Final_State_Shower ("<<number<<") :\n"
		 <<"   Momentum conservation violation in Knot :"<<mo->kn_no<<std::endl
		 <<"          is: "<<ptest<<"  "<<p.Abs2()<<std::endl
		 <<"   should be:"<<mo->part->Momentum()<<"  "<<mo->part->Momentum().Abs2()
		 <<std::endl<<std::endl
		 <<GetMomentum(mo->left,number)<<" + "<<GetMomentum(mo->right,number)<<" vs. "<<std::endl
		 <<mo->left->part->Momentum()<<" + "<<mo->right->part->Momentum()<<std::endl
		 <<"   from : "<<mo->left->part->Flav()<<" / "<<mo->right->part->Flav()<<std::endl;
    }
    return p;
  }
  number++;
  return mo->part->Momentum();
}


Particle * Final_State_Shower::FindAuntParton(Knot * mo) 
{
  Knot * au(mo->prev->left);
  if (au==mo) au = mo->prev->right;
  for (int k1=1;k1<=2;++k1) {
    for (int k2=1;k2<=2;++k2) {
      if ((mo->part->GetFlow(k1)>0) &&
	  (au->part->GetFlow(3-k1)==mo->part->GetFlow(k1))) return au->part;  
    }
  }
  Blob * bl(mo->part->ProductionBlob());
  if (bl) {
    Particle * aup(NULL);
    for (int i=0; i<bl->NInP();++i) {
      aup=bl->InParticle(i);
      for (int k1=1;k1<=2;++k1) {
	for (int k2=1;k2<=2;++k2) {
	  if ((mo->part->GetFlow(k1)>0) &&
	      (aup->GetFlow(k2)==mo->part->GetFlow(k1))) return aup;
	}
      }
    }
    for (int i=0; i<bl->NOutP();++i) {
      aup=bl->OutParticle(i);
      for (int k1=1;k1<=2;++k1)
	for (int k2=1;k2<=2;++k2)
	  if ((mo->part->GetFlow(k1) > 0 ) &&
	    (aup->GetFlow(3-k1)==mo->part->GetFlow(k1))) return aup;
    }
  }
  msg_Tracking()<<METHOD<<"(..): No blob for mother "<<mo->kn_no
		<<"\n   Return normal aunt."<<std::endl;
  return au->part;
}

int Final_State_Shower::TestKinematics(Tree *const tree,const int mode)
{
  msg_Debugging()<<METHOD<<"(..): {"<<std::endl;
  msg_Indent();
  size_t hard(0);
  std::vector<Vec4D> jets;
  if (!CollectMomenta(tree->GetRoot(),jets,hard)) return 0;
  int jfmode(p_jf->Type());
  p_jf->SetType(1);
  bool cluster(true);
  while (cluster) {
    cluster=false;
    double pt2min(p_jf->ShowerPt2());
    std::vector<Vec4D>::iterator jwin(jets.end()), kwin(jets.end());
    for (std::vector<Vec4D>::iterator 
	   jit(jets.begin());jit!=jets.end();++jit) {
      for (std::vector<Vec4D>::iterator kit(jit);kit!=jets.end();++kit) {
	if (kit==jit) ++kit;
	if (kit==jets.end()) break;
	double pt2jk(p_jf->MTij2(*jit,*kit));
	if (pt2jk<pt2min) {
	  pt2min=pt2jk;
	  jwin=jit;
	  kwin=kit;
	}
      }
    }
    if (jwin!=jets.end() && kwin!=jets.end()) {
      *jwin+=*kwin;
      jets.erase(kwin);
      cluster=true;
    }
  }
  if (m_jetmode==2) {
    for (std::vector<Vec4D>::iterator 
	   jit(jets.begin());jit!=jets.end();++jit) {
      double pt2j(p_jf->MTij2(Vec4D(1.,0.,0.,1.),*jit));
      if (pt2j<p_jf->ShowerPt2()) jit=--jets.erase(jit);
    }
  }
  p_jf->SetType(jfmode);
  msg_Debugging()<<"produced "<<jets.size()
		 <<" jets out of "<<hard<<" losejetmode = "
		 <<m_losejetmode<<std::endl;
  if (jets.size()>hard) {
    if (hard<=m_maxjets) {
      msg_Debugging()<<"produced "<<(jets.size()-hard)
		     <<" additional jets"<<std::endl;
      msg_Debugging()<<"}\n";
      if (mode) return -1;
      return 0;
    }
  }
  else if (jets.size()<hard && mode && m_losejetmode>1) {
    msg_Debugging()<<"lost "<<(hard-jets.size())
		   <<" jets"<<std::endl;
    msg_Debugging()<<"}\n";
    return -1;
  }
  msg_Debugging()<<"}\n";
  return 1;
}

int Final_State_Shower::CollectMomenta(Knot *knot,std::vector<Vec4D> &vecs,
				       size_t &hard)
{
  msg_Debugging()<<"Final_State_Shower::CollectMomenta("
		 <<knot->kn_no<<"): {\n";
  msg_Indent();
  if (knot->left) p_kin->DoSingleKinematics(knot);
  if ((knot->left==NULL || knot->left->stat==3)
      &&knot->part->Flav().Strong()) {
    msg_Debugging()<<"take mom "<<knot->kn_no<<" -> "
		   <<knot->part->Momentum()<<"\n";
    vecs.push_back(knot->part->Momentum());
  }
  int dtest(1);
  size_t rhard(hard);
  if (knot->left && knot->left->stat!=3) {
    dtest=CollectMomenta(knot->left,vecs,hard);
    if (dtest==1) 
      dtest=CollectMomenta(knot->right,vecs,hard);
  }
  if (rhard==hard && knot->part->Info()=='H'
      && knot->part->Flav().Strong()) {
    msg_Debugging()<<"hard "<<knot->kn_no<<std::endl;
    ++hard;
  }
  msg_Debugging()<<"}\n";
  return dtest;
}

