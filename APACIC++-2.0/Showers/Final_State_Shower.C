#include "Final_State_Shower.H"
#include "Initial_State_Shower.H"
#include "Run_Parameter.H"
#include "Random.H"
#include "MyStrStream.H"
#include "Veto_Info.H"

#ifdef PROFILE__all
#include "prof.hh"
#else 
#define PROFILE_HERE
#define PROFILE_LOCAL(LOCALNAME)
#endif

using namespace APACIC;
using namespace ATOOLS;

Final_State_Shower::Final_State_Shower(MODEL::Model_Base *const model,
				       ATOOLS::Jet_Finder *const jf,
				       Data_Read *const dataread):
  p_kin(new Timelike_Kinematics(jf)), 
  p_sud(new Timelike_Sudakov(p_kin,model)), 
  p_jv(NULL)
{ 
  p_sud->SetScaleFactor
    (ToType<double>(rpa.gen.Variable("FS_CPL_SCALE_FACTOR","1.0")));
  p_sud->SetEvolutionScheme
    (dataread->GetValue<int>("FS_EVOLUTION_SCHEME",0));  
  p_kin->SetZScheme(dataread->GetValue<int>("FS_Z_SCHEME",1));      
  p_sud->SetOrderingScheme(dataread->GetValue<int>("FS_ORDERING_SCHEME",1));  
  p_sud->SetCouplingScheme(dataread->GetValue<int>("FS_COUPLING_SCHEME",1));
  p_sud->SetMassScheme(dataread->GetValue<int>("FS_MASS_SCHEME",2));   
  p_sud->SetWidthScheme(dataread->GetValue<int>("FS_WIDTH_SCHEME",1));   
  p_sud->SetMECorrectionScheme(dataread->GetValue<int>("FS_ME_SCHEME",0)); 
  p_sud->SetQEDMECorrectionScheme(dataread->GetValue<int>("FS_QED_ME_SCHEME",0)); 
  p_sud->SetCorrelationScheme(dataread->GetValue<int>("FS_CORR_SCHEME",0));
  p_sud->SetKTScheme(dataread->GetValue<int>("FS_KT_SCHEME",1));
  p_sud->SetQEDScheme(dataread->GetValue<int>("FS_QED_SCHEME",0));        
  p_sud->SetPT2Min(dataread->GetValue<double>("FS_PT2MIN",1.0));
  p_sud->SetPT2MinQED(dataread->GetValue<double>("FS_PT2MIN_QED",.0025));
  p_sud->Init(dataread->GetValue<double>("F_MEDIUM",0.0));
  m_showermode=ToType<int>(rpa.gen.Variable("SHOWER_MODE"));
  msg_Debugging()<<METHOD<<"(): Intermediate showering is "
		 <<m_showermode<<std::endl;
}

Final_State_Shower::~Final_State_Shower() 
{
  delete p_sud;
  delete p_kin;
}

int Final_State_Shower::PerformShower(Tree *tree,int jetveto) 
{
  PROFILE_HERE;
#ifdef USING__Veto_Info
  p_sud->ClearVetos();
  p_sud->SetMode(0);
#endif
  tree->GetRoot()->Store();
  if (InitializeJets(tree,tree->GetRoot())) {
    if (p_kin->DoKinematics(tree->GetRoot())) return 1;
    msg.Error()<<METHOD<<"("<<jetveto<<"): "
	       <<"Kinematics failed."<<std::endl;
    return 0;
  }
  else {
    msg.Error()<<METHOD<<"("<<jetveto<<"): "
	       <<"Shower evolution failed."<<std::endl;
    return 0;
  }
  return 0;
}

int Final_State_Shower::
FillISBranch(Initial_State_Shower *const ini,
	     Tree *tree,Knot *mo,const bool jetveto,
	     const double &sprime,const double &z,Knot *partner)
{
#ifdef USING__Veto_Info
  p_sud->AddVeto();
#endif
  mo->t=Max(mo->t,mo->tout);
  while (true) {
    ResetDaughters(mo);
    mo->Restore(1);
    if (!p_sud->Dice(mo)) {
      if ((m_showermode&2) && (mo->stat!=0 && mo->shower==2)) {
	msg_Debugging()<<"decay shower for "<<mo->kn_no<<"\n";
	mo->tmo=mo->t=mo->tout;
	mo->tout=sqr(mo->left->part->Momentum().Mass()+
		     mo->right->part->Momentum().Mass());
	mo->shower=3;
	mo->stat=3;
	mo->qjv=mo->decay->qjv;
	mo->qljv=mo->decay->qljv;
	mo->maxjets=mo->decay->maxjets;
	continue;
      }
      else {
	Reset(mo);
	break;
      }
    }
    else {
      mo->E2 = sqr(((1.0/z-1.0)*sprime-mo->t)/(2.0*sqrt(sprime)));
      mo->part->SetMomentum(Vec4D(sqrt(mo->E2),0.,0.,sqrt(mo->E2-mo->t)));
      msg_Debugging()<<"tlfsl test emission at t = "<<mo->t
		     <<", z = "<<mo->z<<", mode = "<<mo->shower<<"\n";
      if ((m_showermode&2) && (mo->stat!=0 && mo->shower>1)) {
	msg_Debugging()<<"resonance decay "<<mo->shower<<"\n";
	InitDaughters(tree,mo,p_sud->GetFlB(),p_sud->GetFlC(),
		      Simple_Polarisation_Info(),Simple_Polarisation_Info(),1);
	if (!p_kin->Shuffle(mo,0)) {
	  msg_Debugging()<<"shuffle failed\n";
	  continue;
	}
      } 
      if (!ini->DoKinematics()) continue;
      mo->Store();
      int stat(jetveto?p_jv->TestISKinematics(mo->prev,partner):1);
      if (stat!=1) continue;
      switch (p_sud->OrderingScheme()) {
      case 1: {
	double th(p_kin->GetOpeningAngle(mo->z,mo->E2,mo->t,
					 sqr(p_sud->GetFlB().PSMass()),
					 sqr(p_sud->GetFlC().PSMass())));
	double thmo(mo->part->Momentum().Theta());
	if (mo->part->Momentum()[3]<0.0) thmo=M_PI-thmo;
	if (th>thmo) continue;
	mo->sthcrit=mo->thcrit=thmo;
	break;
      }
      case 2: {
	double kt2(p_kin->GetRelativeKT2(mo->z,mo->E2,mo->t,
					 sqr(p_sud->GetFlB().PSMass()),
					 sqr(p_sud->GetFlC().PSMass())));
	double kt2mo(mo->part->Momentum().PPerp2());
	if (kt2>kt2mo) continue;
	mo->smaxpt2=mo->maxpt2=kt2mo;
	break;
      }
      default:
	break;
      }
      mo->stat=1;
      if (!((m_showermode&2) && (mo->stat!=0 && mo->shower>1))) {
	InitDaughters(tree,mo,p_sud->GetFlB(),p_sud->GetFlC(),
		      Simple_Polarisation_Info(),Simple_Polarisation_Info(),1);
      }
      return 1;
    }
  }
  msg_Debugging()<<"reset knot "<<mo->kn_no<<"\n";
  int stat(ini->DoKinematics());
  if (stat!=1) return stat;
  stat=jetveto?p_jv->TestISKinematics(mo->prev,partner):1;
  if (stat!=1) return stat;
  msg_Debugging()<<"}\n";
  return 1;
}

int Final_State_Shower::
TimelikeFromSpacelike(Initial_State_Shower *const ini,Tree *const tree,
		      Knot *const mo,const bool jetveto,
		      const double &sprime,const double &z,Knot *partner)
{
  msg_Debugging()<<METHOD<<"(["<<mo->kn_no<<","<<mo->part->Info()<<"],"
		 <<jetveto<<","<<sprime<<","<<z<<"): {\n";
  msg_Indent();
#ifdef USING__Veto_Info
  p_sud->SetMode(1);
#endif
  mo->sthcrit=mo->thcrit=M_PI;
  mo->smaxpt2=mo->maxpt2=1.0e10;
  if (mo->part->Info()!='H' || mo->left==NULL || mo->right==NULL) {
    int stat(FillISBranch(ini,tree,mo,jetveto,sprime,z,partner));
    if (stat!=1) return stat;
    stat=EvolveJet(tree,mo,-1);
    msg_Debugging()<<"tlfsl stat = "<<stat<<"\n";
    if (stat==1) {
      if (!p_kin->DoKinematics(mo)) return -1;
      msg_Debugging()<<"}\n";
      return 1;
    }
    Reset(mo);
    mo->sthcrit=mo->thcrit=M_PI;
    mo->smaxpt2=mo->maxpt2=1.0e10;
    return 0;
  }
  else {
    if (mo->decay!=NULL) {
      mo->Restore(1);
      mo->shower=2;
      mo->stat=3;
    }
    if (mo->shower==2 && mo->decay==NULL) {
      mo->decay = tree->NewKnot();
      mo->decay->Copy(mo);
      mo->minpt2=sqr(mo->part->Flav().Width());
      mo->Store(1);
    }
    if (mo->qjv==1.0e10) {
      mo->qjv=mo->prev->qjv;
      mo->qljv=mo->prev->qljv;
      mo->maxjets=mo->prev->maxjets;
      mo->Store(1);
    }
    if (mo->stat==0) {
      msg_Debugging()<<"init jets knot "<<mo->kn_no<<"\n";
      mo->Store();
      EstablishRelations(mo,mo->left,mo->right);
      return InitializeJets(tree,mo,1);
    }
    else {
      msg_Debugging()<<"intermediate shower knot "<<mo->kn_no<<"\n";
      while (true) {
	int stat(1);
	if ((stat=FillISBranch(ini,tree,mo,jetveto,sprime,z,partner))==-1)
	  return stat;
	msg_Debugging()<<"fillisbranch returned "<<stat<<"\n";
	if (stat==0) {
	  if (mo->stat==0) return 0;
	  continue;
	}
	if (mo->shower>=2 || mo->decay!=NULL) {
	  mo->Store();
	  EstablishRelations(mo,mo->left,mo->right);
	  return InitializeJets(tree,mo);
	}
	return EvolveJet(tree,mo);
      }
    }
    msg.Error()<<METHOD<<"(): Internal error."<<std::endl;
    return -1;
  }
  msg.Error()<<METHOD<<"(): Internal error."<<std::endl;
  return -1;
}

int Final_State_Shower::InitializeJets(Tree *tree,Knot *mo,int init)
{
  if (mo==NULL || mo->left==NULL || mo->right==NULL) {
    THROW(fatal_error,"No daughters in knot "+ToString(mo->kn_no));
  }
  msg_Debugging()<<METHOD<<"("<<mo->kn_no<<","<<init<<"): {\n";
  msg_Indent();
  Knot *d1(mo->left), *d2(mo->right);
  if (d1->shower==2 && d1->decay==NULL) {
    d1->decay = tree->NewKnot();
    d1->decay->Copy(d1);
    d1->minpt2=sqr(d1->part->Flav().Width());
  }
  if (d2->shower==2 && d2->decay==NULL) {
    d2->decay = tree->NewKnot();
    d2->decay->Copy(d2);
    d2->minpt2=sqr(d2->part->Flav().Width());
  }
  if (d1->qjv==1.0e10) {
    d1->qjv=mo->qjv;
    d1->qljv=mo->qljv;
    d1->maxjets=mo->maxjets;
  }
  if (d2->qjv==1.0e10) {
    d2->qjv=mo->qjv;
    d2->qljv=mo->qljv;
    d2->maxjets=mo->maxjets;
  }
  int first(1);
  bool dice1(d1->stat>0), dice2(d2->stat>0);
  if (dice1 || dice2) {
    SmearDaughters(mo);
    if (mo==tree->GetRoot() && dice1 && dice2) first=2; 
    while (true) {
      int accept1(1), accept2(1);
      if (FillBranch(tree,mo,first)==1) {
	if (d1->t>d2->t) {
	  if (dice1) {
	    if (d1->shower>=2 || d1->decay!=NULL) {
	      EstablishRelations(d1,d1->left,d1->right);
	      accept1=InitializeJets(tree,d1);
	    }
	    else accept1=EvolveJet(tree,d1);
	  }
	  if (dice2 && accept1==1) {
	    if (d2->shower>=2 || d2->decay!=NULL) {
	      EstablishRelations(d2,d2->left,d2->right);
	      accept2=InitializeJets(tree,d2); 
	    }
	    else accept2=EvolveJet(tree,d2);
	  }
	}
	else {
	  if (dice2) {
	    if (d2->shower>=2 || d2->decay!=NULL) {
	      EstablishRelations(d2,d2->left,d2->right);
	      accept2=InitializeJets(tree,d2); 
	    }
	    else accept2=EvolveJet(tree,d2);
	  }
	  if (dice1 && accept2==1) {
	    if (d1->shower>=2 || d1->decay!=NULL) {
	      EstablishRelations(d1,d1->left,d1->right);
	      accept1=InitializeJets(tree,d1); 
	    }
	    else accept1=EvolveJet(tree,d1);
	  }
	}
	if (accept1==1 && accept2==1) {
	  msg_Debugging()<<"accept jets\n";
	  break;
	} 
      }
      if (d1->stat==0 && d2->stat==0) {
	msg_Debugging()<<"reset knot "<<mo->kn_no<<"\n";
	break;
      }
    }
  }
  msg_Debugging()<<"out ("<<d1->kn_no<<","<<d1->stat<<") "
		 <<dice1<<", ("<<d2->kn_no<<","<<d2->stat<<") "
		 <<dice2<<std::endl;
  if (!dice1) {
    msg_Debugging()<<"init jets knot "<<d1->kn_no<<"\n";
    if (init) EstablishRelations(d1,d1->left,d1->right);
    if (!InitializeJets(tree,d1)) {
      msg_Debugging()<<"}\n";
      return 0;
    }
  }
  if (!dice2) {
    msg_Debugging()<<"init jets knot "<<d2->kn_no<<"\n";
    if (init) EstablishRelations(d2,d2->left,d2->right);
    if (!InitializeJets(tree,d2)) {
      msg_Debugging()<<"}\n";
      return 0;
    }
  }
  msg_Debugging()<<"end initialize jets, knot "<<mo->kn_no<<"\n";
  return 1;
}

Knot *Final_State_Shower::ChooseDaughter(Knot * mo)
{
  Knot *d1(mo->left), *d2(mo->right);
  if (d1->stat==3^d2->stat==3) return d1->stat==3?d1:d2;
  if ((d1->t>d1->tout && d1->t>=d2->t) || d2->t<=d2->tout) {
    if (d1->stat!=0) return d1;
    return d2;
  }
  if ((d2->t>d2->tout && d2->t>=d1->t) || d1->t<=d1->tout) {
    if (d2->stat!=0) return d2;
    return d1;
  }
  THROW(fatal_error,"Selection failed.");
  return NULL;
}

int Final_State_Shower::FillBranch(Tree *tree,Knot *mo,int first)
{
  msg_Debugging()<<METHOD<<"(["<<mo->kn_no<<","<<mo->stat<<","
		 <<mo->part->Info()<<"],"<<first<<"): p_mo = "
		 <<mo->part->Momentum()<<", t/t_out-1.0 = "
		 <<mo->t/mo->tout-1.0<<" {\n";
  msg_Indent();
  if (first<1 && mo->t<=mo->tout) {
    msg_Debugging()<<"!first & "<<mo->t<<"<"<<mo->tout<<std::endl; 
    return 0; 
  }
  Knot *d1(mo->left), *d2(mo->right);
  msg_Debugging()<<"p_d1 = "<<d1->part->Momentum()
		 <<", t_d1 = "<<d1->t<<", ("
		 <<d1->kn_no<<","<<d1->stat<<")"<<std::endl;
  msg_Debugging()<<"p_d2 = "<<d2->part->Momentum()
		 <<", t_d2 = "<<d2->t<<", ("
		 <<d2->kn_no<<","<<d2->stat<<")"<<std::endl;
#ifdef USING__Veto_Info
  p_sud->AddVeto();
#endif
  while (true) {
    Knot *d(ChooseDaughter(mo));
    msg_Debugging()<<"selected daughter "<<d->kn_no<<", "
		   <<d->part->Flav()<<", t="<<d->t<<", stats "
		   <<d1->stat<<" "<<d2->stat<<std::endl;
    ResetDaughters(d);
    if (p_sud->Dice(d,first==2?mo:NULL)) { 
      msg_Debugging()<<"test emission at ("<<d->t
		     <<","<<d->z<<") knot "<<d->kn_no<<"\n";
      InitDaughters(tree,d,p_sud->GetFlB(),p_sud->GetFlC(),
		    p_sud->GetPolB(),p_sud->GetPolC(),true);
      d->stat=1;
    }   
    else if ((m_showermode&2) && (d->stat!=0 && d->shower==2)) {
      msg_Debugging()<<"decay shower for "<<d->kn_no<<"\n";
      d->tmo=d->t=d->tout;
      d->tout=sqr(d->left->part->Momentum().Mass()+
		  d->right->part->Momentum().Mass());
      d->shower=3;
      d->stat=3;
      d->qjv=d->decay->qjv;
      d->qljv=d->decay->qljv;
      d->maxjets=d->decay->maxjets;
      continue;
    }
    else {
      Reset(d);
    }
    if (d1->stat!=3 && d2->stat!=3) {
      mo->Restore();
      if (p_kin->Shuffle(mo,Max(0,first))) {
   	if (p_jv->TestFSKinematics(mo)==1) {
	  p_sud->AcceptBranch(mo);
	  msg_Debugging()<<"kinematics check passed"<<std::endl;
	  mo->stat=0;
	  mo->Store();
	  msg_Debugging()<<"}\n";
	  return 1;
	}
	msg_Debugging()<<"kinematics vetoed\n";
#ifdef USING__Veto_Info
	p_sud->SetVeto(svc::jet_veto);
#endif
      }
      else msg_Debugging()<<"shuffle failed\n";
    }
    if (d1->stat==0 && d2->stat==0) break;
  }
  msg_Debugging()<<"found no possible branch\n";
  if (d1->left && d1->left->part->Info()!='H') {
    ResetDaughters(d1);
    d1->stat=3;
  }
  if (d2->left && d2->left->part->Info()!='H') {
    ResetDaughters(d2);
    d2->stat=3;
  }
  msg_Debugging()<<"}\n";
  if (first>0) {
    mo->Restore();
    if ((mo->left->left==NULL || mo->right->left==NULL) && 
	!p_kin->Shuffle(mo,0)) {
      msg_Debugging()<<"shuffle failed\n";
      msg_Debugging()<<"}\n";
      return -1;
    }
  }
  else {
    if (!p_kin->Shuffle(mo,0)) {
      msg_Debugging()<<"shuffle failed\n";
      msg_Debugging()<<"}\n";
      return 0;
    }
//     if (p_jv->TestFSKinematics(mo)!=1) {
//       msg_Debugging()<<"kinematics vetoed\n";
//       msg_Debugging()<<"}\n";
//       return 0;
//     }
  }
  if (d1->stat==0 && d2->stat==0) return 1;
  return 0;
}

int Final_State_Shower::EvolveJet(Tree *tree,Knot *mo,int first)
{
  msg_Debugging()<<METHOD<<"(["<<mo->kn_no<<","<<mo->part->Flav()<<"]): "
		 <<"p_mo = "<<mo->part->Momentum()<<" "
		 <<mo->part->Momentum().Abs2()<<", t_mo = "<<mo->t<<" {\n";
  msg_Indent();
  if (mo->stat==0 && mo->decay==NULL) {
    msg_Debugging()<<"stat = 0, reset daughters\n";
    ResetDaughters(mo);
    mo->left=mo->right=NULL;
    msg_Debugging()<<"}\n";
    return true;
  }
  if (mo->decay!=NULL) mo->zs=mo->z;
  Knot *d1(mo->left), *d2(mo->right);
  int evolve1(1), evolve2(1), stat;
  while (true) {
    stat=FillBranch(tree,mo,first);
    msg_Debugging()<<"status' "<<d1->stat<<" "<<d2->stat<<"\n";
    if (d1->stat==0 && d2->stat==0) break;
    if (d1->t>d2->t) { 
      evolve1=EvolveJet(tree,d1);
      if (evolve1==1) evolve2=EvolveJet(tree,d2);
    }
    else { 
      evolve2=EvolveJet(tree,d2);
      if (evolve2==1) evolve1=EvolveJet(tree,d1);
    }
    if (evolve1==1 && evolve2==1) {
      msg_Debugging()<<"}\n";
      return 1; 
    }
  }
  if (stat==0) {
    Reset(mo);
    msg_Debugging()<<"reset knot "<<mo->kn_no<<", t = "<<mo->t<<"\n";
    msg_Debugging()<<"}\n";
    return 0;
  }
  msg_Debugging()<<"}\n";
  return 1;
}

bool Final_State_Shower::SetAllColours(Knot *mo) 
{ 
  return SetColours(mo,p_kin); 
}

bool Final_State_Shower::
SetColors(Knot *mo,unsigned int oldr,unsigned int newr,
	  unsigned int olda,unsigned int newa)
{
  msg_Debugging()<<METHOD<<"("<<oldr<<"->"<<newr<<","<<olda<<"->"<<newa
		 <<"): sync colors in "<<mo->kn_no
		 <<"("<<mo->part->GetFlow(1)<<","
		 <<mo->part->GetFlow(2)<<")->"
		 <<(mo->left?mo->left->kn_no:0)<<","
		 <<(mo->right?mo->right->kn_no:0)<<"\n";
  msg_Indent();
  for (short unsigned int i(0);i<2;++i) {
    Knot *k(i==0?mo->left:mo->right);
    if (k!=NULL) {
      unsigned int r(k->part->GetFlow(1)), a(k->part->GetFlow(2));
      if (r==oldr || a==olda) {
	if (!SetColors(k,oldr,newr,olda,newa)) {
	  msg.Error()<<METHOD<<"(): Colour adjustment failed."<<std::endl;
	  return false;
	}
 	if (r==oldr) k->part->SetFlow(1,newr);
 	if (a==olda) k->part->SetFlow(2,newa);
      }
    }  
  }
  return true;
}

bool Final_State_Shower::SetColours(Knot *mo,Timelike_Kinematics *kin)
{
  if (mo==NULL) {
    msg.Error()<<METHOD<<"(..): Error. Void mother knot."<<std::endl;
    return false;
  }
  if (mo->left==NULL) {
    if (mo->right==NULL) return true;
    mo->right->part->SetFlow(1,mo->part->GetFlow(1));
    mo->right->part->SetFlow(2,mo->part->GetFlow(2));
    return true;
  }
  if (mo->right==NULL) {
    mo->left->part->SetFlow(1,mo->part->GetFlow(1));
    mo->left->part->SetFlow(2,mo->part->GetFlow(2));
    return true;
  }
  // check if already enough colours
  Knot *test(NULL);
  int all_colors_known(1);
  std::set<int> cnt;
  for (int i=0;i<3;++i) {
    if (i==0) test = mo;
    if (i==1) test = mo->left;
    if (i==2) test = mo->right;
    if (test->part->Flav().Strong()) {
      int nc(0);
      for (short unsigned int i(1);i<=2;++i) {
	unsigned int c(test->part->GetFlow(i));
	if (c) {
	  ++nc;
	  std::set<int>::iterator cit(cnt.find(c));
	  if (cit!=cnt.end()) cnt.erase(cit);
	  else cnt.insert(c);
	}
      }
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
  if (cnt.size()!=0 && mo->decay==NULL) all_colors_known=0;
  msg_Debugging()<<METHOD<<"(): "<<mo->kn_no<<","
		 <<mo->part->Flav()<<"("<<mo->part->GetFlow(1)
		 <<","<<mo->part->GetFlow(2)<<") -> "<<mo->left->kn_no<<","
		 <<mo->left->part->Flav()<<"("<<mo->left->part->GetFlow(1)
		 <<","<<mo->left->part->GetFlow(2)<<");"<<mo->right->kn_no<<","
		 <<mo->right->part->Flav()<<"("<<mo->right->part->GetFlow(1)
		 <<","<<mo->right->part->GetFlow(2)<<") => "
		 <<all_colors_known<<"\n";
  msg_Indent();
  Knot *d1(mo->left), *d2(mo->right);
  if (all_colors_known) {
    if (mo->decay && mo->left->decay!=mo->decay) {
      // sync colors
      SetColors(mo->left,mo->decay->part->GetFlow(1),mo->part->GetFlow(1),
  		mo->decay->part->GetFlow(2),mo->part->GetFlow(2));
      SetColors(mo->right,mo->decay->part->GetFlow(1),mo->part->GetFlow(1),
  		mo->decay->part->GetFlow(2),mo->part->GetFlow(2));
    }
  }
  else if (mo->part->Flav().Kfcode()!=kf::none) {
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
  if (bl) {
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
    else for (int i=0;i<bl_meps->NInP();++i) bl_meps->InParticle(i)->SetStatus(part_status::decayed);
  }
  // deactivate in partons!
  if (!kn) return;
  Particle * p(NULL);
  if (kn->part->Info()=='H') {
    /* 
       New jet : kn = hard parton from ME info = 'HF'
                 and kn outgoing
                 or kn->left or kn->right not from ME
    */
    if (!kn->left) {
      if (bl||bl_meps) {
	p = new Particle(*kn->part);
	p->SetFinalMass(p->Momentum().Mass());
      }
      if (bl) {
	jet = new Blob();
#ifdef USING__Veto_Info
	jet->AddData("FS_VS",new Blob_Data<std::vector<int> >(p_sud->Vetos(0)));
	jet->AddData("IFS_VS",new Blob_Data<std::vector<int> >(p_sud->Vetos(1)));
#endif
	jet->AddToInParticles(p);
      }
      if (bl_meps) {
	bl_meps->AddToOutParticles(p);
	bl_meps->SetStatus(blob_status::inactive);
      }

      p = new Particle(*kn->part);
      p->SetFinalMass(p->Momentum().Mass());
      if (pl) p->SetNumber(pl->size());
         else p->SetNumber(0);
      kn->part->SetNumber(p->Number());
      if (bl) {
	jet->AddToOutParticles(p);
	jet->SetId();
	jet->SetType(btp::FS_Shower);
	jet->SetTypeSpec("APACIC++2.0");
	jet->SetStatus(blob_status::needs_harddecays |
		       blob_status::needs_hadronization);
	bl->push_back(jet);
      }
      return;
    }
    else {
      if ((kn->left->part->Info() != 'H') || (kn->right->part->Info() != 'H')) {
	p = new Particle(*kn->part);
	p->SetFinalMass(p->Momentum().Mass());
      	p->SetStatus(part_status::decayed);
	if (pl) {
	  pl->push_back(p);
	  p->SetNumber(pl->size());
	}
	else p->SetNumber(0);
        kn->part->SetNumber(p->Number());
	if (bl_meps) {
	  bl_meps->AddToOutParticles(p);
	  bl_meps->SetStatus(blob_status::inactive);
	}
	if (bl) {
	  jet = new Blob();
#ifdef USING__Veto_Info
	  jet->AddData("FS_VS",new Blob_Data<std::vector<int> >(p_sud->Vetos(0)));
	  jet->AddData("IFS_VS",new Blob_Data<std::vector<int> >(p_sud->Vetos(1)));
#endif
	  jet->AddToInParticles(p);
	  jet->SetId();
	  jet->SetType(btp::FS_Shower);
	  jet->SetTypeSpec("APACIC++2.0");
	  jet->SetStatus(blob_status::needs_harddecays |
			 blob_status::needs_hadronization);
	  bl->push_back(jet);
	}
      }
    }
  }
  else {
    if (!kn->left) {
      if (bl && !jet) {
	msg.Error()<<"ERROR in Final_State_Shower ::ExtractPartons :\n"
		   <<"    No jet for Parton : "<<kn->part->Number()<<std::endl;
	abort();
      }
      if (pl) kn->part->SetNumber(pl->size());
	 else kn->part->SetNumber(0);
      kn->part->SetStatus(part_status::active);
      if (pl) pl->push_back(kn->part);
      p = new Particle(*kn->part);
      p->SetFinalMass(p->Momentum().Mass());
      if (bl) jet->AddToOutParticles(p);
    }
  }
  ExtractPartons(kn->left,jet,bl,pl); 
  ExtractPartons(kn->right,jet,bl,pl);
}

void Final_State_Shower::ExtractFinalState(Blob *blob,Knot *kn)
{
  if (!kn || !blob) return;
  Particle * p(NULL);
  if (!kn->left) {
    p = new Particle(*kn->part);
    p->SetNumber(0);
    blob->AddToOutParticles(p);
    return;
  }
  ExtractFinalState(blob,kn->left);
  ExtractFinalState(blob,kn->right);
}

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

  if (d1->decay!=mo->decay) {
    d1->part->SetFlow(1,0);
    d1->part->SetFlow(2,0);
    d2->part->SetFlow(1,0);
    d2->part->SetFlow(2,0);
  }  

  // set color connections (if not yet known)
  APACIC::Final_State_Shower::SetColours(mo,0);

  double t_mo(mo->part->Momentum().Abs2()), st_mo(mo->t);
  double tb(d1->part->Momentum().Abs2()), tc(d2->part->Momentum().Abs2());
  double E_mo(mo->part->Momentum()[0]), z_mo(d1->part->Momentum()[0]/E_mo); 
  double th(p_kin->GetOpeningAngle(z_mo,sqr(E_mo),t_mo,tb,tc));
  double maxpt2(p_kin->GetRelativeKT2(z_mo,sqr(E_mo),t_mo,tb,tc));
  if (mo->part->Flav().IsQuark() && 
      d1->part->Flav().Strong() && d2->part->Flav().Strong()) {
    if (d1->part->Flav().IsQuark()) {
      d1->t      = mo->t;
      d1->thcrit = mo->thcrit;
      d2->t      = st_mo;
      d2->thcrit = th;
    }
    else {
      d1->t      = st_mo;
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
	d2->t      = st_mo;
	d2->thcrit = th;
      }
      else {
	d1->t      = st_mo;
	d1->thcrit = th;
	d2->t      = mo->t;
	d2->thcrit = mo->thcrit;
      }
      d1->maxpt2 = d2->maxpt2 = maxpt2;
    }
    else if (!(d1->part->Flav().Strong()) && (d2->part->Flav().Strong())) {
      d1->t      = st_mo;
      d1->thcrit = M_PI;
      d2->t      = mo->t;
      d2->thcrit = mo->thcrit;
      d1->maxpt2 = d2->maxpt2 = mo->maxpt2;
    }
    else if ((d1->part->Flav().Strong()) && !(d2->part->Flav().Strong())) {
      d1->t      = mo->t;
      d1->thcrit = mo->thcrit;
      d2->t      = st_mo;
      d2->thcrit = M_PI;
      d1->maxpt2 = d2->maxpt2 = mo->maxpt2;
    }
  }
  else {
    if ((d1->part->Flav().Strong()) && (d2->part->Flav().Strong())) {
      d1->t      = st_mo;
      d1->thcrit = th;
      d2->t      = st_mo;
      d2->thcrit = th;
      d1->maxpt2 = d2->maxpt2 = st_mo;
    }
    else {
      d1->t      = st_mo;
      d1->thcrit = M_PI;
      d2->t      = st_mo;
      d2->thcrit = M_PI;
      d1->maxpt2 = d2->maxpt2 = st_mo;
    }
  }
  mo->t = t_mo;
  if (d1->part->Info()=='H' && d1->left && 
      d1->left->part->Info()=='H') 
    d1->t=d1->part->Momentum().Abs2();
  if (d2->part->Info()=='H' && d2->left && 
      d2->left->part->Info()=='H') 
    d2->t=d2->part->Momentum().Abs2();
  if (d1->decay!=mo->decay) {
    msg_Debugging()<<"restore saved info in "
		   <<d1->kn_no<<" & "<<d2->kn_no<<"\n";
    d1->thcrit=d1->sthcrit;
    d1->smaxpt2=d1->smaxpt2;
    d2->thcrit=d2->sthcrit;
    d2->smaxpt2=d2->smaxpt2;
  }
}

void Final_State_Shower::
InitDaughters(Tree * tree,Knot * mo,ATOOLS::Flavour flb,ATOOLS::Flavour flc,
	      ATOOLS::Simple_Polarisation_Info polb,
	      ATOOLS::Simple_Polarisation_Info polc,bool diced)
{ 
  if (!mo->left) {
    mo->left          = tree->NewKnot();
    mo->right         = tree->NewKnot();
  }
  else if (mo->shower>=2) {
    Knot *newd1(tree->NewKnot());
    mo->left->prev=newd1;
    mo->right->prev=newd1;
    newd1->left=mo->left;
    newd1->right=mo->right;
//     msg_Debugging()<<"copied mom set to "<<mo->part->Momentum()<<"\n";
//     newd1->part->SetMomentum(mo->part->Momentum());
    mo->left=newd1;
    mo->right = tree->NewKnot();
  }
  if (diced) {
    mo->left->prev     = mo;
    mo->left->polinfo  = polb;
    mo->left->tout     = sqr(flb.PSMass());
    mo->left->stat     = 3;  
    mo->left->part->SetFlav(flb);
    mo->left->part->SetInfo('F');
    mo->left->part->SetStatus(part_status::active);
    mo->left->didkin   = false;
    if (mo->shower>=2) mo->left->tout=mo->tout;
    mo->left->shower   = mo->shower;  
    mo->left->decay    = mo->decay;
    if (mo->shower==3) mo->left->tmo=mo->t;
    mo->left->qjv      = mo->qjv;
    mo->left->qljv     = mo->qljv;
    mo->left->maxjets  = mo->maxjets;

    mo->right->prev    = mo;
    mo->right->polinfo = polc;
    mo->right->tout    = sqr(flc.PSMass());
    mo->right->stat    = 3;  
    mo->right->part->SetFlav(flc);
    mo->right->part->SetInfo('F');
    mo->right->part->SetStatus(part_status::active);
    mo->right->didkin  = false;
    mo->right->qjv     = mo->qjv;
    mo->right->qljv    = mo->qljv;
    mo->right->maxjets = mo->maxjets;

    mo->left->minpt2   = mo->minpt2;
    mo->right->minpt2  = mo->minpt2;

    if (mo->part->Info()!='H') mo->part->SetInfo('f');
    mo->part->SetStatus(part_status::decayed);
  }

  // Reset kinematics
  double th, maxpt2;
  if ((mo->left->part->Flav().Strong()) && 
      (mo->right->part->Flav().Strong())) {
    th     = p_kin->GetOpeningAngle(mo);
    maxpt2 = p_kin->GetRelativeKT2(mo->z,mo->E2,mo->t,
				   sqr(mo->left->part->Flav().PSMass()),
				   sqr(mo->right->part->Flav().PSMass()));
  }
  else {
    if (mo->prev) {
      th     = mo->prev->thcrit;
      maxpt2 = mo->prev->maxpt2;
    }
    else {
      th     = M_PI;
      maxpt2 = mo->maxpt2;
    }
  }

  mo->left->t       = mo->t; 
  mo->left->pt2lcm  = mo->pt2lcm; 
  mo->left->E2      = (mo->z)*(mo->z)*(mo->E2); 
  mo->left->thcrit  = th;
  mo->left->maxpt2  = maxpt2;

  mo->right->t      = mo->t; 
  mo->right->pt2lcm = mo->pt2lcm; 
  mo->right->E2     = (1. - mo->z)*(1. - mo->z)*(mo->E2);
  mo->right->thcrit = th;
  mo->right->maxpt2 = maxpt2; 
}

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
  mo->part->SetStatus(part_status::decayed);
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
  mo->left->part->SetStatus(part_status::active);
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
  mo->right->part->SetStatus(part_status::active);
  mo->right->part->SetInfo('H');
  mo->right->part->SetFlow(2,mo->left->part->GetFlow(1));
  mo->right->t         = mo->t;
  mo->right->tout      = sqr(mo_flavs[1].PSMass());
  mo->right->maxpt2    = 0.;
  mo->right->E2        = E*E;
  mo->right->thcrit    = start_th;
}

bool Final_State_Shower::ResetDaughters(Knot * mo)
{
  msg_Debugging()<<METHOD<<"("<<mo->kn_no<<"): shower = "<<mo->shower<<"\n";
  if (mo->left!=NULL) {
    if (mo->decay!=NULL) {
      mo->left=mo->decay->left;
      mo->right=mo->decay->right;
      mo->decay->left->prev=mo;
      mo->decay->right->prev=mo;
      return true;
    }
    else {
      Reset(mo->left);
      Reset(mo->right);
    }  
  }
  if (mo->part->Info()!='H') mo->part->SetInfo('F');
  mo->part->SetStatus(part_status::active);
  return false;
} 

void Final_State_Shower::Reset(Knot * mo) 
{ 
  msg_Debugging()<<METHOD<<"("<<mo->kn_no<<"):\n";
  if (mo==NULL) return;
  if (mo->shower==3) {
    mo->shower=0;
    mo->tout=mo->tmo;
  }
  mo->t=mo->tout;
  mo->stat=0;
  mo->part->SetStatus(part_status::active);
  if (mo->left!=NULL && mo->left->part->Info()!='H') mo->left=NULL;
  if (mo->right!=NULL && mo->right->part->Info()!='H') mo->right=NULL;  
  if (mo->part->Info()!='H') {
    mo->part->SetInfo('F');
    mo->part->SetMomentum(Vec4D());
  }
}

void Final_State_Shower::OutputTree(Tree * tree) 
{
  if (tree->GetRoot()==0) {
    msg.Out()<<"Empty Tree\n";
  }
  else {
    int number(0);
    msg.Out()<<"Final Tree:\n"<<*tree<<std::endl
	     <<"Total 4 Mom = "<<GetMomentum(tree->GetRoot(),number);
    msg.Out()<<"   for "<<number<<" FS particles.\n";
  }
}


Vec4D  Final_State_Shower::GetMomentum(Knot * mo, int & number) 
{
  if (mo->left) {
    Vec4D p(GetMomentum(mo->left,number)+GetMomentum(mo->right,number));
    Vec4D ptest(mo->left->part->Momentum()+mo->right->part->Momentum());
    if (!(ptest==mo->part->Momentum())) {
      number-=10000;
      msg.Error()<<METHOD<<"(..):  Four momentum not conserved "
		 <<"in knot "<<mo->kn_no<<"\n"
		 <<"  p_miss = "<<(ptest-mo->part->Momentum())<<"\n"
		 <<"  p_old  = "<<mo->part->Momentum()<<"  "
		 <<mo->part->Momentum().Abs2()<<"\n"
		 <<"  p_new  = "<<ptest<<"  "<<p.Abs2()<<std::endl;
      msg_Debugging()<<GetMomentum(mo->left,number)<<" + "
		     <<GetMomentum(mo->right,number)<<" vs. "<<std::endl
		     <<mo->left->part->Momentum()<<" + "
		     <<mo->right->part->Momentum()<<std::endl;
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

