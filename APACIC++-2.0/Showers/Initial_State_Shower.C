#include "Initial_State_Shower.H"

#include "PDF_Handler.H"
#include "Data_Read.H"
#include "MyStrStream.H"
#include "Veto_Info.H"
#include <iomanip>

#ifdef PROFILE__all
#include "prof.hh"
#else 
#define PROFILE_HERE 
#define PROFILE_LOCAL(LOCALNAME) 
#endif

using namespace APACIC;
using namespace ATOOLS;

Initial_State_Shower::Initial_State_Shower(PDF::ISR_Handler *const isr, 
					   ATOOLS::Jet_Finder *const jf,
					   Final_State_Shower *const fin,
					   MODEL::Model_Base *const model,
					   Data_Read *const dataread) : 
  p_fin(fin),
  p_tools(new Sudakov_Tools(model)),
  p_kin(new Spacelike_Kinematics(jf)),
  p_suds(new Spacelike_Sudakov*[2]),
  m_allowed(200)
{
  double cplscalefac(0.25*ToType<double>
		     (rpa.gen.Variable("IS_CPL_SCALE_FACTOR","1.0")));
  m_t0=dabs(dataread->GetValue<double>("IS_PT2MIN",4.0));
  double shadron(dataread->GetValue<double>("IS_MAX_SCALE",
					    sqr(rpa.gen.Ecms())));
  double emin(dataread->GetValue<double>("IS_MINIMAL_E",0.5));
  int cplscheme(dataread->GetValue<int>("IS_COUPLING_SCHEME",1));
  int pdfscheme(dataread->GetValue<int>("IS_PDF_SCALE_SCHEME",1));
  int orderingscheme(dataread->GetValue<int>("IS_ORDERING_SCHEME",0));
  for (short unsigned int i(0);i<2;++i) {
    if (isr->PDF(i)->Q2Min()>m_t0*cplscalefac) {
      msg.Error()<<METHOD<<"(..):\n   IS_PT2MIN("<<m_t0
		 <<")*IS_CPL_SCALE_FACTOR("<<cplscalefac<<") "
		 <<"smaller than minimum scale given by PDF ("
		 <<isr->PDF(i)->Q2Min()<<").\n"
		 <<"   Please change your settings in "
		 <<dataread->FileName()<<"\n";
      THROW(fatal_error,"Minimal PDF scale too low.");
    }
    p_suds[i] = new Spacelike_Sudakov
      (isr->PDF(i),p_tools,p_kin,m_t0,dataread,i);
    p_suds[i]->SetRemnant(isr->GetRemnant(i));
    p_suds[i]->SetMaxScale(shadron);
    p_suds[i]->SetMinEnergy(emin);
    p_suds[i]->SetPDFScaleFactor(cplscalefac);
    p_suds[i]->SetCouplingScaleFactor(cplscalefac);
    p_suds[i]->SetCouplingScheme(cplscheme);
    p_suds[i]->SetPDFScheme(pdfscheme);
    p_suds[i]->SetOrderingScheme(orderingscheme);
    m_extra_pdf[i]=m_to_be_diced[i]=1;
    if (!p_suds[i]->Initialize()) 
      THROW(fatal_error,"failed to intialize is sudakov");
  }
  m_t0 = - dabs(m_t0);
}

Initial_State_Shower::~Initial_State_Shower()
{
  delete p_tools;
  delete p_kin;
  for (short unsigned int i(0);i<2;i++) delete p_suds[i];
  delete p_suds;
}

int Initial_State_Shower::PerformShower(Tree **const trees,Tree *const fintree,
					const int jetvetoflag) 
{
  PROFILE_HERE;
  p_fstree=fintree;
  p_istrees=trees;
  m_jetveto=jetvetoflag>0;
#ifdef USING__Veto_Info
  p_suds[0]->ClearVetos();
  p_suds[1]->ClearVetos();
#endif
  if (jetvetoflag<0 || jetvetoflag>1) m_extra_pdf[0]=m_extra_pdf[1]=0;
  else m_extra_pdf[0]=m_extra_pdf[1]=1;
  if (InitializeSystem(trees,trees[0]->GetRoot(),trees[1]->GetRoot())) {
    p_kin->BoostInCMS(trees,trees[0]->GetInitiator(),trees[1]->GetInitiator());
    m_lab=p_kin->BoostInLab(trees);
    return 1;
  }
  msg.Error()<<METHOD<<"(..): Shower failure."<<std::endl;
  return 0;
}

bool Initial_State_Shower::InitializeSystem(Tree ** trees,Knot * k1,Knot * k2)
{
  msg_Debugging()<<METHOD<<"("<<k1->kn_no<<","<<k2->kn_no<<"): {\n";
  msg_Indent();
  m_to_be_diced[1]=m_to_be_diced[0]=1;
  if (k1==NULL || k2==NULL) THROW(fatal_error,"No trees.");
  bool decay1(k1->stat>=1), decay2(k2->stat>=1);
  int first(0);
  if (decay1^decay2) first=1;
  if (!decay1 && !decay2) first=2;
  int mismatch(0), caught_jetveto(0);
  bool accepted(true); 
  while (true) {
    accepted=true;
    m_sprime=(k1->part->Momentum()+k2->part->Momentum()).Abs2();
    msg_Debugging()<<"decay1 = "<<decay1<<", decay2 = "<<decay2
		   <<", jetveto = "<<caught_jetveto<<"\n";  
    // Parton 1/Tree 1 is the one to decay.
    if (decay1 && caught_jetveto!=3) {
      msg_Debugging()<<"evolve 1\n";
      m_to_be_diced[0]=0;
      if (FillBranch(trees,k1,k2,0)) {
	if (k1->z>0.0) m_sprime=m_sprime/k1->z;
      }
      else {
	accepted=false;
      }
    }
    // Parton 2/Tree 2 is the one to decay.    
    if (decay2 && caught_jetveto!=2) {
      msg_Debugging()<<"evolve 2\n";
      m_to_be_diced[1]=0;
      if (FillBranch(trees,k2,k1,1)) {
	if (k2->z>0.0) m_sprime=m_sprime/k2->z;
      }
      else {
	accepted=false;
      }
    }
    if (accepted) {
      p_kin->InitKinematics(trees,k1,k2,first);
      if (!decay1) SetColours(k1);
      if (!decay2) SetColours(k2);
      p_suds[0]->AcceptBranch(k1);
      p_suds[1]->AcceptBranch(k2);
      int stat(EvolveSystem(trees,k1,k2));
      if (stat==1) {
	msg_Debugging()<<"}\n";
	return true;
      }
      else if (stat==2 || stat==3) {
	msg_Debugging()<<"init system veto "<<stat<<"\n";
	if (stat==2) m_sprime=m_sprime*k1->z;
	else m_sprime=m_sprime*k2->z;
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
      if (mismatch>m_allowed) {
	msg.Error()<<METHOD<<"(..): Shower failure, "
		   <<mismatch<<" trials."<<std::endl;
	msg_Debugging()<<"}\n";
	return false;
      }
      trees[0]->Restore();
      trees[1]->Restore();
    }
  }
  msg_Debugging()<<"iss reset";
  msg_Debugging()<<"}\n";
  return false;
}

int Initial_State_Shower::DoKinematics()
{
  if (!p_kin->DoKinematics(p_istrees,p_k1,p_k2,m_ntree0,false)) return -1;
  return 1-p_suds[m_ntree0]->PTVeto(p_k1->prev);
}

int Initial_State_Shower::EvolveSystem(Tree **const trees,Knot *k1,Knot *k2)
{
  msg_Debugging()<<METHOD<<"("<<k1->kn_no<<","<<k2->kn_no<<"): {\n";
  msg_Indent();
  if (k1->t==k1->tout && k2->t==k2->tout) return 1;  
  int ntree0(0), ntree1(1);
  ChooseMother(ntree0,ntree1,k1,k2);
  int caught_jetveto(0);
  while (true) {
    if (k1->stat>0 || caught_jetveto) {
      k1->stat           = 0;
      k1->E2             = sqr(k1->part->Momentum()[0]);
      k1->prev->E2       = k1->E2/sqr(k1->z);                      
      k1->prev->left->E2 = k1->prev->E2*sqr(1.-k1->z);
    }  
    else {
      double sprime_a((k1->part->Momentum()+k2->part->Momentum()).Abs2());
      double sprime_b((k1->prev->part->Momentum()+k2->part->Momentum()).Abs2());
      k1->z=sprime_a/sprime_b;
      k1->prev->x=k1->x/k1->z;

      if (k1->prev->x>1.) return 0;

      m_sprime/=k1->z;
    }
    Knot *mo(k1->prev);
    if (mo->stat>0) {
      m_to_be_diced[ntree0]=0;
      // do not dice t for me knot
      if (mo->part->Info()!='H' ||
	  mo->prev==NULL || mo->prev->part->Info()!='H') 
	if (!FillBranch(trees,k1->prev,k2,ntree0)) return 0; 
    }

    double maxt(p_kin->CalculateMaxT(k1,k2));
    if (maxt<mo->left->tout) return 0;
    else {
      if (caught_jetveto==0 && mo->left->stat>0) {
	if (mo->left->part->Info()!='H') k1->prev->left->t=dabs(k1->t);
	mo->left->tmax=maxt;
      }
    } 
  
    if (caught_jetveto==0) {
      if (p_fin) {
	double sprime_a = (k1->part->Momentum()+k2->part->Momentum()).Abs2();

	if (!p_kin->DoKinematics(trees,k1,k2,ntree0,true)) return 0;
	p_k1=k1; 
	p_k2=k2; 
	m_ntree0=ntree0;
	switch (p_fin->TimelikeFromSpacelike(this,trees[ntree0],k1->prev->left,
					     m_jetveto,sprime_a,k1->z,k2)) {
	case -1:
	  msg.Error()<<METHOD<<"(..): FS Shower failure."<<std::endl;
	  return 0;
	case 0: 
	  caught_jetveto=-10;
#ifdef USING__Veto_Info
	  p_suds[ntree0]->SetVeto(svc::jet_veto);
#endif
	  break;
	case 1:
	  break;
	}
      }
      else {
	Knot * mo = k1->prev->left;
	mo->t     = mo->tout;
	mo->stat  = 0;
	mo->part->SetStatus(part_status::active);
	if (!p_kin->DoKinematics(trees,k1,k2,ntree0,false)) return 0;
      }
    }
    else {
      if (!p_kin->DoKinematics(trees,k1,k2,ntree0,false)) return 0;
    }

    if (caught_jetveto==-10) {
      m_to_be_diced[ntree0]=1;
      p_kin->ResetMomenta(k1,trees[ntree0]);
      msg_Debugging()<<"caught jetveto\n";
      return (ntree0+2);
    }
    
    if (p_fin) p_fin->SetAllColours(k1->prev->left);
    p_suds[ntree0]->AcceptBranch(k1->prev);
  
    if (k1->prev->z>0.) m_sprime = m_sprime/k1->prev->z;
    
    if (ntree0==0) {
      int stat(EvolveSystem(trees,k1->prev,k2));
      if (stat==0 || stat ==1) return stat;
      caught_jetveto=stat;
    }
    else {  
      int stat(EvolveSystem(trees,k2,k1->prev));
      if (stat==0 || stat ==1) return stat;
      caught_jetveto=stat;
    }
    
    if (caught_jetveto && k1->prev->z>0.) m_sprime = m_sprime*k1->prev->z;

    if (caught_jetveto!=ntree0+2) {
      m_to_be_diced[ntree0]=1;
      k1->prev->t=k1->prev->tmax;
      p_kin->ResetMomenta(k1,trees[ntree0]);
      msg_Debugging()<<"jet veto not on tree "<<ntree0<<"\n";
      msg_Debugging()<<"}\n";
      return caught_jetveto;
    }
  }
  return 0;
}  

int Initial_State_Shower::FillBranch(Tree ** trees,Knot * active,
				     Knot * partner,int leg) {
  msg_Debugging()<<METHOD<<"("<<active->kn_no<<","<<partner->kn_no<<","
 		 <<leg<<"): {\n";
  msg_Indent();
  Flavour flavs[2];
#ifdef USING__Veto_Info
  p_suds[leg]->AddVeto();
#endif
  if (active->right) {
    active->qjv=active->right->qjv;
    active->qljv=active->right->qljv;
    active->left->qjv=active->qjv;
    active->left->qljv=active->qljv;
  }
  if (p_suds[leg]->Dice(active,m_sprime,m_extra_pdf[leg])) {
    flavs[0] = p_suds[leg]->GetFlA();
    flavs[1] = p_suds[leg]->GetFlC();    
    FillMotherAndSister(trees[leg],active,flavs,leg);
    msg_Debugging()<<"test emission at t = "<<active->t
 		   <<", z = "<<active->z<<"\n";
    msg_Debugging()<<"}\n";
    return 1;
  }
  active->prev   = 0;
  active->stat   = 0;
  active->t      = active->tout;
  active->part->SetStatus(part_status::active);
  msg_Debugging()<<"no branch";
  msg_Debugging()<<"}\n";
  return 1;
}

void Initial_State_Shower::ChooseMother(int &ntree0,int &ntree1, 
					Knot *&k1, Knot *&k2)
{
  bool swap((k1->t > k2->t) && (k2->t != k2->tout));
  bool known1(k1->stat==0 && k1->t != k1->tout), known2(k2->stat==0 && k2->t != k2->tout);

  if (known1 || known2) {
    bool save_swap(swap);
    if (known1 && !known2) swap=false;
    else if (!known1 && known2) swap=true;

    if (known1 && known2) {
      if (k1->prev->kn_no<k2->prev->kn_no) swap=false;
      else swap=true;
    }

    if (swap!=save_swap) {
      msg_Tracking()<<"Initial_State_Shower::ChooseMother changed swap according to known history \n";
    }
  }

  if (swap) {
    Knot * kh(k1);
    k1 =k2;
    k2 =kh;
    ntree0=1;
    ntree1=0;
  }
}

void Initial_State_Shower::FillMotherAndSister(Tree * tree,Knot * k,
					       Flavour * k_flavs,int leg)
{
  double pt2max(sqr(rpa.gen.Ecms())), th(4.*k->z*k->z*k->t/(4.*k->z*k->z*k->t-(1.-k->z)*k->x*k->x*pt2max));
  if (!k->part->Flav().Strong() || 
      !k_flavs[0].Strong() || !k_flavs[1].Strong() ) th=k->thcrit;
  Knot * mother(NULL);
  if (k->prev) {
    mother = k->prev;
    // delete old color informations if flav content changed
    if (mother->part->Flav()!=k_flavs[0]) {
      mother->part->SetFlow(1,0);
      mother->part->SetFlow(2,0);
    }
    if (mother->prev) mother->prev=NULL;
  }
  else {
    mother  = tree->NewKnot();
    k->prev = mother;
  }

  mother->right  = k;
  mother->part->SetFlav(k_flavs[0]);
  mother->part->SetInfo('I');
  mother->part->SetStatus(part_status::active);
  mother->t      = k->t;
  mother->tout   = sqr(k_flavs[0].PSMass()); 
  mother->x      = k->x/k->z;
  mother->stat   = 1;
  mother->E2     = 0.;
  mother->thcrit = th;
  mother->pt2lcm = p_suds[leg]->PT2();
  mother->qjv    = k->qjv;
  mother->qljv   = k->qljv;
  mother->maxjets= k->maxjets;

  Knot * sister = 0;
  if (mother->left) {
    sister = mother->left;
    // delete old color informations if flav content changed
    if (sister->part->Flav()!=k_flavs[1]) {
      sister->part->SetFlow(1,0);
      sister->part->SetFlow(2,0);
    }
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
  sister->part->SetStatus(part_status::active);
  sister->t      = 0.;
  sister->tout   = sqr(k_flavs[1].PSMass());
  sister->x      = (mother->x)*(1.-k->z);
  sister->stat   = 1;
  sister->E2     = 0.;
  sister->thcrit = th; 
  sister->qjv    = k->qjv;
  sister->qljv   = k->qljv;
  sister->maxjets= k->maxjets;

  if (k->part->Info() != 'G' && k->part->Info() != 'H') k->part->SetInfo('i');
  k->part->SetStatus(part_status::decayed);
  SetColours(k);
}

void Initial_State_Shower::SetColours(Knot * k)
{
  if (k==NULL || k->prev==NULL || k->prev->left==NULL) return;
  Knot * mother(k->prev), * sister(mother->left), * test(NULL);
  
  int all_colors_known(1), nquark(0), ngluon(0);
  for (int i=0;i<3;++i) {
    if (i==0) test = k;
    if (i==1) test = mother;
    if (i==2) test = sister;
    
    if (test->part->Flav().Strong()) {
      int nc(0);
      if (test->part->GetFlow(1)) ++nc;
      if (test->part->GetFlow(2)) ++nc;
      if (test->part->Flav().IsQuark()) {
	++nquark;
	if (nc!=1) all_colors_known=0;
      }
      else if (test->part->Flav().IsGluon()) {
	++ngluon;
	if (nc!=2) all_colors_known=0;
      }
      else {
	msg.Error()<<"ERROR in Initial_State_Shower::SetColours:"<<std::endl
		   <<"   Strongly interacting particle "<<test->part->Flav()
		   <<" not covered by SetColours "<<std::endl;
      }
    }      
  }

  if (all_colors_known) return;

  if (nquark+ngluon==2) {
    if (ngluon==2) {
      // a) g -> g + X
      // b) g -> X + g
      // c) X -> g + g
      int col[2];
      bool swap(k->part->Flav().Strong());
      if (swap) for (size_t i=0; i<2; ++i) col[i]=k->part->GetFlow(i+1);
      else for (size_t i=0; i<2; ++i) col[i]=Flow::Counter();
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


void Initial_State_Shower::InitTwoTrees(Tree ** trees,double E2) 
{
  double x1(0.005+ran.Get()*0.295), x2(0.005+ran.Get()*0.295);

  double scale(x1*x2*E2), E(0.5 * sqrt(E2));
  m_pt2_1 = m_pt2_2 = scale;
  m_th_1  = m_th_2  = M_PI;

  Knot * d1   = trees[0]->NewKnot();
  *(d1->part) = Particle(1,Flavour(kf::u),x1*E*Vec4D(1.,0.,0.,1.));
  d1->part->SetStatus(part_status::active);
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
  d2->part->SetStatus(part_status::active);
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

ATOOLS::Vec4D Initial_State_Shower::GetMomentum(Knot *const mo,int &number) 
{
  if (mo->left) 
    return GetMomentum(mo->left,number)+GetMomentum(mo->right,number);
  ++number;
  return mo->part->Momentum();
}

void Initial_State_Shower::SetStartingConditions(double k1,double k2,int flag) 
{
  switch (flag) {
  case  1 :  m_th_1  = k1; m_th_2  = k2; break;
  default :  m_pt2_1 = k1; m_pt2_2 = k2; break;
  }
}

void Initial_State_Shower::BoostFS() 
{
  if (p_fstree==NULL) return;
  Vec4D mom1(p_istrees[0]->GetRoot()->part->Momentum());
  Vec4D mom2(p_istrees[1]->GetRoot()->part->Momentum());
  Vec4D vl(mom1[0]+mom2[0],-1.*Vec3D(mom1+mom2));
  m_labboost=Poincare(vl);
  m_labboost.BoostBack(mom1);
  m_cmsrot=Poincare(Vec4D::ZVEC,mom1);
  p_fstree->BoRo(m_cmsrot);
  p_fstree->BoRo(m_labboost);
}

void Initial_State_Shower::BoostBackFS() 
{
  if (p_fstree==NULL) return;
  m_labboost.Invert();
  m_cmsrot.Invert();
  p_fstree->BoRo(m_labboost);
  p_fstree->BoRo(m_cmsrot);
}

void Initial_State_Shower::ExtractPartons(Knot *const kn,const int &beam,
					  Blob *const jet,Blob_List *const bl,
					  Particle_List *const pl) 
{
  if (kn==NULL) return;
  // fetch last PSME blobs
  m_bl_meps_is=0;
  m_bl_meps_fs=0;
  for (Blob_List::iterator blit=bl->begin();blit!=bl->end();++blit) {
    if ((*blit)->Type()==btp::ME_PS_Interface_IS) {
      m_bl_meps_is=(*blit);
    }
    if ((*blit)->Type()==btp::ME_PS_Interface_FS) {
      m_bl_meps_fs=(*blit);
    }
  }
  if (m_bl_meps_is==NULL) {
    msg.Error()<<METHOD<<"(..): No Interface found. Abort."<<std::endl;
    abort();
  }
  m_bl_meps_is->SetStatus(blob_status::inactive);
  int nr(1000);
  SingleExtract(kn,beam,jet,bl,nr);
}


void Initial_State_Shower::SingleExtract(Knot *const kn,const int &beam,
					 Blob *jet,Blob_List *const bl,
					 int &nr) 
{
  if (kn==NULL) {
    return;
  }
  Particle *p(NULL);
  bool newblob(false), lastknot(false), is_is(true), ignore(false);

  if (kn->prev==NULL) {
    newblob=true;
  }
  else {
    for (Knot * k=kn;is_is && k->prev; k=k->prev) {
      if (k==k->prev->left) is_is=false;
    }
    if (is_is && (kn->prev->part->Info()=='G' || kn->prev->part->Info()=='H'))
      ignore = true;
  }
  if (kn->left==NULL) lastknot=true;

  if (!is_is && !lastknot) {
    if (kn->left->part->Info()=='H') ignore=true;
  }
  if (!ignore && !is_is && kn->part->Info()=='H') newblob=true;

  // --- create new blob ---
  if (newblob) {
    jet = new Blob();
    jet->SetStatus(blob_status::needs_harddecays |
		   blob_status::needs_beams |
		   blob_status::needs_hadronization);
    jet->SetId();
#ifdef USING__Veto_Info
    jet->AddData("IS_VS",new Blob_Data<std::vector<int> >(p_suds[beam]->Vetos()));
#endif
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
    p = new Particle(*kn->part);
    p->SetFinalMass(sqrt(kn->tout));
    p->SetStatus(part_status::decayed);
    jet->AddToInParticles(p);
  }

  // --- add to MEPS blob ---
  if (!ignore && (kn->part->Info()=='H' || kn->part->Info()=='G')) {
    if  (is_is && m_bl_meps_is) {
      p = new Particle(*kn->part);
      p->SetFinalMass(sqrt(kn->tout));
      jet->AddToOutParticles(p);
      p->SetStatus(part_status::decayed);
      m_bl_meps_is->AddToInParticles(p);
    }
    else if (!is_is && m_bl_meps_fs) {
      if (!p) {
	p = new Particle(*kn->part);
	p->SetFinalMass(sqrt(kn->tout));
	jet->AddToInParticles(p);
      }
      else {
      }
      p->SetStatus(part_status::decayed); 
      m_bl_meps_fs->AddToOutParticles(p); 
    }
  }

  // --- add final state particle ---
  if (lastknot && !is_is) {
    p = new Particle(*kn->part);
    p->SetFinalMass(sqrt(kn->tout));
    p->SetStatus(part_status::active);
    jet->AddToOutParticles(p);
  }

  SingleExtract(kn->left,beam,jet,bl,nr); 
  SingleExtract(kn->right,beam,jet,bl,nr); 
}

bool Initial_State_Shower::TestShower(Tree ** trees) 
{
  double E2(sqr(rpa.gen.Ecms()));

  ran.ReadInStatus("RandomA.dat",2735);

  msg.Out()<<" Starting Test IS Shower :"<<std::endl;
  for (long int n=1;n<=rpa.gen.NumberOfEvents();n++) { 
    if (n%2500==0) msg.Out()<<" "<<n<<" events"<<std::endl;

    for (int i=0;i<2;i++) trees[i]->Reset();
    InitTwoTrees(trees,E2);
    if (!PerformShower(trees,NULL,0)) return 0;
  }
  msg_Events()<<"Initial_State_Shower::TestShower : "
	      <<"Terminated loops over events successfully."<<std::endl;
  return 1;
}

void Initial_State_Shower::OutputTree(Tree * tree) 
{
  if (tree->GetInitiator()==NULL) {
    msg.Out()<<"empty Tree"<<std::endl;
  }
  else {
    int number(0);
    msg.Out()<<"final Tree:"<<std::endl<<*tree<<std::endl
	     <<"Total 4 Mom = "<<GetMomentum(tree->GetInitiator(),number)
	     <<" for "<<number<<" FS particles."<<std::endl;
  }
}

