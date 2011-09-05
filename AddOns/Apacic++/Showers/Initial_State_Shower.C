#include "AddOns/Apacic++/Showers/Initial_State_Shower.H"

#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/MyStrStream.H"
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
					   Final_State_Shower *const fin,
					   MODEL::Model_Base *const model,
					   Data_Reader *const dataread) : 
  p_fin(fin),
  p_tools(new Sudakov_Tools(model)),
  p_kin(new Spacelike_Kinematics()),
  p_suds(new Spacelike_Sudakov*[2]),
  m_allowed(100)
{
  double cplscalefac(0.25*ToType<double>
		     (rpa->gen.Variable("IS_CPL_SCALE_FACTOR","1.0")));
  m_t0=dabs(dataread->GetValue<double>("IS_PT2MIN",1.0));
  double shadron(dataread->GetValue<double>("IS_MAX_SCALE",
					    sqr(rpa->gen.Ecms())));
  double emin(dataread->GetValue<double>("IS_MINIMAL_E",0.5));
  int cplscheme(dataread->GetValue<int>("IS_COUPLING_SCHEME",1));
  int pdfscheme(dataread->GetValue<int>("IS_PDF_SCALE_SCHEME",1));
  int orderingscheme(dataread->GetValue<int>("IS_ORDERING_SCHEME",0));
  for (short unsigned int i(0);i<2;++i) {
    if (isr->PDF(i)->Q2Min()>m_t0*
	ToType<double>(rpa->gen.Variable("FACTORIZATION_SCALE_FACTOR"))) {
      msg_Error()<<METHOD<<"(..):\n   IS_PT2MIN("<<m_t0
		 <<")*FACTORIZATION_SCALE_FACTOR("
		 <<rpa->gen.Variable("FACTORIZATION_SCALE_FACTOR")<<") "
		 <<"smaller than minimum scale given by PDF ("
		 <<isr->PDF(i)->Q2Min()<<").\n"
		 <<"   Please change your settings in "
		 <<dataread->InputFile()<<"\n";
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
    m_extra_pdf[i]=0;
    if (!p_suds[i]->Initialize()) 
      THROW(fatal_error,"failed to intialize is sudakov");
  }
  int maxem(-1);
  Data_Reader read(" ",";","!","=");
  if (read.ReadFromFile(maxem,"PS_MAX_EMISSIONS")) {
    msg_Info()<<METHOD<<"(): Set maximum emission number "<<maxem<<".\n";
    Spacelike_Sudakov::SetMaxEmissions(maxem);
  }
  m_t0 = - dabs(m_t0);
}

Initial_State_Shower::~Initial_State_Shower()
{
  delete p_tools;
  delete p_kin;
  for (short unsigned int i(0);i<2;i++) delete p_suds[i];
  delete [] p_suds;
}

int Initial_State_Shower::PerformShower
(Tree **const trees,Tree *const fintree) 
{
  PROFILE_HERE;
  p_fstree=fintree;
  p_istrees=trees;
  for (int i(0);i<m_allowed;++i) {
    m_extra_pdf[0]=m_extra_pdf[1]=1;
    if (InitializeSystem(trees,0,trees[0]->GetRoot(),trees[1]->GetRoot())) {
      p_kin->BoostInCMS(trees,p_fstree,0,trees[0]->GetInitiator(),
			trees[1]->GetInitiator());
      m_lab=p_kin->BoostInLab(trees,p_fstree);
      return 1;
    }
  }
  return 1;
}

void Initial_State_Shower::SetDirection(Knot *const k)
{
  if (k->right) {
    if (k->dir!=0) k->right->dir=k->dir;
    SetDirection(k->right);
  }
}

void Initial_State_Shower::SetVetoScales(Knot *const d)
{
  if (d==NULL || d->prev==NULL) return;
  if (d->prev->qjv==1.0e10) {
    d->prev->qjv=d->qjv;
  }
  if (d->prev->left->qjv==1.0e10) {
    d->prev->left->qjv=d->qjv;
  }
}

bool Initial_State_Shower::InitializeSystem(Tree **const trees,int tree1,
					    Knot *k1,Knot *k2)
{
  msg_Debugging()<<METHOD<<"("<<k1->kn_no<<","<<k2->kn_no<<"): {\n";
  msg_Indent();
  Spacelike_Sudakov::SetEmissions(0);
  if (k1==NULL || k2==NULL) THROW(fatal_error,"No trees.");
  int decay1(k1->shower>0), decay2(k2->shower>0);
  SetDirection(trees[0]->GetInitiator());
  SetDirection(trees[1]->GetInitiator());
  m_sprime=(k1->part->Momentum()+k2->part->Momentum()).Abs2();
  if (decay1 && decay2) FillBranch(trees,tree1,k1,k2);
  while (true) {
    m_sprime=(k1->part->Momentum()+k2->part->Momentum()).Abs2();
    if (decay1 || decay2) {
      int tree(tree1);
      ChooseMother(tree1,k1,k2,true);
      if (tree!=tree1) std::swap<int>(decay1,decay2);
      if (!FillBranch(trees,tree1,k1,k2)) decay1=false;
    }
    msg_Debugging()<<"t_{"<<k1->kn_no<<"} = "<<k1->t
		   <<" ("<<decay1<<"), t_{"<<k2->kn_no<<"} "
		   <<k2->t<<" ("<<decay2<<")\n";
    if (decay1 && k1->t<k1->tout) m_sprime/=k1->z;
    if (decay2 && k2->t<k2->tout) m_sprime/=k2->z;
    int stat(p_kin->InitKinematics(trees,p_fstree,tree1,k1,k2));
    if (stat==-1) break;
    if (stat==1) { 
      if (!decay1) SetColours(k1);
      if (!decay2) SetColours(k2);
      p_suds[0]->AcceptBranch(k1);
      p_suds[1]->AcceptBranch(k2);
      bool t1(k1->shower==0 || k1->t<k1->tout);
      bool t2(k2->shower==0 || k2->t<k2->tout);
      if (t1 || t2)
	if ((stat=EvolveSystem(trees,tree1,t1?k1->prev:k1,
			       t2?k2->prev:k2))==-1) 
	  break;
      if (stat==1) return true;
    }
    if (!(decay1 || decay2)) break;
  }
  trees[0]->Restore();
  trees[1]->Restore();
  if (p_fstree!=NULL) p_fstree->Restore();
  msg_Debugging()<<"iss reset";
  msg_Debugging()<<"}\n";
  return false;
}

int Initial_State_Shower::DoKinematics()
{
  if (!p_kin->DoKinematics(p_istrees,p_fstree,
			   m_tree1,p_k1,p_k2,false)) return 0;
  return 1-p_suds[m_tree1]->PTVeto(p_k1->prev);
}

int Initial_State_Shower::EvolveSystem(Tree **const trees,int tree1,
				       Knot *k3,Knot *k5)
{
  msg_Debugging()<<METHOD<<"("<<tree1<<";"<<k3->kn_no<<","<<k5->kn_no
		 <<"): { m_sprime = "<<m_sprime<<"\n";
  m_extra_pdf[tree1]=0;
  msg_Indent();
  k3->tmo=k3->t;
  k5->tmo=k5->t;
  int veto(0);
  double oldsprime(m_sprime);
  // select parton / tree to evolve according either to history or virtuality
  ChooseMother(tree1,k3,k5);
  Knot *k1(k3->right), *k4(k3->left);
  Knot *k2((k5->t<0.0&&k5->right!=NULL)?k5->right:k5); 
  while (true) {
    m_sprime=oldsprime;
    if ((k3->shower==0 || k3->t<k3->tout) && k1!=NULL) 
      p_kin->ResetMomenta(k1,trees[tree1]);
    if (k3->shower && k3->t==k3->tout) {
      if (k5->shower && k5->t==k5->tout) {
	msg_Debugging()<<"evolution finished, veto = "<<veto<<"\n";
	if (veto) {
	  if (k3->part->Info()!='H' || k3->prev==NULL ||
	      k3->prev->part->Info()!='H') {
	    k3->t=k3->tmo;
	    k3->stat=1;
	  }
	  return 0;
	}
	return 1;
      }
      return EvolveSystem(trees,tree1,k3->prev!=NULL?k3->prev:k3,k5);
    }
    if (veto && !k3->shower) return 0;
    msg_Debugging()<<"start evolution of "<<k3->kn_no<<"->"<<k1->kn_no
		   <<" at "<<k3->t<<", veto = "<<veto<<"\n";
    ++veto;
    if (k1->shower>0) {
      // set properties of t-channel particle
      k1->E2       = sqr(k1->part->Momentum()[0]);
      // set properties of mother
      k3->E2       = k1->E2/sqr(k1->z);                      
      k3->left->E2 = k3->E2*sqr(1.-k1->z);
    }  
    else {
      double s12((k1->part->Momentum()+k2->part->Momentum()).Abs2());
      double s32((k3->part->Momentum()+k2->part->Momentum()).Abs2());
      msg_Debugging()<<"x_{"<<k1->kn_no<<"} = "<<k1->x
		     <<" -> x_{"<<k3->kn_no<<"} = "
		     <<(k1->x*s32/s12)<<"\n";
      // x>1 can happen as outermost partons are allowed to radiate
      // if left radiates and right remains to be unfolded, momenta might be 
      // shifted such that x>1. in that case, retry left emission
      if (k1->x>s12/s32) return 0;
      k1->z=s12/s32;
      k3->x=k1->x/k1->z;
      m_sprime/=k1->z;
    }
    SetVetoScales(k1);
    if (k3->stat>0) {
      // dice t for non-me knot
      if (k3->part->Info()!='H' ||
	  k3->prev==NULL || k3->prev->part->Info()!='H')
	if (!FillBranch(trees,tree1,k3,k2)) k3->stat=0;
    }
    double maxt(p_kin->CalculateMaxT(k1,k2));
    if (k4->part->Info()!='H') {
      if (maxt<k4->tout) continue;
      k4->t=dabs(k1->t);
      k4->stat=3;
    }
    k4->tmax=maxt;
    if (p_fin) {
      if (!p_kin->DoKinematics(trees,p_fstree,tree1,k1,k2,true)) continue;
      double s12((k1->part->Momentum()+k2->part->Momentum()).Abs2());
      p_k1=k1;
      p_k2=k2;
      m_tree1=tree1;
      switch (p_fin->TimelikeFromSpacelike
	      (this,trees[tree1],k4,s12,k1->z,k2)) {
      case -1:
	return -1;
      case 0: 
	msg_Debugging()<<"caught jetveto\n";
	return 0;
      case 1:
	p_fin->SetAllColours(k4);
	break;
      }
    }
    else {
      k4->t     = k4->tout;
      k4->stat  = 0;
      k4->part->SetStatus(part_status::active);
      if (!p_kin->DoKinematics(trees,p_fstree,tree1,k1,k2,false)) continue;
    }
    if (k3->shower==1 && k3->t<k3->tout) {
      double mmaxt(p_kin->CalculateMaxT(k3,k2,false));
      if (mmaxt<k3->prev->left->tout) continue;
    }
    p_suds[tree1]->AcceptBranch(k3);
    k1->stat=0;
    if (k3->z>0.) m_sprime/=k3->z;
    int stat(-1);
    if (abs(stat=EvolveSystem(trees,tree1,k3->prev!=NULL?k3->prev:k3,k5))==1)
      return stat;
    msg_Debugging()<<"back in ("<<k3->kn_no<<","<<k5->kn_no<<") with "
		   <<k3->t<<" "<<k5->t<<"\n";
  }
  return 0;
}  

int Initial_State_Shower::FillBranch(Tree **const trees,const int tree1,
				     Knot *const active,Knot *const partner) {
  msg_Debugging()<<METHOD<<"("<<active->kn_no<<","<<partner->kn_no<<","
 		 <<tree1<<"): {\n";
  msg_Indent();
  Flavour flavs[2];
  SetVetoScales(active->right);
  if (active->prev && active->prev->part->Info()=='H') return true;
  if (p_suds[tree1]->Dice(active,m_sprime,
			  sqr(p_kin->MS()->Mass(partner->part->Flav())),
			  m_extra_pdf[tree1])) {
    flavs[0] = p_suds[tree1]->GetFlA();
    flavs[1] = p_suds[tree1]->GetFlC();    
    FillMotherAndSister(trees[tree1],active,flavs);
    active->prev->pt2lcm = p_suds[tree1]->PT2();
    msg_Debugging()<<"test emission at t = "<<active->t
 		   <<", z = "<<active->z<<", pt = "
		   <<sqrt(active->prev->pt2lcm)<<"\n";
    msg_Debugging()<<"}\n";
    p_suds[tree1]->AddEmission();
    return 1;
  }
  active->prev   = 0;
  active->stat   = 0;
  active->t      = active->tout;
  active->part->SetStatus(part_status::active);
  msg_Debugging()<<"no branch\n";
  msg_Debugging()<<"}\n";
  return 0;
}

void Initial_State_Shower::
ChooseMother(int &ntree1,Knot *&k1, Knot *&k2,bool init)
{
  bool swap(k1->t>k2->t);
  if (init && (k1->prev==NULL || k2->prev==NULL)) swap=k1->prev!=NULL;
  bool known1(k1->part->Info()=='H' || k1->shower==0);
  bool known2(k2->part->Info()=='H' || k2->shower==0);
  msg_Debugging()<<METHOD<<"("<<k1->kn_no<<","<<k2->kn_no
		 <<","<<init<<"): t1 = "<<k1->t<<", t2 = "<<k2->t
		 <<", known1 = "<<known1<<", known2 = "<<known2<<"\n";
  if (known1 || known2) {
    if (known1 && !known2) swap=init;
    else if (!known1 && known2) swap=!init;
    if (known1 && known2) {
      if (k1->kn_no<k2->kn_no) swap=init;
      else swap=!init;
    }
  }
  if (k1->t>=0.0 && k2->t<0.0) swap=true;
  if (k2->t>=0.0 && k1->t<0.0) swap=false;
  if (swap) {
    std::swap<Knot*>(k1,k2);
    ntree1=1-ntree1;
  }
  msg_Debugging()<<"  swap = "<<swap<<" -> k1 = "<<k1->kn_no
		 <<", k2 = "<<k2->kn_no<<"\n";
}

void Initial_State_Shower::FillMotherAndSister(Tree * tree,Knot * k,
					       Flavour * k_flavs)
{
  double pt2max(sqr(rpa->gen.Ecms()));
  double th(4.*k->z*k->z*k->t/(4.*k->z*k->z*k->t-(1.-k->z)*k->x*k->x*pt2max));
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
  mother->tout   = sqr(p_kin->MS()->Mass(k_flavs[0])); 
  mother->x      = k->x/k->z;
  mother->stat   = 1;
  mother->E2     = 0.;
  mother->thcrit = th;
  mother->qjv    = k->qjv;
  mother->dir    = k->dir;

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
  sister->tout   = sqr(p_kin->MS()->Mass(k_flavs[1]));
  sister->x      = (mother->x)*(1.-k->z);
  sister->stat   = 1;
  sister->E2     = 0.;
  sister->thcrit = th; 
  sister->qjv    = k->qjv;

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
      if (test->part->Flav().IsQuark() || test->part->Flav().IsSquark()) {
	++nquark;
	if (nc!=1) all_colors_known=0;
      }
      else if (test->part->Flav().IsGluon() || test->part->Flav().IsGluino()) {
	++ngluon;
	if (nc!=2) all_colors_known=0;
      }
      else {
	msg_Error()<<"ERROR in Initial_State_Shower::SetColours:"<<std::endl
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

  if (mother->part->Flav().IsQuark() || mother->part->Flav().IsSquark()) {
    if (mother->part->Flav().IsAnti()) {
      if (k->part->Flav().IsQuark() || k->part->Flav().IsSquark()) {
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
      if (k->part->Flav().IsQuark() || k->part->Flav().IsSquark()) {
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
    if (k->part->Flav().IsQuark() || k->part->Flav().IsSquark()) {
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
  *(d1->part) = Particle(1,Flavour(kf_u),x1*E*Vec4D(1.,0.,0.,1.));
  d1->part->SetStatus(part_status::active);
  d1->part->SetInfo('G');
  d1->part->SetFlow(1,500);
  d1->part->SetFlow(2,501);
  d1->t       = -scale;
  d1->tout    = sqr(p_kin->MS()->Mass(Flavour(kf_u))); 
  d1->x       = x1;
  d1->E2      = sqr(x1*E);
  d1->maxpt2  = scale;
  d1->costh   = -1.; 
  d1->thcrit  = M_PI;
  d1->stat    = 1;
  
  Knot * d2   = trees[1]->NewKnot();
  *(d2->part) = Particle(2,Flavour(kf_u).Bar(),x2*E*Vec4D(1.,0.,0.,-1.));
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

void Initial_State_Shower::SingleExtract(Blob *jet,Knot *const kn) 
{
  if (kn==NULL) return;
  if (kn->prev==NULL) {
    Particle *p(new Particle(*kn->part));
    p->SetInfo('I');
    p->SetFinalMass(sqrt(kn->tout));
    p->SetStatus(part_status::active);
    p->SetNumber();
    jet->AddToInParticles(p);
  }
  p_fin->ExtractFinalState(jet,kn->left); 
  SingleExtract(jet,kn->right); 
}

void Initial_State_Shower::OutputTree(Tree * tree) 
{
  if (tree->GetInitiator()==NULL) {
    msg_Out()<<"empty Tree"<<std::endl;
  }
  else {
    int number(0);
    msg_Out()<<"final Tree:"<<std::endl<<*tree<<std::endl
	     <<"Total 4 Mom = "<<GetMomentum(tree->GetInitiator(),number)
	     <<" for "<<number<<" FS particles."<<std::endl;
  }
}

void Initial_State_Shower::SetMS(ATOOLS::Mass_Selector *const ms)
{
  p_suds[0]->SetMS(ms);
  p_suds[1]->SetMS(ms);
  p_kin->SetMS(ms);
}
