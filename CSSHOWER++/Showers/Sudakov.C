#include "CSSHOWER++/Showers/Sudakov.H"

#include "CSSHOWER++/Showers/Splitting_Function_Base.H"
#include "CSSHOWER++/Tools/Singlet.H"
#include "CSSHOWER++/Showers/Shower.H"
#include "MODEL/Main/Single_Vertex.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/My_Limits.H"
#include "ATOOLS/Org/My_MPI.H"

using namespace CSSHOWER;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;

bool CSSHOWER::Sudakov::s_init=false;

Sudakov::Sudakov(PDF::ISR_Handler *isr,const int qcd,const int qed) :
  m_qcdmode(qcd), m_ewmode(qed), m_pdfmin(1.0e-4, 1.0e-2),
  m_reweightscalecutoff {0.0}, m_keeprewinfo {false}
{
  p_pdf = new PDF::PDF_Base*[2];
  for (int i=0;i<2; i++) p_pdf[i] = isr->PDF(i);
}

Sudakov::~Sudakov()
{
  delete [] p_pdf;
  for (size_t i(0);i<m_addsplittings.size();++i) delete m_addsplittings[i];
  for (size_t i(0);i<m_cgets.size();++i) delete m_cgets[i];
  s_init=false;
}

struct FTrip {
public:
  Flavour m_a, m_b, m_c;
public:
  inline FTrip(const Flavour &a,const Flavour &b,const Flavour &c):
    m_a(a), m_b(b), m_c(c) {}
  bool operator<(const FTrip &f) const
  {
    if (m_a<f.m_a) return true;
    if (m_a==f.m_a) {
      if (m_b<f.m_b) return true;
      if (m_b==f.m_b) {
	return m_c<f.m_c;
      }
    }
    return false;
  }
};

void Sudakov::InitSplittingFunctions(MODEL::Model_Base *md,const int kfmode)
{
  if (!s_init) {
    s_init=true;
    SFC_Filler_Getter::Getter_List flist(SFC_Filler_Getter::GetGetters());
    for (SFC_Filler_Getter::Getter_List::const_iterator git(flist.begin());
         git!=flist.end();++git) (*git)->GetObject(SFC_Filler_Key(md,&m_cgets));
    if (msg_LevelIsDebugging()) {
      msg_Out()<<METHOD<<"(): {\n\n"
               <<"   // available coupling calcs\n\n";
      SFC_Getter::PrintGetterInfo(msg->Out(),25);
      msg_Out()<<"\n   // available lorentz calcs\n\n";
      SFL_Getter::PrintGetterInfo(msg->Out(),25);
      msg_Out()<<"\n}"<<std::endl;
    }
  }
  msg_Debugging()<<METHOD<<"(): Init splitting functions {\n";
  msg_Indent();
  std::set<FTrip> sfs;
  const Vertex_Table *vtab(md->VertexTable());
  for (Vertex_Table::const_iterator
	 vlit=vtab->begin();vlit!=vtab->end();++vlit) {
    for (Vertex_List::const_iterator
	   vit=vlit->second.begin();vit!=vlit->second.end();++vit) {
      Single_Vertex *v(*vit);
      if (v->NLegs()>3) continue;
      if (sfs.find(FTrip(v->in[0].Bar(),v->in[1],v->in[2]))!=sfs.end()) continue;
      sfs.insert(FTrip(v->in[0].Bar(),v->in[1],v->in[2]));
      sfs.insert(FTrip(v->in[0].Bar(),v->in[2],v->in[1]));
      msg_Debugging()<<"Add "<<v->in[0].Bar()<<" -> "<<v->in[1]<<" "<<v->in[2]<<" {\n";
      {
	msg_Indent();
	int dmode(0);
	if (v->in[2]==v->in[0].Bar()) dmode=1;
	else if (v->in[1]!=v->in[0].Bar() &&
		 v->in[1].IsAnti() && !v->in[2].IsAnti()) dmode=1;
	Add(new Splitting_Function_Base(SF_Key(v,dmode,cstp::FF,kfmode,m_qcdmode,m_ewmode,1,m_pdfmin)));
	Add(new Splitting_Function_Base(SF_Key(v,dmode,cstp::FF,kfmode,m_qcdmode,m_ewmode,-1,m_pdfmin)));
	Add(new Splitting_Function_Base(SF_Key(v,dmode,cstp::FI,kfmode,m_qcdmode,m_ewmode,1,m_pdfmin)));
	Add(new Splitting_Function_Base(SF_Key(v,dmode,cstp::FI,kfmode,m_qcdmode,m_ewmode,-1,m_pdfmin)));
	if (v->in[0].Bar().Mass()<100.0 && v->in[1].Mass()<100.0 && v->in[2].Mass()<100.0) {
	  Add(new Splitting_Function_Base(SF_Key(v,dmode,cstp::IF,kfmode,m_qcdmode,m_ewmode,1,m_pdfmin)));
	  Add(new Splitting_Function_Base(SF_Key(v,dmode,cstp::IF,kfmode,m_qcdmode,m_ewmode,-1,m_pdfmin)));
	  Add(new Splitting_Function_Base(SF_Key(v,dmode,cstp::II,kfmode,m_qcdmode,m_ewmode,1,m_pdfmin)));
	  Add(new Splitting_Function_Base(SF_Key(v,dmode,cstp::II,kfmode,m_qcdmode,m_ewmode,-1,m_pdfmin)));
	}
	if (v->in[1]!=v->in[2]) {
	  if (v->in[0].Bar().Mass()<100.0 && v->in[1].Mass()<100.0 && v->in[2].Mass()<100.0) {
	    Add(new Splitting_Function_Base(SF_Key(v,1-dmode,cstp::IF,kfmode,m_qcdmode,m_ewmode,1,m_pdfmin)));
	    Add(new Splitting_Function_Base(SF_Key(v,1-dmode,cstp::IF,kfmode,m_qcdmode,m_ewmode,-1,m_pdfmin)));
	    Add(new Splitting_Function_Base(SF_Key(v,1-dmode,cstp::II,kfmode,m_qcdmode,m_ewmode,1,m_pdfmin)));
	    Add(new Splitting_Function_Base(SF_Key(v,1-dmode,cstp::II,kfmode,m_qcdmode,m_ewmode,-1,m_pdfmin)));
	  }
	}
      }
      msg_Debugging()<<"}\n";
    }
  }
  msg_Debugging()<<"}\n";
}

void Sudakov::SetCoupling(MODEL::Model_Base *md,
			  const double &k0sqi,const double &k0sqf,
			  const double &isfac,const double &fsfac,
			  const double &k0sq_gsplit_fac)
{
  m_k0sqi=k0sqi;
  m_k0sqf=k0sqf;
  m_k0sq_gsplit_fac=k0sq_gsplit_fac;
  for (std::vector<Splitting_Function_Base*>::iterator
	 sit(m_splittings.begin());sit!=m_splittings.end();)
    if (!(*sit)->Coupling()->SetCoupling(md,m_k0sqi,m_k0sqf,isfac,fsfac)) {
      delete *sit;
      sit=m_splittings.erase(sit);
    }
    else {
      ++sit;
    }
  for (std::vector<Splitting_Function_Base*>::iterator
	 sit(m_addsplittings.begin());sit!=m_addsplittings.end();)
    if (!(*sit)->Coupling()->SetCoupling(md,m_k0sqi,m_k0sqf,isfac,fsfac)) {
      delete *sit;
      sit=m_addsplittings.erase(sit);
    }
    else {
      ++sit;
    }
}

void Sudakov::Add(Splitting_Function_Base * split)
{
  if (split->On()<0) {
    delete split;
    return;
  }
  if (split->On()) {
    split->SetFacScaleFactor(m_facscalefactor);
    Splitting_Function_Group::Add(split);
    msg_Debugging()<<" -> add\n";
  }
  AddToMaps(split,!split->On());
}

void Sudakov::AddToMaps(Splitting_Function_Base * split,const int mode)
{
  if (split->On()<0) {
    delete split;
    return;
  }
  split->SetMassThreshold(m_mth);
  split->SetScaleScheme(m_scs);
  split->SetPDF(p_pdf);
  split->SetEFac(p_shower);
  if (mode) {
    m_addsplittings.push_back(split);
    msg_Debugging()<<"\n";
  }
  Flavour flavA = split->GetFlavourA();
  Flavour flavB = split->GetFlavourB();
  Flavour flavC = split->GetFlavourC();
  if (split->GetCol()<0) {
    switch(split->GetType()) {
    case cstp::IF: m_fifmap[flavA][flavC][flavB]=split; break;
    case cstp::II: m_fiimap[flavA][flavC][flavB]=split; break;
    default: break;
    }
    return;
  }
  switch(split->GetType()) {
  case cstp::FF:
    m_fffmap[flavB][flavC][flavA]=split;
    if (split->On()) m_sffmap[flavB][flavC][flavA]=split;
    break;
  case cstp::FI:
    m_ffimap[flavB][flavC][flavA]=split;
    if (split->On()) m_sfimap[flavB][flavC][flavA]=split;
    break;
  case cstp::IF:
    m_iffmap[flavA][flavC][flavB]=split;
    if (split->On()) m_sifmap[flavA][flavC][flavB]=split;
    break;
  case cstp::II:
    m_ifimap[flavA][flavC][flavB]=split;
    if (split->On()) m_siimap[flavA][flavC][flavB]=split;
    break;
  case cstp::none: break;
  }
}

bool Sudakov::Generate(Parton * split, double kt2win)
{
  p_split  = split;
  m_cfl    = p_split->GetFlavour();
  std::list<Parton*> slist;
  // No spectators means that the parton cannot be split.
  if (!MakeSpectatorList(slist)) return false;
  // Initialise trial run.
  p_split->SetForcedDecay(false);
  m_weight = 1.0;
  m_type   = cstp::none;
  double t, z, t0 = Min(m_k0sqi,m_k0sqf), y, phi;
  Parton * spect  = NULL;
  Splitting_Function_Base * selected = NULL;
  for (list<Parton *>::iterator pit=slist.begin();pit!=slist.end();pit++) {
    spect = (*pit);
    //msg_Out()<<"--- "<<METHOD<<" checks colours ("<<spect<<"/ size ="<<slist.size()<<"): "
    //	     <<"["<<p_split->GetFlow(1)<<" "<<p_split->GetFlow(2)<<"] "
    //	     <<"["<<spect->GetFlow(1)<<" "<<spect->GetFlow(2)<<"]\n";			 
    if (p_split->ForcedDecay()) break;
    int success = Generate(spect,t0,kt2win,t,y,z,phi);
    //if (success==-1)
    //msg_Out()<<"--- "<<METHOD<<" results in forced splitting for "<<p_split<<", t = "<<t<<")\n";
    if (success!=0) {
      //msg_IODebugging()<<"shrink evolution window "<<t0<<" -> "<<t<<"\n";
      //msg_Out()<<"--- "<<METHOD<<" shrink evolution window "<<t0<<" -> "<<t<<", "
      //	       <<"spect = "<<spect<<", selected = "<<p_selected<<".\n";
      m_sy     = y;
      m_sz     = z;
      m_sphi   = phi;
      m_st     = t0 = t;
      selected = p_selected;
      //if (success==-1) msg_Out()<<"--- "<<METHOD<<" forced successful. "
      //			<<"Spectator = "<<spect<<", "
      //			<<"split = "<<p_split<<" ["<<p_split->ForcedDecay()<<"]\n";
    }
  }
  //if (p_split->ForcedDecay()) msg_Out()<<"--> forced, split = "<<p_split<<", spect = "<<spect<<"\n";
  p_spect    = NULL;
  p_selected = NULL;
  if (selected) {
    p_split->SetSpect(spect);
    p_spect    = spect;
    p_selected = selected;
    //msg_Out()<<"--- "<<METHOD<<" selected "<<p_spect<<", t = "<<m_st
    //	     <<", y = "<<m_sy<<", z = "<<m_sz<<", phi = "<<m_sphi<<", "
    //	     <<"forced = "<<p_split->ForcedDecay()<<"\n";
  }
  ClearSpecs();
  ResetLastInt();
  return (p_spect!=NULL);
}

bool Sudakov::MakeSpectatorList(list<Parton*> & slist) {
  int colourcharge = m_cfl.StrongCharge();
  if (((colourcharge==8 || (p_split->GetType()==pst::FS?colourcharge:-colourcharge)==3) &&
       p_split->GetLeft()==NULL) ||
      ((colourcharge==8 || (p_split->GetType()==pst::FS?colourcharge:-colourcharge)==-3) &&
       p_split->GetRight()==NULL)) {
    msg_Error()<<METHOD<<" throws an error exception:\n"<<(*p_split)<<".\n";
    THROW(fatal_error,"Invalid color flow.");
  }
  msg_IODebugging()<<"--- "<<METHOD<<":\n"
		   <<"   adds spectators for [type = "<<p_split->GetType()<<"]"
		   <<" ("<<p_split->GetFlow(1)<<", "<<p_split->GetFlow(2)<<").\n";
  // Adding colour spectators under the assumption that the left and right partons are
  // correctly colour-connected to the splitter.
  if (p_split->GetLeft() && p_split!=p_split->GetLeft())  {
    slist.push_back(p_split->GetLeft());
    msg_IODebugging()<<"    --> add left: "<<p_split->GetLeft()->GetType()
		     <<" ("<<p_split->GetLeft()->GetFlow(1)<<", "
		     <<p_split->GetLeft()->GetFlow(2)<<").\n";
  }
  if (p_split->GetRight() && p_split!=p_split->GetRight()) {
    slist.push_back(p_split->GetRight());
    msg_IODebugging()<<"   --> add right: "<<p_split->GetRight()->GetType()
		     <<" ("<<p_split->GetRight()->GetFlow(1)<<", "
		     <<p_split->GetRight()->GetFlow(2)<<").\n";
  }
  msg_IODebugging()<<"    ===> found "<<slist.size()<<" spectator(s). for "<<m_cfl<<"\n";
  // Adding additional (EW) charge spectators under the assumption that they may be
  // different from the left and right partons, which are correctly colour-connected to the
  // splitter and already in the spectator list.  Condition is that the splitter and spectator
  // either have opposite charges or that the splitter is a photon or Z and the spectator is
  // charged.
  int splitcharge = (p_split->GetType()==pst::IS?-1:1) * m_cfl.IntCharge();
  if (splitcharge!=0 || m_cfl.IsPhoton() || m_cfl.Kfcode()==kf_Z) {
    Singlet *sing(p_split->GetSing());
    for (PLiter pit(sing->begin());pit!=sing->end();++pit) {
      int spectcharge = ((*pit)->GetType()==pst::IS?-1:1) * (*pit)->GetFlavour().IntCharge();
      if (*pit!=p_split->GetLeft() && *pit!=p_split->GetRight() &&
	  spectcharge!=0 && (splitcharge==0 || splitcharge*spectcharge<0))
	slist.push_back(*pit);
    }
  }
  return (slist.size()>0);
}

int Sudakov::Generate(Parton *spect,double t0,double kt2win,double &t,double &y,double &z,double &phi)
{
  ClearSpecs();
  ResetLastInt();
  p_spect   = spect;
  double Q2 = 0.;
  int beam  = -1;
  m_flspec  = spect->GetFlavour();
  if (!DefineStartingPoint(spect,Q2,beam)) return false;
  double last=0.0;
  for (size_t i(0);i<m_splittings.size();++i)
    m_partint[i] = last += m_splittings[i]->Last();
  if (!IsEqual(m_partint.back(),m_lastint))
    msg_Error()<<METHOD<<"(): Error, last = "<<m_lastint
	       <<" vs. "<<m_partint.back()<<"."<<std::endl;
  m_lastint = m_partint.back();

  t = p_split->KtStart();
  msg_IODebugging()<<"starting scale "<<t<<", cutoff scale "<<t0<<"\n";
  if (p_spect==p_split->GetLeft() && t>p_split->KtSoft(0)) {
    msg_IODebugging()<<"reset starting scale "<<t<<" -> "<<p_split->KtSoft(0)<<"\n";
    t=p_split->KtSoft(0);
  }
  if (p_spect==p_split->GetRight() && t>p_split->KtSoft(1)) {
    msg_IODebugging()<<"reset starting scale "<<t<<" -> "<<p_split->KtSoft(1)<<"\n";
    t=p_split->KtSoft(1);
  }
  double x = 0.0;
  //msg_Out()<<"--- "<<METHOD<<" starting trials for "<<p_split<<" ["<<p_split->ForcedDecay()<<"] & "
  //	   <<spect<<", t = "<<t<<"\n";
  while (t>=Max(t0,kt2win)) {
    t = ProduceT(t);
    //if (p_split->GetType()==pst::IS && p_split->GetFlavour().HadMass()>3.)
    //msg_Out()<<"------ new t = "<<t<<" vs. "
    //	       <<sqr(p_split->GetFlavour().HadMass())<<" and "<<Max(t0,kt2win)<<"\n";
    if (m_forced_splittings &&
	p_split->GetType()==pst::IS &&
	t<sqr(p_split->GetFlavour().HadMass()) && t>t0) {
      if (!((p_split->GetFlow(1)!=0 && p_split->GetFlow(1)==spect->GetFlow(2)) ||
	    (p_split->GetFlow(2)!=0 && p_split->GetFlow(2)==spect->GetFlow(1)))) {
	t = t0;
	return 0;
      }
      if (MakeForcedSplitting(spect,t,z,phi)) return -1;
    }
    else p_split->SetForcedDecay(false);
    SelectOne();
    if (!ValidT(t,t0)) return false;
    p_split->SetSpect(p_spect=p_selected->SelectSpec());
    z         = Z();
    m_flspec   = p_spect->GetFlavour();
    switch (FixKinematics(spect,t,z,Q2,x,y)) {
    case -1: return false;
    case 0:  continue;
    case 1:
    default: break;
    }
    const bool veto(Veto(Q2,x,t,y,z));
    if (m_keeprewinfo) {
      const double accwgt(Selected()->LastAcceptanceWeight());
      const double lastscale(Selected()->LastScale());
      if (accwgt < 1.0 && accwgt > 0.0 && lastscale > m_reweightscalecutoff) {
        Sudakov_Reweighting_Info info;
        info.accepted = veto;
        info.scale    = lastscale;
        info.accwgt   = accwgt;
        info.lastj    = Selected()->Lorentz()->LastJ();
        info.lastcpl  = Selected()->Coupling()->Last();
        info.sf       = Selected();
        info.x        = x;
        info.y        = y;
        info.z        = z;
        info.flspec   = Selected()->Lorentz()->FlSpec();
        p_split->SudakovReweightingInfos().push_back(info);
      }
    }
    if (veto) {
      phi = 2.0*M_PI*ran->Get();
      //msg_Out()<<"--- "<<METHOD<<" yields normal splitting with t = "<<t<<", z = "<<z<<", "
      //       <<"colours = ["<<p_spect->GetFlow(1)<<" "<<p_spect->GetFlow(2)<<"]\n";
      return 1;
    }
  }
  return 0;
}

bool Sudakov::MakeForcedSplitting(Parton * spect,double & t,double & z,double & phi) {
  if (!FixOne(Flavour(kf_gluon),p_split->GetFlavour(),
	      spect->GetType()==pst::IS ? cstp::II : cstp::IF)) return false;
  t      = sqr(p_split->GetFlavour().HadMass());
  do { z = Z(); } while (pow(z/m_zmax,m_gluon_xscaling_in_forced_splittings)<ran->Get());
  phi    = 2.*M_PI*ran->Get();
  p_split->SetSpect(spect);
  p_split->SetForcedDecay(true);
  //msg_Out()<<"--- "<<METHOD<<" tries a forced decay t = "<<t<<" "
  //	   <<"[forced = "<<p_split->ForcedDecay()<<"], "
  //	   <<"spect = "<<p_split->GetSpect()<<"\n";
  return true;
}

int Sudakov::FixKinematics(Parton * spect,const double & t,const double & z,
			   double & Q2,double & x,double & y) {
  m_type=p_split->GetType()==pst::FS?
    (p_split->GetSpect()->GetType()==pst::FS?cstp::FF:cstp::FI):
    (p_split->GetSpect()->GetType()==pst::FS?cstp::IF:cstp::II);
  const Mass_Selector *ms(p_selected->Lorentz()->MS());
  switch (m_type) {
  case (cstp::FF) : {
    double mi2 = sqr(ms->Mass(((*m_splitter)->GetFlavourB())));
    double mj2 = sqr(ms->Mass(((*m_splitter)->GetFlavourC())));
    double mk2 = sqr(ms->Mass(m_flspec));
    Q2 = (p_split->Momentum()+p_split->GetSpect()->Momentum()).Abs2();
    if (Q2<=mi2+mj2+mk2) return -1;
    y  = p_shower->KinFF()->GetY(Q2,t,z,mi2,mj2,mk2,
				 (*m_splitter)->GetFlavourA(),
				 (*m_splitter)->GetFlavourC());
    x  = 0.;
    if (y<0.0 || y>1.0) return 0;
    break;
  }
  case (cstp::FI) : {
    double mi2  = sqr(ms->Mass(((*m_splitter)->GetFlavourB())));
    double mj2  = sqr(ms->Mass(((*m_splitter)->GetFlavourC())));
    double ma2  = sqr(ms->Mass(m_flspec));
    double mij2 = sqr(ms->Mass(((*m_splitter)->GetFlavourA())));
    Q2 = -(p_split->Momentum()-p_split->GetSpect()->Momentum()).Abs2();
    y  = p_shower->KinFI()->GetY(-Q2,t,z,mi2,mj2,ma2,
				 (*m_splitter)->GetFlavourA(),
				 (*m_splitter)->GetFlavourC());
    y  = 1.0-y*(-Q2-mij2-ma2)/(-Q2-mi2-mj2-ma2);
    x  = p_split->GetSpect()->Xbj();
    if (y<0.0 || y>1.0-x) return 0;
    break;
  }
  case (cstp::IF) : {
    double ma2 = sqr(ms->Mass(((*m_splitter)->GetFlavourA())));
    double mi2 = sqr(ms->Mass(((*m_splitter)->GetFlavourC())));
    double mk2 = sqr(ms->Mass(m_flspec));
    Q2 = -(p_split->Momentum()-p_split->GetSpect()->Momentum()).Abs2();
    y  = p_shower->KinIF()->GetY(-Q2,t,z,ma2,mi2,mk2,
				 (*m_splitter)->GetFlavourB(),
				 (*m_splitter)->GetFlavourC());
    x  = p_split->Xbj();
    if (y<0.0 || y>1.0 || z<x) return 0;
    break;
  }
  case (cstp::II) : {
    double ma2 = sqr(ms->Mass(((*m_splitter)->GetFlavourA())));
    double mi2 = sqr(ms->Mass(((*m_splitter)->GetFlavourC())));
    double mb2 = sqr(ms->Mass(m_flspec));
    Q2 = (p_split->Momentum()+p_split->GetSpect()->Momentum()).Abs2();
    y  = p_shower->KinII()->GetY(Q2,t,z,ma2,mi2,mb2,
				 (*m_splitter)->GetFlavourB(),
				 (*m_splitter)->GetFlavourC());
    x  = p_split->Xbj();
    if (y<0.0 || y>1.0-z || z<x) return 0;
    break;
  }
  default:
    THROW(fatal_error,"Undefined kinematics type");
  }
  return 1;
}

 bool Sudakov::ValidT(const double & t,const double & t0) {
  double k0sq(p_split->GetType()==pst::IS?m_k0sqi:m_k0sqf);
  if (t<Max(t0,k0sq) ||
      (p_selected->GetFlavourA().IsGluon() &&
       p_selected->GetFlavourB().IsQuark() &&
       t<k0sq*m_k0sq_gsplit_fac)) return false;
  return true;
}

bool Sudakov::DefineStartingPoint(Parton * spect,double & Q2, int & beam) {
  switch (p_split->GetType()) {
  case pst::FS:
    if (spect->GetType()==pst::FS) {
      Q2   = (p_split->Momentum()+spect->Momentum()).Abs2();
      beam = -1;
      if (!DefineFFBoundaries(Q2,1.)) return false;
    }
    if (spect->GetType()==pst::IS) {
      Q2   = -(p_split->Momentum()-spect->Momentum()).Abs2();
      beam = spect->Beam();
      if (!DefineFIBoundaries(Q2,spect->Xbj(),beam)) return false;
    }
    break;
  case pst::IS:
    if (spect->GetType()==pst::FS) {
      Q2   = -(p_split->Momentum()-spect->Momentum()).Abs2();
      beam = p_split->Beam();
      if (IsNan(Q2)) {
	msg_Out()<<METHOD<<" has no meaningful Q2 from:\n"
		 <<"   "<<p_split->Momentum()<<" - "<<spect->Momentum()<<"\n"<<(*p_split)<<(*spect)
		 <<"   will return false and hope for the best.\n";
	exit(1);
	return false;
      }
      if (!DefineIFBoundaries(Q2,p_split->Xbj(),beam)) return false;
    }
    if (spect->GetType()==pst::IS) {
      Q2   = (p_split->Momentum()+spect->Momentum()).Abs2();
      beam = p_split->Beam();
      if (!DefineIIBoundaries(Q2,p_split->Xbj(),beam)) return false;
    }
    break;
  case pst::none:
    msg_Out()<<"Error in "<<METHOD<<": No pst-type for splitter.\n"
	     <<(*p_split);
    //return false;
    Abort();
  }
  if (m_type==cstp::none) {
    msg_Out()<<"Error in "<<METHOD<<": No type for splitter.\n"<<(*p_split);
    //return false;
    Abort();
  }
  if (m_lastint<=0.0 || IsBad(m_lastint)) return false;
  return true;
}

bool Sudakov::DefineFFBoundaries(double Q2,double x)
{
  if (4.*m_k0sqf>Q2) return false;

  m_type   = cstp::FF;
  double deltaz(sqrt(1.-4.*m_k0sqf/Q2));
  m_zmin   = 0.5*(1.-deltaz);
  m_zmax   = 0.5*(1.+deltaz);
  m_scale  = p_split->KtStart();
  double over(OverIntegrated(m_zmin,m_zmax,m_scale,Q2));
  if (over<0. || IsNan(over)) {
    msg_Error()<<"Error in "<<METHOD<<"\n"
    	       <<"   Integral for SF's<0 :"
	       <<"{"<<m_zmin<<","<<m_zmax<<","<<m_scale<<"}"<<endl;
    return false;
  }
  return true;
}

bool Sudakov::DefineFIBoundaries(double Q2,double x,int beam)
{
  if (p_pdf[beam]==NULL) return false;
  double xmax = Min(0.999999,p_pdf[beam]->XMax());
  double xmin = Max(1.e-6,p_pdf[beam]->XMin());
  if (x>=xmax || x<=xmin)                                   return false;
  if (m_k0sqf*x>Q2*(1.-x))                                   return false;
  if (Q2<=p_pdf[beam]->Q2Min() || Q2>=p_pdf[beam]->Q2Max()) return false;

  m_type = cstp::FI;
  double deltaz(1.0-4.0*Min(1.0,x/(1.0-x))*(m_k0sqf/Q2));
  if (deltaz<0.0) return false;
  deltaz=sqrt(deltaz);
  m_zmin   = 0.5*(1.0-deltaz);
  m_zmax   = 0.5*(1.0+deltaz);
  m_scale  = p_split->KtStart();
  double over(OverIntegrated(m_zmin,m_zmax,m_scale,x,beam));
  if (over<0. || IsNan(over)) {
    msg_Error()<<"Error in "<<METHOD<<"\n"
    	       <<"   Integral for SF's<0 :"
	       <<"{"<<m_zmin<<","<<m_zmax<<","<<m_scale<<"}"<<endl;
    return false;
  }
  return true;
}

bool Sudakov::DefineIFBoundaries(double Q2,double x,int beam)
{
  if (p_pdf[beam]==NULL) return false;
  double xmax = Min(0.999999,p_pdf[beam]->XMax());
  double xmin = Max(1.e-6,p_pdf[beam]->XMin());
  if (x>=xmax || x<=xmin)                                   return false;
  if (m_k0sqi>Q2)                                           return false;
  if (Q2<=p_pdf[beam]->Q2Min() || Q2>=p_pdf[beam]->Q2Max()) return false;

  m_type   = cstp::IF;
  m_zmin   = x/xmax;
  m_zmax   = Q2/(Q2+m_k0sqi);
  if (m_zmin>m_zmax) return false;
  m_scale  = p_split->KtStart();
  double over(OverIntegrated(m_zmin,m_zmax,m_scale,x,beam));
  if (over<0. || IsNan(over)) {
    msg_Error()<<"Error in "<<METHOD<<"\n"
    	       <<"   Integral for SF's<0 :"
	       <<"{"<<m_zmin<<","<<m_zmax<<","<<m_scale<<"}"<<endl;
    return false;
  }
  return true;
}

bool Sudakov::DefineIIBoundaries(double Q2,double x,int beam)
{
  if (p_pdf[beam]==NULL) return false;
  double xmax = Min(0.999999,p_pdf[beam]->XMax());
  double xmin = Max(1.e-6,p_pdf[beam]->XMin());
  if (x>=xmax || x<=xmin)                                   return false;
  if (m_k0sqi>Q2)                                           return false;
  if (Q2<=p_pdf[beam]->Q2Min() || Q2>=p_pdf[beam]->Q2Max()) return false;

  m_type=cstp::II;
  m_zmin   = x/xmax;
  m_zmax   = Q2/(Q2+m_k0sqi);
  if (m_zmin>m_zmax) return false;
  m_scale  = p_split->KtStart();
  double over(OverIntegrated(m_zmin,m_zmax,m_scale,x,beam));
  if (over<0. || IsNan(over)) {
    msg_Error()<<"Error in "<<METHOD<<"\n"
    	       <<"   Integral for SF's<0 :"
	       <<"{"<<m_zmin<<","<<m_zmax<<","<<m_scale<<"}"<<endl;
    return false;
  }
  return true;
}


double Sudakov::OverIntegrated(const double zmin,const double zmax,
			       const double scale,const double xbj,int beam) {
  for (m_splitter=m_splittings.begin();m_splitter!=m_splittings.end();m_splitter++) {
    if ((*m_splitter)->GetType()==m_type &&
	(p_split->GetLeft()==p_spect || p_split->GetRight()==p_spect ||
	 (*m_splitter)->Coupling()->AllowSpec(m_flspec))) {
      if ((*m_splitter)->PureQCD() &&
	  !(p_split->GetLeft()==p_spect || p_split->GetRight()==p_spect)) continue;
      bool match=false;
      switch (m_type) {
      case cstp::FF:
      case cstp::FI:
	if ((*m_splitter)->GetFlavourA()==m_cfl) match=true;
	break;
      case cstp::IF:
      case cstp::II:
	if ((*m_splitter)->GetFlavourB()==m_cfl) match=true;
	break;
      case cstp::none:
	THROW(fatal_error,"Internal error");
      }
      if (match) {
	(*m_splitter)->AddSpec(p_spect);
	(*m_splitter)->SetSpec(p_spect);
	if (beam!=-1) (*m_splitter)->Lorentz()->SetBeam(beam);
	m_lastint += (*m_splitter)->OverIntegrated(zmin,zmax,scale,xbj);
	if (m_lastint>0. && m_lastint <0.) cout<<(*this);
      }
    }
  }
  return m_lastint;
}

double Sudakov::ProduceT(double t)
{
  double ne = 2.*M_PI/m_lastint;
  return t * exp(log(ran->Get())*Max(ne,1.0e-3));
}

bool Sudakov::Veto(double Q2,double x,double t,double y,double z) {
  return Splitting(Q2,x,t,y,z);
}

bool Sudakov::Splitting(double Q2,double x,double t,double y,double z) {
  double wt = RejectionWeight(z,y,x,t,Q2), efac = p_selected->EFac();
  if (ran->Get()>wt) {
    if (efac!=1.0) {
      m_weight*=(1.0-wt/efac)/(1.0-wt);
      p_split->Weights().push_back(std::make_pair(t,(1.0-wt/efac)/(1.0-wt)));
    }
    return false;
  }
  m_weight*=1.0/efac;
  return true;
}

const SF_E_Map *Sudakov::HasKernel(const ATOOLS::Flavour &fli,
				   const ATOOLS::Flavour &flj,
				   const cstp::code type) const
{
  const SF_EEE_Map *cmap(&m_sffmap);
  if (type==cstp::FI) cmap=&m_sfimap;
  else if (type==cstp::IF) cmap=&m_sifmap;
  else if (type==cstp::II) cmap=&m_siimap;
  SF_EEE_Map::const_iterator eees(cmap->find(fli));
  if (eees==cmap->end()) return NULL;
  SF_EE_Map::const_iterator ees(eees->second.find(flj));
  if (ees==eees->second.end()) return NULL;
  return &ees->second;
}

int Sudakov::HasKernel(const ATOOLS::Flavour &fli,
                       const ATOOLS::Flavour &flj,
                       const ATOOLS::Flavour &flk,
                       const cstp::code type) const
{
  // find whether a kernel for ij k -> i j k exists and which coupling type
  // it has, i.e. doesn't exist (=0), qcd (=1), ew (=2), qcd and ew (=3)
  // the latter can happen e.g. if  i=q  j=qbar  k=q'  where ij = {G | P}
  int cpl=0;
  const SF_E_Map * sfmap = HasKernel(fli, flj, type);
  if (sfmap==NULL) return 0;
  for (SF_E_Map::const_iterator es=sfmap->begin(); es!=sfmap->end(); ++es) {
    Splitting_Function_Base* sf = es->second;
    if (sf->Coupling()->AllowSpec(flk)) {
      if (sf->PureQCD()) cpl|=1;
      else cpl|=2;
    }
  }
  return cpl;
}

double Sudakov::CplFac(const ATOOLS::Flavour &fli,const ATOOLS::Flavour &flj,
                       const ATOOLS::Flavour &flk,const cstp::code type,
		       const int cpl,const double &mu2) const
{
  const SF_E_Map * sfmap = HasKernel(fli, flj, type);
  if (sfmap==NULL) return 0;
  for (SF_E_Map::const_iterator es=sfmap->begin(); es!=sfmap->end(); ++es) {
    Splitting_Function_Base* sf = es->second;
    if (sf->Coupling()->AllowSpec(flk)) {
      if (cpl==1 && sf->PureQCD()) return sf->Coupling()->CplFac(mu2);
      if (cpl==2 && !sf->PureQCD()) return sf->Coupling()->CplFac(mu2);
    }
  }
  return -1.0;
}
