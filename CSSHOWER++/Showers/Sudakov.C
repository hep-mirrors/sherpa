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
			  const double &isfac,const double &fsfac)
{
  m_k0sqi=k0sqi;
  m_k0sqf=k0sqf;
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
  if (split->GetCol()<0) {
    switch(split->GetType()) {
    case cstp::IF:
      m_fifmap[split->GetFlavourA()]
	[split->GetFlavourC()]
	[split->GetFlavourB()]=split;
      break;
    case cstp::II:
      m_fiimap[split->GetFlavourA()]
	[split->GetFlavourC()]
	[split->GetFlavourB()]=split;
      break;
    default: break;
    }
    if (split->GetCol()<0) return;
  }
  switch(split->GetType()) {
  case cstp::FF:
    m_fffmap[split->GetFlavourB()]
      [split->GetFlavourC()]
      [split->GetFlavourA()]=split;
    if (split->On())
      m_sffmap[split->GetFlavourB()]
	[split->GetFlavourC()]
	[split->GetFlavourA()]=split;
    break;
  case cstp::FI:
    m_ffimap[split->GetFlavourB()]
      [split->GetFlavourC()]
      [split->GetFlavourA()]=split;
    if (split->On())
      m_sfimap[split->GetFlavourB()]
	[split->GetFlavourC()]
	[split->GetFlavourA()]=split;
    break;
  case cstp::IF:
    m_iffmap[split->GetFlavourA()]
      [split->GetFlavourC()]
      [split->GetFlavourB()]=split;
    if (split->On())
      m_sifmap[split->GetFlavourA()]
	[split->GetFlavourC()]
	[split->GetFlavourB()]=split;
    break;
  case cstp::II:
    m_ifimap[split->GetFlavourA()]
      [split->GetFlavourC()]
      [split->GetFlavourB()]=split;
    if (split->On())
      m_siimap[split->GetFlavourA()]
	[split->GetFlavourC()]
	[split->GetFlavourB()]=split;
    break;
  case cstp::none: break;
  }
}

Parton* Sudakov::SelectKinSpect(Parton *split, Parton *spect){
  std::list<Parton *>::const_iterator it_w = p_shower->GetSinglet()[0].begin();
  if      (split->GetFlavour().IsbQuark())    std::advance(it_w, 2); // move to W+
  else if (split->GetFlavour().IsbbarQuark()) std::advance(it_w, 3); // move to W-
  else    THROW(fatal_error, "Invalid assigment.");
  return *it_w;
}

void Sudakov::SetTrialVariables(Parton *split, Parton *spect, Parton * kinspect){
  const double mw2              = sqr(Flavour(24).Mass());
  const ATOOLS::Vec4D pai_tilde = split->Momentum();
  const ATOOLS::Vec4D pb        = spect->Momentum();
  const ATOOLS::Vec4D pw_tilde  = kinspect->Momentum();
  m_trialvariables.m_Q2         = (pai_tilde + pw_tilde).Abs2();
  const ATOOLS::Vec4D p_minus   = pw_tilde-mw2/(m_trialvariables.m_Q2-mw2)*pai_tilde;
  m_trialvariables.m_paipb      = pai_tilde*pb;
  m_trialvariables.m_alpha      = pb*p_minus / (pai_tilde*p_minus);
}

double Sudakov::GetViMax(Parton *split, const double &z) const {
  const Vec4D pai_tilde = split->Momentum();
  const Vec4D pw_tilde  = split->GetKinSpect()->Momentum();
  const Vec4D pa        = z*pai_tilde;
  const Vec4D n         = pw_tilde - pa + pai_tilde;
  return 2.*pa*n/(n*n)*(1.-z)/z;
}

double Sudakov::KT2(const double &vi, const double &z, const double &phi) const {
  const double mw2     = sqr(Flavour(24).Mass());
  const double alpha   = m_trialvariables.m_alpha;
  const double Qprime2 = m_trialvariables.m_Q2-mw2;
  const double paipb   = m_trialvariables.m_paipb;
  const double A       = alpha*vi*Qprime2/2.;
  const double B       = paipb*(vi*(z-1.-mw2/Qprime2)+1-z);
  const double cosphi  = cos(phi);
  const double kt2     = vi*Qprime2*(A+B-2.*sqrt(A*B)*cosphi)/(2.*paipb);
  return kt2;
}

double Sudakov::GetVi(const double &z, const double &phi, const double &t,
                      const double &vimax) const {
  const double Q2    = m_trialvariables.m_Q2;
  const double mw2   = sqr(Flavour(24).Mass());
  if(p_shower->KinFF()->m_evolscheme == 2) return t/(z/2.*(Q2-mw2));
  else if(p_shower->KinFF()->m_evolscheme == 1){
    const double paipb = m_trialvariables.m_paipb;
    const double alpha = m_trialvariables.m_alpha;

    p_shower->KinFF()->m_calcV.calculate_helpers(Q2,t,z,paipb,alpha,sqr(cos(phi)));

    /* select any vi randomly */
    const double r = ran->Get();
    double vi;

    /* first check */
    if     (r<=0.25)   vi = p_shower->KinFF()->m_calcV.GetV1();
    else if(r<=0.50)   vi = p_shower->KinFF()->m_calcV.GetV2();
    else if(r<=0.75)   vi = p_shower->KinFF()->m_calcV.GetV3();
    else if(r<=1.00)   vi = p_shower->KinFF()->m_calcV.GetV4();

    if(!IsBad(vi) && (vi>0. && vi<vimax) && IsEqual(KT2(vi,z,phi),t,1e-2)) return vi;
    return -1.;
  }
  else THROW(fatal_error, "Invalid evolution scheme for res-aware kinematics!");
}

double Sudakov::GetViab(Parton *split, const double &z, const double &phi, const double &t,
                        const double &vi) const {
  Vec4D pai_tilde = split->Momentum();
  Vec4D pw_tilde  = split->GetKinSpect()->Momentum();
  Vec4D pb        = split->GetSpect()->Momentum();
  const double mw2      = sqr(Flavour(24).Mass());
  const double y_wia    = 1.-z; 
  const double ztilde_w = 1-vi;
  PHASIC::Kin_Args ff(y_wia,ztilde_w,phi);

  /* boost into pai+p_ rest-frame in order to be able to use phi_ib and construct additional
     momentum */
        ATOOLS::Vec4D pminus = pw_tilde - mw2/(m_trialvariables.m_Q2-mw2)*pai_tilde;
  const ATOOLS::Vec4D pboost = pai_tilde+pminus;
  Poincare bst(pboost);
  bst.Boost(pai_tilde);
  bst.Boost(pw_tilde);
  bst.Boost(pminus);
  bst.Boost(pb);
  ff.m_res    = true;
  ff.m_pb     = pb;
  ff.m_pminus = pminus;

  if (PHASIC::ConstructFFDipole(mw2,0,mw2,0,pw_tilde,pai_tilde,ff)<0) return -1.;
  Vec4D pi = ff.m_pj;
  Vec4D pa = ff.m_pk;

  if(p_shower->KinFF()->m_evolscheme == 1){
    const double kt2 = (pa*pi)*(pb*pi)/(pa*pb);
    if(IsEqual(t, kt2, 1e-3)) return pa*pb / (pi*(pa+pb));
  }
  else if(p_shower->KinFF()->m_evolscheme == 2){
    const double virtuality = pa*pi;
    if(IsEqual(t, virtuality, 1e-3)) return pa*pb / (pi*(pa+pb));
  }
  return -1.;
}

bool Sudakov::Generate(Parton * split) 
{
  m_weight=1.0;
  int cc(split->GetFlavour().StrongCharge());
  if (((cc==8 || (split->GetType()==pst::FS?cc:-cc)==3) &&
       split->GetLeft()==NULL) ||
      ((cc==8 || (split->GetType()==pst::FS?cc:-cc)==-3) &&
       split->GetRight()==NULL)) {
    msg_Out()<<METHOD<<":\n"<<(*split)<<".\n";
    THROW(fatal_error,"Invalid color flow.");
  }
  m_cfl  = split->GetFlavour();
  m_type = cstp::none;
  std::vector<Parton*> slist;
  msg_IODebugging()<<"---- "<<METHOD<<":\n"
	   <<"   Check spectators for [type = "<<split->GetType()<<"]"
	   <<" ("<<split->GetFlow(1)<<", "<<split->GetFlow(2)<<").\n";
  if (split->GetLeft() && split!=split->GetLeft())  {
    slist.push_back(split->GetLeft());
    msg_IODebugging()<<"   --> add left: "<<split->GetLeft()->GetType()
	     <<" ("<<split->GetLeft()->GetFlow(1)<<", "
	     <<split->GetLeft()->GetFlow(2)<<").\n";
  }
  if (split->GetRight() && split!=split->GetRight()) {
    slist.push_back(split->GetRight());
    msg_IODebugging()<<"   --> add right: "<<split->GetRight()->GetType()
	     <<" ("<<split->GetRight()->GetFlow(1)<<", "
	     <<split->GetRight()->GetFlow(2)<<").\n";
  }
  msg_IODebugging()<<"   ===> found "<<slist.size()<<" spectator(s).\n";
  int sc=split->GetFlavour().IntCharge();
  if (split->GetType()==pst::IS) sc=-sc;
  if (sc!=0 || split->GetFlavour().IsPhoton() ||
      split->GetFlavour().Kfcode()==kf_Z) {
    Singlet *sing(split->GetSing());
    for (PLiter pit(sing->begin());pit!=sing->end();++pit) {
      int cc=(*pit)->GetFlavour().IntCharge();
      if ((*pit)->GetType()==pst::IS) cc=-cc;
      if (*pit!=split->GetLeft() && *pit!=split->GetRight() &&
	  cc!=0 && (sc==0 || sc*cc<0)) 
	slist.push_back(*pit);
    }
  }
  double t0(m_k0sqf), t, y, z, phi;
  Parton *spect(NULL);
  Splitting_Function_Base *selected(NULL);
  for (size_t i(0);i<slist.size();++i) {
    int success(Generate(split,slist[i],t0,t,y,z,phi));
    if (success) {
      msg_IODebugging()<<"shrink evolution window "<<t0<<" -> "<<t<<"\n";
      m_sy=y;
      m_sz=z;
      m_sphi=phi;
      m_st=t0=t;
      spect=slist[i];
      selected=p_selected;
    }
  }
  p_spect=NULL;
  p_selected=NULL;
  if (spect) {
    p_spect=spect;
    p_split->SetSpect(p_spect);
    p_selected=selected;
    msg_IODebugging()<<"selected "<<p_spect<<", t = "<<m_st
		     <<", y = "<<m_sy<<", z = "<<m_sz<<", phi = "<<m_sphi<<"\n";
  }
  ClearSpecs();
  ResetLastInt();
  return p_spect!=NULL;
}

bool Sudakov::Generate(Parton *split,Parton *spect,double t0,double &t,double &y,double &z,double &phi)
{
  ClearSpecs();
  ResetLastInt();
  p_split=split;
  p_spect=spect;
  double Q2(0.);
  int beam = -1;
  m_flspec = spect->GetFlavour();
  switch (split->GetType()) {
  case pst::FS: 
    if (spect->GetType()==pst::FS) {
      switch(p_shower->KinFF()->m_dipole_case){
        case EXTAMP::CS:
          m_trialvariables.m_Q2 = (split->Momentum()+spect->Momentum()).Abs2();
          break;
        case EXTAMP::IDa:{
          if(!(split->GetFlavour().IsbQuark() || split->GetFlavour().IsbbarQuark())) return false;
          SetTrialVariables(split, spect, SelectKinSpect(split, spect));
          break;}
        default:
          THROW(fatal_error, "Unsupported dipole-case for matching.");
      }
      if (!DefineFFBoundaries(1.)) return false;
      break;
    }
    if (spect->GetType()==pst::IS) {
      Q2       = -(split->Momentum()-spect->Momentum()).Abs2();
      beam     = spect->Beam();
      if (!DefineFIBoundaries(Q2,spect->Xbj(),beam)) return false;
      break;
    }
  case pst::IS:
    if (spect->GetType()==pst::FS) {
      Q2       = -(split->Momentum()-spect->Momentum()).Abs2();
      beam     = split->Beam();
      if (!DefineIFBoundaries(Q2,split->Xbj(),beam)) return false;
      break;
    }
    if (spect->GetType()==pst::IS) {
      Q2 = (split->Momentum()+spect->Momentum()).Abs2();
      beam = split->Beam();
      if (!DefineIIBoundaries(Q2,split->Xbj(),beam)) return false;
      break;
    }
  case pst::none: 
    msg_Error()<<"Error in Sudakov::Generate : No pst-type for splitter.\n"
	       <<(*split);
    return false;
    Abort();
  }
  if (m_type==cstp::none) {
    msg_Error()<<"Error in Sudakov::Generate : No type for splitter.\n"
	       <<(*split);
    return false;
    Abort();
  }
  if (m_lastint<=0.0 || IsBad(m_lastint)) return false;
  double last=0.0;
  for (size_t i(0);i<m_splittings.size();++i)
    m_partint[i]=last+=m_splittings[i]->Last();
  if (!IsEqual(m_partint.back(),m_lastint))
    msg_Error()<<METHOD<<"(): Error, last = "<<m_lastint
	       <<" vs. "<<m_partint.back()<<"."<<std::endl;
  m_lastint=m_partint.back();

  t = split->KtStart();
  msg_IODebugging()<<"starting scale "<<t<<", cutoff scale "<<t0<<"\n";
  if (p_spect==p_split->GetLeft() && t>p_split->KtSoft(0)) {
    msg_IODebugging()<<"reset starting scale "<<t<<" -> "<<split->KtSoft(0)<<"\n";
    t=split->KtSoft(0);
  }
  if (p_spect==p_split->GetRight() && t>p_split->KtSoft(1)) {
    msg_IODebugging()<<"reset starting scale "<<t<<" -> "<<split->KtSoft(1)<<"\n";
    t=split->KtSoft(1);
  }
  double x = 0.0;
  
  bool success(false);
  while (t>=t0) {
    t=ProduceT(t);
    SelectOne();
    split->SetSpect(p_spect=p_selected->SelectSpec());
    if(p_shower->KinFF()->m_dipole_case==EXTAMP::IDa)
      split->SetKinSpect(SelectKinSpect(split,split->GetSpect()));
    m_flspec = p_spect->GetFlavour();
    z = Z();
    phi = Phi();
    double k0sq(p_split->GetType()==pst::IS?m_k0sqi:m_k0sqf);
    if (t<Max(t0,k0sq))  return false;
    double Q2 = 0.;
    m_type=split->GetType()==pst::FS?
      (split->GetSpect()->GetType()==pst::FS?cstp::FF:cstp::FI):
      (split->GetSpect()->GetType()==pst::FS?cstp::IF:cstp::II);
    const Mass_Selector *ms(p_selected->Lorentz()->MS());
    switch (m_type) {
    case (cstp::FF) : {
      switch(p_shower->KinFF()->m_dipole_case){
        case EXTAMP::CS:
        {
        double mi2 = sqr(ms->Mass(((*m_splitter)->GetFlavourB())));
        double mj2 = sqr(ms->Mass(((*m_splitter)->GetFlavourC())));
        double mk2 = sqr(ms->Mass(m_flspec));
        Q2 = (split->Momentum()+split->GetSpect()->Momentum()).Abs2();
        if (Q2<=mi2+mj2+mk2) return false;
        y = p_shower->KinFF()->GetY(Q2,t,z,mi2,mj2,mk2,
				    (*m_splitter)->GetFlavourA(),
				    (*m_splitter)->GetFlavourC());
        x = 0.;
        if (y<0.0 || y>1.0) continue;
        break;
        }
        case EXTAMP::IDa:
        {
        Q2 = m_trialvariables.m_Q2;
        const double vimax = GetViMax(split,z);
        SetTrialVariables(split, split->GetSpect(), split->GetKinSpect());
        y      = GetVi(z,phi,t,vimax);         if(y < 0.)   continue; // watch out y=^vi
        m_viab = GetViab(split,z,phi,t,y);     if(m_viab < 0.) continue;
        break;
        }
      }
    }
      break; 
    case (cstp::FI) : {
      double mi2 = sqr(ms->Mass(((*m_splitter)->GetFlavourB())));
      double mj2 = sqr(ms->Mass(((*m_splitter)->GetFlavourC())));
      double ma2 = sqr(ms->Mass(m_flspec));
      double mij2= sqr(ms->Mass(((*m_splitter)->GetFlavourA())));
      Q2 = -(split->Momentum()-split->GetSpect()->Momentum()).Abs2();
      y = p_shower->KinFI()->GetY(-Q2,t,z,mi2,mj2,ma2,
				    (*m_splitter)->GetFlavourA(),
				    (*m_splitter)->GetFlavourC());
      y = 1.0-y*(-Q2-mij2-ma2)/(-Q2-mi2-mj2-ma2);
      x = split->GetSpect()->Xbj();
      if (y<0.0 || y>1.0-x) continue;
    }
      break; 
    case (cstp::IF) : {
      double ma2 = sqr(ms->Mass(((*m_splitter)->GetFlavourA())));
      double mi2 = sqr(ms->Mass(((*m_splitter)->GetFlavourC())));
      double mk2 = sqr(ms->Mass(m_flspec));
      Q2 = -(split->Momentum()-split->GetSpect()->Momentum()).Abs2();
      y = p_shower->KinIF()->GetY(-Q2,t,z,ma2,mi2,mk2,
				    (*m_splitter)->GetFlavourB(),
				    (*m_splitter)->GetFlavourC());
      x = split->Xbj();
      if (y<0.0 || y>1.0 || z<x) continue;
    }
      break;
    case (cstp::II) : {
      double ma2 = sqr(ms->Mass(((*m_splitter)->GetFlavourA())));
      double mi2 = sqr(ms->Mass(((*m_splitter)->GetFlavourC())));
      double mb2 = sqr(ms->Mass(m_flspec));
      Q2 = (split->Momentum()+split->GetSpect()->Momentum()).Abs2();
      y = p_shower->KinII()->GetY(Q2,t,z,ma2,mi2,mb2,
				    (*m_splitter)->GetFlavourB(),
				    (*m_splitter)->GetFlavourC());
      x   = split->Xbj();
      if (y<0.0 || y>1.0-z || z<x) continue;
    }
      break;
  default:
      msg_Error()<<"Error in Sudakov::Generate!"<<std::endl;
      Abort();
    }
    double kt2=0.;
    if(p_shower->KinFF()->m_dipole_case==EXTAMP::IDa)
      kt2 = GetKt2(split,z,phi,y);
    const bool veto(Veto(Q2, x,t,y,z,phi,kt2));
    if (m_keeprewinfo) {
      const double accwgt(Selected()->LastAcceptanceWeight());
      const double lastscale(Selected()->LastScale());
      if (accwgt < 1.0 && accwgt > 0.0 && lastscale > m_reweightscalecutoff) {
        Sudakov_Reweighting_Info info;
        info.accepted = veto;
        info.scale = lastscale;
        info.accwgt = accwgt;
        info.lastj = Selected()->Lorentz()->LastJ();
        info.lastcpl = Selected()->Coupling()->Last();
        info.sf = Selected();
        info.x = x;
        info.y = y;
        info.z = z;
        info.flspec = Selected()->Lorentz()->FlSpec();
        p_split->SudakovReweightingInfos().push_back(info);
      }
    }
    if (veto) {
      success = true;
      break;
    }
  }
  msg_IODebugging()<<"trial "<<p_spect<<", t = "<<t
		   <<", y = "<<y<<", z = "<<z<<", phi = "<<phi<<"\n";
  return success;
}

bool Sudakov::DefineFFBoundaries(double x)
{
  const double Q2 = m_trialvariables.m_Q2;
  if (4.*m_k0sqf>Q2) return false;
  
  m_type=cstp::FF;
  switch(p_shower->KinFF()->m_dipole_case){
    case EXTAMP::CS:
      {
      double deltaz(sqrt(1.-4.*m_k0sqf/Q2));
      m_zmin   = 0.5*(1.-deltaz);
      m_zmax   = 0.5*(1.+deltaz);
      break;
      }
    case EXTAMP::IDa:
      {
      const double     mw2 = sqr(Flavour(24).Mass());
      const double Qprime2 = Q2-mw2;
      m_zmin = 0;
      if(p_shower->KinFF()->m_evolscheme == 1){
        const double alpha   = m_trialvariables.m_alpha;
        const double paipb   = m_trialvariables.m_paipb;
        m_zmax = (-2*m_k0sqf*paipb*Qprime2 + paipb*pow(Qprime2,2) - alpha*pow(Qprime2,3) +
                  sqrt(alpha*pow(Qprime2,4)*(4*m_k0sqf*paipb +
                  alpha*pow(Qprime2,2))))/(paipb*pow(Qprime2,2));
      }
      else if(p_shower->KinFF()->m_evolscheme == 2){
        const double radicand = sqr(1.+m_k0sqf/Qprime2)/4.-m_k0sqf/Qprime2*(1+mw2/Qprime2);
        if(radicand >= 0.) m_zmax = (1.+m_k0sqf/Qprime2)/2 + sqrt(radicand);
        else               m_zmax = (1.+m_k0sqf/Qprime2)/2;
        /* following line is rarely relevant: only if Qprime2 is extremely small */
        if(m_zmax>1.) m_zmax = 0.99999;        // TODO
      }
      /* in case zmax is slightly above 1, due to numerics */
      if(m_zmax>1. && IsEqual(m_zmax,1.,1.e-8))           m_zmax=0.99999;
      if(!(m_zmax>0.) || !(m_zmax<=1.))    THROW(fatal_error, "zmax wrong");
      }
  }
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
  
  m_type=cstp::FI;
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
  
  m_type=cstp::IF;
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
  double ne=2.*M_PI/m_lastint;
  return t * exp(log(ran->Get())*Max(ne,1.0e-3));
}

bool Sudakov::Veto(double Q2,double x,double t,double y,double z, double phi, double kt2) {
  if (!Splitting(Q2,x,t,y,z,phi,kt2))          return false;
  return true;
}

bool Sudakov::Splitting(double Q2,double x,double t,double y,double z, double phi, double kt2) {
  double wt;
  switch(p_shower->KinFF()->m_dipole_case){
    case EXTAMP::CS:
      wt = RejectionWeight(z,y,x,t,Q2);
      break;
    case EXTAMP::IDa:
      wt = RejectionWeight(z,m_viab,x,kt2,Q2,phi,
      y,m_trialvariables.m_alpha,m_trialvariables.m_paipb);
      break;
  }

  double efac=p_selected->EFac();
  if (ran->Get()>wt) {
    if (efac!=1.0) {
      m_weight*=(1.0-wt/efac)/(1.0-wt);
      p_split->Weights().push_back
	(std::make_pair(t,(1.0-wt/efac)/(1.0-wt)));
    }
    return false;  
  }
  else {
    m_weight*=1.0/efac;
  }
  return true;
}

double Sudakov::GetKt2(Parton * split, const double &z, const double &phi, const double &vi) const {
  Vec4D pai_tilde = split->Momentum();
  Vec4D pw_tilde  = split->GetKinSpect()->Momentum();
  Vec4D pb        = split->GetSpect()->Momentum();
  const double mw2      = sqr(Flavour(24).Mass());
  const double y_wia    = 1.-z;
  const double ztilde_w = 1-vi;
  PHASIC::Kin_Args ff(y_wia,ztilde_w,phi);

  /* boost into pai+p_ rest-frame in order to be able to use phi_ib and construct additional
     momentum */
        ATOOLS::Vec4D pminus = pw_tilde - mw2/(m_trialvariables.m_Q2-mw2)*pai_tilde;
  const ATOOLS::Vec4D pboost = pai_tilde+pminus;
  Poincare bst(pboost);
  bst.Boost(pai_tilde);
  bst.Boost(pw_tilde);
  bst.Boost(pminus);
  bst.Boost(pb);
  ff.m_res    = true;
  ff.m_pb     = pb;
  ff.m_pminus = pminus;

  if (PHASIC::ConstructFFDipole(mw2,0,mw2,0,pw_tilde,pai_tilde,ff)<0)
    THROW(fatal_error, "Must not happend!")

  Vec4D pi = ff.m_pj;
  Vec4D pa = ff.m_pk;

  const Vec4D p_aib = pa+pi+pb;
  Poincare bst_aib(p_aib);
  bst_aib.Boost(pa);
  bst_aib.Boost(pi);

  const double theta = pa.Theta(pi);
  const double pi_perp = Vec3D(pi).Abs()*sin(theta);
  return sqr(pi_perp);
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
