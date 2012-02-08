#include "CSSHOWER++/Showers/Sudakov.H"

#include "CSSHOWER++/Showers/Splitting_Function_Base.H"
#include "CSSHOWER++/Tools/Singlet.H"
#include "CSSHOWER++/Showers/Shower.H"
#include "MODEL/Interaction_Models/Single_Vertex.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/My_Limits.H"

using namespace CSSHOWER;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;

bool CSSHOWER::Sudakov::s_init=false;

Sudakov::Sudakov(PDF::ISR_Handler *isr,const int qed) : 
  p_rms(NULL)
{
  m_ewmode=qed;
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
  Vertex_Table *vtab(md->GetVertexTable());
  for (Vertex_Table::const_iterator
	 vlit=vtab->begin();vlit!=vtab->end();++vlit) {
    for (Vertex_List::const_iterator 
	   vit=vlit->second.begin();vit!=vlit->second.end();++vit) {
      Single_Vertex *v(*vit);
      if (v->nleg>3 || !v->on) continue;
      if (sfs.find(FTrip(v->in[0],v->in[1],v->in[2]))!=sfs.end()) continue;
      sfs.insert(FTrip(v->in[0],v->in[1],v->in[2]));
      sfs.insert(FTrip(v->in[0],v->in[2],v->in[1]));
      msg_Debugging()<<"Add "<<v->in[0]<<" -> "<<v->in[1]<<" "<<v->in[2]<<" {\n";
      {
	msg_Indent();
	int dmode(0);
	if (v->in[2]==v->in[0]) dmode=1;
	else if (v->in[1]!=v->in[0] && 
		 v->in[1].IsAnti() && !v->in[2].IsAnti()) dmode=1;
	Add(new Splitting_Function_Base(SF_Key(p_rms,v,dmode,cstp::FF,kfmode,m_ewmode,1)));
	Add(new Splitting_Function_Base(SF_Key(p_rms,v,dmode,cstp::FF,kfmode,m_ewmode,-1)));
	Add(new Splitting_Function_Base(SF_Key(p_rms,v,dmode,cstp::FI,kfmode,m_ewmode,1)));
	Add(new Splitting_Function_Base(SF_Key(p_rms,v,dmode,cstp::FI,kfmode,m_ewmode,-1)));
	if (v->in[0].Mass()<100.0) {
  	  Add(new Splitting_Function_Base(SF_Key(p_rms,v,dmode,cstp::IF,kfmode,m_ewmode,1)));
  	  Add(new Splitting_Function_Base(SF_Key(p_rms,v,dmode,cstp::IF,kfmode,m_ewmode,-1)));
 	  Add(new Splitting_Function_Base(SF_Key(p_rms,v,dmode,cstp::II,kfmode,m_ewmode,1)));
 	  Add(new Splitting_Function_Base(SF_Key(p_rms,v,dmode,cstp::II,kfmode,m_ewmode,-1)));
	}
	if (v->in[1]!=v->in[2]) {
	  AddToMaps(new Splitting_Function_Base(SF_Key(p_rms,v,1-dmode,cstp::FF,kfmode,m_ewmode,1)));
	  AddToMaps(new Splitting_Function_Base(SF_Key(p_rms,v,1-dmode,cstp::FF,kfmode,m_ewmode,-1)));
	  AddToMaps(new Splitting_Function_Base(SF_Key(p_rms,v,1-dmode,cstp::FI,kfmode,m_ewmode,1)));
	  AddToMaps(new Splitting_Function_Base(SF_Key(p_rms,v,1-dmode,cstp::FI,kfmode,m_ewmode,-1)));
	  if (v->in[0].Mass()<100.0) {
  	    Add(new Splitting_Function_Base(SF_Key(p_rms,v,1-dmode,cstp::IF,kfmode,m_ewmode,1)));
  	    Add(new Splitting_Function_Base(SF_Key(p_rms,v,1-dmode,cstp::IF,kfmode,m_ewmode,-1)));
 	    Add(new Splitting_Function_Base(SF_Key(p_rms,v,1-dmode,cstp::II,kfmode,m_ewmode,1)));
 	    Add(new Splitting_Function_Base(SF_Key(p_rms,v,1-dmode,cstp::II,kfmode,m_ewmode,-1)));
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
  split->SetPDF(p_pdf);
  split->SetEFac(p_shower);
  if (mode) {
    m_addsplittings.push_back(split);
    msg_Debugging()<<"\n";
  }
  if (split->GetCol()<0) {
    switch(split->GetType()) {
    case cstp::IF:
      m_fifmap[split->GetFlavourA().Bar()]
	[split->GetFlavourC()]
	[split->GetFlavourB().Bar()]=split;
      break;
    case cstp::II:
      m_fiimap[split->GetFlavourA().Bar()]
	[split->GetFlavourC()]
	[split->GetFlavourB().Bar()]=split;
      break;
    default: break;
    }
    return;
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
    m_iffmap[split->GetFlavourA().Bar()]
      [split->GetFlavourC()]
      [split->GetFlavourB().Bar()]=split;
    if (split->On())
      m_sifmap[split->GetFlavourA()]
	[split->GetFlavourC()]
	[split->GetFlavourB()]=split;
    break;
  case cstp::II:
    m_ifimap[split->GetFlavourA().Bar()]
      [split->GetFlavourC()]
      [split->GetFlavourB().Bar()]=split;
    if (split->On())
      m_siimap[split->GetFlavourA()]
	[split->GetFlavourC()]
	[split->GetFlavourB()]=split;
    break;
  case cstp::none: break;
  }
}

bool Sudakov::Generate(Parton * split) 
{
  m_weight=1.0;
  ClearSpecs();
  ResetLastInt();
  m_cfl = split->GetFlavour();
  
  m_type     = cstp::none;
  std::vector<Parton*> slist;
  if (split->GetLeft()) slist.push_back(split->GetLeft());
  if (split->GetRight()) slist.push_back(split->GetRight());
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
  p_split=split;
  for (size_t i(0);i<slist.size();++i) {
    Parton *spect(slist[i]);
    p_spect=spect;
  double Q2(0.);
  int beam = -1;
  m_flspec = spect->GetFlavour();
  switch (split->GetType()) {
  case pst::FS: 
    if (spect->GetType()==pst::FS) {
      Q2    = (split->Momentum()+spect->Momentum()).Abs2();
      if (!DefineFFBoundaries(Q2,1.)) continue; 
      break;
    }
    if (spect->GetType()==pst::IS) {
      Q2       = -(split->Momentum()-spect->Momentum()).Abs2();
      beam     = spect->Beam();
      if (!DefineFIBoundaries(Q2,spect->Xbj(),beam)) continue;
      break;
    }
  case pst::IS:
    if (spect->GetType()==pst::FS) {
      Q2       = -(split->Momentum()-spect->Momentum()).Abs2();
      beam     = split->Beam();
      if (!DefineIFBoundaries(Q2,split->Xbj(),beam)) continue;
      break;
    }
    if (spect->GetType()==pst::IS) {
      Q2 = (split->Momentum()+spect->Momentum()).Abs2();
      beam = split->Beam();
      if (!DefineIIBoundaries(Q2,split->Xbj(),beam)) continue; 
      break;
    }
  case pst::none: 
    msg_Error()<<"Error in Sudakov::Generate : No pst-type for splitter. "<<endl
	       <<(*split)<<(*spect);
    abort();
  }
  if (m_type==cstp::none) {
    msg_Error()<<"Error in Sudakov::Generate : No type for splitter. "<<endl
	       <<(*split)<<(*spect);
    abort();
  }
  }  
  if (m_lastint<=0.0 || IsBad(m_lastint)) return false;
  double last=0.0;
  for (size_t i(0);i<m_splittings.size();++i)
    m_partint[i]=last+=m_splittings[i]->Last();
  if (!IsEqual(m_partint.back(),m_lastint))
    msg_Error()<<METHOD<<"(): Error, last = "<<m_lastint
	       <<" vs. "<<m_partint.back()<<"."<<std::endl;
  m_lastint=m_partint.back();

  m_kperp2       = split->KtStart();
  double x(0.); 
  
  bool success(false);
  while (m_kperp2>=m_k0sqf) {
    ProduceT();
    SelectOne();
    split->SetSpect(p_selected->SelectSpec());
    m_z = Z();
    double k0sq(p_split->GetType()==pst::IS?m_k0sqi:m_k0sqf);
    if (m_kperp2<k0sq)  return false;
    double Q2 = 0.;
    m_type=split->GetType()==pst::FS?
      (split->GetSpect()->GetType()==pst::FS?cstp::FF:cstp::FI):
      (split->GetSpect()->GetType()==pst::FS?cstp::IF:cstp::II);
    switch (m_type) {
    case (cstp::FF) : {
      double mi2 = sqr(p_rms->Mass(((*m_splitter)->GetFlavourB())));
      double mj2 = sqr(p_rms->Mass(((*m_splitter)->GetFlavourC())));
      double mk2 = sqr(p_rms->Mass(m_flspec));
      Q2 = (split->Momentum()+split->GetSpect()->Momentum()).Abs2();
      if (Q2<=mi2+mj2+mk2) return false;
      m_y = p_shower->KinFF()->GetY(Q2,m_kperp2,m_z,mi2,mj2,mk2);
      if (m_y<0.0) continue;
      x   = 0.;
    }    
      break; 
    case (cstp::FI) : {
      double mi2 = sqr(p_rms->Mass(((*m_splitter)->GetFlavourB())));
      double mj2 = sqr(p_rms->Mass(((*m_splitter)->GetFlavourC())));
      double ma2 = sqr(p_rms->Mass(m_flspec));
      Q2 = -(split->Momentum()-split->GetSpect()->Momentum()).Abs2();
      m_y = 1.0-p_shower->KinFI()->GetY(-Q2,m_kperp2,m_z,mi2,mj2,ma2);
      if (m_y>1.0) continue;
      x   = split->GetSpect()->Xbj();
    }
      break; 
    case (cstp::IF) : {
      double ma2 = sqr(p_rms->Mass(((*m_splitter)->GetFlavourA())));
      double mi2 = sqr(p_rms->Mass(((*m_splitter)->GetFlavourC())));
      double mk2 = sqr(p_rms->Mass(m_flspec));
      Q2 = -(split->Momentum()-split->GetSpect()->Momentum()).Abs2();
      m_y = p_shower->KinIF()->GetY(-Q2,m_kperp2,m_z,ma2,mi2,mk2);
      if (m_y<0.0) continue;
      x   = split->Xbj();
    }
      break;
    case (cstp::II) : {
      double ma2 = sqr(p_rms->Mass(((*m_splitter)->GetFlavourA())));
      double mi2 = sqr(p_rms->Mass(((*m_splitter)->GetFlavourC())));
      double mb2 = sqr(p_rms->Mass(m_flspec));
      Q2 = (split->Momentum()+split->GetSpect()->Momentum()).Abs2();
      m_y = p_shower->KinII()->GetY(Q2,m_kperp2,m_z,ma2,mi2,mb2);
      if (m_y<0.0) continue;
      x   = split->Xbj();
    }
      break;
  default:
      msg_Error()<<"Error in Sudakov::Generate!"<<std::endl;
      abort();
    }
    if (Veto(Q2,x)) { 
      success=true; 
      break; 
    } 
  }
  /*
  if (success) {
    std::cout<<"----------------------------------"<<std::endl;
    std::cout<<m_type<<" SUCCEED split "<<split->GetFlavour()
    <<"("<<(*m_splitter)->GetFlavourA()<<" -> "
    <<(*m_splitter)->GetFlavourB()<<"+"<<(*m_splitter)->GetFlavourC()<<") "
    <<" + spect : "<<spect->GetFlavour()
    <<" with kt2 = "<<m_kperp2<<", y = "<<m_y<<", z = "<<m_z<<endl;
    }
  */
  m_phi = 2.0*M_PI*ran->Get();
  return success;
}

bool Sudakov::DefineFFBoundaries(double Q2,double x)
{
  if (4.*m_k0sqf>Q2) return false;
  
  m_type=cstp::FF;
  double deltaz(sqrt(1.-4.*m_k0sqf/Q2));
  m_zmin   = 0.5*(1.-deltaz);
  m_zmax   = 0.5*(1.+deltaz);
  m_scale  = p_split->KtStart();
  if (OverIntegrated(m_zmin,m_zmax,m_scale,Q2)<0.) {
    msg_Error()<<"Error in Sudakov::DefineFFBoundaries : "<<endl
    	       <<"   Integral for SF's<0 : {"<<m_zmin<<","<<m_zmax<<","<<m_scale<<"}"<<endl;
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
  if (OverIntegrated(m_zmin,m_zmax,m_scale,x,beam)<0.) {
    msg_Error()<<"Error in Sudakov::DefineFIBoundaries : "<<endl
	       <<"   Integral for SF's<0 : {"<<m_zmin<<","<<m_zmax<<","<<m_scale<<"}"<<endl;
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
  if (OverIntegrated(m_zmin,m_zmax,m_scale,x,beam)<0.) {
    msg_Error()<<"Error in Sudakov::DefineIFBoundaries : "<<endl
	       <<"   Integral for SF's<0 : {"<<m_zmin<<","<<m_zmax<<","<<m_scale<<"}"<<endl;
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
  if (OverIntegrated(m_zmin,m_zmax,m_scale,x,beam)<0.) {
    msg_Error()<<"Error in Sudakov::DefineIIBoundaries : "<<endl
	       <<"   Integral for SF's<0 : {"<<m_zmin<<","<<m_zmax<<","<<m_scale<<"}"<<endl;
    return false;
  }
  return true;
}


double Sudakov::OverIntegrated(const double zmin,const double zmax,
			       const double scale,const double xbj,int beam) {
  for (m_splitter=m_splittings.begin();m_splitter!=m_splittings.end();m_splitter++) {
    if ((*m_splitter)->GetType()==m_type && 
	(*m_splitter)->Coupling()->AllowSpec(m_flspec)) {
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

void Sudakov::ProduceT()
{
  double ne=2.*M_PI/m_lastint;
  m_kperp2 *= exp(log(ran->Get())*Max(ne,1.0e-3));
}

bool Sudakov::Veto(double Q2,double x) {
  if (!KinCheck(Q2,x)) return false;
  if (!Splitting(Q2,x))          return false;
  return true;
}

bool Sudakov::KinCheck(double Q2,double x) {
  switch (m_type) {
  case cstp::FF: {
    if (m_y<0. || m_y>1.)      return false;
    double mui2 = sqr(p_rms->Mass((*m_splitter)->GetFlavourB()))/Q2;
    double muj2 = sqr(p_rms->Mass((*m_splitter)->GetFlavourC()))/Q2;;
    double muk2 = sqr(p_rms->Mass((*m_splitter)->GetFlavourSpec()))/Q2; 
    double muk  = sqrt(muk2);
    if (mui2!=0. || muj2!=0. || muk2!=0.) {
      //massive kinematics
      double fac  = 1.-mui2-muj2-muk2;
      double ym   = 2.*sqrt(mui2*muj2)/fac;
      double yp   = 1.- 2.*muk*(1.-muk)/fac;
      
      if (m_y<ym || m_y>yp) return false;
      
      double viji  = sqrt(sqr(fac*m_y)-4.*mui2*muj2)/(fac*m_y+2.*mui2);
      double vijk  = sqrt(sqr(2.*muk2+fac*(1.-m_y))-4.*muk2)/(fac*(1.-m_y));
      double frac  = (2.*mui2+fac*m_y)/(2.*(mui2+muj2+fac*m_y));
      double delta = frac*viji*vijk;
      double zm    = frac - delta;  
      double zp    = frac + delta;
      if (m_z<zm || m_z>zp) return false;
    }
  }
    break;
  case cstp::FI: {
    if (m_y<0. || m_y>1.-x) return false;
    double Qt2 = (Q2-sqr(p_rms->Mass((*m_splitter)->GetFlavourA())))/(1.-m_y);
    double muij2 = sqr(p_rms->Mass((*m_splitter)->GetFlavourA()))/Qt2;
    double mui2  = sqr(p_rms->Mass((*m_splitter)->GetFlavourB()))/Qt2;
    double muj2  = sqr(p_rms->Mass((*m_splitter)->GetFlavourC()))/Qt2;;
    if (m_y < sqr(mui2+muj2)-muij2) return false;
    if (muij2!=0. || mui2!=0. || muj2!=0.) {
      //massive kinematics
      double div   = m_y+muij2;
      double term1 = div+mui2-muj2; 
      double term2 = sqrt(sqr(div-mui2-muj2)-4.*mui2*muj2); 
      double zp = (term1 + term2)/(2.*div);
      double zm = (term1 - term2)/(2.*div);
      if (m_z<zm || m_z>zp) return false;
    }
  }
    break;
  case cstp::IF: { 
    if (m_y<0. || m_y>1.)      return false;
    if (x>m_z  || 0.>x)        return false;
    double mk2 = sqr(p_rms->Mass(m_flspec));
    if (mk2!=0) {
      //massive spectator
      double muk2 = mk2*m_z/(Q2+mk2);
      double yp = (1.-m_z)/(1.-m_z+muk2);
      if (m_y>yp)               return false;
    }
  }
    break;
  case cstp::II: 
    if (m_y<0.       || m_y>1.)  return false;
    if (m_y>(1.-m_z) || x>m_z)   return false;
    break;
  case cstp::none:             return false;
  }
  return true;
}

bool Sudakov::Splitting(double Q2,double x) {
  if (m_kperp2>p_split->KtPrev()) return false;
  double wt(RejectionWeight(m_z,m_y,x,m_kperp2,Q2));
  double efac=p_selected->EFac();
  if (ran->Get()>wt) {
    if (efac!=1.0) {
      m_weight*=(1.0-wt/efac)/(1.0-wt);
      p_split->Weights().push_back
	(std::make_pair(m_kperp2,(1.0-wt/efac)/(1.0-wt)));
    }
    return false;  
  }
  else {
    m_weight*=1.0/efac;
    return true;
  }
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
