#include "CSSHOWER++/Showers/Sudakov.H"

#include "CSSHOWER++/Showers/Splitting_Function_Base.H"
#include "CSSHOWER++/Tools/Singlet.H"
#include "CSSHOWER++/Showers/Shower.H"
#include "MODEL/Interaction_Models/Single_Vertex.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Math/Random.H"

using namespace CSSHOWER;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;

Sudakov::Sudakov(PDF::ISR_Handler *isr,const int qed) : 
  p_rms(NULL)
{
  //int hadron = rpa.gen.Beam1().Strong()==1?0:1;
  for (int i=0;i<2; i++) p_pdf[i] = isr->PDF(i);

}

Sudakov::~Sudakov() 
{
  for (size_t i(0);i<m_addsplittings.size();++i) delete m_addsplittings[i];
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
    if (m_a>f.m_a) return false;
    if (m_b<f.m_b) return true;
    if (m_b>f.m_b) return false;
    if (m_c<f.m_c) return true;
    if (m_c>f.m_c) return false;
    return false;
  }
};

void Sudakov::InitSplittingFunctions(MODEL::Model_Base *md)
{
  SFC_Filler_Getter::Getter_List flist(SFC_Filler_Getter::GetGetters());
  for (SFC_Filler_Getter::Getter_List::const_iterator git(flist.begin());
       git!=flist.end();++git) (*git)->GetObject(md);
  if (msg_LevelIsDebugging()) {
    msg_Out()<<METHOD<<"(): {\n\n"
	     <<"   // available coupling calcs\n\n";
    SFC_Getter::PrintGetterInfo(msg->Out(),25);
    msg_Out()<<"\n   // available lorentz calcs\n\n";
    SFL_Getter::PrintGetterInfo(msg->Out(),25);
    msg_Out()<<"\n}"<<std::endl;
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
	Add(new Splitting_Function_Base(SF_Key(p_rms,v,dmode,cstp::FF)));
	Add(new Splitting_Function_Base(SF_Key(p_rms,v,dmode,cstp::FI)));
	if (v->in[0].Mass()<100.0) {
  	  Add(new Splitting_Function_Base(SF_Key(p_rms,v,dmode,cstp::IF)));
 	  Add(new Splitting_Function_Base(SF_Key(p_rms,v,dmode,cstp::II)));
	}
	if (v->in[1]!=v->in[2]) {
	  AddToMaps(new Splitting_Function_Base(SF_Key(p_rms,v,1-dmode,cstp::FF)));
	  AddToMaps(new Splitting_Function_Base(SF_Key(p_rms,v,1-dmode,cstp::FI)));
	  if (v->in[0].Mass()<100.0) {
  	    Add(new Splitting_Function_Base(SF_Key(p_rms,v,1-dmode,cstp::IF)));
 	    Add(new Splitting_Function_Base(SF_Key(p_rms,v,1-dmode,cstp::II)));
	  }
	}
      }
      msg_Debugging()<<"}\n";
    }
  }
  msg_Debugging()<<"}\n";
}

void Sudakov::SetCoupling(MODEL::Model_Base *md,const double &k0sq,
			  const double &isfac,const double &fsfac)
{
  m_k0sq=k0sq;
  m_as_is_fac=isfac;
  m_as_fs_fac=fsfac;
  for (std::vector<Splitting_Function_Base*>::iterator
	 sit(m_splittings.begin());sit!=m_splittings.end();)
    if (!(*sit)->Coupling()->SetCoupling(md,m_k0sq,m_as_is_fac,m_as_fs_fac)) {
      delete *sit;
      sit=m_splittings.erase(sit);
    }
    else {
      ++sit;
    }
  for (std::vector<Splitting_Function_Base*>::iterator
	 sit(m_addsplittings.begin());sit!=m_addsplittings.end();)
    if (!(*sit)->Coupling()->SetCoupling(md,m_k0sq,m_as_is_fac,m_as_fs_fac)) {
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
  if (mode) {
    m_addsplittings.push_back(split);
    msg_Debugging()<<"\n";
  }
  switch(split->GetType()) {
  case cstp::FF:
    m_ffmap[split->GetFlavourB()]
      [split->GetFlavourC()]
      [split->GetFlavourA()]=split;
    break;
  case cstp::FI:
    m_fimap[split->GetFlavourB()]
      [split->GetFlavourC()]
      [split->GetFlavourA()]=split;
    break;
  case cstp::IF:
    m_ifmap[split->GetFlavourA().Bar()]
      [split->GetFlavourC()]
      [split->GetFlavourB().Bar()]=split;
    break;
  case cstp::II:
    m_iimap[split->GetFlavourA().Bar()]
      [split->GetFlavourC()]
      [split->GetFlavourB().Bar()]=split;
    break;
  case cstp::none: break;
  }
}

bool Sudakov::Dice(Parton * split) 
{
  ClearSpecs();
  ResetLastInt();
  m_cfl = split->GetFlavour();
  if (m_cfl.StrongCharge()==8) m_nspect = 2.;
                               else m_nspect = 1.;  
  
  m_type     = cstp::none;
  std::vector<Parton*> slist;
  if (split->GetLeft()) slist.push_back(split->GetLeft());
  if (split->GetRight()) slist.push_back(split->GetRight());
  int sc=split->GetFlavour().IntCharge();
  if (sc!=0 || split->GetFlavour().IsPhoton() ||
      split->GetFlavour().Kfcode()==kf_Z) {
    Singlet *sing(split->GetSing());
    for (PLiter pit(sing->begin());pit!=sing->end();++pit)
      if (*pit!=split->GetLeft() && *pit!=split->GetRight() &&
	  (*pit)->GetFlavour().IntCharge()!=0 &&
	  (sc==0 || sc*(*pit)->GetFlavour().IntCharge()<0)) 
	slist.push_back(*pit);
  }
  p_split=split;
  for (size_t i(0);i<slist.size();++i) {
    Parton *spect(slist[i]);
    p_spec=spect;
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
    msg_Error()<<"Error in Sudakov::Dice : No pst-type for splitter. "<<endl
	       <<(*split)<<(*spect);
    abort();
  }
  if (m_type==cstp::none) {
    msg_Error()<<"Error in Sudakov::Dice : No type for splitter. "<<endl
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
  
  double cf(p_split->GetType()==pst::IS?m_as_is_fac:m_as_fs_fac);
  bool success(false);
  while (m_kperp2>=m_k0sq/cf) {
    ProduceT();
    SelectOne();
    split->SetSpect(p_selected->SelectSpec());
    m_z = Z();
    if (m_kperp2<m_k0sq/cf)  return false;
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
      msg_Error()<<"Error in Sudakov::Dice!"<<std::endl;
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
  m_phi = 2.0*M_PI*ran.Get();
  return success;
}

bool Sudakov::DefineFFBoundaries(double Q2,double x)
{
  if (4.*m_k0sq>Q2) return false;
  
  m_type=cstp::FF;
  double deltaz(sqrt(1.-4.*m_k0sq/Q2));
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
  if (m_k0sq*x>Q2*(1.-x))                                   return false;
  if (Q2<=p_pdf[beam]->Q2Min() || Q2>=p_pdf[beam]->Q2Max()) return false;
  
  m_type=cstp::FI;
  m_zmin   = Min(0.999,m_k0sq/Q2*x/(1.0-x));
  m_zmax   = Min(0.999,m_k0sq/Q2*xmax/(1.0-xmax));
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
  if (m_k0sq>Q2)                                            return false;
  if (Q2<=p_pdf[beam]->Q2Min() || Q2>=p_pdf[beam]->Q2Max()) return false;
  
  m_type=cstp::IF;
  m_zmin   = x/xmax;
  m_zmax   = Q2/(Q2+m_k0sq);
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
  if (m_k0sq>Q2)                                            return false;
  if (Q2<=p_pdf[beam]->Q2Min() || Q2>=p_pdf[beam]->Q2Max()) return false;
 
  m_type=cstp::II;
  m_zmin   = x/xmax;
  m_zmax   = Q2/(Q2+m_k0sq);
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
	  !(p_split->GetLeft()==p_spec || p_split->GetRight()==p_spec)) continue;
      (*m_splitter)->AddSpec(p_spec);
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
	(*m_splitter)->SetFlavourSpec(m_flspec);
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
  double ne=m_nspect*2.*M_PI/(m_lastint*m_rbmax);
  m_kperp2 *= exp(log(ran.Get())*Max(ne,1.0e-3));
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
  double wt(RejectionWeight(m_z,m_y,x,m_kperp2,Q2));
  if (ran.Get()>wt) return false;  
  return true;
}

