#include "PHASIC++/Scales/KFactor_Setter_Base.H"

#include "PHASIC++/Scales/MINLO_Scale_Setter.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Org/Default_Reader.H"
#include "ATOOLS/Org/Run_Parameter.H"

#include <map>

namespace PHASIC {

  class Sudakov: public ATOOLS::Function_Base {
  private:

    ATOOLS::Flavour m_fl;
    ATOOLS::Gauss_Integrator m_gauss;
    MODEL::Running_AlphaS *p_as;

    double m_Q2, m_mur2, m_prec;
    int m_fo, m_mode, m_nfgs;

    double K(const double &nf) const;

    double Ggq(const double &e,const double &q2,const double &m=0.0);

  public:

    Sudakov(const ATOOLS::Flavour &fl,const int mode,
	    const int nfgs,const double prec);

    double Delta(const double &q2,const double &Q2);
    double Delta1(const double &q2,const double &Q2,const double &mur2);

    double operator()(double q2);

  };// end of class Sudakov

  class MINLO_KFactor_Setter: public KFactor_Setter_Base {
  private:

    MINLO_Scale_Setter *p_minlo;

    std::map<ATOOLS::Flavour,Sudakov*> m_suds;

    double m_sudweight, m_lastmuR2, m_lastq02[2];
    int    m_vmode, m_rsfvar, m_ordonly;

  public:

    MINLO_KFactor_Setter(const KFactor_Setter_Arguments &args);

    ~MINLO_KFactor_Setter();

    double KFactor(const int mode);

    bool UpdateKFactor(const ATOOLS::Variation_Parameters &var);

  };// end of class MINLO_KFactor_Setter

}// end of namespace PHASIC

using namespace PHASIC;
using namespace ATOOLS;

DECLARE_GETTER(MINLO_KFactor_Setter,"MINLO",
	       KFactor_Setter_Base,KFactor_Setter_Arguments);

KFactor_Setter_Base *ATOOLS::Getter
<KFactor_Setter_Base,KFactor_Setter_Arguments,MINLO_KFactor_Setter>::
operator()(const KFactor_Setter_Arguments &args) const
{
  return new MINLO_KFactor_Setter(args);
}

void ATOOLS::Getter<KFactor_Setter_Base,KFactor_Setter_Arguments,
		    MINLO_KFactor_Setter>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"MINLO kfactor scheme\n";
}

MINLO_KFactor_Setter::MINLO_KFactor_Setter
(const KFactor_Setter_Arguments &args):
  KFactor_Setter_Base(args), m_rsfvar(0)
{
  p_minlo=dynamic_cast<MINLO_Scale_Setter*>(p_proc->ScaleSetter());
  if (p_minlo==NULL) THROW(fatal_error,"Must use MINLO scale");
  Default_Reader reader;
  int mode(reader.Get<int>("MINLO_SUDAKOV_MODE",3));
  m_ordonly=reader.Get<int>("MINLO_ORDERED_ONLY",0);
  int nfgs(reader.Get<int>("MINLO_SUDAKOV_NF_GSPLIT",6));
  double prec(reader.Get<double>("MINLO_SUDAKOV_PRECISION",1.0e-4));
  m_suds[Flavour(kf_gluon)] = new Sudakov(Flavour(kf_gluon),mode,nfgs,prec);
  for (size_t i(0);i<=6;++i) {
    m_suds[Flavour(i,0)] = new Sudakov(Flavour(i,0),mode,nfgs,prec);
    m_suds[Flavour(i,1)] = new Sudakov(Flavour(i,1),mode,nfgs,prec);
  }
  if (reader.Get<int>("MINLO_SELF_TEST",0)) {
    int fl(reader.Get<int>("MINLO_SELF_TEST_FLAV",21));
    double ecm2(reader.Get<double>("MINLO_SELF_TEST_ECM",91.2));
    std::ofstream q(("R2_"+ToString(fl)+"_"
		     +ToString(ecm2)+".dat").c_str(),std::ios::trunc);
    ecm2=sqr(ecm2);
    for (double y=0.0;y>=-4.0;y-=0.01)
      q<<y<<" "<<sqr(m_suds[Flavour(fl)]->Delta(pow(10,y)*ecm2,ecm2))<<"\n";
    q.close();
    THROW(normal_exit,"Done");
  }
}

MINLO_KFactor_Setter::~MINLO_KFactor_Setter()
{
  for (std::map<ATOOLS::Flavour,Sudakov*>::const_iterator
	 sit(m_suds.begin());sit!=m_suds.end();++sit)
    delete sit->second;
}

bool MINLO_KFactor_Setter::UpdateKFactor(const Variation_Parameters &var)
{
  DEBUG_FUNC("K = "<<m_sudweight<<" * ( 1 + "<<m_weight/m_sudweight-1.0<<" )");
  if (p_minlo->Q02(0)==m_lastq02[0] &&
      p_minlo->Q02(1)==m_lastq02[1] &&
      p_proc->ScaleSetter()->Scale(stp::ren)==m_lastmuR2) return false;
  if (p_minlo->Q02(0)==m_lastq02[0] &&
      p_minlo->Q02(1)==m_lastq02[1]) m_rsfvar=1;
  KFactor(m_vmode);
  m_rsfvar=0;
  return true;
}

double MINLO_KFactor_Setter::KFactor(const int mode) 
{
  m_vmode=mode;
  if (!m_on) return 1.0;
  DEBUG_FUNC(p_proc->Name()<<", mode = "<<mode);
#ifdef DEBUG__MINLO
  std::cout.precision(12);
  std::string name(p_proc->Name().substr(5));
  for (size_t pos(name.find('_'));pos!=std::string::npos;pos=name.find('_')) name.replace(pos,1," ");
  msg_Debugging()<<"DEBUG MINLO Event "<<rpa->gen.NumberOfGeneratedEvents()+1
		 <<" { "<<name<<"\n";
#endif
  m_weight=1.0;
  double muR2(p_minlo->MuRAvg(1));
  double Q02[2]={p_minlo->Q02(0),p_minlo->Q02(1)}, sub(0.0);
  m_lastq02[0]=p_minlo->Q02(0);
  m_lastq02[1]=p_minlo->Q02(1);
#ifdef DEBUG__MINLO
  msg_Debugging()<<"DEBUG MINLO   \\mu_{NLO} = "<<sqrt(muR2)<<" -> "<<(*MODEL::as)(muR2)<<"\n";
#endif
  m_lastmuR2=muR2;
  Cluster_Amplitude *ampl=p_proc->Info().Has(nlo_type::real)&&
    p_minlo->Ampl()->Next()?p_minlo->Ampl()->Next():p_minlo->Ampl();
#ifdef DEBUG__MINLO
  bool ord(!(ampl->Flag()&1));
  int step(0);
#endif
  for (;ampl->Next();ampl=ampl->Next()) {
    Cluster_Amplitude *next(ampl->Next());
    msg_Debugging()<<*ampl<<"\n";
#ifdef DEBUG__MINLO
    msg_Debugging()<<"DEBUG MINLO   Step "<<++step<<" {\n";
#endif
    for (size_t i(0);i<next->Legs().size();++i) {
      Cluster_Leg *l(next->Leg(i));
#ifdef DEBUG__MINLO
      if (l->K()) {
	Cluster_Leg *li(NULL), *lj(NULL);
	for (size_t j(0);j<ampl->Legs().size();++j)
	  if (ampl->Leg(j)->Id()&l->Id())
	    if (li==NULL) li=ampl->Leg(j);
	    else lj=ampl->Leg(j);
	msg_Debugging()<<"DEBUG MINLO     Clustering "<<ID(li->Id())
		       <<ID(lj->Id())<<"->"<<ID(l->Id())<<" "
		       <<l->Flav()<<" "<<ampl->KT2()<<"\n";
      }
#endif
      std::map<ATOOLS::Flavour,Sudakov*>::iterator sit(m_suds.find(l->Flav()));
      if (sit==m_suds.end()) continue;
      int is((l->Id()&3)?1:0);
      double gamma[2]={0.0,0.0};
      if (mode==1 && p_proc->Info().m_fi.m_nlotype!=nlo_type::lo) {
	gamma[0]=sit->second->Delta1(Q02[is],next->KT2(),muR2);
	gamma[1]=sit->second->Delta1(Q02[is],ampl->KT2(),muR2);
      }
      double delta[2]={m_rsfvar?0.0:sit->second->Delta(Q02[is],next->KT2()),
		       m_rsfvar?0.0:sit->second->Delta(Q02[is],ampl->KT2())};
      msg_Debugging()<<"Sudakov for "<<ID(l->Id())<<"["<<is
		     <<"] -> \\Delta_{"<<l->Flav()<<"}("<<sqrt(Q02[is])
		     <<","<<sqrt(next->KT2())<<") / \\Delta_{"<<l->Flav()
		     <<"}("<<sqrt(Q02[is])<<","<<sqrt(ampl->KT2())
		     <<") = "<<delta[0]<<" / "<<delta[1]<<" = "
		     <<delta[0]/delta[1]<<" ("<<gamma[0]-gamma[1]<<")\n";
#ifdef DEBUG__MINLO
      msg_Debugging()<<"DEBUG MINLO     Sudakovs "<<ID(l->Id())
		     <<" -> \\Delta_{"<<l->Flav()<<"}("<<sqrt(Q02[is])
		     <<","<<sqrt(next->KT2())<<") / \\Delta_{"<<l->Flav()
		     <<"}("<<sqrt(Q02[is])<<","<<sqrt(ampl->KT2())
		     <<") = "<<delta[0]<<" / "<<delta[1]<<" = "
		     <<delta[0]/delta[1]<<"\n";
      msg_Debugging()<<"DEBUG MINLO     Sudakov subtractions "<<ID(l->Id())
		     <<" -> \\Delta^1_{"<<l->Flav()<<"}("<<sqrt(Q02[is])
		     <<","<<sqrt(next->KT2())<<") - \\Delta^1_{"<<l->Flav()
		     <<"}("<<sqrt(Q02[is])<<","<<sqrt(ampl->KT2())
		     <<") = "<<-gamma[0]<<" - "<<-gamma[1]<<" = "
		     <<(-gamma[0]+gamma[1])<<"\n";
#endif
      m_weight*=delta[0]/delta[1];
      sub+=gamma[0]-gamma[1];
    }
#ifdef DEBUG__MINLO
    msg_Debugging()<<"DEBUG MINLO   }\n";
#endif
    if (next->Next()==NULL) {
      msg_Debugging()<<*next<<"\n";
#ifdef DEBUG__MINLO
      msg_Debugging()<<"DEBUG MINLO   Core "<<next->KT2()<<"\n";
#endif
      if ((next->Flag()&1) && m_ordonly) {
	msg_Debugging()<<"Unordered configuration\n";
	if (m_ordonly&2) m_weight=0.0;
#ifdef DEBUG__MINLO
	ord=false;
#endif
      }
    }
  }
#ifdef DEBUG__MINLO
  msg_Debugging()<<"DEBUG MINLO } "<<(ord?"ordered":"unordered")
		 <<" -> weight = "<<m_weight<<" * ( 1 - "<<-sub
		 <<" ) = "<<m_weight*(1.0+sub)<<"\n";
#endif
  if (m_rsfvar) {
    msg_Debugging()<<"w = "<<m_sudweight<<" * ( 1 + "<<sub<<" )\n";
    return m_weight=m_sudweight*(1.0+sub);
  }
  msg_Debugging()<<"w = "<<m_weight<<" * ( 1 + "<<sub<<" )\n";
  m_sudweight=m_weight;
  return m_weight*=(1.0+sub);
}

Sudakov::Sudakov(const ATOOLS::Flavour &fl,const int mode,
		 const int nfgs,const double prec):
  m_fl(fl), m_gauss(this), p_as(MODEL::as),
  m_mode(mode), m_nfgs(nfgs), m_prec(prec)
{
}

double Sudakov::Delta(const double &q2,const double &Q2)
{
  if (q2>=Q2) return 1.0;
  m_fo=0;
  m_Q2=Q2;
  return exp(-m_gauss.Integrate(q2,Q2,m_prec));
}

double Sudakov::Delta1(const double &q2,const double &Q2,const double &mur2)
{
  if (q2>=Q2) return 0.0;
  m_fo=1;
  m_Q2=Q2;
  m_mur2=mur2;
  return m_gauss.Integrate(q2,Q2,m_prec);
}

double Sudakov::K(const double &nf) const
{
  if (m_fo || !(m_mode&2)) return 0.0;
  return 3.*(67./18.-sqr(M_PI)/6.)-10./9.*nf/2.;
}

double Sudakov::operator()(double q2)
{
  double eps(sqrt(q2/m_Q2));
  double e((m_mode&1)?eps:0.0);
  double nf(Min(m_nfgs,p_as->Nf(m_fo?m_mur2:q2)));
  double as2pi((*p_as)(m_fo?m_mur2:q2)/(2.0*M_PI));
  if (m_fl.IsQuark()) {
    double gam=as2pi/q2*4.0/3.0*
      (2.0*log(1.0/eps)*(1.0+as2pi*K(nf))
       -3.0/2.0*sqr(1.0-e));
    if (!m_fl.IsMassive()) return gam;
    double k=sqrt(q2)/m_fl.Mass();
    gam+=as2pi/q2*4.0/3.0*as2pi*K(nf)*
      (-(1.0-e)*k*k/(e*e+k*k)
       +k*atan((1-e)*k/(e+k*k))+log((e*e+k*k)/(1.0+k*k)));
    gam+=as2pi/q2*4.0/3.0*
      ((1.0-e*e)/2.0+e*sqr(1.0-e)*(e/(k*k+e*e)-(1.0-e)/(k*k+sqr(1.0-e)))
       -k*(atan(1.0/k)+(1.0-k*k)*atan(e*k/(1.0-e+k*k)))
       -(1.0-k*k/2.0)*(log((k*k+1.0)/(k*k+e*e))-2.0*k*atan(e/k)));
    return gam;
  }
  if (m_fl.IsGluon()) {
    double gam0=Ggq(e,q2), gam=3.0*gam0;
    for (long int i(4);i<=m_nfgs;++i) {
      if (Flavour(i).Mass()) gam+=Ggq(e,q2,Flavour(i).Mass());
      else if (nf>=i) gam+=gam0;
    }
    return as2pi/q2*3.0*
      (2.0*log(1.0/eps)*(1.0+as2pi*K(nf))
       -sqr(1.0-e)/6.0*(11.0-e*(2.0-3.0*e)))
      +as2pi*gam;
  }
  return 0.0;
}

double Sudakov::Ggq(const double &e,const double &q2,const double &m)
{
  if (m*m>q2) return 0.0;
  return 0.5/(q2+m*m)*sqr(1.0-e)*(1.0-(1.0-e)*(1.0+3.0*e)/3.0*q2/(q2+m*m));
}
