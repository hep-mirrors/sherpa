#include "Integrable_Base.H"
#include "Beam_Spectra_Handler.H"
#include "ISR_Handler.H"

#include "Run_Parameter.H"
#include "Combined_Selector.H"
#include "Standard_Selector.H"
#include "Phase_Space_Handler.H"
#include "Regulator_Base.H"
#include "Running_AlphaS.H"
#include "Message.H"
#include "MyStrStream.H"
#include "Data_Reader.H"
#include <algorithm>

using namespace PHASIC;
using namespace MODEL;
using namespace ATOOLS;

struct TDouble: public Term {
  double m_value;
};// end of struct Double

struct TVec4D: public Term {
  Vec4D m_value;
  TVec4D(const Vec4D &value): m_value(value) {}
};// end of struct Vec4D

#define CA 3.0
#define TR 0.5
#define Nf 5.0

class Order_Y {
public:
  bool operator()(const Vec4D &a,const Vec4D &b)
  {
    return a.Y()>b.Y();
  }
  bool operator()(const std::pair<Vec4D,Flavour> &a,
		  const std::pair<Vec4D,Flavour> &b)
  {
    return a.first.Y()>b.first.Y();
  }
};

std::ostream &PHASIC::operator<<(std::ostream &str,const cls::scheme &s)
{
  switch (s) {
  case cls::unknown: return str<<"<unknown>";
  case cls::sum: return str<<"sum";
  case cls::sample: return str<<"sample";
  }
  return str<<"<error>";
}

std::ostream &PHASIC::operator<<(std::ostream &str,const hls::scheme &s)
{
  switch (s) {
  case hls::unknown: return str<<"<unknown>";
  case hls::sum: return str<<"sum";
  case hls::sample: return str<<"sample";
  }
  return str<<"<error>";
}

Integrable_Base::Integrable_Base(const size_t nin,const size_t nout,
				 const scl::scheme scalescheme,const int kfactorscheme,
				 BEAM::Beam_Spectra_Handler *const beamhandler,
				 PDF::ISR_Handler *const isrhandler,
				 Selector_Data *const selectordata,
				 const cls::scheme &clsc,const hls::scheme &hlsc):
  m_name(""), m_nin(nin), m_nout(nout),
  m_nvector(ATOOLS::Max(nin+nout,(size_t)1)), m_corenout(nout),
  p_flavours(NULL), p_momenta(new Vec4D[ATOOLS::Max(nin+nout,(size_t)1)]), 
  m_scalescheme(scalescheme), m_colorscheme(clsc), m_helicityscheme(hlsc), m_kfactorscheme(kfactorscheme), 
  m_maxjetnumber(99), m_coremaxjetnumber(99),
  m_nstrong(0), m_neweak(0), m_orderQCD(-1), m_orderEW(-1), m_usepi(0),
  m_threshold(0.), m_overflow(0.), m_enhancefac(1.0), m_rfactor(1.0), m_xinfo(std::vector<double>(4)),
  m_n(0), m_expevents(1), m_dicedevents(0), m_accevents(0), m_last(0.), m_lastlumi(0.), m_lastdxs(0.), 
  m_max(0.), m_totalxs(0.),m_totalsum (0.), m_totalsumsqr(0.), m_totalerr(0.), 
  m_ssum(0.), m_ssumsqr(0.), m_smax(0.), m_ssigma2(0.), m_wmin(0.), m_me_as_factor(1.0), m_ycut(0.0), 
  m_sn(0), m_son(1), m_swaped(false), 
  p_selected(this), p_parent(this), 
  p_regulator(Regulator_Base::GetRegulator(this,"Identity",std::vector<double>())),
  p_beamhandler(beamhandler), p_isrhandler(isrhandler), 
  p_pshandler(NULL), p_activepshandler(NULL), p_selector(NULL), 
  p_cuts(NULL), p_whisto(NULL), p_jf(NULL), m_ownselector(true), m_efunc("1"), m_muf2tag(""),
  p_muf2calc(NULL), p_mur2calc(NULL), m_muf2tagset(this), m_mur2tagset(this)
{
  m_anasum=m_validanasum=0.0;
  m_expevents=m_dicedevents=m_accevents=0;
  m_gmin=-1.0;
  int kfs(ToType<int>(rpa.gen.Variable("S_KFACTOR_SCHEME","1"))&2);        
  double fssf(ToType<double>(rpa.gen.Variable("FS_CPL_SCALE_FACTOR","1.0")));
  double issf(ToType<double>(rpa.gen.Variable("IS_CPL_SCALE_FACTOR","1.0")));
  m_ps_kfactor=kfs?CA*(67./18.-M_PI*M_PI/6.)-10./9.*TR*Nf:0.0;
  m_ps_cpl_factor=Min(fssf,issf);
}

Integrable_Base::~Integrable_Base()
{
  if (p_mur2calc!=NULL) delete p_mur2calc;
  if (p_muf2calc!=NULL) delete p_muf2calc;
  if (p_selector!=NULL && m_ownselector) delete p_selector;
  if (p_flavours!=NULL) delete [] p_flavours;
  if (p_momenta!=NULL) delete [] p_momenta;
  if (p_whisto!=NULL) delete p_whisto;
  if (p_pshandler) delete p_pshandler;
  delete p_regulator;
}

Integrable_Base *const Integrable_Base::Selected()
{ 
  if (p_selected!=this && p_selected!=NULL) return p_selected->Selected();
  return this; 
}

Integrable_Base *const Integrable_Base::Parent()
{ 
  if (p_parent!=this && p_parent!=NULL) return p_parent->Parent();
  return this; 
}

double Integrable_Base::TotalResult()
{ 
//   if (m_ssigma2>0. && m_sn<1000) return m_totalsum/m_ssigma2; 
//   if (m_sn<1000) return m_ssum/m_sn;
//   double ssigma2 = (m_ssumsqr/m_sn - sqr(m_ssum/m_sn))/(m_sn-1);
//   return (m_totalsum+m_ssum/ssigma2/m_sn)/(m_ssigma2+1./ssigma2);
  if (m_ssigma2>0. && m_sn<1000) return m_totalsum/m_ssigma2; 
  if (m_sn<1000) return m_ssum/m_sn;
  double ssigma2 = sqr(m_ssum/m_sn)/((m_ssumsqr/m_sn - sqr(m_ssum/m_sn))/(m_sn-1));
  return (m_totalsum+m_ssum*ssigma2/m_sn)/(m_ssigma2+ssigma2);
}

double Integrable_Base::TotalVar() 
{
  if (m_nin==1 && m_nout==2) return 0.;
  if (m_sn<1000) {
    if (m_ssigma2>0.) return m_totalsum/m_ssigma2/sqrt(m_ssigma2); 
    else return TotalResult(); 
  }

  double disc = sqr(m_ssum/m_sn)/((m_ssumsqr/m_sn - sqr(m_ssum/m_sn))/(m_sn-1));
  if (disc>0.) return TotalResult()/sqrt(m_ssigma2+disc);
  
  return m_totalsum/m_ssigma2/sqrt(m_ssigma2);
}

double Integrable_Base::RemainTimeFactor(double maxerr) 
{
  if (m_sn<1000) return 0.;
  double disc = sqr(m_ssum/m_sn)/((m_ssumsqr/m_sn - sqr(m_ssum/m_sn))/(m_sn-1));
  return (sqr(1./maxerr)-m_ssigma2)/disc;
}

void Integrable_Base::SetMomenta(const Vec4D *momenta) 
{ 
  if (!p_momenta) {
    msg_Error()<<"Integrable_Base::SetMomenta("<<momenta<<"): "
		       <<"p_momenta = NULL. Abort."<<std::endl;
    abort();
  }
  for (size_t i=0;i<NVector();++i) p_momenta[i]=momenta[i];
  if (Selected()!=this) 
    for (size_t i=0;i<NVector();++i) 
      Selected()->p_momenta[i]=momenta[i];
}

void Integrable_Base::InitWeightHistogram() 
{
  if (p_whisto!=NULL) {
    delete p_whisto; };
  double av=TotalResult();
  if (!av>0.) {
    msg_Error()<<"Integrable_Base::InitWeightHistogram(): "
		       <<"No valid result: "<<av<<std::endl;
    return;
  }
  if (av<.3) av/=10.;
  av = exp(log(10.)*int(log(av)/log(10.)+0.5));
  p_whisto = new Histogram(10,av*1.e-4,av*1.e6,100);
}

void Integrable_Base::ReadInHistogram(std::string dir)
{
  std::string filename = dir+"/"+m_name;
  std::ifstream from;
  bool     hit = 0;
  from.open(filename.c_str());
  if (from) hit = 1;
  from.close();
  if (!hit) return;
  if (p_whisto) delete p_whisto; 
  p_whisto = new Histogram(filename);	
}

void Integrable_Base::WriteOutHistogram(std::string dir)
{
  if (!p_whisto) return;
  p_whisto->Output(dir+"/"+m_name);	
}

double Integrable_Base::GetMaxEps(double epsilon)
{
  if (!p_whisto) return m_max;
  int ovn=int(p_whisto->Fills()*epsilon);
  if (ovn<1) return m_max;
  double maxeps = m_max;
  double min = m_max*epsilon; //upper boundary for weight reduction
  for (int i=p_whisto->Nbin()+1;i>0;i--) {
    if (ovn>=p_whisto->Value(i)) {
      ovn-=(int)p_whisto->Value(i);
      maxeps=exp(log(10.)*(p_whisto->Xmin()+(i-1)*p_whisto->BinSize()));
    }
    else return ATOOLS::Max(Min(maxeps,m_max),min);
  }
  return m_max;
}

void Integrable_Base::AddPoint(const double value) 
{
  m_n++;
  m_sn++;
  m_ssum    += value;
  m_ssumsqr += value*value;
  if (value>m_max)  m_max  = value;
  if (value>m_smax) m_smax = value;
  if (p_whisto) p_whisto->Insert(value);
}

void Integrable_Base::SetMax(const double max, int depth) 
{
  if (max!=0.) m_max=max;
} 

void Integrable_Base::SetMax() 
{
  msg_Error()<<"Integrable_Base::SetMax(): Virtual function called !"<<std::endl;
} 

void Integrable_Base::ResetMax(int) 
{
  msg_Error()<<"Integrable_Base::ResetMax(): Virtual function called !"<<std::endl;
} 

Blob_Data_Base *Integrable_Base::OneEvent() 
{
  msg_Error()<<"Integrable_Base::OneEvent(): Virtual function called !"<<std::endl;
  return false;
} 

bool Integrable_Base::Trigger(const Vec4D *const momenta)
{
  return p_selector->Trigger(momenta);
} 

Blob_Data_Base *Integrable_Base::OneEvent(const double mass,const int mode) 
{
  return p_activepshandler->OneEvent(mass,mode);
} 

Blob_Data_Base *Integrable_Base::SameEvent() 
{
  return p_activepshandler->SameEvent();
  msg_Error()<<"Integrable_Base::SameEvent(): Virtual function called !"<<std::endl;
  return false;
} 

Blob_Data_Base *Integrable_Base::WeightedEvent(const int mode) 
{
  msg_Error()<<"Integrable_Base::WeightedEvent(): Virtual function called !"<<std::endl;
  return NULL;
} 

Blob_Data_Base *Integrable_Base::SameWeightedEvent() 
{
  msg_Error()<<"Integrable_Base::SameWeightedEvent(): Virtual function called !"<<std::endl;
  return NULL;
} 

void Integrable_Base::SetPSHandler(Phase_Space_Handler *const pshandler) 
{
  p_activepshandler=pshandler;
} 

void Integrable_Base::OptimizeResult()
{
  msg_Error()<<"Integrable_Base::OptimizeResult(): Virtual function called !"<<std::endl;
} 

void Integrable_Base::SetMomenta()   
{ 
  const Vec4D *p=p_activepshandler->Point();
  for (size_t i=0;i<m_nvector;++i) p_momenta[i]=p[i];
}

double Integrable_Base::CalculateScale(const Vec4D *momenta) 
{
  if (!m_kfkey.Assigned()) {
    std::string kfname(p_activepshandler->Process()->Name());
    std::string kfinfo("O(QCD)="+ToString(m_orderQCD));
    msg_Debugging()<<METHOD<<"(): Assign '"<<Name()
		   <<"' to '"<<kfname<<"','"<<kfinfo<<"'\n";
    m_kfkey.Assign(kfname,3,0,p_activepshandler->GetInfo());
    m_kfkey.SetInfo(kfinfo);
  }
  scl::scheme scheme(m_scalescheme);
  if (scheme==scl::unknown) 
    THROW(fatal_error,"Unknown scale scheme: "+ToString(m_scalescheme));
  SetMomenta(momenta);
  if (m_nin==1) {
    m_scale[stp::ren]=momenta[0].Abs2();
    m_kfkey[0]=m_scale[stp::ren];
    return momenta[0].Abs2();
  }
  if (m_nin!=2) THROW(fatal_error,"Too many incoming particles.");
  if (scheme==scl::ckkw) {
    double S(sqr(rpa.gen.Ecms()));
    m_scale[stp::ren]=m_ycut*sqr(p_jf->DeltaR())*S;
    m_scale[stp::fac]=m_cycut*S;
    double pt2(p_jf->ActualValue()*S);
    if ((int)m_corenout==m_coremaxjetnumber) {
      // highest multiplicity treatment
      if (m_nstrong>2) m_scale[stp::fac]=pt2;
    }
    else {
      // two scale treatment
      m_scale[stp::fac]=Min(m_scale[stp::fac],pt2);
    }
    SetRenormalizationScale();
    SetFactorizationScale();
    m_kfkey[0]=m_scale[stp::ren];
    m_kfkey[1]=m_scale[stp::fac];
    return m_scale[stp::fac]*m_ps_cpl_factor;
  }
  double pt2(-1.0);
  if (scheme & scl::fixed) {
    pt2=m_scale[stp::fac];
    scheme-=scl::fixed;
  }
  else if (scheme & scl::shat) {
    pt2=(momenta[0]+momenta[1]).Abs2();
    scheme-=scl::shat;
  }
  else if (scheme & scl::gmeanpt) {
    if (m_nin+m_nout==4) {
      double s = (momenta[0]+momenta[1]).Abs2();
      double t = (momenta[0]-momenta[2]).Abs2()-
	(sqr(p_flavours[2].PSMass())+sqr(p_flavours[3].PSMass()))/2.;
      double u = (momenta[0]-momenta[3]).Abs2()-
	(sqr(p_flavours[2].PSMass())+sqr(p_flavours[3].PSMass()))/2.;
      pt2 = 2.*s*t*u/(s*s+t*t+u*u);
    }
    else {
      pt2 = 1.;
      int count(0);
      for (size_t i=m_nin;i<m_nin+m_nout;i++) {
	if ((scheme & scl::strong) && !p_flavours[i].Strong()) continue;
	if (scheme & scl::mass) pt2 *= momenta[i].MPerp2();
	else pt2 *= momenta[i].PPerp2();
	count++;
      }
      if (count==0) pt2=(momenta[0]+momenta[1]).Abs2();
      else pt2 = pow(pt2,1./count);
    }
    scheme-=scl::gmeanpt;
    scheme-=scl::strong;
    scheme-=scl::mass;
  }
  else if (scheme & scl::ameanpt) {
    pt2 = 0.;
    int count(0);
    for (size_t i=m_nin;i<m_nin+m_nout;i++) {
      if ((scheme & scl::strong) && !p_flavours[i].Strong()) continue;
      if (scheme & scl::mass) pt2 += momenta[i].MPerp2();
      else pt2 += momenta[i].PPerp2();
      count++;
    }
    if (count==0) pt2=(momenta[0]+momenta[1]).Abs2();
    else pt2/=count;
    scheme-=scl::ameanpt;
    scheme-=scl::strong;
    scheme-=scl::mass;
  }
  else if (scheme & scl::minpt) {
    pt2 = (momenta[0]+momenta[1]).Abs2();
    for (size_t i=m_nin;i<m_nin+m_nout;i++) {
      if ((scheme & scl::strong) && !p_flavours[i].Strong()) continue;
      if (scheme & scl::mass) pt2 = Min(pt2,momenta[i].MPerp2());
      else pt2 = Min(pt2,momenta[i].PPerp2());
    }
    scheme-=scl::minpt;
    scheme-=scl::strong;
    scheme-=scl::mass;
  }
  else if (scheme & scl::sumpt) {
    pt2 = 0.;
    for (size_t i=m_nin;i<m_nin+m_nout;i++) {
      if ((scheme & scl::strong) && !p_flavours[i].Strong()) continue;
      if (scheme & scl::mass) pt2 += momenta[i].MPerp2();
      else pt2 += momenta[i].PPerp2();
    }
    scheme-=scl::sumpt;
    scheme-=scl::strong;
    scheme-=scl::mass;
  }
  else if (scheme & scl::updf) {
    double M2=0.;
    if (m_resonances.size()>0) {
      M2=ATOOLS::sqr(m_resonances[0].Mass());
    }
    ATOOLS::Vec4D *p=p_momenta;
    double S2=p[4]*p[5], x1=p[5]*p[0]/S2, x2=p[4]*p[1]/S2;
    double xi=(p[0]+p[1]).PMinus()/(p[0]+p[1]).PPlus();
    m_scale[PHASIC::stp::kp21]=x1*x1*2.0*S2*xi;
    m_scale[PHASIC::stp::kp22]=x2*x2*2.0*S2/xi;
    pt2=m_scale[PHASIC::stp::ren]=
      sqrt(m_scale[PHASIC::stp::kp21]*m_scale[PHASIC::stp::kp22]);
    scheme-=scl::updf;
    scheme-=scl::bfkl;
  }
  else if (scheme & scl::bfkl) {
    std::vector<Vec4D> moms(m_nout);
    for (size_t i(0);i<m_nout;++i) moms[i]=p_momenta[m_nin+i];
    std::sort(moms.begin(),moms.end(),Order_Y());
    m_kfkey[0]=m_scale[stp::kp21]=moms.front().PPerp2();
    m_kfkey[1]=m_scale[stp::kp22]=moms.back().PPerp2();
    msg_Debugging()<<"set scale "<<m_scale[stp::kp21]<<" "<<m_scale[stp::kp21]<<"\n";
    scheme-=scl::bfkl;
    if (scheme!=scl::unknown && scheme!=scl::reggeise) 
      THROW(fatal_error,"Unknown scale scheme: "+ToString(m_scalescheme));
    return m_scale[stp::fac]=-1.0;
  }
  if (scheme & scl::div_by_2) { 
    pt2/=2.; 
    scheme-=scl::div_by_2; 
  }
  if (scheme & scl::mult_by_2) { 
    pt2*=2.; 
    scheme-=scl::mult_by_2;  
  }
  if (pt2<0. || scheme!=scl::unknown) 
    THROW(fatal_error,"Unknown scale scheme: "+ToString(m_scalescheme));
  m_scale[stp::fac]=m_scale[stp::ren]=pt2;
  SetRenormalizationScale();
  SetFactorizationScale();
  m_kfkey[0]=m_scale[stp::ren];
  m_kfkey[1]=m_scale[stp::fac];
  msg_Debugging()<<"scale: "<<m_scale[stp::ren]<<" "<<m_scale[stp::ren]<<"\n";
  if (Selected()==NULL) return m_scale[stp::fac]*m_ps_cpl_factor;
  return (*Selected()->p_regulator)[m_scale[stp::fac]*m_ps_cpl_factor]; 
}

double Integrable_Base::KFactor() 
{
  if (!m_kfkey.Assigned()) {
    std::string kfname(p_activepshandler->Process()->Name());
    std::string kfinfo("O(QCD)="+ToString(m_orderQCD));
    msg_Debugging()<<METHOD<<"(): Assign '"<<Name()
		   <<"' to '"<<kfname<<"','"<<kfinfo<<"'\n";
    m_kfkey.Assign(kfname,3,0,p_activepshandler->GetInfo());
    m_kfkey.SetInfo(kfinfo);
  }
  if (m_scalescheme&scl::ckkw) {
    if (m_kfkey.Weight()!=ATOOLS::UNDEFINED_WEIGHT) return m_kfkey.Weight();
    if (m_orderQCD<0 || m_orderEW<0) {
      THROW(fatal_error,"Couplings not set for process '"+Name()+"'");
    }
    if (m_nstrong<=2) return m_rfactor;
    m_scale[stp::ren]=m_kfkey[0];
    m_scale[stp::fac]=m_kfkey[1];
    double asn(as->AlphaS(m_me_as_factor*m_ps_cpl_factor*m_scale[stp::ren]));
    m_kfkey[2]=asn;
    if (m_ps_kfactor!=0.0) asn*=1.+asn/(2.0*M_PI)*m_ps_kfactor;
    m_kfkey<<m_rfactor*pow(asn/as->AlphaS(rpa.gen.CplScale()),m_orderQCD);
    msg_Debugging()<<METHOD<<"(): "<<Name()<<" ("<<m_nstrong<<","
		   <<m_orderQCD<<") {\n"
		   <<"  \\mu_{fac}   = "<<sqrt(m_scale[stp::fac])<<"\n"
		   <<"  \\mu_{ren}   = "<<sqrt(m_scale[stp::ren])<<"\n"
		   <<"  me scalefac = "<<m_me_as_factor<<"\n"
		   <<"  r scalefac  = "<<m_rfactor<<"\n"
		   <<"  ps scalefac = "<<m_ps_cpl_factor<<"\n"
		   <<"  ps k factor = "<<m_ps_kfactor
		   <<"\n} -> as = "<<asn<<" => K = "<<m_kfkey.Weight()<<"\n";
    return m_kfkey.Weight();
  }
  else if (m_scalescheme&scl::bfkl) {
    static int nf(-1);
    if (nf<0) {
      Data_Reader read(" ",";","!","=");
      int sm;
      if (!read.ReadFromFile(sm,"BFKL_SPLIT_MODE")) sm=(1<<6)-1;
      else msg_Info()<<METHOD<<"(): Set splitting mode "<<sm<<".\n";
      nf=0;
      for (size_t i(1);i<=5;++i) if (sm&(1<<i)) ++nf;
      msg_Info()<<METHOD<<"(): N_f = "<<nf<<"\n";
    }
    if (m_orderQCD<0 || m_orderEW<0) {
      THROW(fatal_error,"Couplings not set for process '"+Name()+"'");
    }
    if (m_nstrong<=2) return 1.0;
    SetMomenta();
    double weight(1.0), ascplscale((*as)(rpa.gen.CplScale()));
    std::vector<std::pair<Vec4D,Flavour> > moms(m_nout);
    for (size_t i(0);i<m_nout;++i)
      moms[i]=std::pair<Vec4D,Flavour>
	(p_momenta[m_nin+i],p_flavours[m_nin+i]);
    std::sort(moms.begin(),moms.end(),Order_Y());
    BFKL_PT_Selector *ptsel
      ((BFKL_PT_Selector*)((Combined_Selector*)p_selector)->
       GetSelector("BFKL_PT_Selector"));
    double ptmin(0.0);
    if (ptsel!=NULL) ptmin=ptsel->PTMin();
    else THROW(critical_error,"'SCALE_SCHEME = BFKL' implies BFKL_PT <min> <max>' in 'Selector.dat'.");
    double asr((*as)(sqr(ptmin)));
    Vec4D q(p_momenta[p_momenta[0][3]>p_momenta[1][3]?0:1]);
    Flavour p(p_flavours[p_momenta[0][3]>p_momenta[1][3]?0:1]);
    q-=moms.front().first;
    msg_Debugging()<<METHOD<<"(): "<<Name()<<" ("<<m_nstrong<<","
		   <<m_orderQCD<<") "<<p<<" {\n";
    weight*=(*as)(moms.front().first.PPerp2())/ascplscale;
    msg_Debugging()<<"  k_{T,0} = "<<moms.front().first.PPerp()
		   <<" -> as = "<<(*as)(moms.front().first.PPerp2())<<"\n";
    for (size_t i(1);i<m_nout;++i) {
      weight*=(*as)(moms[i].first.PPerp2())/ascplscale;
      if (moms[i].second.IsQuark()) 
	p=p.IsGluon()?moms[i].second.Bar():Flavour(kf_gluon);
      if (m_scalescheme&scl::reggeise) {
	double qt2(q.PPerp2()), asc((*as)(qt2));
	double b0(as->Beta0(qt2)/M_PI);
	double z((q-moms[i].first).PPlus()/q.PPlus());
	if (z>1) z=q.PMinus()/(q-moms[i].first).PMinus();
 	double sud(exp(-(p.IsQuark()?4.0/3.0*(2.0+z):2.0*3.0+nf*0.5*z)
 		       /(2.0*M_PI*b0)*log(asr/asc)
 		       *(moms[i-1].first.Y()-moms[i].first.Y())));
	weight*=sud;
	msg_Debugging()<<"  q_{T,"<<i<<"} = "<<q.PPerp()<<" ("<<p
		       <<"<-"<<moms[i].second<<"), dy = "
		       <<moms[i-1].first.Y()-moms[i].first.Y()
		       <<" -> \\Delta = "<<sud<<"\n";
      }
      msg_Debugging()<<"  k_{T,"<<i<<"} = "
		     <<moms[i].first.PPerp()<<" -> as = "
		     <<(*as)(moms[i].first.PPerp2())<<"\n";
      q-=moms[i].first;
    }
    m_kfkey<<weight;
    msg_Debugging()<<"} -> K = "<<m_kfkey.Weight()<<"\n";
    return m_kfkey.Weight();
  }
  if (m_kfactorscheme==1) {
    if (m_kfkey.Weight()!=ATOOLS::UNDEFINED_WEIGHT) return m_kfkey.Weight();
    if (m_orderQCD<0 || m_orderEW<0) {
      THROW(fatal_error,"Couplings not set for process '"+Name()+"'");
    }
    if (m_orderQCD>0) {
      m_scale[stp::ren]=m_kfkey[0];
      m_kfkey<<m_rfactor*pow(as->AlphaS(m_scale[stp::ren])/
			     as->AlphaS(rpa.gen.CplScale()),m_orderQCD);
      if (m_kfkey.Weight()!=ATOOLS::UNDEFINED_WEIGHT)
	return m_kfkey.Weight();
    } 
    else 
      return m_rfactor;
  }
  return 1.;
}

double Integrable_Base::TriggerEfficiency()
{
  if (p_pshandler!=NULL) {
    Multi_Channel *fsr=p_pshandler->FSRIntegrator();
    return fsr->ValidN()/(double)fsr->N();
  }
  return 1.0;
}

void Integrable_Base::GetGMin(double &g, double &meff)
{
  if (m_validanasum==0.) return;
  double q = m_anasum/m_validanasum;
  double add = q*(double)m_dicedevents;
  add = Min(add,sqrt((double)m_expevents)*(double)m_dicedevents);
  meff += add;
  if (m_expevents<100) return;
  if (q<1.) if (m_expevents<sqrt(Parent()->m_expevents)) q=1.;
  g = Min(g,q);
}

void Integrable_Base::SetISRThreshold(const double threshold) 
{
  m_threshold=threshold;
}

void Integrable_Base::CreateMomenta(const size_t n)
{ 
  delete [] p_momenta; 
  p_momenta = new Vec4D[n]; 
  m_nvector=n; 
}

Spin_Correlation_Tensor* Integrable_Base::GetSpinCorrelations()
{
  return NULL;
}

void Integrable_Base::FillAmplitudes(HELICITIES::Amplitude_Tensor*,double)
{
  return;
}

std::map<std::string,std::string> Integrable_Base::ScaleTags()
{
  std::map<std::string,std::string> tags;
  tags["UNKNOWN"]=ToString(scl::unknown);
  tags["CKKW"]=ToString(scl::ckkw);
  tags["FIX_SCALE"]=ToString(scl::fixed);
  tags["S_HAT"]=ToString(scl::shat);
  tags["A_MEAN_PT2"]=ToString(scl::ameanpt);
  tags["G_MEAN_PT2"]=ToString(scl::gmeanpt);
  tags["MIN_PT2"]=ToString(scl::minpt);
  tags["SUM_PT2"]=ToString(scl::sumpt);
  tags["STRONG"]=ToString(scl::strong);
  tags["MASSIVE"]=ToString(scl::mass);
  tags["HALVE"]=ToString(scl::div_by_2);
  tags["DOUBLE"]=ToString(scl::mult_by_2);
  tags["UPDF"]=ToString(scl::updf);
  tags["BFKL"]=ToString(scl::bfkl);
  tags["REGGEISATION"]=ToString(scl::reggeise);
  return tags;
}
 
std::map<std::string,std::string> Integrable_Base::ColorSchemeTags()
{
  std::map<std::string,std::string> tags;
  tags["UNKNOWN"]=ToString((int)cls::unknown);
  tags["SUM"]=ToString((int)cls::sum);
  tags["SAMPLE"]=ToString((int)cls::sample);
  return tags;
}
 
std::map<std::string,std::string> Integrable_Base::HelicitySchemeTags()
{
  std::map<std::string,std::string> tags;
  tags["UNKNOWN"]=ToString((int)cls::unknown);
  tags["SUM"]=ToString((int)cls::sum);
  tags["SAMPLE"]=ToString((int)cls::sample);
  return tags;
}
 
void Integrable_Base::SetScaleScheme(const scl::scheme s)
{ 
  m_scalescheme=s;   
  if ((m_scalescheme&scl::ckkw) && m_ycut<=0.0 && p_selector!=NULL) {
    if (p_jf==NULL) {
      if (p_selector->Name()!="Combined_Selector")
	THROW(critical_error,"'SCALE_SCHEME = CKKW' implies JetFinder <ycut> <deltar>' in 'Selector.dat'.");
      p_jf=(Jet_Finder_Base *)
	((Combined_Selector*)p_selector)->GetSelector("Jetfinder");
      if (p_jf==NULL) {
	p_jf=(Jet_Finder_Base *)
	  ((Combined_Selector*)p_selector)->GetSelector("Dipole_Jetfinder");
	if (p_jf==NULL) 
	  THROW(critical_error,"'SCALE_SCHEME = CKKW' implies JetFinder <ycut> <deltar>' in 'Selector.dat'.");
      }
      m_me_as_factor=p_jf->Type()>1?1.0:0.25;
      p_jf->FillCombinations();
    }
    m_ycut=p_jf->Ycut();
    m_cycut=p_jf->GlobalCoreYcut();
    m_scale[stp::ren]=m_ycut*sqr(rpa.gen.Ecms()*p_jf->DeltaR());
    m_scale[stp::fac]=m_cycut*sqr(rpa.gen.Ecms());
  }
}

void Integrable_Base::SetFactorizationScale(const std::string &muf2)
{ 
  m_muf2tag=muf2;   
  if (m_muf2tag=="") return;
  msg_Debugging()<<METHOD<<"(): Set scale '"<<muf2<<"' in '"<<Name()<<"' {\n";
  msg_Indent();
  if (p_muf2calc==NULL) {
    p_muf2calc = new Algebra_Interpreter();
    p_muf2calc->SetTagReplacer(&m_muf2tagset);
    m_muf2tagset.SetCalculator(p_muf2calc);
  }
  p_muf2calc->AddTag("MU_F","1.0");
  p_muf2calc->AddTag("MU_R","1.0");
  p_muf2calc->AddTag("H_T","1.0");
  p_muf2calc->AddTag("Q_MIN","1.0");
  for (size_t i=0;i<NIn()+NOut();++i) 
    p_muf2calc->AddTag("p["+ToString(i)+"]",ToString(p_momenta[i]));
  p_muf2calc->Interprete(m_muf2tag);
  msg_Debugging()<<"}\n";
}

void Integrable_Base::SetRenormalizationScale(const std::string &mur2)
{ 
  m_mur2tag=mur2;   
  if (m_mur2tag=="") return;
  msg_Debugging()<<METHOD<<"(): Set scale '"<<mur2<<"' in '"<<Name()<<"' {\n";
  msg_Indent();
  if (p_mur2calc==NULL) {
    p_mur2calc = new Algebra_Interpreter();
    p_mur2calc->SetTagReplacer(&m_mur2tagset);
    m_mur2tagset.SetCalculator(p_mur2calc);
  }
  p_mur2calc->AddTag("MU_F","1.0");
  p_mur2calc->AddTag("MU_R","1.0");
  p_mur2calc->AddTag("H_T","1.0");
  p_mur2calc->AddTag("Q_MIN","1.0");
  for (size_t i=0;i<NIn()+NOut();++i) 
    p_mur2calc->AddTag("p["+ToString(i)+"]",ToString(p_momenta[i]));
  p_mur2calc->Interprete(m_mur2tag);
  msg_Debugging()<<"}\n";
}

void Integrable_Base::SetFactorizationScale()
{ 
  if (m_muf2tag!="") {
    m_scale[stp::fac]=((TDouble*)p_muf2calc->Calculate())->m_value;
    msg_Debugging()<<METHOD<<"(): Set \\mu_f = "
		   <<sqrt(m_scale[stp::fac])<<"\n";
  }
}

void Integrable_Base::SetRenormalizationScale()
{ 
  if (m_mur2tag!="") {
    m_scale[stp::ren]=((TDouble*)p_mur2calc->Calculate())->m_value;
    msg_Debugging()<<METHOD<<"(): Set \\mu_r = "
		   <<sqrt(m_scale[stp::ren])<<"\n";
  }
}

bool Integrable_Base::FillSIntegrator(Multi_Channel *&mc)
{
  return false;
}

void Integrable_Base::UpdateIntegrator(Multi_Channel *&mc)
{
}

std::string Tag_Setter::ReplaceTags(std::string &expr) const
{
  return p_calc->ReplaceTags(expr);
}

Term *Tag_Setter::ReplaceTags(Term *term) const
{
  if (term->m_tag=="MU_F") {
    ((TDouble*)term)->m_value=p_ib->Scale(stp::fac);
    return term;
  }
  if (term->m_tag=="MU_R") {
    ((TDouble*)term)->m_value=p_ib->Scale(stp::ren);
    return term;
  }
  if (term->m_tag=="H_T") {
    double ht(0.0);
    for (size_t i(0);i<p_ib->NOut();++i) ht+=p_ib->Momenta()[i+p_ib->NIn()].PPerp();
    ((TDouble*)term)->m_value=sqr(ht);
    return term;
  }
  if (term->m_tag=="Q_MIN") {
    if (p_ib->JetFinder()==NULL) 
      THROW(fatal_error,"Q_MIN cannot be used without jet finder");
    ((TDouble*)term)->m_value=p_ib->JetFinder()->ActualValue()*sqr(rpa.gen.Ecms());
    return term;
  }
  int i=ATOOLS::ToType<int>(term->m_tag.substr(2,term->m_tag.length()-3));
  ((TVec4D*)term)->m_value=p_ib->Momenta()[i];
  return term;
}

