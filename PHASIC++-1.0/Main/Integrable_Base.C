#include "Integrable_Base.H"
#include "Beam_Spectra_Handler.H"
#include "ISR_Handler.H"

#include "Run_Parameter.H"
#include "Combined_Selector.H"
#include "Jet_Finder.H"
#include "Phase_Space_Handler.H"
#include "Regulator_Base.H"
#include "Running_AlphaS.H"
#include "Message.H"

using namespace PHASIC;

Integrable_Base::Integrable_Base(const size_t nin,const size_t nout,
				 const int scalescheme,const int kfactorscheme,
				 BEAM::Beam_Spectra_Handler *const beamhandler,
				 PDF::ISR_Handler *const isrhandler,
				 ATOOLS::Selector_Data *const selectordata):
  m_name(""), m_nin(nin), m_nout(nout), m_naddin(0), m_naddout(0), 
  m_nvector(ATOOLS::Max(nin+nout,(size_t)1)), p_flavours(NULL), p_addflavours(NULL), 
  p_momenta(new ATOOLS::Vec4D[ATOOLS::Max(nin+nout,(size_t)1)]), p_addmomenta(NULL), 
  m_scalescheme(scalescheme), m_kfactorscheme(kfactorscheme), 
  m_nstrong(0), m_neweak(0), m_usepi(0),
  m_threshold(0.), m_overflow(0.), m_rfactor(1.0), m_enhancefac(1.0), m_xinfo(std::vector<double>(4)),
  m_n(0), m_expevents(1), m_dicedevents(0), m_accevents(0), m_last(0.), m_lastlumi(0.), m_lastdxs(0.), 
  m_max(0.), m_totalxs(0.),m_totalsum (0.), m_totalsumsqr(0.), m_totalerr(0.), 
  m_ssum(0.), m_ssumsqr(0.), m_smax(0.), m_ssigma2(0.), m_wmin(0.), 
  m_sn(0), m_son(1), m_swaped(false), 
  p_selected(this), p_parent(this), 
  p_regulator(Regulator_Base::GetRegulator(this,"Identity",std::vector<double>())),
  p_beamhandler(beamhandler), p_isrhandler(isrhandler), 
  p_pshandler(NULL), p_activepshandler(NULL), p_selector(NULL), 
  p_cuts(NULL), p_whisto(NULL), m_ownselector(true), m_efunc("1") {}

Integrable_Base::~Integrable_Base()
{
  if (p_selector!=NULL && m_ownselector) delete p_selector;
  if (p_flavours!=NULL) delete [] p_flavours;
  if (p_momenta!=NULL) delete [] p_momenta;
  if (p_addflavours!=NULL) delete [] p_addflavours;
  if (p_addmomenta!=NULL) delete [] p_addmomenta;
  if (p_whisto!=NULL) delete p_whisto;
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
//   double ssigma2 = (m_ssumsqr/m_sn - ATOOLS::sqr(m_ssum/m_sn))/(m_sn-1);
//   return (m_totalsum+m_ssum/ssigma2/m_sn)/(m_ssigma2+1./ssigma2);
  if (m_ssigma2>0. && m_sn<1000) return m_totalsum/m_ssigma2; 
  if (m_sn<1000) return m_ssum/m_sn;
  double ssigma2 = ATOOLS::sqr(m_ssum/m_sn)/((m_ssumsqr/m_sn - ATOOLS::sqr(m_ssum/m_sn))/(m_sn-1));
  return (m_totalsum+m_ssum*ssigma2/m_sn)/(m_ssigma2+ssigma2);
}

double Integrable_Base::TotalVar() 
{
  if (m_nin==1 && m_nout==2) return 0.;
  if (m_sn<1000) {
    if (m_ssigma2>0.) return m_totalsum/m_ssigma2/sqrt(m_ssigma2); 
    else return TotalResult(); 
  }

  double disc = ATOOLS::sqr(m_ssum/m_sn)/((m_ssumsqr/m_sn - ATOOLS::sqr(m_ssum/m_sn))/(m_sn-1));
  if (disc>0.) return TotalResult()/sqrt(m_ssigma2+disc);
  
  return m_totalsum/m_ssigma2/sqrt(m_ssigma2);
}

double Integrable_Base::RemainTimeFactor(double maxerr) 
{
  if (m_sn<1000) return 0.;
  double disc = ATOOLS::sqr(m_ssum/m_sn)/((m_ssumsqr/m_sn - ATOOLS::sqr(m_ssum/m_sn))/(m_sn-1));
  return (ATOOLS::sqr(1./maxerr)-m_ssigma2)/disc;
}

void Integrable_Base::SetMomenta(const ATOOLS::Vec4D *momenta) 
{ 
  if (!p_momenta) {
    ATOOLS::msg.Error()<<"Integrable_Base::SetMomenta("<<momenta<<"): "
		       <<"p_momenta = NULL. Abort."<<std::endl;
    abort();
  }
  for (size_t i=0;i<NVector();++i) p_momenta[i]=momenta[i];
  if (Selected()!=this) 
    for (size_t i=0;i<NVector();++i) 
      Selected()->p_momenta[i]=momenta[i];
}

void Integrable_Base::SetAddMomenta(const ATOOLS::Vec4D *momenta) 
{ 
  if (!p_addmomenta) {
    ATOOLS::msg.Error()<<"Integrable_Base::SetAddMomenta("<<momenta<<"): "
		       <<"p_addmomenta = NULL. Abort."<<std::endl;
    abort();
  }
  for (size_t i=0;i<NAddIn()+NAddOut();++i) p_addmomenta[i]=momenta[i];
  if (Selected()!=this) 
    for (size_t i=0;i<NAddIn()+NAddOut();++i) 
      Selected()->p_addmomenta[i]=momenta[i];
}

void Integrable_Base::SetAddFlavours(const ATOOLS::Flavour *flavours) 
{ 
  if (!p_addflavours) {
    ATOOLS::msg.Error()<<"Integrable_Base::SetAddFlavours("<<flavours<<"): "
		       <<"p_addflavours = NULL. Abort."<<std::endl;
    abort();
  }
  for (size_t i=0;i<NAddIn()+NAddOut();++i) p_addflavours[i]=flavours[i];
}

void Integrable_Base::InitWeightHistogram() 
{
  if (p_whisto!=NULL) {
    delete p_whisto; };
  double av=TotalResult();
  if (!av>0.) {
    ATOOLS::msg.Error()<<"Integrable_Base::InitWeightHistogram(): "
		       <<"No valid result: "<<av<<std::endl;
    return;
  }
  if (av<.3) av/=10.;
  av = exp(log(10.)*int(log(av)/log(10.)+0.5));
  p_whisto = new ATOOLS::Histogram(10,av*1.e-4,av*1.e6,100);
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
  p_whisto = new ATOOLS::Histogram(filename);	
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
    else return ATOOLS::Max(ATOOLS::Min(maxeps,m_max),min);
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

void Integrable_Base::SetScale(const double scale) 
{ 
  ATOOLS::msg.Error()<<"Integrable_Base::SetScale("<<scale
		     <<"): Virtual function called."<<std::endl;
}

void Integrable_Base::SetMax(const double max, int depth) 
{
  if (max!=0.) m_max=max;
} 

void Integrable_Base::SetMax() 
{
  ATOOLS::msg.Error()<<"Integrable_Base::SetMax(): Virtual function called !"<<std::endl;
} 

void Integrable_Base::ResetMax(int) 
{
  ATOOLS::msg.Error()<<"Integrable_Base::ResetMax(): Virtual function called !"<<std::endl;
} 

ATOOLS::Blob_Data_Base *Integrable_Base::OneEvent() 
{
  ATOOLS::msg.Error()<<"Integrable_Base::OneEvent(): Virtual function called !"<<std::endl;
  return false;
} 

#ifdef EVENT_STATISTICS 
#include <map>
#include <iomanip>

#define DFORMAT std::setw(15)
#define LIFORMAT std::setw(8)

typedef std::map<std::string,std::vector<double> > Mem_Map;

struct Memory_Map: public Mem_Map {

  // destructor
  ~Memory_Map();

};// end of struct Memory_Map

Memory_Map::~Memory_Map() 
{
  msg_Info()<<"Memory_Map::~Memory_Map(): Weighted analysis statistics {\n\n";
  {
    msg_Indent();
    std::cout.precision(6);
    double expectedxs=0.0, dicedxs=0.0;
    double expectedevents=0.0, dicedevents=0.0, untriggered=0.0;
    msg_Info()<<std::setw(25)<<std::setiosflags(std::ios::left)<<"Process"
	      <<std::resetiosflags(std::ios::left)<<std::setiosflags(std::ios::right)
	      <<DFORMAT<<"exp. xs"<<DFORMAT<<"diced xs"
	      <<LIFORMAT<<"exp. n"<<LIFORMAT<<"diced n"
	      <<LIFORMAT<<"acc. n"<<LIFORMAT<<"extra n"
	      <<std::resetiosflags(std::ios::right)<<std::endl;
    for (Mem_Map::const_iterator mit=begin();mit!=end();++mit) {
      if (mit->second.size()==0) return;
      expectedxs+=mit->second[0];
      expectedevents+=mit->second[1];
      dicedxs+=mit->second[2];
      dicedevents+=mit->second[3];
      untriggered+=mit->second[4];
      msg_Info()<<std::setw(25)<<std::setiosflags(std::ios::left)<<mit->first
		<<std::resetiosflags(std::ios::left)<<std::setiosflags(std::ios::right)
		<<DFORMAT<<(mit->second[0]*ATOOLS::rpa.Picobarn())
		<<DFORMAT<<(mit->second[2]/mit->second[4]*ATOOLS::rpa.Picobarn())
		<<LIFORMAT<<mit->second[1]<<LIFORMAT<<mit->second[4]
		<<LIFORMAT<<mit->second[3]<<LIFORMAT<<mit->second[5]
		<<std::resetiosflags(std::ios::right)<<std::endl;
    }
    msg_Info()<<"\n";
    msg_Info()<<"exp. n = "<<expectedevents<<", diced n = "<<untriggered
	      <<", acc. n = "<<dicedevents<<"\n";
    msg_Info()<<"exp. xs = "<<(expectedxs/untriggered*ATOOLS::rpa.Picobarn())
	      <<", diced xs = "<<(dicedxs/untriggered*ATOOLS::rpa.Picobarn())<<std::endl;
  }
  msg_Info()<<"\n}"<<std::endl;
}

static Memory_Map mem;
#endif

bool Integrable_Base::Trigger(const ATOOLS::Vec4D *const momenta)
{
  return p_selector->Trigger(momenta);
} 

ATOOLS::Blob_Data_Base *Integrable_Base::OneEvent(const double mass,const int mode) 
{
  return p_activepshandler->OneEvent(mass,mode);
} 

ATOOLS::Blob_Data_Base *Integrable_Base::SameEvent() 
{
  return p_activepshandler->SameEvent();
  ATOOLS::msg.Error()<<"Integrable_Base::SameEvent(): Virtual function called !"<<std::endl;
  return false;
} 

ATOOLS::Blob_Data_Base *Integrable_Base::WeightedEvent(const int mode) 
{
  ATOOLS::msg.Error()<<"Integrable_Base::WeightedEvent(): Virtual function called !"<<std::endl;
  return NULL;
} 

ATOOLS::Blob_Data_Base *Integrable_Base::WeightedEventNS(const int mode) 
{
  ATOOLS::msg.Error()<<"Integrable_Base::WeightedEventNS(): Virtual function called !"<<std::endl;
  return NULL;
} 

ATOOLS::Blob_Data_Base *Integrable_Base::SameWeightedEvent() 
{
  ATOOLS::msg.Error()<<"Integrable_Base::SameWeightedEvent(): Virtual function called !"<<std::endl;
  return NULL;
} 

void Integrable_Base::SetPSHandler(Phase_Space_Handler *const pshandler) 
{
  p_activepshandler=pshandler;
} 

void Integrable_Base::OptimizeResult()
{
  ATOOLS::msg.Error()<<"Integrable_Base::OptimizeResult(): Virtual function called !"<<std::endl;
} 

void Integrable_Base::SetMomenta()   
{ 
  const ATOOLS::Vec4D *p=p_activepshandler->Point();
  for (size_t i=0;i<m_nvector;++i) p_momenta[i]=p[i];
}

double Integrable_Base::CalculateScale(const ATOOLS::Vec4D *momenta) 
{
  SetMomenta(momenta);
  if (m_nin==1) return momenta[0].Abs2();
  if (m_nin!=2) {
    ATOOLS::msg.Error()<<"Integrable_Base::CalculateScale(..): "
		       <<"Too many particles."<<std::endl;
    abort();
  }
  double s;
  if (m_nin==2) s = (momenta[0]+momenta[1]).Abs2();
  if (m_nin==1) s = (momenta[0]).Abs2();
  double pt2 = 0.;

  //new
  switch (m_scalescheme) {
  case 1 :   
    pt2 = m_scale[stp::as];
    break;
  case 2  :
    if (m_nin+m_nout==4) {
      double t = (momenta[0]-momenta[2]).Abs2()-
	(ATOOLS::sqr(p_flavours[2].PSMass())+ATOOLS::sqr(p_flavours[3].PSMass()))/2.;
      double u = (momenta[0]-momenta[3]).Abs2()-
	(ATOOLS::sqr(p_flavours[2].PSMass())+ATOOLS::sqr(p_flavours[3].PSMass()))/2.;
      pt2 = 4.*s*t*u/(s*s+t*t+u*u);
    }
    else {
      pt2 = 0.;
      for (size_t i=m_nin;i<m_nin+m_nout;i++) {
	pt2 += ATOOLS::sqr(momenta[i][1])+ATOOLS::sqr(momenta[i][2]);
      }
    }
    break;
  case 3  :
    pt2 = s;
    double pt2i;
    for (size_t i=m_nin;i<m_nin+m_nout;i++) {
      pt2i = ATOOLS::sqr(momenta[i][1])+ATOOLS::sqr(momenta[i][2]);
      if (pt2i<pt2) pt2 = pt2i;
    }
    break;
  case 4  :
    pt2=1.0;
    for (size_t i=m_nin;i<m_nin+m_nout;i++) pt2*=momenta[i].PPerp2();
    pt2=pow(pt2,1.0/m_nout);
    break;
  case 21 :
    if (m_nin+m_nout==4) {
      double t = (momenta[0]-momenta[2]).Abs2();
      double u = (momenta[0]-momenta[3]).Abs2();
      pt2 = 2.*s*t*u/(s*s+t*t+u*u);
    }
    else {
      pt2 = 0.;
      for (size_t i=m_nin;i<m_nin+m_nout;i++) {
	pt2 += ATOOLS::sqr(momenta[i][1])+ATOOLS::sqr(momenta[i][2]);
      }
    }
    break;
  case 63 :
    if ((int)m_nout!=m_maxjetnumber) {
      pt2 = m_scale[stp::as];
    }
    else {
      pt2 = ATOOLS::sqr(momenta[m_nin][1])+ATOOLS::sqr(momenta[m_nin][2]);
      for (size_t i=m_nin+1;i<m_nin+m_nout;++i) {
	if (p_flavours[i].Strong())
	  pt2 =  ATOOLS::Min(pt2,ATOOLS::sqr(momenta[i][1])+ATOOLS::sqr(momenta[i][2]));
      }
    }
    break;
  case 64 :
    pt2 = ATOOLS::sqr(momenta[m_nin][1])+ATOOLS::sqr(momenta[m_nin][2]);
    for (size_t i=m_nin+1;i<m_nin+m_nout;++i) {
      pt2 =  ATOOLS::Min(pt2,ATOOLS::sqr(momenta[i][1])+ATOOLS::sqr(momenta[i][2]));
    }
    break;
  case 65: {
    double ycut = ATOOLS::rpa.gen.Ycut();
    m_scale[stp::fac] = ATOOLS::rpa.gen.FactorizationScaleFactor()*ycut*ATOOLS::sqr(ATOOLS::rpa.gen.Ecms());
    pt2 = m_scale[stp::fac];
      
    double y=2.;
    if (p_selector->Name()=="Combined_Selector") {
      ATOOLS::Jet_Finder * jf = (ATOOLS::Jet_Finder *)
	((ATOOLS::Combined_Selector*)p_selector)->GetSelector("Jetfinder");
      if (jf) y=jf->ActualValue();
      else {
	ATOOLS::msg.Out()<<"WARNING in Process_Base::Scale : "<<std::endl
			 <<"   No jetfinder found, cannot use SCALESCHEME=="<<m_scalescheme<<"."
			 <<" Return s as scale."<<std::endl;
	pt2 = ATOOLS::rpa.gen.FactorizationScaleFactor()*s;
      }
    }
    
    // if highest number of jets
    if ((int)m_nout==m_maxjetnumber) {
      if (y==2.) pt2 = ATOOLS::rpa.gen.FactorizationScaleFactor()*s;
      else pt2 = ATOOLS::rpa.gen.FactorizationScaleFactor()*y*ATOOLS::sqr(ATOOLS::rpa.gen.Ecms());
    }
    else {
      // tanju two scale treatment
      pt2 = ATOOLS::rpa.gen.FactorizationScaleFactor()*y*ATOOLS::sqr(ATOOLS::rpa.gen.Ecms());
      if (m_scale[stp::fac]<pt2) pt2 = m_scale[stp::fac];
    }
    break;
  }
  case 101: {// pp->V scheme
    double M2=0.;
    if (m_resonances.size()>0) {
      M2=ATOOLS::sqr(m_resonances[0].Mass());
    }
    ATOOLS::Vec4D *p=p_momenta;
    double S2=p[4]*p[5], x1=p[5]*p[0]/S2, x2=p[4]*p[1]/S2;
    double xi=(p[0]+p[1]).PMinus()/(p[0]+p[1]).PPlus();
    m_scale[PHASIC::stp::kp21]=x1*x1*2.0*S2*xi;
    m_scale[PHASIC::stp::kp22]=x2*x2*2.0*S2/xi;
    // ew scale a la watt
    double sc=(p[0]+p[1]).PPerp2();
    pt2=m_scale[PHASIC::stp::as]=pow(sc,2./3.)*pow(M2,1./3.);
    break;
  }
  case 102: {// g*g*->qqb scheme
    const ATOOLS::Vec4D *p=momenta;
    ATOOLS::Vec4D p2=p[m_nin], p3=p[m_nin+m_nout-1];
    double S2=p[m_nin+m_nout]*p[m_nin+m_nout+1];
    double a1=p[m_nin+m_nout+1]*p[0]/S2;
    double b2=p[m_nin+m_nout]*p[1]/S2;
    m_scale[PHASIC::stp::kp21]=a1*a1*2.*S2*p2.PMinus()/p2.PPlus();
    m_scale[PHASIC::stp::kp22]=b2*b2*2.*S2*p3.PPlus()/p3.PMinus();
    // qcd scale
    double s=(p[0]+p[1]).Abs2();
    double t=(p[0]-p[2]).Abs2();
    double u=(p[0]-p[3]).Abs2();
    pt2=m_scale[PHASIC::stp::as]=2.*s*t*u/(s*s+t*t+u*u);
    break;
  }
  case 103: {// g*g*->gg scheme
    const ATOOLS::Vec4D *p=momenta;
    ATOOLS::Vec4D p2=p[m_nin], p3=p[m_nin+m_nout-1];
    if (p2.PMinus()/p2.PPlus()>p3.PMinus()/p3.PPlus())
      std::swap<ATOOLS::Vec4D>(p2,p3);
    double S2=p[m_nin+m_nout]*p[m_nin+m_nout+1];
    double a1=p[m_nin+m_nout+1]*p[0]/S2;
    double b2=p[m_nin+m_nout]*p[1]/S2;
    m_scale[stp::kp21]=a1*a1*2.*S2*p2.PMinus()/p2.PPlus();
    m_scale[stp::kp22]=b2*b2*2.*S2*p3.PPlus()/p3.PMinus();
    // mean of factorisation scales
    pt2=m_scale[stp::as]=
      ATOOLS::sqr(sqrt(m_scale[stp::kp21])+sqrt(m_scale[stp::kp22]))/4.0;
    break;
  }
  case 104: {// g*g*->gg scheme
    const ATOOLS::Vec4D *p=momenta;
    ATOOLS::Vec4D p2=p[m_nin], p3=p[m_nin+m_nout-1];
    if (p2.PMinus()/p2.PPlus()>p3.PMinus()/p3.PPlus())
      std::swap<ATOOLS::Vec4D>(p2,p3);
    double S2=p[m_nin+m_nout]*p[m_nin+m_nout+1];
    double a1=p[m_nin+m_nout+1]*p[0]/S2;
    double b2=p[m_nin+m_nout]*p[1]/S2;
    double fac=ATOOLS::rpa.gen.FactorizationScaleFactor();
    m_scale[stp::kp21]=fac*a1*a1*2.*S2*p2.PMinus()/p2.PPlus();
    m_scale[stp::kp22]=fac*b2*b2*2.*S2*p3.PPlus()/p3.PMinus();
    // mean of factorisation scales
    pt2=m_scale[stp::as]=
      sqrt(m_scale[stp::kp21]*m_scale[stp::kp22]);
    break;
  }
  default :
    pt2 = s;
  }
  m_scale[stp::fac]=pt2;
  if (Selected()==NULL) return m_scale[stp::fac];
  return (*Selected()->p_regulator)[m_scale[PHASIC::stp::fac]]; 
}

double Integrable_Base::KFactor(const double scale) 
{
  switch (m_kfactorscheme) {
  case 1  :
    if (m_nstrong>2) {
      return m_rfactor*pow(MODEL::as->AlphaS(scale)/
			   MODEL::as->AlphaS(ATOOLS::sqr(ATOOLS::rpa.gen.Ecms())),m_nstrong-2);
    } 
    else 
      return m_rfactor;
  case 65:
    m_scale[stp::fac]=scale;
    msg_Debugging()<<Name()<<" : "<<m_nstrong<<std::endl;
    msg_Debugging()<<"  Q_F^2 = "<<m_scale[stp::fac]<<std::endl;
    msg_Debugging()<<"  Q_R^2 = "<<m_scale[stp::as]<<std::endl;
    if (m_nstrong>2) {
      return m_rfactor*pow(MODEL::as->AlphaS(m_scale[stp::as])/
			   MODEL::as->AlphaS(ATOOLS::sqr(ATOOLS::rpa.gen.Ecms())),m_nstrong-2);
    } 
    else 
      return m_rfactor;
  case 101:{// pp -> V scheme
    const double CF=4./3.;
    return exp(CF*MODEL::as->AlphaS(scale)*M_PI/2.);
  }
  case 102  :
    return MODEL::as->AlphaS(m_scale[stp::kp21])
      *MODEL::as->AlphaS(m_scale[stp::kp22])
      /ATOOLS::sqr(MODEL::as->AlphaS(ATOOLS::sqr(ATOOLS::rpa.gen.Ecms())));
  default :
    return 1.;
  }
}

void Integrable_Base::SetEvents(const double number)
{
  m_anasum=m_validanasum=0.0;
  m_gmin=1.;
  m_expevents=(long int)(0.6*number*TotalXS()/TriggerEfficiency()+1.0);
  m_dicedevents=m_accevents=0;
#ifdef EVENT_STATISTICS  
  std::vector<double> smem(6,0.0);
  smem[0]=TotalResult();
  smem[1]=m_expevents;
  mem[Name()]=smem;
#endif
}

void Integrable_Base::AddEvent(const double xs,const double validxs,
			       const int ncounts)
{
  m_anasum+=xs;
  m_validanasum+=validxs;
  m_dicedevents+=ncounts; 
  if (xs!=0.0) m_accevents++;
#ifdef EVENT_STATISTICS  
  if (mem[Name()].size()==0) {
    std::vector<double> smem(6,0.0);
    smem[0]=TotalResult();
    smem[1]=m_expevents;
    mem[Name()]=smem;
  }
  mem[Name()][4]+=ncounts;
  mem[Name()][2]+=xs;
  if (xs!=0.0) mem[Name()][3]++;
#endif
}

double Integrable_Base::TriggerEfficiency()
{
  if (p_pshandler!=NULL) {
    Multi_Channel *fsr=p_pshandler->FSRIntegrator();
    return fsr->ValidN()/(double)fsr->N();
  }
  return 1.0;
}

void Integrable_Base::ResetEvents(double gmin)
{
  m_gmin = gmin;
#ifdef EVENT_STATISTICS  
  mem[Name()][5]=m_expevents-mem[Name()][1];
#endif
  if (m_validanasum==0.) return;
  long int expevents=m_expevents;
  m_expevents=(long int)(m_anasum/m_validanasum/gmin*(double)m_expevents+.5);
  m_expevents=(long int)ATOOLS::Min((double)m_expevents,expevents*sqrt(expevents)/gmin);
  if (!(m_expevents>=expevents)) m_expevents=expevents;
#ifdef EVENT_STATISTICS  
  mem[Name()][5]=m_expevents-expevents;
  msg_Tracking()<<Name()<<" -> extra events "<<expevents<<" "
		<<m_expevents<<" "<<mem[Name()][5]<<" "
		<<m_anasum/m_validanasum<<" "<<gmin<<std::endl;
#endif
}

void Integrable_Base::GetGMin(double &g, double &meff)
{
  if (m_validanasum==0.) return;
  double q = m_anasum/m_validanasum;
  double add = q*(double)m_dicedevents;
  add = ATOOLS::Min(add,sqrt((double)m_expevents)*(double)m_dicedevents);
  meff += add;
  if (m_expevents<100) return;
  if (q<1.) if (m_expevents<sqrt(Parent()->m_expevents)) q=1.;
  g = ATOOLS::Min(g,q);
}

void Integrable_Base::SetISRThreshold(const double threshold) 
{
  m_threshold=threshold;
}

void Integrable_Base::CreateMomenta(const size_t n)
{ 
  delete [] p_momenta; 
  p_momenta = new ATOOLS::Vec4D[n]; 
  m_nvector=n; 
}

