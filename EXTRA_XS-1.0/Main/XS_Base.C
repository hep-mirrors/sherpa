#include "XS_Base.H"

#include "Phase_Space_Handler.H"
#include "Running_AlphaS.H"
#include "Combined_Selector.H"
#include "Standard_Selector.H"
#include "Run_Parameter.H"
#include "Message.H"

#include <stdio.h>

using namespace EXTRAXS;

XS_Base::XS_Base()
{
  m_name="Empty XS";
} 

XS_Base::XS_Base(const size_t nin,const size_t nout,const ATOOLS::Flavour *flavours,
		 const int scalescheme,const int kfactorscheme,const double scalefactor,
		 BEAM::Beam_Spectra_Handler *const beamhandler,PDF::ISR_Handler *const isrhandler,
		 ATOOLS::Selector_Data *const selectordata):
  Integrable_Base(nin,nout,flavours,scalescheme,kfactorscheme,scalefactor,
		  beamhandler,isrhandler,selectordata),
  p_regulator(Regulator_Base::GetRegulator(this,"Identity",std::vector<double>())),
  p_colours(NULL)
{
  Init(flavours);
  ResetSelector(selectordata);
  p_pshandler = new PHASIC::Phase_Space_Handler(this,isrhandler,beamhandler);
  p_activepshandler=p_pshandler;
}

XS_Base::XS_Base(const size_t nin,const size_t nout,const ATOOLS::Flavour *flavours):
  Integrable_Base(nin,nout,flavours),
  p_regulator(Regulator_Base::GetRegulator(this,"Identity",std::vector<double>())),
  p_colours(NULL)
{
  Init(flavours);
  p_selector = new ATOOLS::No_Selector();
}

XS_Base::~XS_Base() 
{
 if (p_colours!=NULL) { 
    for (size_t i=0;i<m_nin+m_nout;++i) delete p_colours[i];
    delete [] p_colours;
  }
}

void XS_Base::Init(const ATOOLS::Flavour *flavours)
{
  p_flavours = new ATOOLS::Flavour[m_nin+m_nout];
  if (flavours!=NULL) {
    for (size_t i=0;i<m_nin+m_nout;++i) p_flavours[i]=flavours[i];
    m_name=GenerateName(m_nin,m_nout,flavours);
  }
  p_colours = new int*[m_nin+m_nout];
  for (size_t i=0;i<m_nin+m_nout;++i) { 
    p_colours[i] = new int[2]; 
    p_colours[i][0]=p_colours[i][1]=0; 
  }  
  double massin=0.0, massout=0.0;
  for (size_t i=0;i<m_nin;++i) massin+=p_flavours[i].Mass();
  for (size_t i=m_nin;i<m_nout;++i) massout+=p_flavours[i].Mass();
  if (massin>massout) m_threshold=ATOOLS::sqr(massin);
  else m_threshold=ATOOLS::sqr(massout);
}

std::string XS_Base::GenerateName(const size_t nin,const size_t nout,
				  const ATOOLS::Flavour *flavours) 
{
  char help[20];
  std::string name;
  sprintf(help,"%i",nin);
  name=std::string(help);
  name+=std::string("_");
  sprintf(help,"%i",nout);
  name+=std::string(help);
  name+=std::string("__");
  for (size_t i=0;i<nin+nout;) {
    name+=std::string(flavours[i].Name());
    if (flavours[i].IsAnti()) {
      if (flavours[i].Kfcode()==ATOOLS::kf::e ||
	  flavours[i].Kfcode()==ATOOLS::kf::mu ||
	  flavours[i].Kfcode()==ATOOLS::kf::tau ||
	  flavours[i].Kfcode()==ATOOLS::kf::Hmin) {
	name.replace(name.length()-1,1,"+");
      }
      else {
	name+=std::string("b"); 
      }
    }
    name+=std::string("_");
    if (++i==nin) name+=std::string("_");
  }
  return name.substr(0,name.length()-1);
}

double XS_Base::Differential(const ATOOLS::Vec4D *momenta) 
{
  SetMomenta(momenta);
  SetSTU(momenta);
  return Differential(m_s,m_t,m_u);
}

bool XS_Base::SetColours(const ATOOLS::Vec4D *momenta) 
{
  SetMomenta(momenta);
  SetSTU(momenta);
  SetScale(momenta[2].PPerp2(),PHASIC::stp::as);  
  return SetColours(m_s,m_t,m_u);
}

double XS_Base::CalculateScale(const ATOOLS::Vec4D *momenta) 
{
  SetMomenta(momenta);
  if (m_nin==1) return momenta[0].Abs2();
  SetSTU(momenta);
  switch (m_scalescheme) {
  case 1:
    m_scale[PHASIC::stp::as]=momenta[2].PPerp2();
    break;
  case 2:
    m_scale[PHASIC::stp::as]=2.*m_s*m_t*m_u/(m_s*m_s+m_t*m_t+m_u*m_u);
    break;
  case 10: {
    double M2=0.;
    if (m_resonances.size()>0) {
      M2=ATOOLS::sqr(m_resonances[0].Mass());
    }
    ATOOLS::Vec4D *p=p_momenta;
    double S2=p[4]*p[5], x1=p[5]*p[0]/S2, x2=p[4]*p[1]/S2;
    double xi=(p[0]+p[1]).PMinus()/(p[0]+p[1]).PPlus();
    m_scale[PHASIC::stp::kp21]=x1*x1*2.*S2*xi;
    m_scale[PHASIC::stp::kp22]=x2*x2*2.*S2/xi;
    double sc=(p[0]+p[1]).PPerp2();
    m_scale[PHASIC::stp::as]=pow(sc,2./3.)*pow(M2,1./3.);
    break;
  }
  case 21: {// hadron scheme
    const ATOOLS::Vec4D *p=momenta;
    double S2=p[4]*p[5];
    double a1=p[5]*p[0]/S2;
    double b2=p[4]*p[1]/S2;
    m_scale[PHASIC::stp::kp21]=a1*a1*2.*S2*p[2].PMinus()/p[2].PPlus();
    m_scale[PHASIC::stp::kp22]=b2*b2*2.*S2*p[3].PPlus()/p[3].PMinus();
    // average of transverse momenta w.r.t. incoming
    m_scale[PHASIC::stp::as]=ATOOLS::sqr((p[2].PPerp(p[0])+p[3].PPerp(p[1]))/2);
    break;
  }
  case 22: {// hadron scheme
    const ATOOLS::Vec4D *p=momenta;
    double S2=p[4]*p[5];
    double a1=p[5]*p[0]/S2;
    double b2=p[4]*p[1]/S2;
    m_scale[PHASIC::stp::kp21]=a1*a1*2.*S2*p[2].PMinus()/p[2].PPlus();
    m_scale[PHASIC::stp::kp22]=b2*b2*2.*S2*p[3].PPlus()/p[3].PMinus();
    // average of transverse momenta 
    m_scale[PHASIC::stp::as]=ATOOLS::sqr((p[2].PPerp()+p[3].PPerp())/2);
    break;
  }
  case 23: {// dis scheme
    const ATOOLS::Vec4D *p=momenta;
    ATOOLS::Vec4D k=p[0]-p[2];
    double z1=p[5]*k/(p[5]*p[0]);
    double z2=p[4]*k/(p[4]*p[1]);
    m_scale[PHASIC::stp::kp21]=p[2].PPerp2()/ATOOLS::sqr(1.-z1);
    m_scale[PHASIC::stp::kp22]=p[3].PPerp2()/ATOOLS::sqr(1.+z2);
    // average of transverse momenta w.r.t. incoming
    m_scale[PHASIC::stp::as]=ATOOLS::sqr((p[2].PPerp(p[0])+p[3].PPerp(p[1]))/2);
    break;
  }
  case 24: {// dis scheme
    const ATOOLS::Vec4D *p=momenta;
    m_scale[PHASIC::stp::kp21]=p[2].PPerp2();
    m_scale[PHASIC::stp::kp22]=p[3].PPerp2();
    // average of transverse momenta w.r.t. incoming
    m_scale[PHASIC::stp::as]=ATOOLS::sqr((p[2].PPerp(p[0])+p[3].PPerp(p[1]))/2);
    break;
  }
  case 31: {// hadron scheme qq
    double M2=0.;
    if (m_resonances.size()>0) {
      M2=ATOOLS::sqr(m_resonances[0].Mass());
    }
    ATOOLS::Vec4D *p=p_momenta;
    double S2=p[4]*p[5], x1=p[5]*p[0]/S2, x2=p[4]*p[1]/S2;
    double xi=(p[0]+p[1]).PMinus()/(p[0]+p[1]).PPlus();
    m_scale[PHASIC::stp::kp21]=x1*x1*2.*S2*xi;
    m_scale[PHASIC::stp::kp22]=x2*x2*2.*S2/xi;
    double sc=(momenta[0]+momenta[1]).PPerp2();
    m_scale[PHASIC::stp::as]=2.*m_s*m_t*m_u/(m_s*m_s+m_t*m_t+m_u*m_u);
    break;
  }
  case 20:
    m_scale[PHASIC::stp::as]=1.;
    break;
  default:
    m_scale[PHASIC::stp::as]=m_s;
    break;
  }
  return (*p_regulator)[m_scale[PHASIC::stp::as]];
}

double XS_Base::KFactor(const double scale) 
{
  switch (m_kfactorscheme) {
  case 1:
    return pow(MODEL::as->AlphaS(scale*m_scalefactor)/
	       MODEL::as->AlphaS(ATOOLS::sqr(ATOOLS::rpa.gen.Ecms())),m_nin+m_nout-2);
  case 10:{
    const double CF=4./3.;
    return exp(CF*MODEL::as->AlphaS(scale)*M_PI/2.);
  }
  default:
    return 1.;
  }
}

void XS_Base::SwapInOrder() 
{
  std::swap(p_flavours[0],p_flavours[1]);
  std::swap(p_momenta[0],p_momenta[1]);
  std::swap(p_colours[0],p_colours[1]);
  m_swaped=true;
}

void XS_Base::RestoreInOrder() 
{
  if (m_swaped) {
    std::swap(p_flavours[0],p_flavours[1]);
    std::swap(p_momenta[0],p_momenta[1]);
    std::swap(p_colours[0],p_colours[1]);
    m_swaped=false;
  }
}

void XS_Base::ResetSelector(ATOOLS::Selector_Data *const selectordata)
{
  if (p_selector!=NULL) delete p_selector;
  if (selectordata!=NULL) {
    p_selector = new ATOOLS::Combined_Selector(m_nin,m_nout,p_flavours,selectordata);
  }
  else {
    ATOOLS::msg.Error()<<"XS_Base::ResetSelector("<<selectordata<<"): "
		       <<"(\""<<m_name<<"\")"<<std::endl
		       <<"   No cuts specified. Initialize 'No_Selector'."<<std::endl;
    p_selector = new ATOOLS::No_Selector();
  }
}

void XS_Base::SetSTU(const ATOOLS::Vec4D *momenta)
{
  m_s=(momenta[0]+momenta[1]).Abs2();
  m_t=(momenta[0]-momenta[2]).Abs2();
  m_u=(momenta[0]-momenta[3]).Abs2();
}

void XS_Base::SetMax(const double max,const int flag)             
{ 
  if (flag==1) m_max=max;
  SetMax();
}

void XS_Base::SetMax()             
{ 
}

XS_Base *const XS_Base::operator[](const size_t i) const 
{
  return NULL;
}

size_t XS_Base::Size() const
{ 
  return 0; 
}

void XS_Base::SelectOne()
{ 
}

void XS_Base::DeSelect()
{ 
}

void XS_Base::Reset()
{
  m_n=0;
  m_last=m_lastlumi=m_lastdxs=0.0;
  m_totalxs=m_totalsum=m_totalsumsqr=m_totalerr=0.0;
}

ATOOLS::Blob_Data_Base *XS_Base::SameWeightedEvent()
{
  return p_activepshandler->SameWeightedEvent();
}

void XS_Base::AssignRegulator(const std::string &regulator,
			      const std::vector<double> &parameters)
{
  Regulator_Base *function=NULL;
  if ((function=Regulator_Base::GetRegulator(this,regulator,parameters))!=NULL) {
    delete p_regulator;
    p_regulator=function;
  }
}
