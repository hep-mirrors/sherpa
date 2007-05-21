#include "XS_Base.H"

#include "Phase_Space_Handler.H"
#include "Running_AlphaS.H"
#include "Combined_Selector.H"
#include "Standard_Selector.H"
#include "Run_Parameter.H"
#include "Regulator_Base.H"
#include "Message.H"
#include "MyStrStream.H"

#include <stdio.h>

using namespace EXTRAXS;

XS_Base::XS_Base():
  p_colours(NULL)
{
  m_name="Empty XS";
  m_scale[PHASIC::stp::sfs]=m_scale[PHASIC::stp::sis]=std::numeric_limits<double>::max();
} 

XS_Base::XS_Base(const size_t nin,const size_t nout,const ATOOLS::Flavour *flavours,
		 const PHASIC::scl::scheme scalescheme,const int kfactorscheme,
		 BEAM::Beam_Spectra_Handler *const beamhandler,
		 PDF::ISR_Handler *const isrhandler,
		 ATOOLS::Selector_Data *const selectordata):
  Integrable_Base(nin,nout,scalescheme,kfactorscheme,
		  beamhandler,isrhandler,selectordata),
  p_colours(NULL)
{
  Init(flavours);
  ResetSelector(selectordata);
  p_pshandler = new PHASIC::Phase_Space_Handler(this,isrhandler,beamhandler);
  p_activepshandler=p_pshandler;
  m_scale[PHASIC::stp::sfs]=m_scale[PHASIC::stp::sis]=std::numeric_limits<double>::max();
}

XS_Base::XS_Base(const size_t nin,const size_t nout,const ATOOLS::Flavour *flavours):
  Integrable_Base(nin,nout),
  p_colours(NULL)
{
  Init(flavours);
  p_selector = new ATOOLS::No_Selector();
  m_scale[PHASIC::stp::sfs]=m_scale[PHASIC::stp::sis]=std::numeric_limits<double>::max();
}

XS_Base::~XS_Base() 
{
  if (p_colours!=NULL) { 
    for (size_t i=0;i<m_nin+m_nout;++i) delete [] p_colours[i];
    delete [] p_colours;
  }
}

void XS_Base::Init(const ATOOLS::Flavour *flavours)
{
  if (m_nin+m_nout==0) return;
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
  std::string name(ATOOLS::ToString(nin)+"_"+ATOOLS::ToString(nout));
  for (size_t i(0);i<nin;++i)
    name+="_"+std::string(flavours[i].IDName());
  name+="_";
  for (size_t i(nin);i<nin+nout;++i)
    name+="_"+std::string(flavours[i].IDName());
  return name;
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
  SetScale(momenta[2].PPerp2(),PHASIC::stp::ren);  
  return SetColours(m_s,m_t,m_u);
}

void XS_Base::SwapInOrder() 
{
  std::swap(p_flavours[0],p_flavours[1]);
  std::swap(p_momenta[0],p_momenta[1]);
  std::swap(p_colours[0],p_colours[1]);
  if (m_naddout>0) {
    std::swap(p_addmomenta[0],p_addmomenta[1]);
  }
  m_swaped=true;
}

void XS_Base::RestoreInOrder() 
{
  if (m_swaped) {
    std::swap(p_flavours[0],p_flavours[1]);
    std::swap(p_momenta[0],p_momenta[1]);
    std::swap(p_colours[0],p_colours[1]);
    if (m_naddout>0) {
      std::swap(p_addmomenta[0],p_addmomenta[1]);
    }
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
  p_selector->SetProcessName(Name());
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

bool XS_Base::SelectOne()
{ 
  return true;
}

bool XS_Base::ReSelect(int)
{ 
  return true;
}

bool XS_Base::SelectOneFromList()
{ 
  return true;
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
  PHASIC::Regulator_Base *function=NULL;
  if ((function=PHASIC::Regulator_Base::GetRegulator(this,regulator,parameters))!=NULL) {
    delete p_regulator;
    p_regulator=function;
  }
}

void XS_Base::Print()
{
  ATOOLS::msg.Out()<<m_name<<" {"<<m_colorscheme<<","<<m_helicityscheme
		   <<"} ("<<m_orderEW<<","<<m_orderQCD
		   <<")  ->  "<<m_totalxs*ATOOLS::rpa.Picobarn()<<" pb\n";
}
