#include "QED_Processes.H"
#include "Single_XS.H"
#include "FSR_Channel.H"

#include "MathTools.H"
#include "Run_Parameter.H"

using namespace EXTRAXS;
using namespace PHASIC;
using namespace AORGTOOLS;
using namespace APHYTOOLS;
using namespace AMATOOLS;

QED_Processes::QED_Processes(ISR::ISR_Handler * _isr,BEAM::Beam_Spectra_Handler * _beam,
			     APHYTOOLS::Flavour * _fl,APHYTOOLS::Selector_Data * _seldata,
			     int _scalescheme,int _kfactorscheme,double _scalefactor) : 
  XS_Group(2,2,_fl,_isr,_beam,_seldata,_scalescheme,_kfactorscheme,_scalefactor)
{
  m_name       = std::string("ee -> qqbar");
  p_xsselector = new XS_Selector();
  for (int ifl=1;ifl<6;++ifl) {
    p_fl[2] = APHYTOOLS::Flavour(kf::code(ifl));
    p_fl[3] = APHYTOOLS::Flavour(kf::code(ifl)).Bar();
    Add(p_xsselector->GetXS(m_nin,m_nout,p_fl));
  }
  p_fl[2] = Flavour(kf::jet);
  p_fl[3] = Flavour(kf::jet);
}

void QED_Processes::CreateFSRChannels() {
  p_ps->FSRIntegrator()->DropAllChannels();
  p_ps->FSRIntegrator()->Add(new S1Channel(2,2,p_fl,Flavour(kf::photon)));
  p_ps->FSRIntegrator()->Add(new S1Channel(2,2,p_fl,Flavour(kf::Z)));
}



