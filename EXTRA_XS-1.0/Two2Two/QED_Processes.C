#include "QED_Processes.H"
#include "Single_XS.H"
#include "FSR_Channel.H"

#include "MathTools.H"
#include "Run_Parameter.H"
#include "XS_Selector.H"
#include "Phase_Space_Handler.H"

using namespace EXTRAXS;
using namespace PHASIC;
using namespace ATOOLS;


QED_Processes::QED_Processes(PDF::ISR_Handler *const isr,BEAM::Beam_Spectra_Handler *const beam,
			     const ATOOLS::Flavour *fl,ATOOLS::Selector_Data *const seldata,
			     const int scalescheme,const int kfactorscheme,const double scalefactor): 
  XS_Group(2,2,fl,scalescheme,kfactorscheme,scalefactor,beam,isr,seldata)
{
  m_name       = std::string("ee -> qqbar");
  for (int ifl=1;ifl<6;++ifl) {
    p_flavours[2] = ATOOLS::Flavour(kf::code(ifl));
    p_flavours[3] = ATOOLS::Flavour(kf::code(ifl)).Bar();
    p_xsselector->SetOffShell(p_isrhandler->KMROn());
    Add(p_xsselector->GetXS(m_nin,m_nout,p_flavours));
  }
  p_flavours[2] = Flavour(kf::jet);
  p_flavours[3] = Flavour(kf::jet);
}

void QED_Processes::CreateFSRChannels() {
  p_pshandler->FSRIntegrator()->DropAllChannels();
  p_pshandler->FSRIntegrator()->Add(new S1Channel(2,2,p_flavours,Flavour(kf::photon)));
  p_pshandler->FSRIntegrator()->Add(new S1Channel(2,2,p_flavours,Flavour(kf::Z)));
}



