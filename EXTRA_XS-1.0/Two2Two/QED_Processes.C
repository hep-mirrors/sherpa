#include "QED_Processes.H"
#include "Running_AlphaS.H"
#include "Running_AlphaQED.H"
#include "MathTools.H"
#include "Single_XS.H"
#include "XS_Selector.H"
#include "Process_Group.H"

#include "Run_Parameter.H"

using namespace AMEGIC;
using namespace EXTRAXS;
using namespace AORGTOOLS;
using namespace APHYTOOLS;
using namespace AMATOOLS;

QED_Processes::QED_Processes(int initflag) : 
  XS_Group(2,2,std::string(" e+ + e- -> q + qbar "))  
{
  xsselector = new XS_Selector();

  if ((rpa.gen.Beam1() == Flavour(kf::e)) &&
      (rpa.gen.Beam2() == Flavour(kf::e).Bar()) ) {
    fl[0] = APHYTOOLS::Flavour(APHYTOOLS::kf::e);
    fl[1] = APHYTOOLS::Flavour(APHYTOOLS::kf::e).Bar();
  }
  else if ((rpa.gen.Beam1() == Flavour(kf::e).Bar()) &&
	   (rpa.gen.Beam2() == Flavour(kf::e)) ) {
    fl[0] = APHYTOOLS::Flavour(APHYTOOLS::kf::e).Bar();
    fl[1] = APHYTOOLS::Flavour(APHYTOOLS::kf::e);
  }
  else {
    msg.Error()<<"ERROR in EXTRAXS::QED_Processes() : "<<std::endl
	       <<"   Mismatch of flavours : "
	       <<rpa.gen.Beam1()<<" and  "<<rpa.gen.Beam2()<<std::endl;
  }
}

void QED_Processes::Initialize(ISR::ISR_Handler * isr, BEAM::Beam_Handler * beam,
			       APHYTOOLS::Selector_Data * _seldata, AMEGIC::Process_Group * _broker)
{
  fl[nin+0] = APHYTOOLS::Flavour(kf::quark);
  fl[nin+1] = APHYTOOLS::Flavour(kf::quark);

  MakeBroker(isr, beam, _seldata, _broker);

  for (int ifl=1;ifl<6;++ifl) {
    fl[nin+0] = APHYTOOLS::Flavour(ifl);
    fl[nin+1] = APHYTOOLS::Flavour(ifl).Bar();
    Add(xsselector->GetXS(nin,nout,fl),true);
  }
  SetISRTypes(fl);
}




