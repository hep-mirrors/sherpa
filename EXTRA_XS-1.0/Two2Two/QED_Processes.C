#include "QED_Processes.H"
#include "Running_AlphaS.H"
#include "Running_AlphaQED.H"
#include "MathTools.H"
#include "Single_XS.H"
#include "XS_Selector.H"

#include "Run_Parameter.H"

using namespace EXTRAXS;
using namespace AORGTOOLS;
using namespace APHYTOOLS;
using namespace AMATOOLS;

QED_Processes::QED_Processes() : 
  XS_Group(2,2,std::string(" e+ + e- -> q + qbar "))  
{
  xsselector = new XS_Selector();

  Init(2,2,0);

  if ((rpa.gen.Beam1() == Flavour(kf::e)) &&
      (rpa.gen.Beam2() == Flavour(kf::e).bar()) ) {
    fl[0] = APHYTOOLS::Flavour(APHYTOOLS::kf::e);
    fl[1] = APHYTOOLS::Flavour(APHYTOOLS::kf::e).bar();
  }
  else if ((rpa.gen.Beam1() == Flavour(kf::e).bar()) &&
	   (rpa.gen.Beam2() == Flavour(kf::e)) ) {
    fl[0] = APHYTOOLS::Flavour(APHYTOOLS::kf::e).bar();
    fl[1] = APHYTOOLS::Flavour(APHYTOOLS::kf::e);
  }
  else {
    msg.Error()<<"ERROR in EXTRAXS::QED_Processes() : "<<std::endl
	       <<"   Mismatch of flavours : "
	       <<rpa.gen.Beam1()<<" and  "<<rpa.gen.Beam2()<<std::endl;
  }

  for (int ifl=1;ifl<6;++ifl) {
    fl[nin+0] = APHYTOOLS::Flavour(ifl);
    fl[nin+1] = APHYTOOLS::Flavour(ifl).bar();
    Add(xsselector->GetXS(nin,nout,fl) );
  }

  fl[0]  = Flavour(kf::e);
  fl[1]  = Flavour(kf::e).bar();
  fl[2]  = Flavour(kf::u);
  fl[3]  = Flavour(kf::u).bar();

  isr_types.push_back(0);
  isr_masses.push_back(0.);
  isr_widths.push_back(0.);

  isr_types.push_back(3);
  isr_masses.push_back(0.);
  isr_widths.push_back(0.);

  if (APHYTOOLS::Flavour(APHYTOOLS::kf::Z).ison()) {
    isr_types.push_back(1);
    isr_masses.push_back(Flavour(APHYTOOLS::kf::Z).mass());
    isr_widths.push_back(Flavour(APHYTOOLS::kf::Z).width());
  }

  CreateSelector();
}

void QED_Processes::CreateSelector() 
{
  msg.Tracking()<<"In QED_Processes::CreateSelector() :"<<std::endl;
  Data_Read dr(rpa.GetPath()+std::string("/ISR.dat"));

  taumin = dr.GetValue<double>("SMIN");
  taumax = dr.GetValue<double>("SMAX");

  sel = new No_Selector();
  msg.Tracking()<<" s-range : "<<taumin<<" ... "<<taumax<<std::endl;
}




