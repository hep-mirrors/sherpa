#include "QCD_Processes.H"
#include "Single_XS.H"
#include "XS_Selector.H"
#include "FSR_Channel.H"

#include "Running_AlphaS.H"

#include "Run_Parameter.H"

#include "MathTools.H"

using namespace EXTRAXS;
using namespace PHASIC;
using namespace MODEL;
using namespace ATOOLS;


QCD_Processes::QCD_Processes(PDF::ISR_Handler * _isr,BEAM::Beam_Spectra_Handler * _beam,
			     ATOOLS::Flavour * _fl,ATOOLS::Selector_Data * _seldata,
			     int _scalescheme,int _kfactorscheme,double _scalefactor,bool fillmodes) : 
  XS_Group(2,2,_fl,_isr,_beam,_seldata,_scalescheme,_kfactorscheme,_scalefactor)
{
  SetFSRInterface(NULL);
  SetFSRMode(0);

  m_name       = std::string("parton parton -> parton parton");
  
  p_xsselector = new XS_Selector();
  aS = (*as)(sqr(rpa.gen.Ecms()));
  
  p_fl[0] = p_fl[1] = p_fl[2] = p_fl[3] = ATOOLS::Flavour(ATOOLS::kf::jet);

  if (fillmodes) FillMode(All);
}

void QCD_Processes::FillMode(Mode mode) 
{
  XS_Group *group;
  switch (mode) {
  case All:
  case gggg:
    group = new XS_Group(2,2,p_fl);
    group->SetName(std::string(" gg -> gg"));
    p_fl[0] = p_fl[1] = p_fl[2] = p_fl[3] = ATOOLS::Flavour(ATOOLS::kf::gluon);
    group->Add(p_xsselector->GetXS(2,2,p_fl));
    Add(group);
    if (mode==gggg) break;
  case qqbgg:
    group   = new XS_Group(2,2,p_fl);
    group->SetName(std::string(" qqb -> gg"));
    p_fl[2] = p_fl[3] = Flavour(kf::gluon);
    for (int i=1;i<6;i++) {
      p_fl[0] = ATOOLS::Flavour(ATOOLS::kf::code(i));
      p_fl[1] = p_fl[0].Bar();
      group->Add(p_xsselector->GetXS(2,2,p_fl));
    }
    Add(group);
    if (mode==qqbgg) break;
  case ggqqb:
    group   = new XS_Group(2,2,p_fl);
    group->SetName(std::string(" gg -> qqb"));
    p_fl[0] = p_fl[1] = Flavour(kf::gluon);
    for (int i=1;i<6;i++) {
      p_fl[2] = ATOOLS::Flavour(ATOOLS::kf::code(i));
      p_fl[3] = p_fl[2].Bar();
      group->Add(p_xsselector->GetXS(2,2,p_fl));
    }
    Add(group);
    if (mode==ggqqb) break;
  case qgqg:
    group = new XS_Group(2,2,p_fl);
    group->SetName(std::string(" qg -> qg "));
    p_fl[0] = p_fl[2] = Flavour(kf::gluon);
    for (int i=1;i<6;i++) {
      p_fl[1] = p_fl[3] = ATOOLS::Flavour(ATOOLS::kf::code(i));
      group->Add(p_xsselector->GetXS(2,2,p_fl));
    }
    Add(group);
    if (mode==qgqg) break;
  case q1q2q1q2:
    group = new XS_Group(2,2,p_fl);
    group->SetName(std::string(" q1q2 -> q1q2 "));
    for (int i=1;i<5;i++) {
      p_fl[0] = p_fl[2] = ATOOLS::Flavour(ATOOLS::kf::code(i));
      for (int j=i+1;j<6;j++) {
	p_fl[1] = p_fl[3] = ATOOLS::Flavour(ATOOLS::kf::code(j));
	group->Add(p_xsselector->GetXS(2,2,p_fl));
      }
    }
  case q1q1q1q1:
    if (mode==q1q1q1q1) group = new XS_Group(2,2,p_fl);
    group->SetName(std::string(" q1q1 -> q1q1 "));
    for (int i=1;i<6;i++) {
      p_fl[0] = p_fl[1] = p_fl[2] = p_fl[3] = ATOOLS::Flavour(ATOOLS::kf::code(i));
      group->Add(p_xsselector->GetXS(2,2,p_fl));
    }
    Add(group);
    if (mode==q1q1q1q1) break;
    if (mode==q1q2q1q2) break;
  case q1q1bq1q1b:
    group = new XS_Group(2,2,p_fl);
    group->SetName(std::string(" q1q1b -> q1q1b "));
    for (int i=1;i<6;i++) {
      p_fl[0] = p_fl[2] = ATOOLS::Flavour(ATOOLS::kf::code(i));
      p_fl[1] = p_fl[3] = p_fl[0].Bar();
      group->Add(p_xsselector->GetXS(2,2,p_fl));
    }
    Add(group);
    if (mode==q1q1bq1q1b) break;
  case q1q2bq1q2b:
    group = new XS_Group(2,2,p_fl);
    group->SetName(std::string(" q1q2b -> q1q2b "));
    for (int i=1;i<5;i++) {
      p_fl[0] = p_fl[2] = ATOOLS::Flavour(ATOOLS::kf::code(i));
      for (int j=i+1;j<6;j++) {
	p_fl[1] = p_fl[3] = ATOOLS::Flavour(ATOOLS::kf::code(j)).Bar();
	group->Add(p_xsselector->GetXS(2,2,p_fl));
      }
    }
    Add(group);
    if (mode==q1q2bq1q2b) break;
  case q1q1bq2q2b:
    group = new XS_Group(2,2,p_fl);
    group->SetName(std::string(" q1q1b -> q2q2b "));
    for (int i=1;i<5;i++) {
      p_fl[0] = ATOOLS::Flavour(ATOOLS::kf::code(i));
      p_fl[1] = p_fl[0].Bar();
      for (int j=i+1;j<6;j++) {
	p_fl[2] = ATOOLS::Flavour(ATOOLS::kf::code(j)).Bar();
	p_fl[3] = p_fl[2].Bar();
	group->Add(p_xsselector->GetXS(2,2,p_fl));
      }
    }
    Add(group);
    if (mode==q1q1bq2q2b) break;
    break;
  case Unknown:
    ATOOLS::msg.Error()<<"QCD_Processes::FillMode("<<mode<<"): No mode specified! Abort."<<std::endl;
  case None:
    break;
  }
}

void QCD_Processes::CreateFSRChannels() 
{
  if ((m_fsrmode==0)||(p_fsrinterface==NULL)) {
    p_ps->FSRIntegrator()->DropAllChannels();
    p_ps->FSRIntegrator()->Add(new S1Channel(2,2,p_fl,Flavour(kf::photon)));
    p_ps->FSRIntegrator()->Add(new T1Channel(2,2,p_fl));
    p_ps->FSRIntegrator()->Add(new U1Channel(2,2,p_fl));
    m_fsrmode=1;
  }
  else {
    if (m_fsrmode==2) {
      p_ps->FSRIntegrator()->DropAllChannels();
      p_ps->FSRIntegrator()->Add(p_fsrinterface);
      m_fsrmode=1;
    }
  }
}
