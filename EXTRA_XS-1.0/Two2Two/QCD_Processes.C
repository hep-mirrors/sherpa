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
using namespace AMATOOLS;
using namespace APHYTOOLS;
using namespace AORGTOOLS;

QCD_Processes::QCD_Processes(ISR::ISR_Handler * _isr,BEAM::Beam_Spectra_Handler * _beam,
			     APHYTOOLS::Flavour * _fl,APHYTOOLS::Selector_Data * _seldata,
			     int _scalescheme,int _kfactorscheme,double _scalefactor) : 
  XS_Group(2,2,_fl,_isr,_beam,_seldata,_scalescheme,_kfactorscheme,_scalefactor)
{
  m_name       = std::string("parton parton -> parton parton");
  
  p_xsselector = new XS_Selector();
  p_fl[0] = p_fl[1] = p_fl[2] = p_fl[3] = APHYTOOLS::Flavour(APHYTOOLS::kf::gluon);
  aS = (*as)(sqr(rpa.gen.Ecms()));
  
  Fill4qmodes();
  Fill2q2gmodes();
  Fill4gmodes();
}

void QCD_Processes::Fill4gmodes() {
  p_fl[0] = p_fl[1] = p_fl[2] = p_fl[3] = APHYTOOLS::Flavour(APHYTOOLS::kf::jet);
  XS_Group * gggg = new XS_Group(2,2,p_fl);
  gggg->SetName(std::string(" gg -> gg"));
  p_fl[0] = p_fl[1] = p_fl[2] = p_fl[3] = APHYTOOLS::Flavour(APHYTOOLS::kf::gluon);
  gggg->Add(p_xsselector->GetXS(2,2,p_fl));
  Add(gggg);
}

void QCD_Processes::Fill2q2gmodes() {
  p_fl[0] = p_fl[1] = p_fl[2] = p_fl[3] = APHYTOOLS::Flavour(APHYTOOLS::kf::jet);
  XS_Group * qqbgg   = new XS_Group(2,2,p_fl);
  qqbgg->SetName(std::string(" qqb -> gg"));
  p_fl[2] = p_fl[3] = Flavour(kf::gluon);
  for (int i=1;i<6;i++) {
    p_fl[0] = APHYTOOLS::Flavour(APHYTOOLS::kf::code(i));
    p_fl[1] = p_fl[0].Bar();
    qqbgg->Add(p_xsselector->GetXS(2,2,p_fl));
  }
  Add(qqbgg);

  p_fl[0] = p_fl[1] = p_fl[2] = p_fl[3] = APHYTOOLS::Flavour(APHYTOOLS::kf::jet);
  XS_Group * ggqqb = new XS_Group(2,2,p_fl);
  ggqqb->SetName(std::string(" gg -> qqb "));
  p_fl[0] = p_fl[1] = Flavour(kf::gluon);
  for (int i=1;i<6;i++) {
    p_fl[2] = APHYTOOLS::Flavour(APHYTOOLS::kf::code(i));
    p_fl[3] = p_fl[2].Bar();
    ggqqb->Add(p_xsselector->GetXS(2,2,p_fl));
  }
  Add(ggqqb);

  p_fl[0] = p_fl[1] = p_fl[2] = p_fl[3] = APHYTOOLS::Flavour(APHYTOOLS::kf::jet);
  XS_Group * gqgq = new XS_Group(2,2,p_fl);
  gqgq->SetName(std::string(" qg -> qg "));
  p_fl[0] = p_fl[2] = Flavour(kf::gluon);
  for (int i=1;i<6;i++) {
    p_fl[1] = p_fl[3] = APHYTOOLS::Flavour(APHYTOOLS::kf::code(i));
    gqgq->Add(p_xsselector->GetXS(2,2,p_fl));
  }
  Add(gqgq);

  p_fl[0] = p_fl[1] = p_fl[2] = p_fl[3] = APHYTOOLS::Flavour(APHYTOOLS::kf::jet);
  XS_Group * gqbgqb = new XS_Group(2,2,p_fl);
  gqbgqb->SetName(std::string(" qbg -> qbg "));
  p_fl[0] = p_fl[2] = Flavour(kf::gluon);
  for (int i=1;i<6;i++) {
    p_fl[1] = p_fl[3] = APHYTOOLS::Flavour(APHYTOOLS::kf::code(i)).Bar();
    gqbgqb->Add(p_xsselector->GetXS(2,2,p_fl));
  }
  Add(gqbgqb);
}

void QCD_Processes::Fill4qmodes() {
  p_fl[0] = p_fl[1] = p_fl[2] = p_fl[3] = APHYTOOLS::Flavour(APHYTOOLS::kf::jet);
  XS_Group * qqqq = new XS_Group(2,2,p_fl);
  qqqq->SetName(std::string(" qq -> qq "));
  for (int i=1;i<6;i++) {
    p_fl[0] = p_fl[1] = p_fl[2] = p_fl[3] = APHYTOOLS::Flavour(APHYTOOLS::kf::code(i));
    qqqq->Add(p_xsselector->GetXS(2,2,p_fl));
  }
  Add(qqqq);


  p_fl[0] = p_fl[1] = p_fl[2] = p_fl[3] = APHYTOOLS::Flavour(APHYTOOLS::kf::jet);
  XS_Group * q1q2q1q2 = new XS_Group(2,2,p_fl);
  q1q2q1q2->SetName(std::string(" q1q2 -> q1q2 "));
  for (int i=1;i<5;i++) {
    p_fl[0] = p_fl[2] = APHYTOOLS::Flavour(APHYTOOLS::kf::code(i));
    for (int j=i+1;j<6;j++) {
      p_fl[1] = p_fl[3] = APHYTOOLS::Flavour(APHYTOOLS::kf::code(j));
      q1q2q1q2->Add(p_xsselector->GetXS(2,2,p_fl));
    }
  }
  Add(q1q2q1q2);

  p_fl[0] = p_fl[1] = p_fl[2] = p_fl[3] = APHYTOOLS::Flavour(APHYTOOLS::kf::jet);
  XS_Group * qqbqqb = new XS_Group(2,2,p_fl);
  qqbqqb->SetName(std::string(" qqb -> qqb "));
  for (int i=1;i<6;i++) {
    p_fl[0] = p_fl[2] = APHYTOOLS::Flavour(APHYTOOLS::kf::code(i));
    p_fl[1] = p_fl[3] = p_fl[0].Bar();
    qqbqqb->Add(p_xsselector->GetXS(2,2,p_fl));
  }
  Add(qqbqqb);

  p_fl[0] = p_fl[1] = p_fl[2] = p_fl[3] = APHYTOOLS::Flavour(APHYTOOLS::kf::jet);
  XS_Group * q1q2bq1q2b = new XS_Group(2,2,p_fl);
  q1q2bq1q2b->SetName(std::string(" q1q2b -> q1q2b "));
  for (int i=1;i<5;i++) {
    p_fl[0] = p_fl[2] = APHYTOOLS::Flavour(APHYTOOLS::kf::code(i));
    for (int j=i+1;j<6;j++) {
      p_fl[1] = p_fl[3] = APHYTOOLS::Flavour(APHYTOOLS::kf::code(j)).Bar();
      q1q2bq1q2b->Add(p_xsselector->GetXS(2,2,p_fl));
    }
  }
  Add(q1q2bq1q2b);

  p_fl[0] = p_fl[1] = p_fl[2] = p_fl[3] = APHYTOOLS::Flavour(APHYTOOLS::kf::jet);
  XS_Group * q1q1bq2q2b = new XS_Group(2,2,p_fl);
  q1q1bq2q2b->SetName(std::string(" q1q1b -> q2q2b "));
  for (int i=1;i<5;i++) {
    p_fl[0] = APHYTOOLS::Flavour(APHYTOOLS::kf::code(i));
    p_fl[1] = p_fl[0].Bar();
    for (int j=i+1;j<6;j++) {
      p_fl[2] = APHYTOOLS::Flavour(APHYTOOLS::kf::code(j)).Bar();
      p_fl[3] = p_fl[2].Bar();
      q1q1bq2q2b->Add(p_xsselector->GetXS(2,2,p_fl));
    }
  }
  Add(q1q1bq2q2b);


  p_fl[0] = p_fl[1] = p_fl[2] = p_fl[3] = APHYTOOLS::Flavour(APHYTOOLS::kf::jet);
  XS_Group * qbqbqbqb = new XS_Group(2,2,p_fl);
  qbqbqbqb->SetName(std::string(" qbqb -> qbqb "));
  for (int i=1;i<6;i++) {
    p_fl[0] = p_fl[1] = p_fl[2] = p_fl[3] = APHYTOOLS::Flavour(APHYTOOLS::kf::code(i)).Bar();
    qbqbqbqb->Add(p_xsselector->GetXS(2,2,p_fl));
  }
  Add(qbqbqbqb);
  
  p_fl[0] = p_fl[1] = p_fl[2] = p_fl[3] = APHYTOOLS::Flavour(APHYTOOLS::kf::jet);
  XS_Group * q1bq2bq1bq2b = new XS_Group(2,2,p_fl);
  q1bq2bq1bq2b->SetName(std::string(" q1bq2b -> q1bq2b "));
  for (int i=1;i<5;i++) {
    p_fl[0] = p_fl[2] = APHYTOOLS::Flavour(APHYTOOLS::kf::code(i)).Bar();
    for (int j=i+1;j<6;j++) {
      p_fl[1] = p_fl[3] = APHYTOOLS::Flavour(APHYTOOLS::kf::code(j)).Bar();
      q1bq2bq1bq2b->Add(p_xsselector->GetXS(2,2,p_fl));
    }
  }
  Add(q1bq2bq1bq2b);
}

void QCD_Processes::CreateFSRChannels() {
  p_ps->FSRIntegrator()->DropAllChannels();
  p_ps->FSRIntegrator()->Add(new S1Channel(2,2,p_fl,Flavour(kf::photon)));
  p_ps->FSRIntegrator()->Add(new T1Channel(2,2,p_fl));
  p_ps->FSRIntegrator()->Add(new U1Channel(2,2,p_fl));
}
