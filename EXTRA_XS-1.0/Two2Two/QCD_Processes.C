#include "QCD_Processes.H"
#include "Single_XS.H"
#include "XS_Selector.H"

#include "Running_AlphaS.H"
#include "Jet_Finder.H"

#include "Run_Parameter.H"

#include "MathTools.H"

using namespace EXTRAXS;
using namespace AMATOOLS;
using namespace APHYTOOLS;
using namespace AORGTOOLS;

QCD_Processes::QCD_Processes() 
  : XS_Group(2,2,std::string(" parton parton -> parton parton"))  
{
  xsselector = new XS_Selector();

  Init(2,2,0);
  fl[0] = fl[1] = fl[2] = fl[3] = APHYTOOLS::Flavour(APHYTOOLS::kf::gluon);
  aS = (*as)(sqr(rpa.gen.Ecms()));

  Fill4qmodes();
  Fill2q2gmodes();
  Fill4gmodes();

  CreateSelector();
}

void QCD_Processes::Fill4gmodes() {
  XS_Group * gggg;
  gggg = new XS_Group(2,2,std::string("gg -> gg"));
  fl[0] = fl[1] = fl[2] = fl[3] = APHYTOOLS::Flavour(APHYTOOLS::kf::gluon);
  gggg->Add(xsselector->GetXS(2,2,fl));
  Add(gggg);
}

void QCD_Processes::Fill2q2gmodes() {
  XS_Group * qqbgg;
  qqbgg = new XS_Group(2,2,std::string("qqb -> gg"));
  fl[2] = fl[3] = APHYTOOLS::Flavour(APHYTOOLS::kf::gluon);
  for (int i=1;i<6;i++) {
    fl[0] = APHYTOOLS::Flavour(APHYTOOLS::kf::code(i));
    fl[1] = fl[0].bar();
    qqbgg->Add(xsselector->GetXS(2,2,fl));
  }
  Add(qqbgg);

  XS_Group * ggqqb;
  ggqqb = new XS_Group(2,2,std::string("gg -> qqb"));
  fl[0] = fl[1] = APHYTOOLS::Flavour(APHYTOOLS::kf::gluon);
  for (int i=1;i<6;i++) {
    fl[2] = APHYTOOLS::Flavour(APHYTOOLS::kf::code(i));
    fl[3] = fl[2].bar();
    ggqqb->Add(xsselector->GetXS(2,2,fl));
  }
  Add(ggqqb);

  XS_Group * gqgq;
  gqgq = new XS_Group(2,2,std::string("gq -> gq"));
  fl[0] = fl[2] = APHYTOOLS::Flavour(APHYTOOLS::kf::gluon);
  for (int i=1;i<6;i++) {
    fl[1] = fl[3] = APHYTOOLS::Flavour(APHYTOOLS::kf::code(i));
    gqgq->Add(xsselector->GetXS(2,2,fl));
  }
  Add(gqgq);

  XS_Group * gqbgqb;
  gqbgqb = new XS_Group(2,2,std::string("gqb -> gqb"));
  fl[0] = fl[2] = APHYTOOLS::Flavour(APHYTOOLS::kf::gluon);
  for (int i=1;i<6;i++) {
    fl[1] = fl[3] = APHYTOOLS::Flavour(APHYTOOLS::kf::code(i)).bar();
    gqbgqb->Add(xsselector->GetXS(2,2,fl));
  }
  Add(gqbgqb);
}

void QCD_Processes::Fill4qmodes() {
  XS_Group * qqqq;
  qqqq = new XS_Group(2,2,std::string("qq -> qq"));
  for (int i=1;i<6;i++) {
    fl[0] = fl[1] = fl[2] = fl[3] = APHYTOOLS::Flavour(APHYTOOLS::kf::code(i));
    qqqq->Add(xsselector->GetXS(2,2,fl));
  }
  Add(qqqq);

  XS_Group * q1q2q1q2;
  q1q2q1q2 = new XS_Group(2,2,std::string("q1q2 -> q1q2"));
  for (int i=1;i<5;i++) {
    fl[0] = fl[2] = APHYTOOLS::Flavour(APHYTOOLS::kf::code(i));
    for (int j=i+1;j<6;j++) {
      fl[1] = fl[3] = APHYTOOLS::Flavour(APHYTOOLS::kf::code(j));
      q1q2q1q2->Add(xsselector->GetXS(2,2,fl));
    }
  }
  Add(q1q2q1q2);

  XS_Group * qqbqqb;
  qqbqqb = new XS_Group(2,2,std::string("qqb -> qqb"));
  for (int i=1;i<6;i++) {
    fl[0] = fl[2] = APHYTOOLS::Flavour(APHYTOOLS::kf::code(i));
    fl[1] = fl[3] = fl[0].bar();
    qqbqqb->Add(xsselector->GetXS(2,2,fl));
  }
  Add(qqbqqb);

  XS_Group * q1q2bq1q2b;
  q1q2bq1q2b = new XS_Group(2,2,std::string("q1q2b -> q1q2b"));
  for (int i=1;i<5;i++) {
    fl[0] = fl[2] = APHYTOOLS::Flavour(APHYTOOLS::kf::code(i));
    for (int j=i+1;j<6;j++) {
      fl[1] = fl[3] = APHYTOOLS::Flavour(APHYTOOLS::kf::code(j)).bar();
      q1q2bq1q2b->Add(xsselector->GetXS(2,2,fl));
    }
  }
  Add(q1q2bq1q2b);

  XS_Group * q1q1bq2q2b;
  q1q2bq1q2b = new XS_Group(2,2,std::string("q1q1b -> q2q2b"));
  for (int i=1;i<5;i++) {
    fl[0] = APHYTOOLS::Flavour(APHYTOOLS::kf::code(i));
    fl[1] = fl[0].bar();
    for (int j=i+1;j<6;j++) {
      fl[2] = APHYTOOLS::Flavour(APHYTOOLS::kf::code(j)).bar();
      fl[3] = fl[2].bar();
      q1q2bq1q2b->Add(xsselector->GetXS(2,2,fl));
    }
  }
  Add(q1q2bq1q2b);


  XS_Group * qbqbqbqb;
  qbqbqbqb = new XS_Group(2,2,std::string("qbqb -> qbqb"));
  for (int i=1;i<6;i++) {
    fl[0] = fl[1] = fl[2] = fl[3] = APHYTOOLS::Flavour(APHYTOOLS::kf::code(i)).bar();
    qbqbqbqb->Add(xsselector->GetXS(2,2,fl));
  }
  Add(qbqbqbqb);

  XS_Group * q1bq2bq1bq2b;
  q1bq2bq1bq2b = new XS_Group(2,2,std::string("q1bq2b -> q1bq2b"));
  for (int i=1;i<5;i++) {
    fl[0] = fl[2] = APHYTOOLS::Flavour(APHYTOOLS::kf::code(i)).bar();
    for (int j=i+1;j<6;j++) {
      fl[1] = fl[3] = APHYTOOLS::Flavour(APHYTOOLS::kf::code(j)).bar();
      q1bq2bq1bq2b->Add(xsselector->GetXS(2,2,fl));
    }
  }
  Add(q1bq2bq1bq2b);
}

void QCD_Processes::CreateSelector() 
{
  msg.Tracking()<<"In QCD_Processes::CreateSelector() :"<<std::endl;
  Data_Read dr(rpa.GetPath()+std::string("/Integration.dat"));

  ycut       = dr.GetValue<double>("YCUT");
  jetfinder  = dr.GetValue<int>("JETFINDER");
  
  // Jets is the only cut on the phase space
  Flavour * dummies;
  dummies    = new Flavour[nin+nout];
  dummies[0] = AORGTOOLS::rpa.gen.Beam1();
  dummies[1] = AORGTOOLS::rpa.gen.Beam2();
  for (int i=nin;i<nin+nout;i++) dummies[i] = Flavour(kf::gluon);

//   sel        = new Jet_Finder(nin+nout,dummies,ycut,jetfinder,4);
//   taumin     = ycut;

  msg.Tracking()<<" tau-range : "<<ycut<<" ... 1."<<std::endl
		<<" pt-range  : "<<sqrt(ycut)*rpa.gen.Ecms()<<" ... "<<rpa.gen.Ecms()<<std::endl
		<<" Initialized jet measure."<<std::endl;
}

double QCD_Processes::Scale(AMATOOLS::vec4d * p)
{
//   s = (p[0]+p[1]).abs2();
//   t = (p[0]-p[2]).abs2();
//   u = (p[0]-p[3]).abs2();
//   return scale = (2.*s*t*u)/(s*s+t*t+u*u);
}


double QCD_Processes::KFactor(double _scale)
{
  return sqr((*as)(_scale)/aS);
}
