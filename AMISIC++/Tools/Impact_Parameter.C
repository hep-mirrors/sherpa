#include "AMISIC++/Perturbative/MI_Processes.H"
#include "AMISIC++/Tools/Impact_Parameter.H"
#include "ATOOLS/Math/Random.H"

using namespace AMISIC;
using namespace ATOOLS;


/////////////////////////////////////////////////////////////////////////////////
// All equations in this file refer to 
// Sjostrand-van der Zijl, PRD 36 (1987) 2019.
/////////////////////////////////////////////////////////////////////////////////

Impact_Parameter::Impact_Parameter() : m_test(false), m_ana(false) {}

Impact_Parameter::~Impact_Parameter() {
  if (m_ana) FinishAnalysis();
}

void Impact_Parameter::Initialize(REMNANTS::Remnant_Handler * remnant_handler,
				  MI_Processes * processes) {
  m_pint.Initialize(remnant_handler,processes);
  p_mo      = m_pint.GetOverlap();
  m_bmax    = p_mo->Bmax();
  p_sudakov = processes->GetSudakov();
  if (m_test) Test();
  if (m_ana)  InitAnalysis();
}

double Impact_Parameter::operator()(const double & s,const double & b) {
  /////////////////////////////////////////////////////////////////////////////////
  // This is f(b), the enhancement factor, Eq. (28)
  /////////////////////////////////////////////////////////////////////////////////
  return (b<m_bmax? m_pint.fb(s,b) : 0.);
}

double Impact_Parameter::CalculateB(const double & s,const double & pt2) {
  /////////////////////////////////////////////////////////////////////////////////
  // If no scale is given, select impact parameter for a minimum bias-type event,
  // essentially given by the Matter_Overlap.
  // If this is for the simulation of genuine MB events, the InitMB method in
  // Amisic will repeat the production of trial b's until one is accepted.
  // If in contrast, it is for a rescatter event, the InitRescatter method in
  // Amisic will give it "one shot".
  /////////////////////////////////////////////////////////////////////////////////
  if (pt2<0. || m_pint.expO(s)<=1.e-12) return (m_b = p_mo->SelectB());
  /////////////////////////////////////////////////////////////////////////////////
  // Select b according to f(b) and accept or reject b with probability given by
  // "factorized Sudakov form factor", Eq. (37).
  // Update the relevant quantities to the current c.m. energy.
  /////////////////////////////////////////////////////////////////////////////////
  double fc       = m_pint.fc(s);
  double hardpart = (*p_sudakov)(s,pt2), softpart, sudakov;
  int    trials   = 1000;
  do {
    m_b      = p_mo->SelectB();
    softpart = fc * (*this)(s,m_b);
    sudakov  = exp(-softpart * hardpart);
  } while (sudakov<ran->Get() && (trials--)>0);
  if (trials<=0) {
    msg_Error()<<METHOD<<" throws warning:\n"
	       <<"   no impact parameter in accordance with Sudakov "
	       <<"from hard = "<<hardpart<<"\n"
	       <<"   Return b = "<<m_b<<" for pt = "<<sqrt(pt2)
	       <<" without Sudakov argument.\n";
    return -1.;
  }
  if (m_ana) BAnalyse(pt2,m_b);
  return m_b;
}

//##########################################################################################
//##########################################################################################
//##########################################################################################
//##########################################################################################
//##########################################################################################

void Impact_Parameter::InitAnalysis() {
  m_histos[std::string("B_tot")]       = new Histogram(0, 0.,  2., 100);
  m_histos[std::string("Hard_tot")]    = new Histogram(0, 0.,  1., 100);
  m_histos[std::string("Soft_tot")]    = new Histogram(0, 0., 10., 100);
  m_histos[std::string("Sud")]         = new Histogram(0, 0.,  1., 100);
  m_histos[std::string("B_25")]        = new Histogram(0, 0.,  5.,  10);
  m_histos[std::string("B_40")]        = new Histogram(0, 0.,  5.,  10);
  m_histos[std::string("B_100")]       = new Histogram(0, 0.,  5.,  10);
  m_histos[std::string("Hard_25")]     = new Histogram(0, 0., .05,  50);
  m_histos[std::string("Hard_40")]     = new Histogram(0, 0., .05,  50);
  m_histos[std::string("Hard_100")]    = new Histogram(0, 0., .002, 10);
  m_histos[std::string("Soft_25")]     = new Histogram(0, 0.,  5., 100);
  m_histos[std::string("Soft_40")]     = new Histogram(0, 0.,  5., 100);
  m_histos[std::string("Soft_100")]    = new Histogram(0, 0.,  5., 100);
  m_histos[std::string("Sud_25")]      = new Histogram(0, 0.,  1., 100);
  m_histos[std::string("Sud_40")]      = new Histogram(0, 0.,  1., 100);
  m_histos[std::string("Sud_100")]     = new Histogram(0, 0.,  1., 100);
}

void Impact_Parameter::FinishAnalysis() {
  Histogram * histo;
  std::string name;
  for (std::map<std::string,Histogram *>::iterator
	 hit=m_histos.begin();hit!=m_histos.end();hit++) {
    histo = hit->second;
    name  = std::string("MPI_Analysis/")+hit->first+std::string(".dat");
    histo->Finalize();
    histo->Output(name);
    delete histo;
  }
  m_histos.clear();
}


void Impact_Parameter::BAnalyse(const double & pt2,const double & b) {
  m_histos[std::string("B_tot")]->Insert(b);
  if (sqrt(pt2)<25.)       m_histos[std::string("B_25")]->Insert(b);
  else if (sqrt(pt2)<40.)  m_histos[std::string("B_40")]->Insert(b);
  else if (sqrt(pt2)<100.) m_histos[std::string("B_100")]->Insert(b);
}

void Impact_Parameter::Analyse(const double & pt2,const double & sudakov,
			       const double & softpart, const double & hardpart) {
  m_histos[std::string("Sud")]->Insert(sudakov);
  m_histos[std::string("Hard_tot")]->Insert(hardpart);
  m_histos[std::string("Soft_tot")]->Insert(softpart);
  if (sqrt(pt2)<25.) {
    m_histos[std::string("Sud_25")]->Insert(sudakov);
    m_histos[std::string("Hard_25")]->Insert(hardpart);
    m_histos[std::string("Soft_25")]->Insert(softpart);
  }
  else if (sqrt(pt2)<40.) {
    m_histos[std::string("Sud_40")]->Insert(sudakov);
    m_histos[std::string("Hard_40")]->Insert(hardpart);
    m_histos[std::string("Soft_40")]->Insert(softpart);
  }
  else if (sqrt(pt2)<100) {
    m_histos[std::string("Sud_100")]->Insert(sudakov);
    m_histos[std::string("Hard_100")]->Insert(hardpart);
    m_histos[std::string("Soft_100")]->Insert(softpart);
  }
}

void Impact_Parameter::Test() {
  msg_Out()<<METHOD<<" starts testing enhancement factor.\n";
  double s = sqr(100.);
  Histogram histoOverlap(0,0.,m_bmax,100);
  double b(0.), bstep(m_bmax/100.);
  while (b<m_bmax) {
    histoOverlap.Insert(b+bstep/2.,(*p_mo)(b+bstep/2.));
    b+= bstep;
  }
  histoOverlap.Output("Overlap.dat");

  Histogram histoPInt(0,0.,m_bmax,100);
  b = 0.;
  while (b<m_bmax) {
    histoPInt.Insert(b+bstep/2.,m_pint(s,b+bstep/2.));
    b+= bstep;
  }
  histoPInt.Output("PInt.dat");

  Histogram histoBWeight(0,0.,m_bmax,100);
  b = 0.;
  while (b<m_bmax) {
    histoBWeight.Insert(b+bstep/2.,(*this)(s,b+bstep/2.));
    b+= bstep;
  }
  histoBWeight.Output("Enhancement_Factor.dat");
  
  msg_Out()<<METHOD<<" starts testing b selection.\n";
  double ntrials = 2.5e7;
  Histogram histoB(0,0.,1.1*m_bmax,100);
  for (long int i=0;double(i)<ntrials;i++) {
    b = p_mo->SelectB();
    histoB.Insert(b);
  }
  histoB.Finalize();
  histoB.Output("B_Distribution.dat");
  exit(1);
}

