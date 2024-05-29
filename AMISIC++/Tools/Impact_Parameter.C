#include "AMISIC++/Perturbative/Sudakov_Argument.H"
#include "AMISIC++/Perturbative/MI_Processes.H"
#include "AMISIC++/Tools/Impact_Parameter.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace AMISIC;
using namespace ATOOLS;


/////////////////////////////////////////////////////////////////////////////
// All equations in this file refer to 
// Sjostrand-van der Zijl, PRD 36 (1987) 2019 or
// Corke-Sjostrand, JHEP 05 (2011) 009.
/////////////////////////////////////////////////////////////////////////////

Impact_Parameter::Impact_Parameter() :
  p_integrator(nullptr), m_test(false), m_ana(true) {}

Impact_Parameter::~Impact_Parameter() {
  if (m_ana) FinishAnalysis();
}

void Impact_Parameter::Initialize(Matter_Overlap * mo,Sudakov_Argument * sud,
				  axis * sbins) {
  p_mo         = mo;
  p_sudarg     = sud;
  p_integrator = ( p_mo->IsDynamic() ?
		   p_sudarg->GetProcesses()->GetIntegrator() : nullptr );
  m_pint.Initialize(p_mo,p_sudarg->GetProcesses(),sbins);
  
  if (!p_mo->IsDynamic() && m_test) Test();
  if (m_ana) InitAnalysis();
}

double Impact_Parameter::operator()(const double & s,const double & scale) {
  //////////////////////////////////////////////////////////////////////////
  // Select impact parameter for either minimum bias-type events - FirstB(s),
  // possibly also setting the scale for the first interaction - or for
  // MPI-type events - CalculateB(s,scale).
  //////////////////////////////////////////////////////////////////////////
  if (scale<0.) return FirstB(s);
  return CalculateB(s,scale);
}

double Impact_Parameter::FirstB(const double & s) {
  ///////////////////////////////////////////////////////////////////////////
  // Select impact parameter for a minimum bias-type event.
  ///////////////////////////////////////////////////////////////////////////
  double wt=-1., b = -1., pt2 = -1.;
  do {
    //////////////////////////////////////////////////////////////////////////
    // For "dynamic" matter overlap, we follow the prescrition of SC.
    // We first select a "dummy" hard kinematics according to the
    // differential cross section at peron level, without overlap, cf.
    // Eq. (SC, 23).  This fixes the kinematics input for the selection
    // of the impact parameter according to the matter overlap.
    //////////////////////////////////////////////////////////////////////////
    if (p_mo->IsDynamic()) {
      pt2 = p_integrator->TrialEvent(s);
      p_mo->FixDynamicRadius(p_integrator->X(0),p_integrator->X(1),pt2,pt2);
      b   = p_mo->SelectB();
      wt  = exp(-(*p_sudarg)(b, pt2));
    }
    //////////////////////////////////////////////////////////////////////////
    // For "static" matter overlap, it is distributed according to the
    // interaction probability, d^2b P_int(b) ~ db^2 P_int(b), Eq (SZ, 24),
    // giving the probability of having at least one (non-diffractive)
    // interaction. 
    //////////////////////////////////////////////////////////////////////////
    else {
      b  = sqrt(ran->Get())*p_mo->Bmax();
      wt = m_pint(s,b);
    }
  } while (wt<ran->Get());
  if (m_ana) BAnalyse(-1.,b);
  return b;
}

double Impact_Parameter::CalculateB(const double & s,const double & scale2) {
  int    trials = 1000;
  double b      = -1.;
  if (p_mo->IsDynamic()) {
  }
  //////////////////////////////////////////////////////////////////////////
  // Select b according to f(b) and accept or reject b with probability
  // given by "factorized Sudakov form factor", Eq. (SZ, 37).
  // Update the relevant quantities to the current c.m. energy.
  //////////////////////////////////////////////////////////////////////////
  else {
    double fc       = m_pint.fc(s);
    double hardpart = (*p_sudarg)(s,scale2), softpart, sudakov;
    do {
      b        = p_mo->SelectB();
      softpart = fc * m_pint.fb(s,b);
      sudakov  = exp(-softpart * hardpart);
    } while ((trials--)>0 && sudakov<ran->Get());
    if (trials<=0) {
      msg_Error()<<METHOD<<" throws warning:\n"
		 <<"   no impact parameter in accordance with Sudakov "
		 <<"from hard = "<<hardpart<<"\n"
		 <<"   Return b = "<<b<<" for pt = "<<sqrt(scale2)
		 <<" without Sudakov argument.\n";
      return -1.;
    }
  }
  if (m_ana) BAnalyse(scale2,b);
  return b;
}

double Impact_Parameter::BEnhancement(const double & s,const double & b) {
  ///////////////////////////////////////////////////////////////////////////
  // For "dynamic" matter overlaps, this is the maximal value it can take
  // for a given impact parameter b, and should be independent of s.
  // In the case of SC-like Gaussian matter overlaps with a kinematics
  // dependent width (i.e. effective radius), this becomes
  // 1/R^2_min exp(-b^2/R^2_max). 
  ///////////////////////////////////////////////////////////////////////////
  if (p_mo->IsDynamic()) return p_mo->MaxValue(b);
  ///////////////////////////////////////////////////////////////////////////
  // For "static" matter overlaps, this is f(b), the enhancement factor,
  // Eq. (SZ, 28), taken from the appropriate look-up table in the
  // Interaction_Probability.
  ///////////////////////////////////////////////////////////////////////////
  return (b<p_mo->Bmax() ? m_pint.fc(s)*m_pint.fb(s,b) : 0.);
}


void Impact_Parameter::InitAnalysis() {
  m_histos[std::string("Hard_tot")]    = new Histogram(0, 0.,  1., 100);
  m_histos[std::string("Soft_tot")]    = new Histogram(0, 0., 10., 100);
  m_histos[std::string("Sud")]         = new Histogram(0, 0.,  1., 100);
  m_histos[std::string("B_init")]      = new Histogram(0, 0., 10., 100);
  m_histos[std::string("B_tot")]       = new Histogram(0, 0., 10., 100);
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

  m_histos[std::string("B_Pint")]      = new Histogram(0, 0.,  10., 100);
  for (size_t i=0;i<100;i++) {
    double b = double(i)*0.1+0.05, s = sqr(rpa->gen.Ecms());
    m_histos[std::string("B_Pint")]->Insert(b,2.*M_PI*b*m_pint(s,b));
  }
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
  if (pt2<0.) m_histos[std::string("B_init")]->Insert(b);
  else {
    m_histos[std::string("B_tot")]->Insert(b);
    if (sqrt(pt2)<25.)       m_histos[std::string("B_25")]->Insert(b);
    else if (sqrt(pt2)<40.)  m_histos[std::string("B_40")]->Insert(b);
    else if (sqrt(pt2)<100.) m_histos[std::string("B_100")]->Insert(b);
  }
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
  Histogram histoOverlap(0,0.,p_mo->Bmax(),100);
  double b(0.), bstep(p_mo->Bmax()/100.);
  while (b<p_mo->Bmax()) {
    histoOverlap.Insert(b+bstep/2.,(*p_mo)(b+bstep/2.));
    b+= bstep;
  }
  histoOverlap.Output("Overlap.dat");

  Histogram histoPInt(0,0.,p_mo->Bmax(),100);
  b = 0.;
  while (b<p_mo->Bmax()) {
    histoPInt.Insert(b+bstep/2.,m_pint(s,b+bstep/2.));
    b+= bstep;
  }
  histoPInt.Output("PInt.dat");

  Histogram histoBWeight(0,0.,p_mo->Bmax(),100);
  b = 0.;
  while (b<p_mo->Bmax()) {
    histoBWeight.Insert(b+bstep/2.,(*this)(s,b+bstep/2.));
    b+= bstep;
  }
  histoBWeight.Output("Enhancement_Factor.dat");
  
  msg_Out()<<METHOD<<" starts testing b selection.\n";
  double ntrials = 2.5e7;
  Histogram histoB(0,0.,1.1*p_mo->Bmax(),100);
  for (long int i=0;double(i)<ntrials;i++) {
    b = p_mo->SelectB();
    histoB.Insert(b);
  }
  histoB.Finalize();
  histoB.Output("B_Distribution.dat");
  THROW(normal_exit,"testing complete");
}

