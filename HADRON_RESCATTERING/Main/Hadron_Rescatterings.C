#include "HADRON_RESCATTERING/Main/Hadron_Rescatterings.H"

using namespace HADRON_RESCATTERING;
using namespace ATOOLS;

Hadron_Rescatterings::Hadron_Rescatterings(const bool & on) :
  m_on(on), m_nfails(0), m_mm2mb(1.e25),
  p_BB(NULL), p_MM(NULL)
{
  hrpars = new HR_Parameters();
}

void Hadron_Rescatterings::Initialize() {
  p_BM = new BaryonMeson();
  p_BB = new BaryonBaryon();
  p_MM = new MesonMeson();
}

Hadron_Rescatterings::~Hadron_Rescatterings() {
  if (p_BB)   { delete p_BB;   p_BB   = NULL; }
  if (p_MM)   { delete p_MM;   p_MM   = NULL; }
  if (hrpars) { delete hrpars; hrpars = NULL; }
}
    
Blob * Hadron_Rescatterings::operator()(Particle * A,Particle * B,
					const double & dist2) {
  const double s    = (A->Momentum()+B->Momentum()).Abs2();
  double xstot      = 0.;
  const Flavour flA = A->Flav(), flB = B->Flav();
  if (s<sqr(flA.HadMass()+flB.HadMass())) return NULL;
  if (flA.IsMeson() && flB.IsMeson()) {
    xstot = p_MM->XStot(flA,flB,s);
  }
  if (flA.IsBaryon() && flB.IsBaryon()) {
    xstot = p_BB->XStot(flA,flB,s);
  }
  msg_Out()<<METHOD<<"("<<flA<<", "<<flB<<"): "
	   <<"xstot = "<<xstot<<" vs. dist^2 = "<<dist2
	   <<" --> "<<(M_PI*dist2*m_mm2mb)<<"\n\n";
  return NULL;
}

void Hadron_Rescatterings::Reset() {}
