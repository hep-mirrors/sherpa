#include "N_Gluon_CDBG.H"

#include "Message.H"
#include "Random.H"
#include "Running_AlphaS.H"
#include "Run_Parameter.H"
#include "Color_Integrator.H"
#include "Helicity_Integrator.H"
#ifdef PROFILE__all
#include "prof.hh"
#else
#define PROFILE_HERE
#define PROFILE_LOCAL(NAME)
#endif

using namespace EXTRAXS;
using namespace ATOOLS;

N_Gluon_CDBG::N_Gluon_CDBG(const size_t nin,const size_t nout,
			   const std::vector<Flavour> &flavs):
  m_nin(nin), m_nout(nout), 
  p_colint(NULL), p_helint(NULL),
  m_mode(0), m_tests(0)
{ 
  std::vector<int> types(m_nin+m_nout,0);
  m_ampl.Construct(flavs,types);
  m_moms.resize(m_nin+m_nout);
  m_as=(*MODEL::as)(sqr(rpa.gen.Ecms()));
  m_sf=Factorial(m_nout);
}

N_Gluon_CDBG::~N_Gluon_CDBG()
{
}

double N_Gluon_CDBG::Factorial(const double &n) 
{
  if (n<=0.0) return 1.0;
  return n*Factorial(n-1.0);
}

double N_Gluon_CDBG::Differential(const std::vector<Vec4D> &momenta)
{
  PROFILE_HERE;
  for (size_t i(0);i<m_moms.size();++i) {
    m_moms[i]=momenta[i];
    if (i<m_nin) m_moms[i]=-1.0*m_moms[i];
  }
  m_ampl.SetColors(p_colint->I(),p_colint->J());
  m_ampl.SetMomenta(m_moms);
  Complex csum(m_ampl.Evaluate(p_helint->Chiralities()));
  csum*=std::conj(csum);
  csum/=4.0*64.0*m_sf;
  csum*=pow(4.0*M_PI*m_as,m_nout);
  return csum.real();
}

bool N_Gluon_CDBG::GaugeTest(std::vector<Vec4D> momenta)
{
  PROFILE_HERE;
  m_ampl.SetMode(m_mode&2);
  momenta[0]=-1.0*momenta[0];
  momenta[1]=-1.0*momenta[1];
  m_ampl.SetColors(p_colint->I(),p_colint->J());
  if (!m_ampl.GaugeTest(momenta)) return false;
  if (m_tests) {
    if (!m_ampl.CyclicityTest()) return false;
    if (!m_ampl.ReflectionTest()) return false;
    if (!m_ampl.DWITest()) return false;
  }
  return true;
}

