#include "N_Parton_CDBG.H"

#include "Message.H"
#include "Random.H"
#include "Color_Integrator.H"
#include "Helicity_Integrator.H"
#include "Data_Reader.H"
#include "STL_Tools.H"
#include <iomanip>

using namespace EXTRAXS;
using namespace ATOOLS;

N_Parton_CDBG::N_Parton_CDBG(const size_t nin,const size_t nout,
			     const std::vector<Flavour> &flavs,
			     const std::vector<std::string> &models):
  m_nin(nin), m_nout(nout), 
  p_colint(NULL), p_helint(NULL),
  m_mode(0), m_tests(0)
{ 
  std::vector<Flavour> nflavs(flavs);
  nflavs[0]=nflavs[0].Bar();
  nflavs[1]=nflavs[1].Bar();
  m_ampl.Construct(nflavs,models);
  m_moms.resize(m_nin+m_nout);
  std::map<Flavour,size_t> fc;
  for (size_t i(2);i<flavs.size();++i) {
    std::map<Flavour,size_t>::iterator fit(fc.find(flavs[i]));
    if (fit==fc.end()) {
      fc[flavs[i]]=0;
      fit=fc.find(flavs[i]);
    }
    ++fit->second;
  }
  m_sf=1.0;
  msg_Debugging()<<METHOD<<"(): Construct symmetry factor {\n";
  msg_Debugging()<<"  Final state:\n";
  for (std::map<Flavour,size_t>::const_iterator fit(fc.begin());
       fit!=fc.end();++fit) {
    msg_Debugging()<<"  "<<std::setw(2)<<fit->second<<" "
		   <<std::setw(15)<<fit->first<<" -> "
		   <<std::setw(12)<<Factorial(fit->second)<<"\n";
    m_sf*=Factorial(fit->second);
  }
  m_fsf=m_sf;
  msg_Debugging()<<"  Initial state:\n";
  for (size_t i(0);i<2;++i) {
    if (flavs[i].IsGluon()) {
      msg_Debugging()<<"     "<<std::setw(15)<<flavs[i]<<" -> 2*8\n";
      m_sf*=2*8;
    }
    else if (flavs[i].IsQuark()) {
      msg_Debugging()<<"     "<<std::setw(15)<<flavs[i]<<" -> 2*3\n";
      m_sf*=2*3;
    }
  }
  msg_Debugging()<<"} -> "<<m_sf<<"\n";
}

N_Parton_CDBG::~N_Parton_CDBG()
{
}

double N_Parton_CDBG::Factorial(const double &n) 
{
  if (n<=0.0) return 1.0;
  return n*Factorial(n-1.0);
}

double N_Parton_CDBG::Differential(const std::vector<Vec4D> &momenta)
{
  for (size_t i(0);i<m_moms.size();++i) {
    m_moms[i]=momenta[i];
    if (i<m_nin) m_moms[i]=-1.0*m_moms[i];
  }
  m_ampl.SetColors(p_colint->I(),p_colint->J());
  m_ampl.SetMomenta(m_moms);
  if (p_helint!=NULL) {
    Complex csum(m_ampl.Evaluate(p_helint->Chiralities()));
    csum*=std::conj(csum)/m_sf;
    return csum.real();
  }
  m_ampl.EvaluateAll();
  double csum(0.0);
  for (size_t j(0);j<m_ampl.Results().size();++j) {
    Complex cres(m_ampl.Results()[j]);
#ifdef DEBUG__BG
    msg_Debugging()<<"A["<<j<<"]"<<m_ampl.Chiralities()[j]
		   <<" = "<<cres<<" -> "<<std::abs(cres)<<"\n";
#endif
    csum+=(cres*std::conj(cres)).real();
  }
  return csum/m_sf;
}

bool N_Parton_CDBG::GaugeTest(std::vector<Vec4D> momenta)
{
  momenta[0]=-1.0*momenta[0];
  momenta[1]=-1.0*momenta[1];
  bool cnt(true);
  while (cnt) {
    m_ampl.SetColors(p_colint->I(),p_colint->J());
    if (!m_ampl.GaugeTest(momenta)) return false;
    for (size_t j(0);j<m_ampl.Results().size();++j)
      if (m_ampl.Results()[j]!=Complex(0.0,0.0)) cnt=false;
    if (cnt) {
      msg_Info()<<METHOD<<"(): Zero result. Redo gauge test."<<std::endl;
      while (!p_colint->GeneratePoint());
    }
  }
  return true;
}

void N_Parton_CDBG::PrintGraphs(const std::string &gpath) const
{
  msg_Info()<<METHOD<<"(): Write diagrams to '"<<gpath<<"'.\n";
  m_ampl.WriteOutGraphs(gpath);
}

