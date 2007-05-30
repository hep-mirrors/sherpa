#include "XS_QCD_CDBG_T.H"

#include "N_Parton_CDBG_T.H"
#include "Phase_Space_Handler.H"

using namespace EXTRAXS;
using namespace PHASIC;
using namespace ATOOLS;

namespace EXTRAXS {

  template <> 
  Single_XS *Single_XS::GetProcess<CDXS_T_pp_np>
  (const size_t nin,const size_t nout,const ATOOLS::Flavour *flavours, 
   const size_t nqed, const size_t nqcd)
  {
    if (nqcd!=nout || nqed!=0) return NULL;
    for (size_t i(0);i<nin+nout;++i)
      if (!flavours[i].IsGluon()&&!flavours[i].IsQuark()) return NULL;
    return new CDXS_T_pp_np(nin,nout,flavours); 
  }
  
}

CDXS_T_pp_np::CDXS_T_pp_np(const size_t &nin,const size_t &nout,
			   const ATOOLS::Flavour *const fl): 
  Single_XS(nin,nout,fl) 
{
  for (short unsigned int i(0);i<4;++i) 
    p_colours[i][0]=p_colours[i][1]=0;
  m_nstrong=m_nin+m_nout;
  std::vector<Flavour> flavs(nin+nout);
  for (size_t i(0);i<nin+nout;++i) flavs[i]=p_flavours[i];
  std::vector<std::string> models(1,"QCD");
  p_bg = new N_Parton_CDBG_T(m_nin,m_nout,flavs,models);
  m_p.resize(m_nvector);
}

CDXS_T_pp_np::~CDXS_T_pp_np()
{
  delete p_bg;
}

bool CDXS_T_pp_np::Tests()
{
  p_activepshandler->InitIntegrators();
  p_activepshandler->TestPoint(p_momenta);
  for (size_t i(0);i<m_nvector;++i) m_p[i]=p_momenta[i];
  return p_bg->GaugeTest(m_p);
}

double CDXS_T_pp_np::operator()(const double s,const double t,const double u)
{
  for (size_t i(0);i<m_nvector;++i) m_p[i]=p_momenta[i];
  m_p[0]=Vec4D(m_p[0][0],0.0,0.0,m_p[0][3]<0.0?-m_p[0][0]:m_p[0][0]);
  m_p[1]=Vec4D(m_p[1][0],0.0,0.0,m_p[1][3]<0.0?-m_p[1][0]:m_p[1][0]);
  return p_bg->Differential(m_p);
}

bool CDXS_T_pp_np::SetColours(const double s,const double t,const double u)
{ 
  return false;
}
