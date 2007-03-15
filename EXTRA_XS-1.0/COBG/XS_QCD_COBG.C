#include "XS_QCD_COBG.H"

#include "N_Gluon_BG.H"
#include "Phase_Space_Handler.H"

using namespace EXTRAXS;
using namespace ATOOLS;

namespace EXTRAXS {

  template <> 
  Single_XS *Single_XS::GetProcess<XS_gg_ng>
  (const size_t nin,const size_t nout,const ATOOLS::Flavour *flavours, 
   const size_t nqed, const size_t nqcd)
  {
    if (nqcd!=nout || nqed!=0) return NULL;
    for (size_t i(0);i<nin+nout;++i)
      if (!flavours[i].IsGluon()) return NULL;
    return new XS_gg_ng(nin,nout,flavours); 
  }
  
}

XS_gg_ng::XS_gg_ng(const size_t &nin,const size_t &nout,
		   const ATOOLS::Flavour *const fl): 
  Single_XS(nin,nout,fl) 
{
  for (short unsigned int i(0);i<4;++i) 
    p_colours[i][0]=p_colours[i][1]=0;
  m_nstrong=m_nin+m_nout;
  std::vector<Flavour> flavs(nin+nout);
  for (size_t i(0);i<nin+nout;++i) flavs[i]=fl[i];
  p_bg = new N_Gluon_BG(m_nin,m_nout,flavs);
  m_p.resize(m_nvector);
}

XS_gg_ng::~XS_gg_ng()
{
  delete p_bg;
}

bool XS_gg_ng::Tests()
{
  p_activepshandler->TestPoint(p_momenta);
  for (size_t i(0);i<m_nvector;++i) m_p[i]=p_momenta[i];
  return p_bg->GaugeTest(m_p);
}

double XS_gg_ng::operator()(const double s,const double t,const double u)
{
  for (size_t i(0);i<m_nvector;++i) m_p[i]=p_momenta[i];
  m_p[0]=Vec4D(m_p[0][0],0.0,0.0,m_p[0][3]<0.0?-m_p[0][0]:m_p[0][0]);
  m_p[1]=Vec4D(m_p[1][0],0.0,0.0,m_p[1][3]<0.0?-m_p[1][0]:m_p[1][0]);
  return p_bg->Differential(m_p);
}

bool XS_gg_ng::SetColours(const double s,const double t,const double u)
{ 
  return false;
}
