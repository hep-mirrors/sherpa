#include "COMIX/Phasespace/PS_Current.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/STL_Tools.H"

using namespace COMIX;
using namespace ATOOLS;

void PS_Current::SetSCC(Current *const scc)
{
  p_scc=scc;
  m_psinfo="";
  m_psinfo=PSInfo();
  if (p_scc) m_psinfo+="_SC"+p_scc->PSInfo();
  if (p_dip) m_psinfo+="_DS"+p_dip->PSInfo();
}

void PS_Current::SetDip(NLO_subevt *const dip)
{
  p_dip=dip;
  m_psinfo="";
  m_psinfo=PSInfo();
  if (p_scc) m_psinfo+="_SC"+p_scc->PSInfo();
  if (p_dip) m_psinfo+="_DS"+p_dip->PSInfo();
}

void PS_Current::ConstructJ(const ATOOLS::Vec4D &p,const int ch,
			    const int cr,const int ca,const int mode)
{
  m_p=p;
  ResetJ();
  PS_Info *i(PS_Info::New(PS_Info(cr,ca)));
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"(): '+' "<<m_id<<" "<<*i<<" "<<m_fl<<"\n";
#endif
  AddJ(i);
}

void PS_Current::SetGauge(const ATOOLS::Vec4D &k)
{
}

void PS_Current::AddPropagator()
{
}

std::string PS_Current::Format(const CObject *c) const
{
  return ToString(*(PS_Info*)c);
}

char PS_Current::Type() const
{
  return 'P';
}   

