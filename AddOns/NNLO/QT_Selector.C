#include "QT_Selector.H"

#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"

#define s_ymax std::numeric_limits<double>::max()

using namespace SHNNLO;
using namespace PHASIC;
using namespace ATOOLS;

QT_Selector::QT_Selector(const Selector_Key &key):
  Selector_Base("NNLOqT_Selector"), p_proc(key.p_proc)
{
  int nnj=0;
  for (size_t i(m_nin);i<m_nin+m_nout;++i)
    if (!Flavour(kf_jet).Includes(p_fl[i])) ++nnj;
  m_qtmin=ToType<double>(key.p_read->Interpreter()->Interprete(key[0][0]));
  m_type=m_nout-(p_proc->Info().Has(nlo_type::real)?nnj+1:nnj);
}

bool QT_Selector::Trigger(Selector_List &sl)
{
  Vec4D q;
  for (size_t i(m_nin);i<sl.size();++i)
    if (Flavour(kf_jet).Includes(p_fl[i])) q+=sl[i].Momentum();
  double qt=q.PPerp();
  m_cqtmin=m_qtmin>0.0?m_qtmin:-m_qtmin*(sl[0].Momentum()+sl[1].Momentum()-q).Mass();
  bool trig=(m_type==0 && qt<m_cqtmin) || (m_type==1 && qt>m_cqtmin);
  return 1-m_sel_log->Hit(1-trig);
}

DECLARE_ND_GETTER(QT_Selector,"NNLOqT",Selector_Base,Selector_Key,true);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,QT_Selector>::
operator()(const Selector_Key &key) const
{
  if (key.empty() || key.front().size()<1) THROW(critical_error,"Invalid syntax");
  const Flavour_Vector &fl(key.p_proc->Flavours());
  return new QT_Selector(key);
}

void ATOOLS::Getter<Selector_Base,Selector_Key,QT_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"NNLO selector";
}
