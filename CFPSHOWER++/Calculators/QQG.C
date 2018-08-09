#include "CFPSHOWER++/Calculators/Gauge_Base.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "MODEL/Main/Single_Vertex.H"
#include "ATOOLS/Phys/Flow.H"
#include "ATOOLS/Org/Message.H"

namespace CFPSHOWER {
  class QQG : public Gauge_Base {
    int m_mode;
  public:
    QQG(const Kernel_Info & info) : Gauge_Base(info) {
      m_charge = m_CF;
      SetName("3-3-8");
    }
    inline const double Charge(const double & scale) const { return m_charge; }
    bool SetColours(Splitting & split);
  };
}


using namespace CFPSHOWER;
using namespace ATOOLS;

bool QQG::SetColours(Splitting & split) {
  m_colors.clear();
  unsigned int newcol  = Flow::Counter();
  Flavour_Vector flavs = split.GetKernel()->GetFlavs();
  if (!flavs[0].IsAnti()) {
    m_colors.push_back(Color(newcol,0));
    m_colors.push_back(Color(split.GetSplitter()->GetColor()[0],newcol));
  }
  else if (flavs[0].IsAnti()) {
    m_colors.push_back(Color(0,newcol));
    m_colors.push_back(Color(newcol,split.GetSplitter()->GetColor()[1]));
  }
  // set colours in splitting
  for (size_t i=0;i<2;i++) split.SetCol(i,m_colors[i]);
  return (m_colors.size()==2); 
}

DECLARE_GETTER(QQG,"QQG",Gauge_Base,Kernel_Info);

Gauge_Base * ATOOLS::Getter<Gauge_Base,Kernel_Info,QQG>::
operator()(const Parameter_Type & info) const
{
  if (abs(info.GetFlavs()[0].StrongCharge())==3 &&
      abs(info.GetFlavs()[1].StrongCharge())==3 &&
      abs(info.GetFlavs()[2].StrongCharge())==8)
    return new QQG(info);
  return NULL;
}

void ATOOLS::Getter<Gauge_Base,Kernel_Info,QQG>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"gauge part for 3-3-8 in colour space";
}
