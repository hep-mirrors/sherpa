#include "DIRE/Tools/Parton.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/STL_Tools.H"
#include "ATOOLS/Org/Message.H"

using namespace DIRE;
using namespace ATOOLS;

size_t Parton::s_cnt(0);

Parton::Parton
(Amplitude *const ampl,
 const ATOOLS::Flavour &f,const ATOOLS::Vec4D &p,
 const Color &c,const int h):
  p_ampl(ampl), m_f(f), m_p(p), m_c(c), m_h(h), m_b(0),
  m_id(0), p_in(NULL)
{
  p_out[1]=p_out[0]=NULL;
  ++s_cnt;
}

Parton::~Parton()
{
  --s_cnt;
}

double Parton::GetXB() const
{
  if (m_b==1) return -m_p.PPlus()/rpa->gen.PBeam(0).PPlus();
  if (m_b==2) return -m_p.PMinus()/rpa->gen.PBeam(1).PMinus();
  return 0.0;
}

void Parton::AddWeight(const Parton *s,const double &t,const double &w)
{
  if (w==1.0) return;
  Weight_Map::iterator wit=m_ws.insert(make_pair(s,Weight_Vector())).first;
  double l(wit->second.empty()?1.0:wit->second.back().m_w);
  wit->second.push_back(Weight(t,l*w));
}

double Parton::GetWeight(const double &t) const
{
  if (m_ws.empty()) return 1.0;
  double wgt(1.0);
  for (Weight_Map::const_iterator
	 wit(m_ws.begin());wit!=m_ws.end();++wit) {
    const Weight_Vector &ws(wit->second);
    size_t l(0), r(ws.size()-1), c((l+r)/2);
    double a(ws[c].m_t);
    while (r-l>1) {
      if (t>a) r=c;
      else l=c;
      c=(l+r)/2;
      a=ws[c].m_t;
    }
    if (t<=ws[r].m_t) wgt*=ws[r].m_w;
    else if (t<=ws[l].m_t) wgt*=ws[l].m_w;
  }
  return wgt;
}

void Parton::SetColor(const Color &c)
{
  if (p_out[0]) {
    Color cc(p_out[0]->Col());
    if (cc.m_i==m_c.m_i) cc.m_i=c.m_i;
    if (cc.m_j==m_c.m_j) cc.m_j=c.m_j;
    p_out[0]->SetColor(cc);
  }
  m_c=c;
}

namespace DIRE {

  std::ostream &operator<<(std::ostream &s,const Parton &p)
  {
    std::string hist;
    if (p.Beam()) hist+=ToString(p.Beam())+" ";
    if (p.In()) hist+=ToString(p.In()->Id())+"->";
    if (p.In()||p.Out(0)) hist+=ToString(p.Id());
    if (p.Out(0)) hist+="->"+ToString(p.Out(0)->Id());
    if (p.Out(1)) hist+=","+ToString(p.Out(1)->Id());
    for (Parton::Weight_Map::const_iterator 
	   wit(p.Weights().begin());wit!=p.Weights().end();++wit) {
      const Parton::Weight_Vector &ws(wit->second);
      hist+=" ["+ToString(wit->first->Id())+"]{"
	+ToString(ws[0].m_t)+":"+ToString(ws[0].m_w);
      for (size_t i(1);i<ws.size();++i)
	hist+=","+ToString(ws[i].m_t)+":"+ToString(ws[i].m_w);
      hist+="}";
    }
    double m(p.Mom().Abs2());
    m=m>=0.0?sqrt(dabs(m)):-sqrt(dabs(m));
    return s<<std::setw(6)<<ToString(p.Id())
	    <<std::right<<std::setw(4)<<p.Flav()
	    <<std::left<<" ["<<p.Hel()<<"]"
	    <<std::setw(10)<<ToString(p.Col())
	    <<p.Mom()<<" "<<m<<" "<<hist;
  }

}
