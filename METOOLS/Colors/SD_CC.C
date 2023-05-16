#include "METOOLS/Explicit/Dipole_Color.H"

#include "METOOLS/Explicit/Vertex.H"
#include "ATOOLS/Org/Exception.H"

using namespace ATOOLS;

namespace METOOLS {

  class SD_Calculator: public Dipole_Color {
  private:

    int m_si, m_sjk;

  public:

    inline SD_Calculator(const Vertex_Key &key): 
      Dipole_Color(key)
    {
      m_cpl=sqrt(dabs(key.p_k->Flav().Charge()/
		      key.p_c->Flav().Charge()));
      m_si=key.p_c->Flav().Charge()<0.0?1:0;
      m_sjk=key.p_k->Flav().Charge()<0.0?1:0;
    }

    std::string Label() const
    {
      return "S-D";
    }

    bool Evaluate(const CObject_Vector &j)
    {
      if (j.size()==2) {
      m_ci.clear();
      m_cjk.clear();
      m_ci.push_back(CInfo((*j[0])(0),(*j[0])(1),m_si,0));
      m_cjk.push_back(CInfo((*j[1])(0),(*j[1])(1),m_sjk,0));
      return true;
      }
      m_ci.clear();
      m_cjk.clear();
      m_ci.push_back(CInfo((*j[0])(0)|(*j[1])(0),(*j[0])(1)|(*j[1])(1),m_si,0));
      m_cjk.push_back(CInfo((*j[2])(0),(*j[2])(1),m_sjk,0));
      return p_cc->Evaluate(j);
    }

  };// end of class SD_Calculator

}// end of namespace METOOLS

using namespace METOOLS;
using namespace ATOOLS;

DECLARE_GETTER(SD_Calculator,"S-D",
	       Color_Calculator,Vertex_Key);

Color_Calculator *ATOOLS::Getter
<Color_Calculator,Vertex_Key,SD_Calculator>::
operator()(const Vertex_Key &key) const
{
  return new SD_Calculator(key);
}

void ATOOLS::Getter<Color_Calculator,Vertex_Key,SD_Calculator>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"delta (subtraction)";
}
