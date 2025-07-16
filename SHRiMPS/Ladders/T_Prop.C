#include "SHRiMPS/Ladders/T_Prop.H"

using namespace SHRIMPS;
using namespace ATOOLS;


T_Prop::T_Prop(const colour_type::code & col,const Vec4D & q,const double & q02) : 
  m_col(col), m_q(q), m_q2(ATOOLS::dabs(m_q.Abs2())), 
  m_qt2(m_q.PPerp2()), m_q02(q02)
{}

std::ostream & SHRIMPS::
operator<<(std::ostream & s, const colour_type::code & colour) {
  if (colour==colour_type::singlet)      s<<" singlet ";
  else if (colour==colour_type::triplet) s<<" triplet ";
  else if (colour==colour_type::octet)   s<<"  octet  ";
  else                                   s<<"   none  ";
  return s;
}

namespace SHRIMPS {

std::ostream & operator<<(std::ostream & s, const T_Prop & tprop) {
  s<<"    | ["<<tprop.Col()<<"]         "
   <<"q = "<<tprop.Q()<<" (qt = "<<sqrt(tprop.QT2())<<", q = "
   <<sqrt(dabs(tprop.Q().Abs2()))<<")"
   <<" and Q0 = "<<sqrt(tprop.Q02())<<" | \n";
  return s;
}

std::ostream & operator<<(std::ostream & s, const TPropList & props) {
  s<<"T propagator list ("<<props.size()<<", "<<(&props)<<"): \n";
  if (props.size()>0) {
    for (TPropList::const_iterator piter=props.begin();
	 piter!=props.end();piter++) s<<(*piter);
  }
  s<<"\n";
  return s;
}

}
