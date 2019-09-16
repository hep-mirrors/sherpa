#include "METOOLS/Explicit/Lorentz_Calculator.H"
#include "METOOLS/Currents/C_Scalar.H"
#include "METOOLS/Currents/C_Vector.H"
#include "METOOLS/Explicit/Vertex.H"
#include "MODEL/Main/Single_Vertex.H"

using namespace ATOOLS;

namespace METOOLS {

  template <typename SType>
  class GGV_Calculator: public Lorentz_Calculator {
  public:
    
    typedef CVec4<SType> CVec4Type;
    typedef CScalar<SType> CScalarType;

    GGV_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {}

    std::string Label() const { return "GGV"; }

    CObject *Evaluate(const CObject_Vector &jj)
    {
      // 0 and 1 are outgoing particle, while 2 is the ingoing.
      // const C...Type &a(*jj[0]->Get<C...Type>()) is the
      // polarization of particle a (in a cyclic order).
      /////////////////////
      //       a/	 //
      //       /	 //
      // ------		 //
      // c     \	 //
      //       b\	 //
      /////////////////////
      
      if (p_v->V()->id.back()==0){
	const CScalarType &a(*jj[0]->Get<CScalarType>());
	const CScalarType &b(*jj[1]->Get<CScalarType>());
	Vec4D p12(p_v->J(0)->P()-p_v->J(1)->P());
	CVec4Type *j(CVec4Type::New(CVec4Type(p12)*(a[0]*b[0])));
	j->SetS(a.S()|b.S());
	return j;
      }
      if (p_v->V()->id.back()==1){
      	const CVec4Type &a(*jj[1]->Get<CVec4Type>());
      	const CScalarType &b(*jj[0]->Get<CScalarType>());
      	Vec4D p12(-p_v->J(1)->P()-p_v->J(0)->P()-p_v->J(0)->P());
      	CScalarType *j(CScalarType::New((p12*a)*b[0]));
      	j->SetS(a.S()|b.S());
      	return j;
      }
      if (p_v->V()->id.back()==2){
      	const CVec4Type &a(*jj[0]->Get<CVec4Type>());
      	const CScalarType &b(*jj[1]->Get<CScalarType>());
      	Vec4D p12(p_v->J(0)->P()+p_v->J(0)->P()+p_v->J(1)->P());
      	CScalarType *j(CScalarType::New((p12*a)*b[0]));//
      	j->SetS(a.S()|b.S());
      	return j;
      }
      return NULL;
    }

  };// end of class GGV_Calculator

  template class GGV_Calculator<double>;

}// end of namespace METOOLS

using namespace METOOLS;

DECLARE_GETTER(GGV_Calculator<double>,"DGGV",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,GGV_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new GGV_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    GGV_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"GGV vertex"; }
