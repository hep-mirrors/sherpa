#include "METOOLS/Explicit/Lorentz_Calculator.H"
#include "METOOLS/Currents/C_Scalar.H"
#include "METOOLS/Currents/C_Vector.H"
#include "METOOLS/Explicit/Vertex.H"
#include "MODEL/Main/Single_Vertex.H"
#include "ATOOLS/Phys/Spinor.H"

using namespace ATOOLS;

namespace METOOLS {

  /*!
    Anomalous vector-vector-pseudoscalar vertex,

      L = eps_{mu,nu,al,be} V1^mu V2^nu p1^al p2^be

    i.e. the Wess-Zumino structure behind omega -> rho pi, gamma -> rho pi and
    pi0 -> gamma gamma.
  */
  template <typename SType>
  class AVVP_Calculator: public Lorentz_Calculator {
  public:

    typedef std::complex<SType> SComplex;
    typedef CVec4<SType>   CVec4Type;
    typedef CScalar<SType> CScalarType;

    AVVP_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {}

    std::string Label() const { return "AVVP"; }

  private:

    static inline unsigned int I0() { return 0; }
    static inline unsigned int I1() { return Spinor<SType>::R1(); }
    static inline unsigned int I2() { return Spinor<SType>::R2(); }
    static inline unsigned int I3() { return Spinor<SType>::R3(); }

    template <class A, class B, class C>
    static SComplex Minor(const A &r0,const B &r1,const C &r2,
			  const unsigned int &c0,const unsigned int &c1,
			  const unsigned int &c2) {
      return ( r0[c0]*(r1[c1]*r2[c2]-r1[c2]*r2[c1])
	      -r0[c1]*(r1[c0]*r2[c2]-r1[c2]*r2[c0])
	      +r0[c2]*(r1[c0]*r2[c1]-r1[c1]*r2[c0]) );
    }

    // Returns J with  J.c = det(c,b,p,q)  for the Minkowski product, i.e.
    // J^mu = eps^{mu,nu,al,be} b_nu p_al q_be up to the global sign above.
    template <class B, class P, class Q>
    static CVec4Type EpsVec(const B &b,const P &p,const Q &q) {
      CVec4Type j;
      j[I0()] =  Minor(b,p,q,I1(),I2(),I3());
      j[I1()] =  Minor(b,p,q,I0(),I2(),I3());
      j[I2()] = -Minor(b,p,q,I0(),I1(),I3());
      j[I3()] =  Minor(b,p,q,I0(),I1(),I2());
      return j;
    }

    // eps_{mu,nu,al,be} a^mu b^nu p^al q^be
    template <class A, class B, class P, class Q>
    static SComplex Eps4(const A &a,const B &b,const P &p,const Q &q) {
      return ( a[I0()]*Minor(b,p,q,I1(),I2(),I3())
	      -a[I1()]*Minor(b,p,q,I0(),I2(),I3())
	      +a[I2()]*Minor(b,p,q,I0(),I1(),I3())
	      -a[I3()]*Minor(b,p,q,I0(),I1(),I2()) );
    }

  public:

    CObject *Evaluate(const CObject_Vector &jj)
    {
      // Outgoing leg 2: the pseudoscalar, built from the two vectors.
      if (p_v->V()->id.back()==2) {
	const CVec4Type &a(*jj[0]->Get<CVec4Type>());
	const CVec4Type &b(*jj[1]->Get<CVec4Type>());
	const Vec4D p1(p_v->J(0)->P()), p2(p_v->J(1)->P());
	CScalarType *j(CScalarType::New(Eps4(a,b,p1,p2)));
	j->SetS(a.S()|b.S());
	return j;
      }
      // Outgoing leg 0: the first vector, from the second vector and the
      // pseudoscalar.  p1 is the momentum of the leg being built.
      if (p_v->V()->id.back()==0) {
	const CVec4Type   &b(*jj[0]->Get<CVec4Type>());
	const CScalarType &s(*jj[1]->Get<CScalarType>());
	const Vec4D p2(p_v->J(0)->P());
	const Vec4D p1(-p2-p_v->J(1)->P());
	CVec4Type *j(CVec4Type::New(EpsVec(b,p1,p2)*s[0]));
	j->SetS(b.S()|s.S());
	return j;
      }
      // Outgoing leg 1: the second vector.  Antisymmetry of the two vector
      // slots supplies the relative minus sign.
      if (p_v->V()->id.back()==1) {
	const CVec4Type   &a(*jj[0]->Get<CVec4Type>());
	const CScalarType &s(*jj[1]->Get<CScalarType>());
	const Vec4D p1(p_v->J(0)->P());
	const Vec4D p2(-p1-p_v->J(1)->P());
	CVec4Type *j(CVec4Type::New(EpsVec(a,p1,p2)*(-s[0])));
	j->SetS(a.S()|s.S());
	return j;
      }
      return NULL;
    }

  };// end of class AVVP_Calculator

  template class AVVP_Calculator<double>;

}// end of namespace METOOLS

using namespace METOOLS;

DECLARE_GETTER(AVVP_Calculator<double>,"DAVVP",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,AVVP_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new AVVP_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    AVVP_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"anomalous VVP vertex"; }
