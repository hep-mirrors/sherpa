#include "METOOLS/Explicit/Lorentz_Calculator.H"
#include "METOOLS/Explicit/Form_Factor.H"
#include "METOOLS/Currents/C_Scalar.H"
#include "METOOLS/Currents/C_Vector.H"
#include "METOOLS/Explicit/Vertex.H"
#include "MODEL/Main/Single_Vertex.H"
#include "ATOOLS/Org/Exception.H"

using namespace ATOOLS;

namespace METOOLS {

  /*!
    gamma* -> n scalars through a current linear in the scalar momenta,

      J^mu = sum_l c_l p_l^mu,

    the c_l per declared hadron leg from the vertex's Linear_Current_FF. Legs:
    0 the photon, 1..n the scalars; momenta outgoing from the vertex (the
    output current's reversed).
  */
  template <typename SType>
  class VPPPP_Calculator: public Lorentz_Calculator {
  public:
    typedef std::complex<SType> SComplex;
    typedef CVec4<SType>   CVec4Type;
    typedef CScalar<SType> CScalarType;

    VPPPP_Calculator(const Vertex_Key &key): Lorentz_Calculator(key), p_ff(NULL) {}
    std::string Label() const { return "VPPPP"; }

    CObject *Evaluate(const CObject_Vector &jj)
    {
      const std::vector<size_t> &id(p_v->V()->id);
      const size_t nl(id.size()), out(id.back());
      if (out>=nl || jj.size()!=nl-1) return NULL;
      if (p_ff==NULL) p_ff=LinearCurrent();
      const std::vector<Complex> &c(p_ff->Coefficients());
      std::vector<size_t> slot(nl);
      std::vector<Vec4D> p(nl);
      for (size_t i(0);i<nl-1;++i) { slot[id[i]]=i; p[id[i]]=p_v->J(i)->P(); }
      p[out]=-p_v->JC()->P();
      size_t s(0);
      for (size_t i(0);i<jj.size();++i) s|=jj[i]->S();
      if (out==0) {
        SComplex sc(1.,0.);
        for (size_t l(1);l<nl;++l) sc*=(*jj[slot[l]]->Get<CScalarType>())[0];
        CVec4Type j(SComplex(0.,0.),SComplex(0.,0.),SComplex(0.,0.),SComplex(0.,0.));
        for (size_t l(1);l<nl;++l) j+=CVec4Type(p[l])*SComplex(c[l]);
        CVec4Type *r(CVec4Type::New(j*sc));
        r->SetS(s);
        return r;
      }
      const CVec4Type &a(*jj[slot[0]]->Get<CVec4Type>());
      SComplex sc(1.,0.), cont(0.,0.);
      for (size_t l(1);l<nl;++l) {
        if (l!=out) sc*=(*jj[slot[l]]->Get<CScalarType>())[0];
        cont+=SComplex(c[l])*(a*CVec4Type(p[l]));
      }
      CScalarType *r(CScalarType::New(cont*sc));
      r->SetS(s);
      return r;
    }
  private:
    Linear_Current_FF *p_ff;

    //! the vertex's Linear_Current_FF (its only form factor)
    Linear_Current_FF *LinearCurrent() const
    {
      Linear_Current_FF *ff(p_v->FormFactor().empty() ? NULL :
                            dynamic_cast<Linear_Current_FF*>(p_v->FormFactor()[0]));
      if (ff==NULL) THROW(fatal_error,"VPPPP needs a Linear_Current_FF form factor.");
      return ff;
    }
  };

  template class VPPPP_Calculator<double>;

}// end of namespace METOOLS

using namespace METOOLS;

DECLARE_GETTER(VPPPP_Calculator<double>,"DVPPPP",Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,VPPPP_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new VPPPP_Calculator<double>(key); }
void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,VPPPP_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"gamma* -> scalars, current linear in the scalar momenta"; }
