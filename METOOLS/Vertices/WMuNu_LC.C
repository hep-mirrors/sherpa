#include "METOOLS/Explicit/Lorentz_Calculator.H"
#include "METOOLS/Currents/C_Spinor.H"
#include "METOOLS/Currents/C_Vector.H"
#include "METOOLS/Explicit/Vertex.H"
#include "MODEL/Main/Single_Vertex.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

using namespace ATOOLS;

namespace METOOLS {

  template <typename SType>
  class WMuNu_Calculator: public Lorentz_Calculator {
  public:

    typedef std::complex<SType> SComplex;

    typedef CSpinor<SType> CSpinorType;
    typedef CVec4<SType> CVec4Type;

    inline SComplex PPlus(const CVec4<SType> &p) const
    { return p[0]+p[ATOOLS::Spinor<SType>::R3()]; }
    inline SComplex PMinus(const CVec4<SType> &p) const
    { return p[0]-p[ATOOLS::Spinor<SType>::R3()]; }

    inline SComplex PT(const CVec4<SType> &p) const
    { return p[ATOOLS::Spinor<SType>::R1()]+
	SComplex(0.0,1.0)*p[ATOOLS::Spinor<SType>::R2()]; }
    inline SComplex PTC(const CVec4<SType> &p) const
    { return p[ATOOLS::Spinor<SType>::R1()]-
	SComplex(0.0,1.0)*p[ATOOLS::Spinor<SType>::R2()]; }

    WMuNu_Calculator(const Vertex_Key &key):
      Lorentz_Calculator(key) {}

    std::string Label() const { return "WMuNu"; }

    CSpinor<SType> *LorentzLeftRight(const CSpinorType &a,const CVec4Type &b)
    {
      switch (a.B()) {
      case -1: {
#ifdef DEBUG__BG
	msg_Debugging()<<"<|g LR "<<a<<"\n";
	msg_Debugging()<<"       "<<b<<"\n";
#endif
	CSpinorType *j(CSpinorType::New(a.R(),a.B(),0,0,0,a.S()|b.S(),3));
	SComplex jp(PPlus(b)), jm(PMinus(b)), jt(PT(b)), jtc(PTC(b));
	(*j)[0]=(a[2]*jp+a[3]*jt);
	(*j)[1]=(a[2]*jtc+a[3]*jm);
	(*j)[2]=(a[0]*jm-a[1]*jt);
	(*j)[3]=(-a[0]*jtc+a[1]*jp);
	return j;
      }
      case 1: {
#ifdef DEBUG__BG
	msg_Debugging()<<"g|> LR "<<a<<"\n";
	msg_Debugging()<<"       "<<b<<"\n";
#endif
	CSpinorType *j(CSpinorType::New(a.R(),a.B(),0,0,0,a.S()|b.S(),3));
	SComplex jp(PPlus(b)), jm(PMinus(b)), jt(PT(b)), jtc(PTC(b));
	(*j)[0]=(a[2]*jm-a[3]*jtc);
	(*j)[1]=(-a[2]*jt+a[3]*jp);
	(*j)[2]=(a[0]*jp+a[1]*jtc);
	(*j)[3]=(a[0]*jt+a[1]*jm);
	return j;
      }
      }
      return NULL;
    }

    CObject *Evaluate(const CObject_Vector &jj)
    {
      if (p_v->V()->id.back()==2) THROW(fatal_error,"Invalid call");
      const CSpinorType &a(*jj[p_v->V()->id.back()]->Get<CSpinorType>());
      const CVec4Type &b(*jj[1-p_v->V()->id.back()]->Get<CVec4Type>());
#ifdef DEBUG__BG
      msg_Debugging()<<"<> LR "<<a<<"\n";
      msg_Debugging()<<"      "<<b<<"\n";
#endif
      CSpinorType *j(LorentzLeftRight(a,b));
      return j;
    }

  };// end of class WMuNu_Calculator

  template class WMuNu_Calculator<double>;

}// end of namespace METOOLS

using namespace METOOLS;

DECLARE_GETTER(WMuNu_Calculator<double>,"DW^{\\mu\\nu}",
	       Lorentz_Calculator,Vertex_Key);
Lorentz_Calculator *ATOOLS::Getter
<Lorentz_Calculator,Vertex_Key,WMuNu_Calculator<double> >::
operator()(const Vertex_Key &key) const
{ return new WMuNu_Calculator<double>(key); }

void ATOOLS::Getter<Lorentz_Calculator,Vertex_Key,
		    WMuNu_Calculator<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"WMuNu vertex"; }
