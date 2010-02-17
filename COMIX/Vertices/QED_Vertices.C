#include "COMIX/Vertices/Vertex.H"

namespace COMIX {

  template <typename SType>
  class QED_FFP: public FFV_Vertex<SType> {
  public:

    typedef std::complex<SType> SComplex;

    typedef ATOOLS::CSpinor<SType> CSpinorType;
    typedef ATOOLS::CVec4<SType>   CVec4Type;

  private:

    SComplex m_cpl;

    void Evaluate(const CSpinorType &a,const CSpinorType &b);
    void Evaluate(const CSpinorType &a,const CVec4Type &b);

  protected:

    std::string CVLabel() const;

  public:

    QED_FFP(const Vertex_Key &key);

  };// end of class QED_FFP

  DECLARE_VTX_GETTER(QED_FFP_Getter);

}// end of namespace COMIX

#include "MODEL/Main/Running_AlphaS.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/STL_Tools.H"
#include "ATOOLS/Org/Message.H"

#define ZERO SComplex(0.0,0.0)
#define M_I SComplex(0.0,1.0)

using namespace COMIX;
using namespace ATOOLS;

template <typename SType>
QED_FFP<SType>::QED_FFP(const Vertex_Key &key):
  FFV_Vertex<SType>(key,1,0),
  m_cpl(-M_I*SComplex(key.p_model->Constant("g_1"))*
	SType((key.p_a->Flav().IsFermion()?key.p_a:key.p_b)->RFlav().Charge()))
{
  msg_Debugging()<<"m_cpl = "<<m_cpl<<" <- "<<key.ID()<<"\n";
}

template <typename SType>
void QED_FFP<SType>::Evaluate(const CSpinorType &a,const CSpinorType &b)
{
  if (a()!=b()) return;
  bool cl(CalcLeft(a,b)), cr(CalcRight(a,b));
  if (!(cl || cr)) return;
  CVec4Type j(ZERO,ZERO,ZERO,ZERO,0,0,a.H(0)+b.H(0),a.H(1)+b.H(1));
  if (cl) j+=LorentzLeft(a,b);
  if (cr) j+=LorentzRight(a,b);
  AddJ(j*m_cpl);
}

template <typename SType>
void QED_FFP<SType>::Evaluate(const CSpinorType &a,const CVec4Type &b)
{
  bool cl(CalcLeft(a)), cr(CalcRight(a));
  if (!(cl || cr)) return;
  CSpinorType j(a.R(),a.B(),a(),a.H(0)+b.H(0),a.H(1)+b.H(1),0);
  if (cl) j+=LorentzLeft(a,b);
  if (cr) j+=LorentzRight(a,b);
  AddJ(j*m_cpl);
}

template <typename SType>
std::string QED_FFP<SType>::CVLabel() const
{
  return "("+ToString(m_cpl.real())+",,"+ToString(m_cpl.imag())+")";
}

template <typename SType,char STag> Vertex_Base *
QED_FFP_Getter<SType,STag>::operator()(const Vertex_Key &key) const
{
  return new QED_FFP<SType>(key);
}

template <typename SType,char STag>
void QED_FFP_Getter<SType,STag>::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"qed ffp vertex "<<STag;
}

namespace COMIX {

  template class QED_FFP_Getter<double,'D'>;
  template class QED_FFP_Getter<long double,'Q'>;

}

DECLARE_TEMPLATE_GETTER(QED_Vertex_Filler,"QED_VFill",
			void,const Model *);

template <typename SType,char STag>
void *QED_Vertex_Filler<SType,STag>::operator()
  (const Model *const &model) const
{
  if (!Flavour(kf_photon).IsOn()) return NULL;
  if (!model->IncludesModel("QED")) return NULL;
  std::string ptag("{"+Flavour(kf_photon).IDName()+"}");
  for (int i(1);i<=16;++i) {
    if (i==7) i=11;
    Flavour f((kf_code)i);
    if (!f.IsOn() || f.IntCharge()==0) continue;
    std::string ftag("{"+f.IDName()+"}");
    std::string fbtag ("{"+f.Bar().IDName()+"}");
    new QED_FFP_Getter<SType,STag>(ftag+fbtag+ptag);
    new QED_FFP_Getter<SType,STag>(fbtag+ptag+fbtag);
    new QED_FFP_Getter<SType,STag>(ptag+ftag+ftag);
  }
  return NULL;
}

template <typename SType,char STag>
void QED_Vertex_Filler<SType,STag>::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"qed vertex filler "<<STag;
}

template class QED_Vertex_Filler<double,'D'>;

template class QED_Vertex_Filler<long double,'Q'>;

