#include "COMIX/Vertices/Vertex.H"

#include "ATOOLS/Org/Message.H"

namespace COMIX {

  template <typename SType>
  class QCD_GGG: public VVV_Vertex<SType> {
  public:

    typedef std::complex<SType> SComplex;

    typedef ATOOLS::CVec4<SType> CVec4Type;

  private:

    SComplex m_cpl;

    void Evaluate(const CVec4Type &a,const CVec4Type &b);

  protected:

    std::string CVLabel() const;

  public:

    QCD_GGG(const Vertex_Key &key);

  };// end of class QCD_GGG

  DECLARE_VTX_GETTER(QCD_GGG_Getter);

  template <typename SType>
  class QCD_GGT: public VVT_Vertex<SType> {
  public:

    typedef std::complex<SType> SComplex;

    typedef ATOOLS::CVec4<SType> CVec4Type;
    typedef ATOOLS::CAsT4<SType> CAsT4Type;

  private:

    SComplex m_cpl;

    void Evaluate(const CVec4Type &a,const CVec4Type &b);
    void Evaluate(const CVec4Type &a,const CAsT4Type &b);

  protected:

    std::string CVLabel() const;

  public:

    QCD_GGT(const Vertex_Key &key);

  };// end of class QCD_GGT

  DECLARE_VTX_GETTER(QCD_GGT_Getter);

  template <typename SType>
  class QCD_QQG: public FFV_Vertex<SType> {
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

    QCD_QQG(const Vertex_Key &key);

  };// end of class QCD_QQG

  DECLARE_VTX_GETTER(QCD_QQG_Getter);

}// end of namespace COMIX

#include "MODEL/Main/Running_AlphaS.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/STL_Tools.H"

#define ZERO SComplex(0.0,0.0)
#define M_I SComplex(0.0,1.0)

using namespace COMIX;
using namespace ATOOLS;

template <typename SType>
QCD_GGG<SType>::QCD_GGG(const Vertex_Key &key):
  VVV_Vertex<SType>(key,0,1),
  m_cpl(M_I*SComplex(key.p_model->Constant("g_3"))) {}

template <typename SType>
void QCD_GGG<SType>::Evaluate(const CVec4Type &a,const CVec4Type &b)
{
  static SType invsqrttwo(1.0/sqrt(2.0));
  if (a(0)==b(1)) {
    if (a(1)==b(0) && a(0)==a(1)) return;
    CVec4Type j(invsqrttwo*Lorentz(a,b));
    j(0)=b(0);
    j(1)=a(1);
#ifdef DEBUG__BG
    msg_Debugging()<<"'+' "<<a<<"\n";
    msg_Debugging()<<"    "<<b<<"\n";
#endif
    AddJ(m_cpl*SType(this->m_cplfac)*j);
  }
  if (a(1)==b(0)) {
    CVec4Type j(invsqrttwo*Lorentz(a,b));
    j(0)=a(0);
    j(1)=b(1);
#ifdef DEBUG__BG
    msg_Debugging()<<"'-' "<<a<<"\n";
    msg_Debugging()<<"    "<<b<<"\n";
#endif
    AddJ(-m_cpl*SType(this->m_cplfac)*j);
  }
}

template <typename SType>
std::string QCD_GGG<SType>::CVLabel() const
{
  return "("+ToString(m_cpl.real())+",,"+ToString(m_cpl.imag())+")";
}

template <typename SType>
QCD_GGT<SType>::QCD_GGT(const Vertex_Key &key):
  VVT_Vertex<SType>(key,0,1),
  m_cpl(M_I*SComplex(key.p_model->Constant("g_3"))) {}

template <typename SType>
void QCD_GGT<SType>::Evaluate(const CVec4Type &a,const CVec4Type &b)
{
  static SType invsqrttwo(1.0/sqrt(2.0));
  if (a(0)==b(1)) {
    if (a(1)==b(0) && a(0)==a(1)) return;
#ifdef DEBUG__BG
    msg_Debugging()<<"'+' "<<a<<"\n";
    msg_Debugging()<<"    "<<b<<"\n";
#endif
    AddJ(m_cpl*SType(this->m_cplfac)*invsqrttwo*Lorentz(a,b));
  }
  if (a(1)==b(0)) {
#ifdef DEBUG__BG
    msg_Debugging()<<"'-' "<<a<<"\n";
    msg_Debugging()<<"    "<<b<<"\n";
#endif
    AddJ(m_cpl*SType(this->m_cplfac)*invsqrttwo*Lorentz(b,a));
  }
}

template <typename SType>
void QCD_GGT<SType>::Evaluate(const CVec4Type &a,const CAsT4Type &b)
{
  static SType invsqrttwo(1.0/sqrt(2.0));
  if (a(0)==b(1) || a(1)==b(0)) {
    if (a(0)==b(1)) {
      if (a(1)==b(0) && a(0)==a(1)) return;
      CVec4Type j(invsqrttwo*Lorentz(a,b));
      j(0)=b(0);
      j(1)=a(1);
#ifdef DEBUG__BG
      msg_Debugging()<<"'+' "<<a<<"\n";
      msg_Debugging()<<"    "<<b<<"\n";
#endif
      AddJ(m_cpl*SType(this->m_cplfac)*j);
    }
    if (a(1)==b(0)) {
      CVec4Type j(invsqrttwo*Lorentz(a,b));
      j(0)=a(0);
      j(1)=b(1);
#ifdef DEBUG__BG
      msg_Debugging()<<"'-' "<<a<<"\n";
      msg_Debugging()<<"    "<<b<<"\n";
#endif
      AddJ(-m_cpl*SType(this->m_cplfac)*j);
    }
  }
}

template <typename SType>
std::string QCD_GGT<SType>::CVLabel() const
{
  return "("+ToString(m_cpl.real())+",,"+ToString(m_cpl.imag())+")";
}

template <typename SType>
QCD_QQG<SType>::QCD_QQG(const Vertex_Key &key):
  FFV_Vertex<SType>(key,0,1),
  m_cpl(-M_I*SComplex(key.p_model->Constant("g_3"))) {}

template <typename SType>
void QCD_QQG<SType>::Evaluate(const CSpinorType &a,const CSpinorType &b)
{
  static SType invsqrttwo(1.0/sqrt(2.0));
  bool cl(CalcLeft(a,b)), cr(CalcRight(a,b));
  if (!(cl || cr)) return;
  CVec4Type j(ZERO,ZERO,ZERO,ZERO,0,0,a.H(0)+b.H(0),a.H(1)+b.H(1));
  if (cl) j+=LorentzLeft(a,b);
  if (cr) j+=LorentzRight(a,b);
  bool singlet(a()==b());
  if (a.B()<0) {
    j(0)=a();
    j(1)=b();
  }
  else {
    j(0)=b();
    j(1)=a();
  }
  AddJ(j*=invsqrttwo*m_cpl*SType(this->m_cplfac));
  if (singlet) {
    j*=-1.0/3.0;
    for (size_t i(this->s_cimin);i<=this->s_cimax;++i) {
      j(0)=j(1)=i;
      AddJ(j);
    }
  }
}

template <typename SType>
void QCD_QQG<SType>::Evaluate(const CSpinorType &a,const CVec4Type &b)
{
  static SType invsqrttwo(1.0/sqrt(2.0));
  bool singlet(b(0)==b(1) && this->s_cimin<=this->s_cimax);
  bool match((a.B()<0 && a()==b(1)) || (a.B()>0 && a()==b(0)));
  if (match || singlet) {
    bool cl(CalcLeft(a)), cr(CalcRight(a));
    if (!(cl || cr)) return;
    CSpinorType j(a.R(),a.B(),a(),a.H(0)+b.H(0),a.H(1)+b.H(1),0);
    if (cl) j+=LorentzLeft(a,b);
    if (cr) j+=LorentzRight(a,b);
    j*=invsqrttwo*m_cpl*SType(this->m_cplfac);
    if (match) {
      j()=a.B()<0?b(0):b(1);
      AddJ(j);
    }
    if (singlet) {
      j()=a();
      AddJ(-j/3.0);
    }
  }
}

template <typename SType>
std::string QCD_QQG<SType>::CVLabel() const
{
  return "("+ToString(m_cpl.real())+",,"+ToString(m_cpl.imag())+")";
}

template <typename SType,char STag> Vertex_Base *
QCD_GGG_Getter<SType,STag>::operator()(const Vertex_Key &key) const
{
  return new QCD_GGG<SType>(key);
}

template <typename SType,char STag>
void QCD_GGG_Getter<SType,STag>::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"qcd ggg vertex "<<STag;
}

template <typename SType,char STag> Vertex_Base *
QCD_GGT_Getter<SType,STag>::operator()(const Vertex_Key &key) const
{
  return new QCD_GGT<SType>(key);
}

template <typename SType,char STag> 
void QCD_GGT_Getter<SType,STag>::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"qcd ggt vertex "<<STag;
}

template <typename SType,char STag> Vertex_Base *
QCD_QQG_Getter<SType,STag>::operator()(const Vertex_Key &key) const
{
  return new QCD_QQG<SType>(key);
}

template <typename SType,char STag>
void QCD_QQG_Getter<SType,STag>::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"qcd qqg vertex "<<STag;
}

namespace COMIX {

  template class QCD_GGG_Getter<double,'D'>;
  template class QCD_GGT_Getter<double,'D'>;
  template class QCD_QQG_Getter<double,'D'>;

  template class QCD_GGG_Getter<long double,'Q'>;
  template class QCD_GGT_Getter<long double,'Q'>;
  template class QCD_QQG_Getter<long double,'Q'>;

}

DECLARE_TEMPLATE_GETTER(QCD_Vertex_Filler,"QCD_VFill",
			void,Vertex_Filler_Key);

template <typename SType,char STag>
void *QCD_Vertex_Filler<SType,STag>::operator()
  (const Vertex_Filler_Key &key) const
{
  if (!Flavour(kf_gluon).IsOn()) return NULL;
  if (!key.p_model->IncludesModel("QCD")) return NULL;
  std::string gtag("{"+Flavour(kf_gluon).IDName()+"}");
  std::string ttag("{"+Flavour(kf_gluon_qgc).IDName()+"}");
  key.p_vtcs->push_back(new QCD_GGG_Getter<SType,STag>(gtag+gtag+gtag));
  key.p_vtcs->push_back(new QCD_GGT_Getter<SType,STag>(gtag+gtag+ttag));
  key.p_vtcs->push_back(new QCD_GGT_Getter<SType,STag>(gtag+ttag+gtag));
  for (int i(1);i<=6;++i) {
    Flavour f((kf_code)i);
    if (!f.IsOn()) continue;
    std::string qtag("{"+f.IDName()+"}");
    std::string qbtag ("{"+f.Bar().IDName()+"}");
    key.p_vtcs->push_back(new QCD_QQG_Getter<SType,STag>(qtag+qbtag+gtag));
    key.p_vtcs->push_back(new QCD_QQG_Getter<SType,STag>(qbtag+gtag+qbtag));
    key.p_vtcs->push_back(new QCD_QQG_Getter<SType,STag>(gtag+qtag+qtag));
  }
  return NULL;
}

template <typename SType,char STag>
void QCD_Vertex_Filler<SType,STag>::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"qcd vertex filler "<<STag;
}

template class QCD_Vertex_Filler<double,'D'>;

template class QCD_Vertex_Filler<long double,'Q'>;
