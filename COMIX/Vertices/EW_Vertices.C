#include "COMIX/Vertices/Vertex.H"

namespace COMIX {

  template <typename SType>
  class EW_HHS: public SSS_Vertex<SType> {
  public:

    typedef std::complex<SType> SComplex;

    typedef ATOOLS::CScalar<SType> CScalarType;

  private:

    SComplex m_cpl;

    void Evaluate(const CScalarType &a,const CScalarType &b);

  protected:

    std::string CVLabel() const;

  public:

    EW_HHS(const Vertex_Key &key);

  };// end of class EW_HHS

  DECLARE_VTX_GETTER(EW_HHS_Getter);

  template <typename SType>
  class EW_FFH: public FFS_Vertex<SType> {
  public:

    typedef std::complex<SType> SComplex;

    typedef ATOOLS::CSpinor<SType> CSpinorType;
    typedef ATOOLS::CScalar<SType> CScalarType;

  private:

    SComplex m_cpl;

    void Evaluate(const CSpinorType &a,const CSpinorType &b);
    void Evaluate(const CSpinorType &a,const CScalarType &b);

  protected:

    std::string CVLabel() const;

  public:

    EW_FFH(const Vertex_Key &key);

  };// end of class EW_FFH

  DECLARE_VTX_GETTER(EW_FFH_Getter);

  template <typename SType>
  class EW_FFZ: public FFV_Vertex<SType> {
  public:

    typedef std::complex<SType> SComplex;

    typedef ATOOLS::CSpinor<SType> CSpinorType;
    typedef ATOOLS::CVec4<SType>   CVec4Type;

  private:

    SComplex m_cpll, m_cplr;

    void Evaluate(const CSpinorType &a,const CSpinorType &b);
    void Evaluate(const CSpinorType &a,const CVec4Type &b);

  protected:

    std::string CVLabel() const;

  public:

    EW_FFZ(const Vertex_Key &key);

  };// end of class EW_FFZ

  DECLARE_VTX_GETTER(EW_FFZ_Getter);

  template <typename SType>
  class EW_FFW: public FFV_Vertex<SType> {
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

    EW_FFW(const Vertex_Key &key);

  };// end of class EW_FFW

  DECLARE_VTX_GETTER(EW_FFW_Getter);

  template <typename SType>
  class EW_VVH: public VVS_Vertex<SType> {
  public:

    typedef std::complex<SType> SComplex;

    typedef ATOOLS::CScalar<SType> CScalarType;
    typedef ATOOLS::CVec4<SType>   CVec4Type;

  private:

    SComplex m_cpl;

    void Evaluate(const CVec4Type &a,const CVec4Type &b);
    void Evaluate(const CVec4Type &a,const CScalarType &b);

  protected:

    std::string CVLabel() const;

  public:

    EW_VVH(const Vertex_Key &key);

  };// end of class EW_VVH

  DECLARE_VTX_GETTER(EW_VVH_Getter);

  template <typename SType>
  class EW_WWV: public VVV_Vertex<SType> {
  public:

    typedef std::complex<SType> SComplex;

    typedef ATOOLS::CVec4<SType> CVec4Type;

  private:

    SComplex m_cpl;

    void Evaluate(const CVec4Type &a,const CVec4Type &b);

  protected:

    std::string CVLabel() const;

  public:

    EW_WWV(const Vertex_Key &key);

  };// end of class EW_WWV

  DECLARE_VTX_GETTER(EW_WWV_Getter);

  template <typename SType>
  class EW_WWT: public VVT_Vertex<SType> {
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

    EW_WWT(const Vertex_Key &key);

  };// end of class EW_WWT

  DECLARE_VTX_GETTER(EW_WWT_Getter);

}// end of namespace COMIX

#include "MODEL/Main/Running_AlphaS.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/STL_Tools.H"

#define ZERO SComplex(0.0,0.0)
#define M_I SComplex(0.0,1.0)

using namespace COMIX;
using namespace ATOOLS;

template <typename SType>
EW_HHS<SType>::EW_HHS(const Vertex_Key &key):
  SSS_Vertex<SType>(key,1,0), m_cpl(0.0) 
{
  if (key.p_c->Flav().Kfcode()==kf_h0_qsc ||
      key.p_b->Flav()==Flavour(kf_h0_qsc).Bar()) {
    m_cpl=M_I;
  }
  else {
    m_cpl=M_I*SType(3.0*sqr(key.p_a->Flav().Yuk()))/
      SComplex(key.p_model->Constant("v_{EW}"));
    if (key.p_b->Flav()==Flavour(kf_h0_qsc))
      m_cpl/=SComplex(3.0*key.p_model->Constant("v_{EW}"));
  }
  msg_Debugging()<<"m_cpl = "<<m_cpl<<" <- "<<key.ID()<<"\n";
}

template <typename SType>
void EW_HHS<SType>::Evaluate(const CScalarType &a,const CScalarType &b)
{
  AddJ(m_cpl*SType(this->m_cplfac)*Lorentz(a,b));
}

template <typename SType>
std::string EW_HHS<SType>::CVLabel() const
{
  return "("+ToString(m_cpl.real())+",,"+ToString(m_cpl.imag())+")";
}

template <typename SType>
EW_FFH<SType>::EW_FFH(const Vertex_Key &key):
  FFS_Vertex<SType>(key,1,0),
  m_cpl(-M_I/SComplex(key.p_model->Constant("v_{EW}"))) 
{
  Flavour f((key.p_a->Flav().IsFermion()?key.p_a:key.p_b)->RFlav());
  std::string tag(f.IDName());
  if (f.IsLepton()) tag=tag.substr(0,tag.length()-1);
  tag="Yukawa_"+tag;
  m_cpl*=key.p_model->ScalarConstant(tag);
  if (IsZero(m_cpl)) this->m_act=false;
  msg_Debugging()<<"m_cpl = "<<m_cpl<<" <- "<<key.ID()<<"\n";
}

template <typename SType>
void EW_FFH<SType>::Evaluate(const CSpinorType &a,const CSpinorType &b)
{
  if (a()!=b()) return;
  AddJ(m_cpl*SType(this->m_cplfac)*Lorentz(a,b));
}

template <typename SType>
void EW_FFH<SType>::Evaluate(const CSpinorType &a,const CScalarType &b)
{
  AddJ(Lorentz(a,b)*m_cpl*SType(this->m_cplfac));
}

template <typename SType>
std::string EW_FFH<SType>::CVLabel() const
{
  return "("+ToString(m_cpl.real())+",,"+ToString(m_cpl.imag())+")";
}

template <typename SType>
EW_FFZ<SType>::EW_FFZ(const Vertex_Key &key):
  FFV_Vertex<SType>(key,1,0),
  m_cpll(-M_I*SComplex(key.p_model->Constant("g_2")/2.0/
	 key.p_model->Constant("\\cos\\theta_W"))),
  m_cplr(m_cpll)
{
  Flavour f((key.p_a->Flav().IsFermion()?key.p_a:key.p_b)->RFlav());
  SComplex vf(2.0*f.Charge()), af(f.IsoWeak());
  vf=af-vf*SComplex(csqr(key.p_model->Constant("\\sin\\theta_W")));
  m_cpll*=vf+af;
  m_cplr*=vf-af;
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"(): "<<key.p_a->Flav()<<" -> A_f = "
		 <<af<<" ( T_f^3 = "<<dabs(key.p_a->Flav().IsoWeak())
		 <<" ), V_f = "<<vf<<" ( Q_f = "
		 <<dabs(key.p_a->Flav().Charge())<<" ) -> ("
		 <<m_cpll<<","<<m_cplr<<")\n";
#endif
}

template <typename SType>
void EW_FFZ<SType>::Evaluate(const CSpinorType &a,const CSpinorType &b)
{
  if (a()!=b()) return;
  bool cl(CalcLeft(a,b)), cr(CalcRight(a,b));
  if (!(cl || cr)) return;
  CVec4Type j(ZERO,ZERO,ZERO,ZERO,0,0,a.H(0)+b.H(0),a.H(1)+b.H(1));
  if (cl) j+=m_cpll*SType(this->m_cplfac)*LorentzLeft(a,b);
  if (cr) j+=m_cplr*SType(this->m_cplfac)*LorentzRight(a,b);
  AddJ(j);
}

template <typename SType>
void EW_FFZ<SType>::Evaluate(const CSpinorType &a,const CVec4Type &b)
{
  bool cl(CalcLeft(a)), cr(CalcRight(a));
  if (!(cl || cr)) return;
  CSpinorType j(a.R(),a.B(),a(),a.H(0)+b.H(0),a.H(1)+b.H(1),0);
  if (cl) j+=LorentzLeft(a,b)*m_cpll*SType(this->m_cplfac);
  if (cr) j+=LorentzRight(a,b)*m_cplr*SType(this->m_cplfac);
  AddJ(j);
}

template <typename SType>
std::string EW_FFZ<SType>::CVLabel() const
{
  return "\\begin{array}{l}\\scriptstyle L("+
    ToString(m_cpll.real())+",,"+ToString(m_cpll.imag())+
    ")\\\\\\scriptstyle R("+ToString(m_cplr.real())
    +",,"+ToString(m_cplr.imag())+")\\end{array}";
}

template <typename SType>
EW_FFW<SType>::EW_FFW(const Vertex_Key &key):
  FFV_Vertex<SType>(key,1,0),
  m_cpl(-M_I*SComplex(key.p_model->Constant("g_2")/sqrt(2.0))) 
{
  Flavour f1(key.p_a->Flav()), f2(key.p_b->Flav());
  if (!key.p_a->Flav().IsFermion()) f1=key.p_c->Flav();
  else if (!key.p_b->Flav().IsFermion()) f2=key.p_c->Flav();
  if (f1.IsQuark()) {
    int i((int)(f1.Kfcode())), j((int)(f2.Kfcode()));
    if (f1.IsDowntype()) std::swap<int>(i,j);
    m_cpl*=key.p_model->ComplexMatrixElement("CKM",i/2-1,(j-1)/2);
  }
  msg_Debugging()<<"m_cpl = "<<m_cpl<<" <- "<<key.ID()<<"\n";
  if (IsZero(m_cpl)) this->m_act=false;
}

template <typename SType>
void EW_FFW<SType>::Evaluate(const CSpinorType &a,const CSpinorType &b)
{
  if (a()!=b()) return;
  if (!CalcLeft(a,b)) return;
  AddJ(m_cpl*SType(this->m_cplfac)*LorentzLeft(a,b));
}

template <typename SType>
void EW_FFW<SType>::Evaluate(const CSpinorType &a,const CVec4Type &b)
{
  if (!CalcLeft(a)) return;
  AddJ(LorentzLeft(a,b)*m_cpl*SType(this->m_cplfac));
}

template <typename SType>
std::string EW_FFW<SType>::CVLabel() const
{
  return "("+ToString(m_cpl.real())+",,"+ToString(m_cpl.imag())+")";
}

template <typename SType>
EW_VVH<SType>::EW_VVH(const Vertex_Key &key):
  VVS_Vertex<SType>(key,1,0),
  m_cpl(key.p_model->Constant("g_2")) 
{
  Flavour f((key.p_a->Flav().IsVector()?key.p_a:key.p_b)->RFlav());
  if (f.Kfcode()==kf_Z) m_cpl/=key.p_model->Constant("\\cos\\theta_W");
  if (key.p_a->Flav().Kfcode()==kf_h0_qsc ||
      key.p_b->Flav().Kfcode()==kf_h0_qsc ||
      key.p_c->Flav().Kfcode()==kf_h0_qsc) {
    m_cpl=-M_I*csqr(m_cpl)/SType(2.0);
  }
  else {
    m_cpl*=-M_I;
    switch (f.Kfcode()) {
    case kf_Wplus: m_cpl*=SComplex(key.p_model->Constant("m_W")); break;
    case kf_Z: m_cpl*=SComplex(key.p_model->Constant("m_Z")); break;
    default:
      THROW(fatal_error,"Internal error");
    }
  }
  msg_Debugging()<<"m_cpl = "<<m_cpl<<" <- "<<key.ID()<<"\n";
}

template <typename SType>
void EW_VVH<SType>::Evaluate(const CVec4Type &a,const CVec4Type &b)
{
  AddJ(m_cpl*SType(this->m_cplfac)*Lorentz(a,b));
}

template <typename SType>
void EW_VVH<SType>::Evaluate(const CVec4Type &a,const CScalarType &b)
{
  AddJ(m_cpl*SType(this->m_cplfac)*Lorentz(a,b));
}

template <typename SType>
std::string EW_VVH<SType>::CVLabel() const
{
  return "("+ToString(m_cpl.real())+",,"+ToString(m_cpl.imag())+")";
}

template <typename SType>
EW_WWV<SType>::EW_WWV(const Vertex_Key &key):
  VVV_Vertex<SType>(key,1,0),
  m_cpl(M_I*SComplex(key.p_model->Constant("g_2"))) 
{
  if (key.p_a->Flav().IsPhoton()||key.p_b->Flav().IsPhoton()||
      key.p_c->Flav().IsPhoton())
    m_cpl*=key.p_model->Constant("\\sin\\theta_W");
  else m_cpl*=key.p_model->Constant("\\cos\\theta_W");
  msg_Debugging()<<"m_cpl = "<<m_cpl<<" <- "<<key.ID()<<"\n";
}

template <typename SType>
void EW_WWV<SType>::Evaluate(const CVec4Type &a,const CVec4Type &b)
{
  CVec4Type j(Lorentz(a,b));
#ifdef DEBUG__BG
  msg_Debugging()<<"'+' "<<a<<"\n";
  msg_Debugging()<<"    "<<b<<"\n";
#endif
  j(1)=j(0)=0;
  AddJ(m_cpl*SType(this->m_cplfac)*j);
}

template <typename SType>
std::string EW_WWV<SType>::CVLabel() const
{
  return "("+ToString(m_cpl.real())+",,"+ToString(m_cpl.imag())+")";
}

template <typename SType>
EW_WWT<SType>::EW_WWT(const Vertex_Key &key):
  VVT_Vertex<SType>(key,1,0),
  m_cpl(M_I*SComplex(key.p_model->Constant("g_2"))) 
{
  if (key.p_a->Flav().Kfcode()==kf_photon || 
      key.p_b->Flav().Kfcode()==kf_photon ||
      key.p_c->Flav().Kfcode()==kf_photon) 
    m_cpl*=key.p_model->Constant("\\sin\\theta_W");
  if (key.p_a->Flav().Kfcode()==kf_Z ||
      key.p_b->Flav().Kfcode()==kf_Z ||
      key.p_c->Flav().Kfcode()==kf_Z) 
    m_cpl*=key.p_model->Constant("\\cos\\theta_W");
  if (key.p_a->Flav().IsTensor()) m_cpl=-m_cpl;
  msg_Debugging()<<"m_cpl = "<<m_cpl<<" <- "<<key.ID()<<"\n";
}

template <typename SType>
void EW_WWT<SType>::Evaluate(const CVec4Type &a,const CVec4Type &b)
{
#ifdef DEBUG__BG
  msg_Debugging()<<"'+' "<<a<<"\n";
  msg_Debugging()<<"    "<<b<<"\n";
#endif
  AddJ(m_cpl*SType(this->m_cplfac)*Lorentz(a,b));
}

template <typename SType>
void EW_WWT<SType>::Evaluate(const CVec4Type &a,const CAsT4Type &b)
{
  CVec4Type j(Lorentz(a,b));
#ifdef DEBUG__BG
  msg_Debugging()<<"'+' "<<a<<"\n";
  msg_Debugging()<<"    "<<b<<"\n";
#endif
  j(1)=j(0)=0;
  AddJ(m_cpl*SType(this->m_cplfac)*j);
}

template <typename SType>
std::string EW_WWT<SType>::CVLabel() const
{
  return "("+ToString(m_cpl.real())+",,"+ToString(m_cpl.imag())+")";
}

template <typename SType,char STag> Vertex_Base *
EW_HHS_Getter<SType,STag>::operator()(const Vertex_Key &key) const
{
  return new EW_HHS<SType>(key);
}

template <typename SType,char STag>
void EW_HHS_Getter<SType,STag>::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"ew hhs vertex "<<STag;
}

template <typename SType,char STag> Vertex_Base *
EW_FFH_Getter<SType,STag>::operator()(const Vertex_Key &key) const
{
  return new EW_FFH<SType>(key);
}

template <typename SType,char STag>
void EW_FFH_Getter<SType,STag>::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"ew ffh vertex "<<STag;
}

template <typename SType,char STag> Vertex_Base *
EW_VVH_Getter<SType,STag>::operator()(const Vertex_Key &key) const
{
  return new EW_VVH<SType>(key);
}

template <typename SType,char STag>
void EW_VVH_Getter<SType,STag>::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"ew vvh vertex "<<STag;
}

template <typename SType,char STag> Vertex_Base *
EW_FFZ_Getter<SType,STag>::operator()(const Vertex_Key &key) const
{
  return new EW_FFZ<SType>(key);
}

template <typename SType,char STag>
void EW_FFZ_Getter<SType,STag>::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"ew ffz vertex "<<STag;
}

template <typename SType,char STag> Vertex_Base *
EW_FFW_Getter<SType,STag>::operator()(const Vertex_Key &key) const
{
  return new EW_FFW<SType>(key);
}

template <typename SType,char STag>
void EW_FFW_Getter<SType,STag>::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"ew ffw vertex "<<STag;
}

template <typename SType,char STag> Vertex_Base *
EW_WWV_Getter<SType,STag>::operator()(const Vertex_Key &key) const
{
  return new EW_WWV<SType>(key);
}

template <typename SType,char STag>
void EW_WWV_Getter<SType,STag>::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"ew wwv vertex "<<STag;
}

template <typename SType,char STag> Vertex_Base *
EW_WWT_Getter<SType,STag>::operator()(const Vertex_Key &key) const
{
  return new EW_WWT<SType>(key);
}

template <typename SType,char STag>
void EW_WWT_Getter<SType,STag>::PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"ew wwt vertex "<<STag;
}

namespace COMIX {

  template class EW_HHS_Getter<double,'D'>;
  template class EW_FFH_Getter<double,'D'>;
  template class EW_VVH_Getter<double,'D'>;
  template class EW_FFZ_Getter<double,'D'>;
  template class EW_FFW_Getter<double,'D'>;
  template class EW_WWV_Getter<double,'D'>;
  template class EW_WWT_Getter<double,'D'>;

  template class EW_HHS_Getter<long double,'Q'>;
  template class EW_FFH_Getter<long double,'Q'>;
  template class EW_VVH_Getter<long double,'Q'>;
  template class EW_FFZ_Getter<long double,'Q'>;
  template class EW_FFW_Getter<long double,'Q'>;
  template class EW_WWV_Getter<long double,'Q'>;
  template class EW_WWT_Getter<long double,'Q'>;

}

DECLARE_TEMPLATE_GETTER(EW_Vertex_Filler,"EW_VFill",void,const Model *);

template <typename SType,char STag>
void *EW_Vertex_Filler<SType,STag>::operator()
  (const Model *const &model) const
{
  if (!model->IncludesModel("EW")) return NULL;
  std::string wptag("{"+Flavour(kf_Wplus).IDName()+"}");
  std::string wmtag("{"+Flavour(kf_Wplus).Bar().IDName()+"}");
  std::string ztag("{"+Flavour(kf_Z).IDName()+"}");
  std::string atag("{"+Flavour(kf_photon).IDName()+"}");
  std::string htag("{"+Flavour(kf_h0).IDName()+"}");
  std::string wpttag("{"+Flavour(kf_Wplus_qgc).IDName()+"}");
  std::string wmttag("{"+Flavour(kf_Wplus_qgc).Bar().IDName()+"}");
  std::string zttag("{"+Flavour(kf_Z_qgc).IDName()+"}");
  std::string hsrtag("{"+Flavour(kf_h0_qsc).IDName()+"}");
  std::string hsatag("{"+Flavour(kf_h0_qsc).Bar().IDName()+"}");
  if (Flavour(kf_h0).IsOn()) {
    // ffh couplings
    for (int i(1);i<=16;++i) {
      if (i==7) i=11;
      Flavour f((kf_code)i);
      if (!f.IsOn()) continue;
      std::string ftag("{"+f.IDName()+"}"), fbtag("{"+f.Bar().IDName()+"}");
      if (i<=6 || i%2==1) {
	new EW_FFH_Getter<SType,STag>(ftag+fbtag+htag);
	new EW_FFH_Getter<SType,STag>(fbtag+htag+fbtag);
	new EW_FFH_Getter<SType,STag>(htag+ftag+ftag);
      }
    }
    // hhs couplings
    new EW_HHS_Getter<SType,STag>(htag+htag+htag);
    new EW_HHS_Getter<SType,STag>(htag+htag+hsrtag);
    new EW_HHS_Getter<SType,STag>(htag+hsatag+htag);
    new EW_HHS_Getter<SType,STag>(htag+hsrtag+htag);
  }
  if (Flavour(kf_Z).IsOn()) {
    // ffz couplings
    for (int i(1);i<=16;++i) {
      if (i==7) i=11;
      Flavour f((kf_code)i);
      if (!f.IsOn()) continue;
      std::string ftag("{"+f.IDName()+"}"), fbtag("{"+f.Bar().IDName()+"}");
      new EW_FFZ_Getter<SType,STag>(ftag+fbtag+ztag);
      new EW_FFZ_Getter<SType,STag>(fbtag+ztag+fbtag);
      new EW_FFZ_Getter<SType,STag>(ztag+ftag+ftag);
    }
    // zzs couplings
    new EW_VVH_Getter<SType,STag>(ztag+ztag+htag);
    new EW_VVH_Getter<SType,STag>(ztag+htag+ztag);
    new EW_VVH_Getter<SType,STag>(ztag+ztag+hsatag);
    new EW_VVH_Getter<SType,STag>(ztag+hsrtag+ztag);
  }
  if (Flavour(kf_Wplus).IsOn()) {
    // qqw couplings
    for (int i(1);i<=5;i+=2) {
      for (int j(2);j<=6;j+=2) {
	Flavour f1((kf_code)i), f2((kf_code)j);
	if (!f1.IsOn() || !f2.IsOn()) continue;
	std::string f1tag("{"+f1.IDName()+"}");
	std::string f2tag("{"+f2.IDName()+"}");
	std::string f1btag("{"+f1.Bar().IDName()+"}");
	std::string f2btag("{"+f2.Bar().IDName()+"}");
	new EW_FFW_Getter<SType,STag>(f1tag+f2btag+wmtag);
	new EW_FFW_Getter<SType,STag>(f2tag+f1btag+wptag);
	new EW_FFW_Getter<SType,STag>(f2btag+wptag+f1btag);
	new EW_FFW_Getter<SType,STag>(f1btag+wmtag+f2btag);
	new EW_FFW_Getter<SType,STag>(wptag+f1tag+f2tag);
	new EW_FFW_Getter<SType,STag>(wmtag+f2tag+f1tag);
      }
    }
    // llw couplings
    for (int i(11);i<=16;i+=2) {
      Flavour f1((kf_code)i), f2((kf_code)(i+1));
      if (!f1.IsOn() || !f2.IsOn()) continue;
      std::string f1tag("{"+f1.IDName()+"}");
      std::string f2tag("{"+f2.IDName()+"}");
      std::string f1btag("{"+f1.Bar().IDName()+"}");
      std::string f2btag("{"+f2.Bar().IDName()+"}");
      new EW_FFW_Getter<SType,STag>(f1tag+f2btag+wmtag);
      new EW_FFW_Getter<SType,STag>(f2tag+f1btag+wptag);
      new EW_FFW_Getter<SType,STag>(f2btag+wptag+f1btag);
      new EW_FFW_Getter<SType,STag>(f1btag+wmtag+f2btag);
      new EW_FFW_Getter<SType,STag>(wptag+f1tag+f2tag);
      new EW_FFW_Getter<SType,STag>(wmtag+f2tag+f1tag);
    }
    // wwt couplings
    new EW_WWT_Getter<SType,STag>(wmtag+wptag+zttag);
    new EW_WWT_Getter<SType,STag>(wptag+zttag+wptag);
    new EW_WWT_Getter<SType,STag>(zttag+wmtag+wmtag);
    if (Flavour(kf_photon).IsOn()) {
      // wwa couplings
      new EW_WWV_Getter<SType,STag>(wmtag+wptag+atag);
      new EW_WWV_Getter<SType,STag>(wptag+atag+wptag);
      new EW_WWV_Getter<SType,STag>(atag+wmtag+wmtag);
      // wat couplings
      new EW_WWT_Getter<SType,STag>(atag+wmtag+wmttag);
      new EW_WWT_Getter<SType,STag>(wmtag+wpttag+atag);
      new EW_WWT_Getter<SType,STag>(wpttag+atag+wptag);
      new EW_WWT_Getter<SType,STag>(wptag+atag+wpttag);
      new EW_WWT_Getter<SType,STag>(atag+wmttag+wmtag);
      new EW_WWT_Getter<SType,STag>(wmttag+wptag+atag);
    }
    if (Flavour(kf_Z).IsOn()) {
      // wwz couplings
      new EW_WWV_Getter<SType,STag>(wmtag+wptag+ztag);
      new EW_WWV_Getter<SType,STag>(wptag+ztag+wptag);
      new EW_WWV_Getter<SType,STag>(ztag+wmtag+wmtag);
      // wzt couplings
      new EW_WWT_Getter<SType,STag>(ztag+wmtag+wmttag);
      new EW_WWT_Getter<SType,STag>(wmtag+wpttag+ztag);
      new EW_WWT_Getter<SType,STag>(wpttag+ztag+wptag);
      new EW_WWT_Getter<SType,STag>(wptag+ztag+wpttag);
      new EW_WWT_Getter<SType,STag>(ztag+wmttag+wmtag);
      new EW_WWT_Getter<SType,STag>(wmttag+wptag+ztag);
    }
    if (Flavour(kf_h0).IsOn()) {
      // wws couplings
      new EW_VVH_Getter<SType,STag>(wmtag+wptag+htag);
      new EW_VVH_Getter<SType,STag>(wptag+htag+wptag);
      new EW_VVH_Getter<SType,STag>(htag+wmtag+wmtag);
      new EW_VVH_Getter<SType,STag>(wptag+wmtag+hsatag);
      new EW_VVH_Getter<SType,STag>(wmtag+hsrtag+wmtag);
      new EW_VVH_Getter<SType,STag>(hsrtag+wptag+wptag);
    }
  }
  return NULL;
}

template <typename SType,char STag>
void EW_Vertex_Filler<SType,STag>::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"ew vertex filler "<<STag;
}

template class EW_Vertex_Filler<double,'D'>;

template class EW_Vertex_Filler<long double,'Q'>;
