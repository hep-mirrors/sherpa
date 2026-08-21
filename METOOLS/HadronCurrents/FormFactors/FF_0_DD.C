#include "METOOLS/HadronCurrents/FormFactors/FF_0_DD.H"
#include "METOOLS/HadronCurrents/Tools.H"

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;


void FF_0_DD_Base::FixMode() {
  if (m_flavs[m_pi[0]].Kfcode()==kf_p_plus &&
      m_flavs[m_pi[1]].Kfcode()==kf_p_plus)
    m_mode = FF_0_DD_mode::pp;
  else if (m_flavs[m_pi[0]].Kfcode()==kf_n &&
	   m_flavs[m_pi[1]].Kfcode()==kf_n)
    m_mode = FF_0_DD_mode::nn;
}

Complex FF_0_DD_Base::operator()(const ATOOLS::Vec4D_Vector& moms) {
  return Complex(0.,0.);
}

class FF_0_DD_None : public FF_0_DD_Base {
public:
  FF_0_DD_None(const FF_Parameters & params) :
    FF_0_DD_Base(params) {}  
  Complex operator()(const ATOOLS::Vec4D_Vector& moms) {
    return Complex(0.,0.);
  }
};

class FF_0_DD_Point : public FF_0_DD_Base {
public:
  FF_0_DD_Point(const FF_Parameters & params) :
    FF_0_DD_Base(params) {}
  Complex operator()(const ATOOLS::Vec4D_Vector& moms) {
    return Complex(1.,0.);
  }
};

class FF_0_DD_Dipole : public FF_0_DD_Base {
private:
  double m_Lambda2;
  void FixParameters(const FF_Parameters & params) {
    msg_Out()<<METHOD<<": ";
    m_norm = 1.;
    if (params.m_parameterMap.find("Lambda2")!=
	params.m_parameterMap.end()) {
      m_Lambda2 = params.m_parameterMap.at("Lambda2");
      msg_Out()<<"Lambda^2 = "<<m_Lambda2<<"\n";
    }
    else msg_Out()<<"no Lambda2 tag in parameter map.\n";
  }
public:
  FF_0_DD_Dipole(const FF_Parameters & params) :
    FF_0_DD_Base(params),
    m_Lambda2(1.)
  {
    FixParameters(params);
  }  
  Complex operator()(const ATOOLS::Vec4D_Vector& moms) {
    double Q2 = dabs((moms[m_pi[0]]+moms[m_pi[1]]).Abs2()); 
    return Complex(1./sqr(1.+Q2/m_Lambda2),0.);
  }
};




DECLARE_FF_GETTER(FF_0_DD_Base,"FF_0_DD")

FormFactor_Base * ATOOLS::Getter<FormFactor_Base,FF_Parameters,
				 FF_0_DD_Base>:: 
operator()(const METOOLS::FF_Parameters &params) const
{
  msg_Out()<<METHOD<<":\n";
  size_t Nbaryons = 0;
  for (size_t i=0;i<params.m_pi.size();i++) {
    if (params.m_flavs[params.m_pi[i]].IsBaryon()) Nbaryons++;
  }
  if (Nbaryons!=2) return NULL;
  if (//   pp
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_p_plus &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_p_plus)    ||
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_n &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_n)         ) {
    if (params.m_ffmodel==ff_model::none)   return new FF_0_DD_None(params);
    if (params.m_ffmodel==ff_model::point)  return new FF_0_DD_Point(params);
    if (params.m_ffmodel==ff_model::dipole) return new FF_0_DD_Dipole(params);
  }
  return NULL;
}

