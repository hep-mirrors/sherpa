#ifndef Fixed_Variable_Channel_C
#define Fixed_Variable_Channel_C

#include "Fixed_Variable_Channel.H"

namespace PHASIC {

  template <class Value_Type>
  Fixed_Variable_Channel<Value_Type>::Fixed_Variable_Channel(int _nin,int _nout,ATOOLS::Flavour *_fl,
							 ATOOLS::Variable _m_variable):
    PHASIC::Channel_Interface(_nin,_nout,_fl),
    m_variable(_m_variable) {}
  
  template <class Value_Type>
  void Fixed_Variable_Channel<Value_Type>::GeneratePoint(ATOOLS::Vec4D *_p,double *_ran)
  {
    switch (m_variable.Type()) {
    case ATOOLS::Variable::p_perp:
      Ehat=sqrt((_p[0]+_p[1]).Abs2());
      pt=(double)m_value;
      if (Ehat/2.0>pt) {
	weight=1.0; 
	_p[2]=ATOOLS::Vec4D(Ehat/2.0,pt*cos(2.0*M_PI*_ran[0]),pt*sin(2.0*M_PI*_ran[0]),sqrt(Ehat*Ehat/4.0-pt*pt));
      }
      else {
	weight=0.0;
	_p[2]=ATOOLS::Vec4D(Ehat/2.0,Ehat/2.0*cos(2.0*M_PI*_ran[0]),Ehat/2.0*sin(2.0*M_PI*_ran[0]),0.0);
      }
      _p[3]=ATOOLS::Vec4D(Ehat/2.0,ATOOLS::Vec3D()-ATOOLS::Vec3D(_p[2]));
      break;
    case ATOOLS::Variable::E_perp:
      Ehat=sqrt((_p[0]+_p[1]).Abs2());
      pt=sqrt((double)(m_value*m_value)-_p[0].Abs2());
      if (Ehat/2.0>pt) {
	weight=1.0; 
	_p[2]=ATOOLS::Vec4D(Ehat/2.0,pt*cos(2.0*M_PI*_ran[0]),pt*sin(2.0*M_PI*_ran[0]),sqrt(Ehat*Ehat/4.0-pt*pt));
      }
      else {
	weight=0.0;
	_p[2]=ATOOLS::Vec4D(Ehat/2.0,Ehat/2.0*cos(2.0*M_PI*_ran[0]),Ehat/2.0*sin(2.0*M_PI*_ran[0]),0.0);
      }
      _p[3]=ATOOLS::Vec4D(Ehat/2.0,ATOOLS::Vec3D()-ATOOLS::Vec3D(_p[2]));
      break;
    default:
      ATOOLS::msg.Error()<<"Fixed_Variable_Channel::GeneratePoint(..): "
			 <<"Cannot handle "<<m_variable.Name()<<"! Setting weight to 0."<<std::endl;
      weight=0.0;
      break;
    }
  }
  
  template <class Value_Type>
  void Fixed_Variable_Channel<Value_Type>::GenerateWeight(ATOOLS::Vec4D *_p)
  {
    weight/=PHASIC::CE.Isotropic2Weight(_p[2],_p[3])*pow(2.0*M_PI,2.0*3.0-4.0);
  }

} // end of namespace PHASIC

#endif
