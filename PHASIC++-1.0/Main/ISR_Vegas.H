#ifndef ISR_Vegas_H
#define ISR_Vegas_H

#include "Single_Channel.H"
#include "Vegas.H"

namespace PHASIC {

  class Threshold_Uniform_V: public Single_Channel {
  protected:

    double m_mass,m_sexp;
    bool   m_zchannel;

    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey, m_sgridkey, m_ygridkey;
    Vegas* p_vegas;
  public:

    // constructor
    Threshold_Uniform_V(const double mass,const double sexp,
		      const std::string cinfo,ATOOLS::Integration_Info *info);

    // member functions
    void GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
		       const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);
    void AddPoint(double);
    void Optimize()  {p_vegas->Optimize();} 
    void EndOptimize()  {p_vegas->EndOptimize();} 
    void WriteOut(std::string pId) { p_vegas->WriteOut(pId); }
    void ReadIn(std::string pId)   { p_vegas->ReadIn(pId); }

  };// end of class Threshold_Uniform

  class Threshold_Forward_V: public Single_Channel {
  protected :

    double m_mass,m_sexp, m_yexponent;
    bool   m_zchannel;

    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey, m_sgridkey, m_ygridkey;
    Vegas* p_vegas;
  public :

    // constructor
    Threshold_Forward_V(const double mass,const double sexp,const double yexponent,
		      const std::string cinfo,ATOOLS::Integration_Info *info);

    // member functions
    void GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
		       const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);
    void AddPoint(double);
    void Optimize()  {p_vegas->Optimize();} 
    void EndOptimize()  {p_vegas->EndOptimize();} 
    void WriteOut(std::string pId) { p_vegas->WriteOut(pId); }
    void ReadIn(std::string pId)   { p_vegas->ReadIn(pId); }

  };// end of class Threshold_Forward

  class Threshold_Backward_V: public Single_Channel {
  protected :

    double m_mass,m_sexp, m_yexponent;
    bool   m_zchannel;

    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey, m_sgridkey, m_ygridkey;
    Vegas* p_vegas;
  public :

    // constructor
    Threshold_Backward_V(const double mass,const double sexp,const double yexponent,
		       const std::string cinfo,ATOOLS::Integration_Info *info);

    // member functions
    void GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
		       const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);
    void AddPoint(double);
    void Optimize()  {p_vegas->Optimize();} 
    void EndOptimize()  {p_vegas->EndOptimize();} 
    void WriteOut(std::string pId) { p_vegas->WriteOut(pId); }
    void ReadIn(std::string pId)   { p_vegas->ReadIn(pId); }

  };// end of class Threshold_Backward

  class Threshold_Central_V: public Single_Channel {
  protected:

    double m_mass,m_sexp;
    bool   m_zchannel;

    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey, m_sgridkey, m_ygridkey;
    Vegas* p_vegas;
  public:

    // constructor
    Threshold_Central_V(const double mass,const double sexp,
		      const std::string cinfo,ATOOLS::Integration_Info *info,const int mode=0);

    // member functions
    void GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
		       const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);
    void AddPoint(double);
    void Optimize()  {p_vegas->Optimize();} 
    void EndOptimize()  {p_vegas->EndOptimize();} 
    void WriteOut(std::string pId) { p_vegas->WriteOut(pId); }
    void ReadIn(std::string pId)   { p_vegas->ReadIn(pId); }

  };// end of class Threshold_Central


  ///////////////////////////////////////////////////////////////////////////////////////////

  class Resonance_Uniform_V: public Single_Channel {
  protected:

    double m_mass, m_width;
    bool   m_zchannel;

    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey, m_sgridkey, m_ygridkey;
    Vegas* p_vegas;
  public:

    // constructor
    Resonance_Uniform_V(const double mass,const double width,
		      const std::string cinfo,ATOOLS::Integration_Info *info);

    // member functions
    void GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
		       const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);
    void AddPoint(double);
    void Optimize()  {p_vegas->Optimize();} 
    void EndOptimize()  {p_vegas->EndOptimize();} 
    void WriteOut(std::string pId) { p_vegas->WriteOut(pId); }
    void ReadIn(std::string pId)   { p_vegas->ReadIn(pId); }

  };// end of class Resonance_Uniform

  class Resonance_Forward_V: public Single_Channel {
  protected :

    double m_mass, m_width, m_yexponent;
    bool   m_zchannel;

    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey, m_sgridkey, m_ygridkey;
    Vegas* p_vegas;
  public :

    // constructor
    Resonance_Forward_V(const double mass,const double width,const double yexponent,
		      const std::string cinfo,ATOOLS::Integration_Info *info);

    // member functions
    void GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
		       const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);
    void AddPoint(double);
    void Optimize()  {p_vegas->Optimize();} 
    void EndOptimize()  {p_vegas->EndOptimize();} 
    void WriteOut(std::string pId) { p_vegas->WriteOut(pId); }
    void ReadIn(std::string pId)   { p_vegas->ReadIn(pId); }

  };// end of class Resonance_Forward

  class Resonance_Backward_V: public Single_Channel {
  protected :

    double m_mass, m_width, m_yexponent;
    bool   m_zchannel;

    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey, m_sgridkey, m_ygridkey;
    Vegas* p_vegas;
  public :

    // constructor
    Resonance_Backward_V(const double mass,const double width,const double yexponent,
		       const std::string cinfo,ATOOLS::Integration_Info *info);

    // member functions
    void GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
		       const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);
    void AddPoint(double);
    void Optimize()  {p_vegas->Optimize();} 
    void EndOptimize()  {p_vegas->EndOptimize();} 
    void WriteOut(std::string pId) { p_vegas->WriteOut(pId); }
    void ReadIn(std::string pId)   { p_vegas->ReadIn(pId); }

  };// end of class Resonance_Backward

  class Resonance_Central_V: public Single_Channel {
  protected:

    double m_mass, m_width;
    bool   m_zchannel;

    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey, m_sgridkey, m_ygridkey;
    Vegas* p_vegas;
  public:

    // constructor
    Resonance_Central_V(const double mass,const double width,
		      const std::string cinfo,ATOOLS::Integration_Info *info,const int mode=0);

    // member functions
    void GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
		       const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);
    void AddPoint(double);
    void Optimize()  {p_vegas->Optimize();} 
    void EndOptimize()  {p_vegas->EndOptimize();} 
    void WriteOut(std::string pId) { p_vegas->WriteOut(pId); }
    void ReadIn(std::string pId)   { p_vegas->ReadIn(pId); }

  };// end of class Resonance_Central


  /////////////////////////////////////////////////////////////////////////////////////

  class Simple_Pole_Uniform_V: public Single_Channel {
  protected:

    double m_exponent;
    bool   m_zchannel;

    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey, m_sgridkey, m_ygridkey;
    Vegas* p_vegas;
  public:

    // constructor
    Simple_Pole_Uniform_V(const double exponent,const std::string cinfo,
			ATOOLS::Integration_Info *info);

    // member functions
    void GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
		       const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);
    void AddPoint(double);
    void Optimize()  {p_vegas->Optimize();} 
    void EndOptimize()  {p_vegas->EndOptimize();} 
    void WriteOut(std::string pId) { p_vegas->WriteOut(pId); }
    void ReadIn(std::string pId)   { p_vegas->ReadIn(pId); }

  };// end of class Simple_Pole_Uniform

  class Simple_Pole_Forward_V: public Single_Channel {
  protected :

    double m_sexponent, m_yexponent;
    bool   m_zchannel;

    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey, m_sgridkey, m_ygridkey;
    Vegas* p_vegas;
  public :

    // constructor
    Simple_Pole_Forward_V(const double sexponent,const double yexponent,
			const std::string cinfo,ATOOLS::Integration_Info *info);

    // member functions
    void GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
		       const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);
    void AddPoint(double);
    void Optimize()  {p_vegas->Optimize();} 
    void EndOptimize()  {p_vegas->EndOptimize();} 
    void WriteOut(std::string pId) { p_vegas->WriteOut(pId); }
    void ReadIn(std::string pId)   { p_vegas->ReadIn(pId); }

  };// end of class Simple_Pole_Forward

  class Simple_Pole_Backward_V: public Single_Channel {
  protected :

    double m_sexponent, m_yexponent;
    bool   m_zchannel;

    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey, m_sgridkey, m_ygridkey;
    Vegas* p_vegas;
  public :

    // constructor
    Simple_Pole_Backward_V(const double sexponent,const double yexponent,
			 const std::string cinfo,ATOOLS::Integration_Info *info);

    // member functions
    void GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
		       const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);
    void AddPoint(double);
    void Optimize()  {p_vegas->Optimize();} 
    void EndOptimize()  {p_vegas->EndOptimize();} 
    void WriteOut(std::string pId) { p_vegas->WriteOut(pId); }
    void ReadIn(std::string pId)   { p_vegas->ReadIn(pId); }

  };// end of class Simple_Pole_Backward

  class Simple_Pole_Central_V: public Single_Channel {
  protected:

    double m_exponent;
    bool   m_zchannel;

    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey, m_sgridkey, m_ygridkey;
    Vegas* p_vegas;
  public:

    // constructor
    Simple_Pole_Central_V(const double exponent,const std::string cinfo,
			ATOOLS::Integration_Info *info,const int mode=0);

    // member functions
    void GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
		       const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);
    void AddPoint(double);
    void Optimize()  {p_vegas->Optimize();} 
    void EndOptimize()  {p_vegas->EndOptimize();} 
    void WriteOut(std::string pId) { p_vegas->WriteOut(pId); }
    void ReadIn(std::string pId)   { p_vegas->ReadIn(pId); }

  };// end of class Simple_Pole_Central


  /////////////////////////////////////////////////////////////////////////////////

  class Leading_Log_Uniform_V: public Single_Channel {
  protected:

    double m_beta, m_factor;
    bool   m_zchannel;

    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey, m_sgridkey, m_ygridkey;
    Vegas* p_vegas;
  public:

    // constructor
    Leading_Log_Uniform_V(const double beta,const double factor,
			const std::string cinfo,ATOOLS::Integration_Info *info);

    // member functions
    void GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
		       const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);
    void AddPoint(double);
    void Optimize()  {p_vegas->Optimize();} 
    void EndOptimize()  {p_vegas->EndOptimize();} 
    void WriteOut(std::string pId) { p_vegas->WriteOut(pId); }
    void ReadIn(std::string pId)   { p_vegas->ReadIn(pId); }

  };// end of class Leading_Log_Uniform

  class Leading_Log_Forward_V: public Single_Channel {
  protected :

    double m_beta, m_factor, m_yexponent;
    bool   m_zchannel;

    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey, m_sgridkey, m_ygridkey;
    Vegas* p_vegas;
  public :

    // constructor
    Leading_Log_Forward_V(const double beta,const double factor,const double yexponent,
			const std::string cinfo,ATOOLS::Integration_Info *info);

    // member functions
    void GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
		       const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);
    void AddPoint(double);
    void Optimize()  {p_vegas->Optimize();} 
    void EndOptimize()  {p_vegas->EndOptimize();} 
    void WriteOut(std::string pId) { p_vegas->WriteOut(pId); }
    void ReadIn(std::string pId)   { p_vegas->ReadIn(pId); }

  };// end of class Leading_Log_Forward

  class Leading_Log_Backward_V: public Single_Channel {
  protected :

    double m_beta, m_factor, m_yexponent;
    bool   m_zchannel;

    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey, m_sgridkey, m_ygridkey;
    Vegas* p_vegas;
  public :

    // constructor
    Leading_Log_Backward_V(const double beta,const double factor,const double yexponent,
			 const std::string cinfo,ATOOLS::Integration_Info *info);

    // member functions
    void GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
		       const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);
    void AddPoint(double);
    void Optimize()  {p_vegas->Optimize();} 
    void EndOptimize()  {p_vegas->EndOptimize();} 
    void WriteOut(std::string pId) { p_vegas->WriteOut(pId); }
    void ReadIn(std::string pId)   { p_vegas->ReadIn(pId); }

  };// end of class Leading_Log_Backward

  class Leading_Log_Central_V: public Single_Channel {
  protected:

    double m_beta, m_factor;
    bool   m_zchannel;

    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey, m_sgridkey, m_ygridkey;
    Vegas* p_vegas;
  public:

    // constructor
    Leading_Log_Central_V(const double beta,const double factor,
			const std::string cinfo,ATOOLS::Integration_Info *info,const int mode=0);

    // member functions
    void GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
		       const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);
    void AddPoint(double);
    void Optimize()  {p_vegas->Optimize();} 
    void EndOptimize()  {p_vegas->EndOptimize();} 
    void WriteOut(std::string pId) { p_vegas->WriteOut(pId); }
    void ReadIn(std::string pId)   { p_vegas->ReadIn(pId); }

  };// end of class Leading_Log_Central


  /////////////////////////////////////////////////////////////////////////

  class LBS_Compton_Peak_Uniform_V: public Single_Channel {
  protected:

    double m_exponent, m_pole;
    bool   m_zchannel;

    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey, m_sgridkey, m_ygridkey;
    Vegas* p_vegas;
  public:

    // constructor
    LBS_Compton_Peak_Uniform_V(const double exponent,const double pole,
			     const std::string cinfo,ATOOLS::Integration_Info *info);
    
    // member functions
    void GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
		       const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);
    void AddPoint(double);
    void Optimize()  {p_vegas->Optimize();} 
    void EndOptimize()  {p_vegas->EndOptimize();} 
    void WriteOut(std::string pId) { p_vegas->WriteOut(pId); }
    void ReadIn(std::string pId)   { p_vegas->ReadIn(pId); }

  };// end of class LBS_Compton_Peak_Uniform

  class LBS_Compton_Peak_Forward_V: public Single_Channel {
  protected :

    double m_exponent, m_pole, m_yexponent;
    bool   m_zchannel;

    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey, m_sgridkey, m_ygridkey;
    Vegas* p_vegas;
  public :

    // constructor
    LBS_Compton_Peak_Forward_V(const double exponent,const double pole,const double yexponent,
			     const std::string cinfo,ATOOLS::Integration_Info *info);

    // member functions
    void GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
		       const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);
    void AddPoint(double);
    void Optimize()  {p_vegas->Optimize();} 
    void EndOptimize()  {p_vegas->EndOptimize();} 
    void WriteOut(std::string pId) { p_vegas->WriteOut(pId); }
    void ReadIn(std::string pId)   { p_vegas->ReadIn(pId); }

  };// end of class LBS_Compton_Peak_Forward

  class LBS_Compton_Peak_Backward_V: public Single_Channel {
  protected :

    double m_exponent, m_pole, m_yexponent;
    bool   m_zchannel;

    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey, m_sgridkey, m_ygridkey;
    Vegas* p_vegas;
  public :

    // constructor
    LBS_Compton_Peak_Backward_V(const double exponent,const double pole,const double yexponent,
			      const std::string cinfo,ATOOLS::Integration_Info *info);

    // member functions
    void GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
		       const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);
    void AddPoint(double);
    void Optimize()  {p_vegas->Optimize();} 
    void EndOptimize()  {p_vegas->EndOptimize();} 
    void WriteOut(std::string pId) { p_vegas->WriteOut(pId); }
    void ReadIn(std::string pId)   { p_vegas->ReadIn(pId); }

  };// end of class LBS_Compton_Peak_Backward

  class LBS_Compton_Peak_Central_V: public Single_Channel {
  protected:

    double m_exponent, m_pole;
    bool   m_zchannel;

    ATOOLS::Info_Key m_spkey, m_xkey, m_ykey, m_sgridkey, m_ygridkey;
    Vegas* p_vegas;
  public:

    // constructor
    LBS_Compton_Peak_Central_V(const double exponent,const double pole,
			     const std::string cinfo,ATOOLS::Integration_Info *info,const int mode=0);

    // member functions
    void GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
		       const double *rans,const int mode=0);
    void GenerateWeight(const int mode=0);
    void AddPoint(double);
    void Optimize()  {p_vegas->Optimize();} 
    void EndOptimize()  {p_vegas->EndOptimize();} 
    void WriteOut(std::string pId) { p_vegas->WriteOut(pId); }
    void ReadIn(std::string pId)   { p_vegas->ReadIn(pId); }

  };// end of class LBS_Compton_Peak_Central

}// end of namespace PHASIC

#endif
