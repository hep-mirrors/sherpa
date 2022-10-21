#include "PHASIC++/Channels/ISR_Channel_Interface.H"

#include "PHASIC++/Channels/Channel_Elements.H"
#include "PHASIC++/Channels/Multi_Channel.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"

using namespace PHASIC;
using namespace ATOOLS;

ISR_Channel_Interface::ISR_Channel_Interface
(Multi_Channel *const fsr,const std::string &cinfo,
 ATOOLS::Integration_Info *info):
  ISR_Channel_Base(info), p_fsr(fsr)
{
  m_name="ISR_Interface";
  m_spkey.SetInfo("ISR_Interface");
  m_ykey.SetInfo("ISR_Interface");
  m_spkey.Assign(cinfo+"::s'",5,0,info);
  m_ykey.Assign(cinfo+"::y",3,0,info);
  m_xkey.Assign(cinfo+"::x",6,0,info);
  m_sgridkey.Assign(m_spkey.Info(),1,0,info);
  m_ygridkey.Assign(m_ykey.Info(),1,0,info);
  p_vegas = new Vegas(0,100,"ISR_Interface");
}

void ISR_Channel_Interface::GeneratePoint(const double *rns)
{
  m_spkey[3]=p_fsr->Channels().front()->SPrime();
  m_ykey[2]=p_fsr->Channels().front()->Y();
}

void ISR_Channel_Interface::GenerateWeight(const int &mode)
{
  m_weight=1.0;
}

void ISR_Channel_Interface::AddPoint(double value)
{
  Single_Channel::AddPoint(value);
}
