#include "AMEGIC++/Amplitude/Amplitude_Base.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Exception.H"

using namespace AMEGIC;
using namespace ATOOLS;

void Amplitude_Base::SetStringOn()   { buildstring = 1; }

void Amplitude_Base::SetStringOff()  { buildstring = 0; }

Point*  Amplitude_Base::GetPointlist() 
{
  msg_Error()<<"Error: Virtual "<<METHOD<<" called!"<<std::endl;
  return NULL;
}

void Amplitude_Base::Add(Amplitude_Base* a, int sign)
{
  msg_Error()<<"Error: Virtual "<<METHOD<<" called!"<<std::endl;
}

int Amplitude_Base::Size()
{
  return 1;
}
 
Zfunc_List* Amplitude_Base::GetZlist()
{
  msg_Error()<<"Error: Virtual "<<METHOD<<" called!"<<std::endl;
  return NULL;
}

Pfunc_List* Amplitude_Base::GetPlist()
{
  msg_Error()<<"Error: Virtual "<<METHOD<<" called!"<<std::endl;
  return NULL;
} 

int Amplitude_Base::GetSign()
{
  msg_Error()<<"Error: Virtual "<<METHOD<<" called!"<<std::endl;
  return 0;
}

void Amplitude_Base::SetSign(int)
{
  msg_Error()<<"Error: Virtual "<<METHOD<<" called!"<<std::endl;
}

void Amplitude_Base::BuildGlobalString(int* i,int j,Basic_Sfuncs* BS,
				       ATOOLS::Flavour* fl,
				       String_Handler* shand)
{
  msg_Error()<<"Error: Virtual "<<METHOD<<" called!"<<std::endl;
}

void Amplitude_Base::DefineOrder(const std::vector<int> &o)
{
  Abort();
}

const std::vector<int> &Amplitude_Base::GetOrder()
{
  THROW(fatal_error, "Error: Amplitude_Base::GetOrder() is not implemented.");
}
