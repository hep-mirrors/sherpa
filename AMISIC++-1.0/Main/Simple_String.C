#include "Simple_String.H"

#include "Data_Reader.H"

#ifdef PROFILE__all
#define PROFILE__Simple_String
#endif
#ifdef PROFILE__Simple_String
#include "prof.hh"
#else
#define PROFILE_HERE
#endif

#ifdef DEBUG__Simple_String
#include "Debugger.H"
#endif

using namespace AMISIC;

Simple_String::Simple_String():
  MI_Base("Simple String",MI_Base::SoftEvent,1,1,0)
{
  SetInputFile("MI.dat");
}

Simple_String::~Simple_String()
{
  CleanUp();
}

void Simple_String::CleanUp() 
{
}

bool Simple_String::Initialize()
{
  PROFILE_HERE;
  if (!CheckInputPath()) return false;
  if (!CheckInputFile()) return false;
  CleanUp();
  ATOOLS::Data_Reader *reader = new ATOOLS::Data_Reader("=",";","!");
  reader->SetInputPath(InputPath());
  reader->SetInputFile(InputFile());
  reader->SetVectorType(reader->VHorizontal);
  return true;
}

bool Simple_String::FillBlob(ATOOLS::Blob *blob)
{
  PROFILE_HERE;
  m_filledblob=false;
  ATOOLS::msg.Error()<<"Simple_String::FillBlob(..): "
		     <<"Cannot create momentum configuration."<<std::endl;
  return false;
}

void Simple_String::Reset()
{
}

void Simple_String::Update(const MI_Base *mibase)
{
  return;
}

void Simple_String::PrepareTerminate() 
{
}
