#include "Simple_String.H"

#include "Data_Reader.H"
#include "Hadron_Remnant.H"
#include "Object.H"

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
  MI_Base("Simple String",MI_Base::SoftEvent,5,4,1)
{
  SetInputFile("MI.dat");
  m_start[0]=1.0;
  m_stop[0]=0.0;
  m_start[3]=m_start[2]=0.0;
  m_stop[3]=m_stop[2]=0.0;
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
  CleanUp();
  //   if (!CheckInputPath()) return false;
  //   if (!CheckInputFile()) return false;
  //   ATOOLS::Data_Reader *reader = new ATOOLS::Data_Reader("=",";","!");
  //   reader->SetInputPath(InputPath());
  //   reader->SetInputFile(InputFile());
  //   reader->SetVectorType(reader->VHorizontal);
  p_remnants[0]=GET_OBJECT(SHERPA::Remnant_Base,"Remnant_Base_0");
  p_remnants[1]=GET_OBJECT(SHERPA::Remnant_Base,"Remnant_Base_1");
  if (p_remnants[0]==NULL || p_remnants[1]==NULL) {
    throw(ATOOLS::Exception(ATOOLS::ex::fatal_error,"No beam remnant handler found.",
			    "Simple_String","Initialize"));
  }
  return true;
}

bool Simple_String::FillBlob(ATOOLS::Blob *blob)
{
  PROFILE_HERE;
  m_filledblob=false;
  if (p_remnants[0]==NULL || p_remnants[1]==NULL) {
    ATOOLS::msg.Error()<<"Simple_String::FillBlob(..): "
		       <<"No remnant found."<<std::endl;
    return false;
  }
  blob->DeleteOwnedParticles();
  const unsigned int flow=ATOOLS::Flow::Counter();
  for (short unsigned int i=0;i<2;++i) {
    SHERPA::Hadron_Remnant *hadron=dynamic_cast<SHERPA::Hadron_Remnant*>(p_remnants[i]);
    if (hadron==NULL) {
      ATOOLS::msg.Error()<<"Simple_String::FillBlob(..): "
			 <<"Incoming particle is no hadron."<<std::endl;
      return false;
    }
    const std::vector<ATOOLS::Flavour> &constit=
      hadron->GetConstituents(ATOOLS::kf::none);
    for (size_t j=0;j<constit.size();++j) {
      if (constit[j].IsQuark() && constit[j].IsAnti()==i) {
	ATOOLS::Particle *particle = new ATOOLS::Particle(0,constit[j]);
	double E=hadron->GetXPDF(constit[j],m_start[0]*m_start[0]);
	double pz=sqrt(E*E-ATOOLS::sqr(constit[j].Mass()));
	if (i==1) pz*=-1.0;
	particle->SetMomentum(ATOOLS::Vec4D(E,0.0,0.0,pz));
	particle->SetFlow(1+constit[j].IsAnti(),flow);
	particle->SetFlow(2-constit[j].IsAnti(),0);
	particle->SetStatus(1);
	blob->AddToInParticles(particle);
	blob->AddToOutParticles(particle);
	break;
      }
    }
  }
  m_filledblob=true;
  return true;
}

bool Simple_String::DiceProcess()
{
  s_stopsoft=true;
  return m_dicedprocess=FillBlob(p_blob);
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

bool Simple_String::VetoProcess(ATOOLS::Blob *blob)
{
  return false;
}
