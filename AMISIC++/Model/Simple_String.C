#include "AMISIC++/Model/Simple_String.H"

#include "ATOOLS/Org/Data_Reader.H"
#include "PDF/Remnant/Hadron_Remnant.H"
#include "AMISIC++/Model/Reggeon_Trajectory.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "BEAM/Main/Beam_Base.H"
#include "ATOOLS/Math/Random.H"

#ifdef PROFILE__all
#define PROFILE__Simple_String
#endif
#ifdef PROFILE__Simple_String
#include "prof.hh"
#else
#define PROFILE_HERE
#endif

using namespace AMISIC;
using namespace ATOOLS;

Simple_String::Simple_String():
  MI_Base("Simple String",MI_Base::SoftEvent,5,4,1)
{
  SetInputFile("MI.dat");
  m_start[0]=1.0;
  m_stop[0]=0.0;
  m_start[3]=m_start[2]=0.0;
  m_stop[3]=m_stop[2]=0.0;
  THROW(fatal_error,"Simple_String needs ISR_Handler");
}

Simple_String::Simple_String(PDF::ISR_Handler *const isr):
  MI_Base("Simple String",MI_Base::SoftEvent,5,4,1),
  p_isr(isr)
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
  while (m_reggeons.size()>0) {
    delete m_reggeons.front();
    m_reggeons.erase(m_reggeons.begin());
  }
}

bool Simple_String::Initialize()
{
  PROFILE_HERE;
  CleanUp();
  if (InputPath()=="" && InputFile()=="") return false;
  if (!rpa.gen.Beam1().IsHadron() ||
      !rpa.gen.Beam2().IsHadron()) return false;
  Data_Reader *reader = new Data_Reader(" ",";","!","=");
  reader->AddComment("#");
  reader->AddWordSeparator("\t");
  reader->SetInputPath(InputPath());
  reader->SetInputFile(InputFile());
  std::vector<std::vector<std::string> > helpsvv;
  if (!reader->MatrixFromFile(helpsvv,"REGGE_TRAJECTORY")) {
    helpsvv.push_back(std::vector<std::string>(3));
    helpsvv.back()[0]="Pomeron";
    helpsvv.back()[1]="1.0808";
    helpsvv.back()[2]="0.25";
  }
  msg_Info()<<"Simple_String::Initialize(): Adding Reggeon {\n";
  for (size_t i=0;i<helpsvv.size();++i) {
    if (helpsvv[i].size()<3) continue;
    m_reggeons.push_back(new Reggeon_Trajectory
			 (ToType<double>(helpsvv[i][1]),
			  ToType<double>(helpsvv[i][2])));
    m_reggeons.back()->SetS(sqr(rpa.gen.Ecms()));
    msg_Info()<<"   "<<std::setw(10)<<helpsvv[i][0]
	      <<" "<<std::setw(8)<<helpsvv[i][1]
	      <<" "<<std::setw(8)<<helpsvv[i][2]<<"\n";
  }
  msg_Info()<<"}"<<std::endl;
  p_remnants[0]=p_isr->GetRemnant(0);
  p_remnants[1]=p_isr->GetRemnant(1);
  return true;
}

bool Simple_String::CreateMomenta()
{
  PROFILE_HERE;
  m_filledblob=false;
  if (p_remnants[0]==NULL || p_remnants[1]==NULL) {
    msg_Error()<<"Simple_String::CreateMomenta(): "
	       <<"No remnant found."<<std::endl;
    return false;
  }
  m_reggeons[0]->Fit(m_start[0]*m_start[0],m_start[2]);
  m_start[1]=sqrt(m_reggeons[0]->GetT(0.0,m_start[0]*m_start[0],ran.Get()));
  const unsigned int flow=Flow::Counter();
  for (short unsigned int i=0;i<2;++i) {
    PDF::Hadron_Remnant *hadron=
      dynamic_cast<PDF::Hadron_Remnant*>(p_remnants[i]);
    if (hadron==NULL) {
      msg_Error()<<"Simple_String::CreateMomenta(): "
		 <<"Incoming particle is no hadron."<<std::endl;
      return false;
    }
    const std::vector<Flavour> &constit=
      hadron->GetConstituents(kf_none);
    double pz=0.0, phi=ran.Get()*2.0*M_PI;
    for (size_t j=0;j<constit.size();++j) {
      if (constit[j].IsQuark() && constit[j].IsAnti()==i) {
	Particle *particle = new Particle(0,constit[j]);
	do {
	  double E=hadron->GetBeam()->Energy()*
	    hadron->GetXPDF(constit[j],m_start[0]*m_start[0]);
	  pz=sqrt(E*E-sqr(constit[j].Mass())-m_start[1]*m_start[1]);
	  if (i==1) pz*=-1.0;
	  particle->SetMomentum(Vec4D(E,(i==0?1.0:-1.0)*m_start[1]*cos(phi),
					      (i==0?1.0:-1.0)*m_start[1]*sin(phi),pz));
	} while (!(dabs(pz)>0.0));
	particle->SetFlow(1+constit[j].IsAnti(),flow);
 	particle->SetFlow(2-constit[j].IsAnti(),0);
	particle->SetStatus(part_status::active);
	m_inparticles.push_back(particle);
	m_outparticles.push_back(particle);
	break;
      }
    }
  }
  m_filledblob=true;
  return true;
}

bool Simple_String::GenerateProcess()
{
  s_stopsoft=true;
  return m_generatedprocess=CreateMomenta();
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
