#include "SimpleXS_Apacic_Interface.H"
#include "XS_Base.H"

using namespace SHERPA;
using namespace EXTRAXS;
using namespace ATOOLS;


SimpleXS_Apacic_Interface::SimpleXS_Apacic_Interface(Matrix_Element_Handler * _me,
						     Shower_Handler * _shower) :
  Perturbative_Interface(_me,_shower), p_tools(NULL)
{
  p_tools = new Interface_Tools(p_shower->GetIniTrees(),p_shower->GetFinTree());
}

SimpleXS_Apacic_Interface::~SimpleXS_Apacic_Interface() 
{
  if (p_tools) delete p_tools;
}

bool SimpleXS_Apacic_Interface::DefineInitialConditions(ATOOLS::Blob * blob) 
{
  if (!blob) return 0;
  blob->BoostInCMS();

  if ((blob->NInP()!=2) || (blob->NOutP()!=2)) {
    msg.Error()<<"Error in SimpleXS_Apacic_Interface::DefineInitialConditions."<<std::endl
	       <<"   No idea how to handle blobs with "
	       <<blob->NInP()<<" -> "<<blob->NOutP()<<" legs."<<std::endl
	       <<"   Abort run. "<<std::endl;
    abort();
  }
  return InitColours(blob);
}

bool SimpleXS_Apacic_Interface::InitColours(Blob * blob) {
  XS_Base * xs = p_me->GetXS();
  if (!(xs->SetColours(p_me->Momenta()))) return 0;
  for (int j=0;j<2;j++) {
    for (int i=0;i<2;i++) {
      blob->InParton(i)->SetFlow(j+1,xs->Colours()[i][j]);
      blob->OutParton(i)->SetFlow(j+1,xs->Colours()[i+2][j]);
    }
  }
 
  double scale = xs->Scale();
  double th1,th2;
  if (m_ini) {
    th1   = p_tools->ColourAngle(blob->InParton(0),blob);
    th2   = p_tools->ColourAngle(blob->InParton(1),blob);
    p_tools->InitializeIncoming(blob,scale,th1,th2);
  }
  if (m_fin) {
    th1   = p_tools->ColourAngle(blob->OutParton(0),blob);
    th2   = p_tools->ColourAngle(blob->OutParton(1),blob);
    p_tools->InitializeOutGoing(blob,scale,th1,th2);
  }
  return 1;
}

bool   SimpleXS_Apacic_Interface::FillBlobs(Blob_List * bl)
{
  return 1;
}






