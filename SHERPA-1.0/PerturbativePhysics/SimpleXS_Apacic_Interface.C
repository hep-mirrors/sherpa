#include "SimpleXS_Apacic_Interface.H"
#include "XS_Base.H"

using namespace SHERPA;
using namespace EXTRAXS;
using namespace ATOOLS;


SimpleXS_Apacic_Interface::SimpleXS_Apacic_Interface(Matrix_Element_Handler *_p_mehandler,
						     Shower_Handler *_p_shower) :
  Perturbative_Interface(_p_mehandler,_p_shower),
  p_tools(NULL),
  p_psme_is(NULL),
  p_psme_fs(NULL)
{
  p_tools = new Interface_Tools(p_shower->GetIniTrees(),p_shower->GetFinTree());
}

SimpleXS_Apacic_Interface::~SimpleXS_Apacic_Interface() 
{
  if (p_tools) delete p_tools;
}

int SimpleXS_Apacic_Interface::DefineInitialConditions(ATOOLS::Blob * blob) 
{
  if (blob==NULL) return false;
  if ((blob->NInP()!=2) || (blob->NOutP()!=2)) {
    msg.Error()<<"Error in SimpleXS_Apacic_Interface::DefineInitialConditions."<<std::endl
	       <<"   Cannot handle blobs with "
	       <<blob->NInP()<<" -> "<<blob->NOutP()<<" legs."<<std::endl
	       <<"   Run cannot continue. "<<std::endl;
    exit(131);
  }
  if (p_shower->ISROn()) {
    p_psme_is = new Blob();
    p_psme_is->SetType(btp::ME_PS_Interface_IS);
    p_psme_is->SetTypeSpec("Sherpa");
    p_psme_is->SetStatus(1);
    for (int i=0;i<blob->NInP();++i) {
      p_psme_is->AddToOutParticles(blob->InParticle(i));
      p_psme_is->SetId(-1);
    }
  }
  if (p_shower->FSROn()) {
    p_psme_fs = new Blob();
    p_psme_fs->SetType(btp::ME_PS_Interface_FS);
    p_psme_fs->SetTypeSpec("Sherpa");
    p_psme_fs->SetStatus(1);
    for (int i=0;i<blob->NOutP();++i) {
      p_psme_fs->AddToInParticles(blob->OutParticle(i));
      p_psme_fs->SetId(-1);
    }
  }
  return InitColours(blob);
}

bool SimpleXS_Apacic_Interface::InitColours(Blob * blob) 
{
  XS_Base * xs = p_mehandler->GetXS();
  if (!(xs->SetColours(p_mehandler->Momenta()))) return false;
  for (int j=0;j<2;j++) {
    for (int i=0;i<blob->NInP();++i) {
      blob->InParticle(i)->SetFlow(j+1,xs->Colours()[i][j]);
    }
    for (int i=0;i<blob->NOutP();++i) blob->OutParticle(i)->SetFlow(j+1,xs->Colours()[i+blob->NInP()][j]);
  }
  double scale=dynamic_cast<PHASIC::Integrable_Base*>(xs)->Scale();
  double E=sqrt(p_mehandler->GetISR_Handler()->Pole())/2.;
  double th1,th2;
  if (m_ini) {
    th1=p_tools->ColourAngle(blob->InParticle(0),blob);
    th2=p_tools->ColourAngle(blob->InParticle(1),blob);
    p_tools->InitializeIncoming(blob,scale,E,th1,th2,p_mehandler->GetISR_Handler()->X1(),p_mehandler->GetISR_Handler()->X2());
  }
  if (m_fin) {
    th1=p_tools->ColourAngle(blob->OutParticle(0),blob);
    th2=p_tools->ColourAngle(blob->OutParticle(1),blob);
    p_tools->InitializeOutGoing(blob,scale,E,th1,th2);
  }
  return 1;
}

bool   SimpleXS_Apacic_Interface::FillBlobs(Blob_List *blobs)
{
  if (p_psme_is!=NULL) {
    p_psme_is->SetId(blobs->size());
    blobs->push_back(p_psme_is);
  }
  if (p_psme_fs!=NULL) {
    p_psme_fs->SetId(blobs->size());
    blobs->push_back(p_psme_fs);
  }
  return true;
}

int SimpleXS_Apacic_Interface::PerformShowers()
{
  return p_shower->PerformShowers(false,p_mehandler->GetISR_Handler()->X1(),
				  p_mehandler->GetISR_Handler()->X2());
}





