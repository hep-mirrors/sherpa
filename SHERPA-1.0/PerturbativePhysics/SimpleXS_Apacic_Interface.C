#include "SimpleXS_Apacic_Interface.H"

#include "Exception.H"
#include "XS_Selector.H"

using namespace SHERPA;

SimpleXS_Apacic_Interface::SimpleXS_Apacic_Interface(Matrix_Element_Handler *mehandler,
						     Shower_Handler *showerhandler):
  Perturbative_Interface(mehandler,showerhandler),
  p_tools(new Interface_Tools(showerhandler->GetIniTrees(),showerhandler->GetFinTree())),
  p_twototwo(new EXTRAXS::XS_Group(2,2,"Core Processes")),
  p_momenta(new ATOOLS::Vec4D[4]), 
  p_flavours(new ATOOLS::Flavour[4]), 
  p_psme_is(NULL),
  p_psme_fs(NULL) {}

SimpleXS_Apacic_Interface::~SimpleXS_Apacic_Interface() 
{
  delete p_twototwo;
  delete [] p_flavours;
  delete [] p_moms;
  delete p_tools;
}

int SimpleXS_Apacic_Interface::DefineInitialConditions(ATOOLS::Blob *blob) 
{
  if (blob==NULL) return false;
  if ((blob->NInP()!=2) || (blob->NOutP()!=2)) {
    ATOOLS::msg.Error()<<*blob;
    throw(ATOOLS::Exception(ATOOLS::ex::fatal_error,"Cannot handle blobs with more than 4 legs.",
			    "SimpleXS_Apacic_Interface","DefineInitialConditions"));
  }
  for (size_t i=0;i<(size_t)blob->NInP();++i) {
    p_flavours[i]=blob->InParticle(i)->Flav();
    p_momenta[i]=blob->InParticle(i)->Momentum();
  }
  for (size_t i=0;i<(size_t)blob->NOutP();++i) {
    p_flavours[i+blob->NInP()]=blob->OutParticle(i)->Flav();
    p_momenta[i+blob->NInP()]=blob->OutParticle(i)->Momentum();
  }
  p_xs=NULL;
  if (p_twototwo->XSSelector()->
      FindInGroup(p_twototwo,p_xs,blob->NInP(),2,p_flavours)==std::string::npos) {
    p_xs = p_twototwo->XSSelector()->GetXS(blob->NInP(),2,p_flavours);
    if (p_xs) p_twototwo->Add(p_xs);
  }
  p_hard=blob;
  return InitColours(blob);
}

int SimpleXS_Apacic_Interface::InitColours(ATOOLS::Blob *blob) 
{
  if (!p_xs->SetColours(p_momenta)) return 0; 
  if (blob->InParticle(0)->Momentum()[3]<blob->InParticle(1)->Momentum()[3]) {
    blob->SwapInParticles(0,1);
    blob->SwapOutParticles(0,1);
  }
  for (int j=0;j<2;j++) {
    for (int i=0;i<blob->NInP();++i) {
      blob->InParticle(i)->SetFlow(j+1,p_xs->Colours()[i][j]);
    }
    for (int i=0;i<blob->NOutP();++i) {
      blob->OutParticle(i)->SetFlow(j+1,p_xs->Colours()[i+blob->NInP()][j]);
    }
  }
  double scale=
    dynamic_cast<PHASIC::Integrable_Base*>(p_xs)->Scale(PHASIC::stp::fac);
  double E=sqrt(p_mehandler->GetISR_Handler()->Pole())/2.;
  double th1,th2;
  if (m_ini) {
    th1=p_tools->ColourAngle(blob->InParticle(0),blob);
    th2=p_tools->ColourAngle(blob->InParticle(1),blob);
    double x1=blob->InParticle(0)->Momentum()[0]/E;
    double x2=blob->InParticle(1)->Momentum()[0]/E;
    p_tools->InitializeIncoming(blob,scale,E,th1,th2,x1,x2);
  }
  if (m_fin) {
    th1=p_tools->ColourAngle(blob->OutParticle(0),blob);
    th2=p_tools->ColourAngle(blob->OutParticle(1),blob);
    p_tools->InitializeOutGoing(blob,scale,E,th1,th2);
  }
  return 1;
}

bool SimpleXS_Apacic_Interface::FillBlobs(ATOOLS::Blob_List *blobs)
{
  if (p_hard==NULL) return false;
  if (p_shower->ISROn()) {
    p_psme_is = new ATOOLS::Blob();
    p_psme_is->SetType(ATOOLS::btp::ME_PS_Interface_IS);
    p_psme_is->SetTypeSpec("Sherpa");
    p_psme_is->SetStatus(1);
    p_psme_is->SetId();
    for (int i=0;i<p_hard->NInP();++i) {
      p_psme_is->AddToOutParticles(p_hard->InParticle(i));
    }
    blobs->push_back(p_psme_is);
  }
  if (p_shower->FSROn()) {
    p_psme_fs = new ATOOLS::Blob();
    p_psme_fs->SetType(ATOOLS::btp::ME_PS_Interface_FS);
    p_psme_fs->SetTypeSpec("Sherpa");
    p_psme_fs->SetStatus(1);
    p_psme_fs->SetId();
    for (int i=0;i<p_hard->NOutP();++i) {
      p_psme_fs->AddToInParticles(p_hard->OutParticle(i));
    }
    blobs->push_back(p_psme_fs);
  }
  return true;
}

int SimpleXS_Apacic_Interface::PerformShowers()
{
  return p_shower->PerformShowers(-1,p_mehandler->GetISR_Handler()->X1(),
				  p_mehandler->GetISR_Handler()->X2());
}

int SimpleXS_Apacic_Interface::PerformDecayShowers()
{ 
  return 1; 
}

void SimpleXS_Apacic_Interface::Reset()
{
  m_rescale[1]=m_rescale[0]=1.;
}

void SimpleXS_Apacic_Interface::UpdateEnergy(const double energy,const size_t i)
{
  m_rescale[i]-=2.*energy/sqrt(p_mehandler->GetISR_Handler()->Pole());
}
