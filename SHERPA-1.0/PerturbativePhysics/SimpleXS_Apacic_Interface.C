#include "SimpleXS_Apacic_Interface.H"

#include "Interface_Tools.H"
#include "Run_Parameter.H"
#include "XS_Selector.H"
#include "Exception.H"

using namespace SHERPA;
using namespace EXTRAXS;
using namespace ATOOLS;

namespace ATOOLS {

  std::ostream &operator<<(std::ostream &str,const XS_Base *&xs)
  {
    return str<<"("<<xs<<")";
  }

}

template XS_Base *&Blob_Data_Base::Get<XS_Base*>();

namespace ATOOLS {

template <> Blob_Data<XS_Base*>::~Blob_Data() {}

template class Blob_Data<XS_Base*>;

}

SimpleXS_Apacic_Interface::SimpleXS_Apacic_Interface(Matrix_Element_Handler *mehandler,
						     Shower_Handler *showerhandler):
  Perturbative_Interface(mehandler,showerhandler),
  p_tools(new Interface_Tools(showerhandler->GetIniTrees(),
			      showerhandler->GetFinTree())),
  p_twototwo(new XS_Group(2,2,"Interface Processes")),
  p_momenta(new Vec4D[4]), 
  p_flavours(new Flavour[4]), 
  p_psme_is(NULL),
  p_psme_fs(NULL) 
{
  p_twototwo->SetScaleScheme(mehandler->GetXS()->ScaleScheme());
}

SimpleXS_Apacic_Interface::~SimpleXS_Apacic_Interface() 
{
  delete p_twototwo;
  delete [] p_flavours;
  delete [] p_moms;
  delete p_tools;
}

Return_Value::code SimpleXS_Apacic_Interface::DefineInitialConditions(ATOOLS::Blob *blob) 
{
  if (blob==NULL) {
    return Return_Value::Error;
    msg.Error()<<"ERROR in "<<METHOD<<" : "<<std::endl
	       <<"   No blob found, will return 'Error' and hope for the best."<<std::endl;
  }
  if ((blob->NInP()!=2) || (blob->NOutP()!=2)) {
    msg.Error()<<"ERROR in "<<METHOD<<" : "<<std::endl
	       <<"   Cannot handle blobs with different than four legs:"<<std::endl
	       <<(*blob)
	       <<"   Will return 'Error' and hope for the best."<<std::endl;
    return Return_Value::Error;
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
  if (p_twototwo->XSSelector()->FindInGroup
      (p_twototwo,p_xs,blob->NInP(),2,p_flavours)==std::string::npos) {
    p_xs = p_twototwo->XSSelector()->GetXS(blob->NInP(),2,p_flavours);
    if (p_xs) p_twototwo->Add(p_xs);
  }
  p_hard=blob;
  p_shower->CleanUp();
  return InitColours(blob);
}

Return_Value::code SimpleXS_Apacic_Interface::InitColours(ATOOLS::Blob *blob) 
{
  if (!p_xs->SetColours(p_momenta)) {
    msg.Error()<<"ERROR in "<<METHOD<<" : "<<std::endl
	       <<"   Could not set initial colour in "<<std::endl
	       <<*blob
	       <<"   Return 'Error' and hope for the best."<<std::endl;
    return Return_Value::Error;
  } 
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
  m_scale=p_xs->Scale(PHASIC::stp::as);
  double E=sqrt(p_mehandler->GetISR_Handler()->Pole())/2.;
  if (m_ini) p_tools->InitializeIncoming(blob,E);
  if (m_fin) p_tools->InitializeOutGoing(blob,E);
  return Return_Value::Success;
}

bool SimpleXS_Apacic_Interface::FillBlobs(ATOOLS::Blob_List *blobs)
{
  if (p_hard==NULL) return false;
  if (p_shower->ISROn()) {
    p_psme_is = new Blob();
    p_psme_is->SetType(btp::ME_PS_Interface_IS);
    p_psme_is->SetTypeSpec("Sherpa");
    p_psme_is->SetStatus(blob_status::needs_showers);
    p_psme_is->SetId();
    for (int i=0;i<p_hard->NInP();++i) {
      p_psme_is->AddToOutParticles(p_hard->InParticle(i));
    }
    blobs->push_back(p_psme_is);
  }
  if (p_shower->FSROn()) {
    p_psme_fs = new Blob();
    p_psme_fs->SetType(btp::ME_PS_Interface_FS);
    p_psme_fs->SetTypeSpec("Sherpa");
    p_psme_fs->SetStatus(blob_status::needs_showers);
    p_psme_fs->SetId();
    for (int i=0;i<p_hard->NOutP();++i) {
      p_psme_fs->AddToInParticles(p_hard->OutParticle(i));
    }
    p_psme_fs->AddData("Core_Process",new Blob_Data<XS_Base*>(p_xs));
    p_psme_fs->AddData("OrderStrong",new 
		       Blob_Data<double>((double)p_xs->OrderStrong()));
    p_psme_fs->AddData("OrderEWeak",new
		       Blob_Data<double>((double)p_xs->OrderEW()));
    blobs->push_back(p_psme_fs);
  }
  p_shower->FillBlobs(blobs); 
  return true;
}

int SimpleXS_Apacic_Interface::PerformShowers()
{
  int jetveto=-1;
  double qmin2i=0., qmin2f=0.; 
  if (p_mehandler->UseSudakovWeight()) {
    p_tools->JetVetoPt2(qmin2i,qmin2f);
    p_shower->SetJetvetoPt2(qmin2i,qmin2f);
    double scale=p_mehandler->FactorisationScale();
    p_shower->SetFactorisationScale(scale);
    jetveto=2;
  }
  msg_Debugging()<<"SimpleXS_Apacic_Interface::PerformShowers(): {\n"
		 <<"   initial = "<<m_ini<<", final = "<<m_fin<<"\n"
		 <<"   sudakov weight = "<<p_mehandler->UseSudakovWeight()<<"\n"
		 <<"   maxpt ini = "<<qmin2i<<" maxpt fin = "<<qmin2f
		 <<" vs. "<<p_hard->OutParticle(0)->Momentum().PPerp2()
		 <<"\n}"<<std::endl;
  return p_shower->PerformShowers(jetveto,0,p_mehandler->GetISR_Handler()->X1(),
				  p_mehandler->GetISR_Handler()->X2(),rpa.gen.Ycut());
}

int SimpleXS_Apacic_Interface::PerformDecayShowers()
{ 
  return 1; 
}

