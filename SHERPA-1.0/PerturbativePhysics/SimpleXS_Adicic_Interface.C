#include "SimpleXS_Adicic_Interface.H"

#include "Adicic.H"
#include "Chain.H"
#include "XS_Base.H"
#include "Exception.H"

using namespace SHERPA;
using namespace ADICIC;
using namespace EXTRAXS;
using namespace ATOOLS;


SimpleXS_Adicic_Interface::SimpleXS_Adicic_Interface(Matrix_Element_Handler *_p_mehandler,
						     Shower_Handler *_p_shower) :
  Perturbative_Interface(_p_mehandler,_p_shower)
{ }

SimpleXS_Adicic_Interface::~SimpleXS_Adicic_Interface() 
{ }

int SimpleXS_Adicic_Interface::DefineInitialConditions(Blob * blob) 
{
  if (blob==NULL) return false;
  if ((blob->NInP()!=2) || (blob->NOutP()!=2)) {
    throw(ATOOLS::Exception(ATOOLS::ex::fatal_error,"Cannot handle blobs with more than 4 legs.",
			    "SimpleXS_Adicic_Interface","DefineInitialConditions"));
  }
  return InitColours(blob);
}

bool SimpleXS_Adicic_Interface::InitColours(Blob * blob) 
{
  XS_Base * xs = p_mehandler->GetXS();
  if (!(xs->SetColours(p_mehandler->Momenta()))) return false;
  for (int j=0;j<2;j++) {
    for (int i=0;i<blob->NInP();++i) {
      blob->InParticle(i)->SetFlow(j+1,xs->Colours()[i][j]);
    }
    for (int i=0;i<blob->NOutP();++i) blob->OutParticle(i)->SetFlow(j+1,xs->Colours()[i+blob->NInP()][j]);
  }
  Chain    * chain;
  if(Valid(blob)) {
    myblob = new Blob();
    myblob->AddToInParticles(blob->OutParticle(0));
    myblob->AddToInParticles(blob->OutParticle(1));
    myblob->SetStatus(1);
    myblob->SetType(btp::FS_Shower);
    myblob->SetTypeSpec("ADICIC++0.0");

    Dipole::Branch     q((*blob->OutParticle(0)));
    Dipole::Antibranch qbar((*blob->OutParticle(1)));
    chain = p_shower->GetChain();
    chain->Initialize(q,qbar);
    return 1;
  }
  return 0;
}

bool SimpleXS_Adicic_Interface::Valid(Blob* blob) {
  for(int i=0; i<blob->NInP(); ++i) {
    for(int j=1; j<3; j++) {
      if(blob->InParticle(i)->GetFlow(j)!=0) return false;
    }
  }
  if(blob->OutParticle(0)->Flav().IsAnti()) blob->SwapOutParticles(0,1);
  if(!(blob->OutParticle(0)->GetFlow(1)!=0 &&
       blob->OutParticle(1)->GetFlow(2)!=0)) return false;
  return true;
}


bool SimpleXS_Adicic_Interface::FillBlobs(Blob_List * blobs)
{
  myblob->SetId();
  blobs->push_back(myblob);
  return true;
}

int SimpleXS_Adicic_Interface::PerformShowers()
{
  return p_shower->PerformShowers(false,1.,1.);
}





