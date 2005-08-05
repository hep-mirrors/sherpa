#include "SimpleXS_CSS_Interface.H"

#include "CS_Shower.H"
#include "Run_Parameter.H"
#include "XS_Base.H"
#include "Exception.H"

using namespace SHERPA;
using namespace CS_SHOWER;
using namespace EXTRAXS;
using namespace ATOOLS;


SimpleXS_CSS_Interface::SimpleXS_CSS_Interface(Matrix_Element_Handler *_p_mehandler,
					       Shower_Handler *_p_shower) :
  Perturbative_Interface(_p_mehandler,_p_shower)
{ }

SimpleXS_CSS_Interface::~SimpleXS_CSS_Interface() 
{ }

int SimpleXS_CSS_Interface::DefineInitialConditions(Blob * blob) 
{
  if (blob==NULL) return false;
  if ((blob->NInP()!=2) || (blob->NOutP()!=2)) {
    THROW(fatal_error,"Cannot handle blobs with more than 4 legs.");
  }
  p_shower->CleanUp();
  return InitColours(blob);
}

bool SimpleXS_CSS_Interface::InitColours(Blob * blob) 
{
  XS_Base * xs = p_mehandler->GetXS();
  if (!(xs->SetColours(p_mehandler->Momenta()))) return false;
  for (int j=0;j<2;j++) {
    for (int i=0;i<blob->NInP();++i) {
      blob->InParticle(i)->SetFlow(j+1,xs->Colours()[i][j]);
    }
    for (int i=0;i<blob->NOutP();++i) blob->OutParticle(i)->SetFlow(j+1,xs->Colours()[i+blob->NInP()][j]);
  }
  if(Valid(blob)) {
    double scale = xs->Scale(PHASIC::stp::as);
    p_blob = new Blob();
    p_blob->AddToInParticles(blob->OutParticle(0));
    p_blob->AddToInParticles(blob->OutParticle(1));
    p_blob->SetStatus(1);
    p_blob->SetType(btp::Shower);
    p_blob->SetTypeSpec("CSS++0.0");
    Singlet * singlet = new Singlet(blob->OutParticle(0)->GetFlow(1),true);
    Parton  * parton  = new Parton(blob->OutParticle(0),pst::FS);
    parton->SetStart(scale/4.);
    parton->SetVeto(scale/4.);
    singlet->push_back(parton);
    parton            = new Parton(blob->OutParticle(1),pst::FS);
    parton->SetStart(scale/4.);
    parton->SetVeto(scale/4.);
    singlet->push_back(parton);
    p_shower->GetAllSinglets()->push_back(singlet);
    return true;
  }
  return false;
}

bool SimpleXS_CSS_Interface::Valid(Blob* blob) {
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


bool SimpleXS_CSS_Interface::FillBlobs(Blob_List * blobs)
{
  p_blob->SetId();
  blobs->push_back(p_blob);
  return true;
}

int SimpleXS_CSS_Interface::PerformShowers()
{
  return p_shower->PerformShowers(false,1,1.,1.,rpa.gen.Ycut());
}









