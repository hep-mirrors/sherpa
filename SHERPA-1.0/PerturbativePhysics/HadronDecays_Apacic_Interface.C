#include "HadronDecays_Apacic_Interface.H"
#include "Tree.H"

using namespace SHERPA;
using namespace EXTRAXS;
using namespace ATOOLS;
using namespace std;

HadronDecays_Apacic_Interface::HadronDecays_Apacic_Interface(Matrix_Element_Handler *mehandler,
							     Shower_Handler *showerhandler):
  Perturbative_Interface(NULL,showerhandler)
{ }

HadronDecays_Apacic_Interface::~HadronDecays_Apacic_Interface() 
{ }

Return_Value::code HadronDecays_Apacic_Interface::DefineInitialConditions(ATOOLS::Blob *blob) 
{
  DEBUG_FUNC("blob->Id()="<<blob->Id());
  p_decayblob=blob;
  p_blob = new Blob();
  p_blob->SetType(btp::FS_Shower);
  p_blob->SetTypeSpec("APACIC++2.0");
  p_blob->SetStatus(blob_status::needs_hadronization);
  p_blob->SetId();
  Particle * part;

  for (int i=0;i<blob->NOutP();i++) {
    part = blob->OutParticle(i);
    p_blob->AddToInParticles(part);
    if (part->Info()=='f')      continue;
    if (!part->Flav().Strong() && !part->Flav().IsDiQuark()) continue;
    if (part->GetFlow(1)!=0 || part->GetFlow(2)!=0) {
      part->SetInfo('f');
      int pos(0);
      if (part->Flav().IsGluon())   pos = 1;
      else if (part->Flav().IsQuark()) {
	if (!part->Flav().IsAnti()) pos = 1;
	                       else pos = 2;
      }
      else if (part->Flav().IsDiQuark()) {
	if (part->Flav().IsAnti())  pos = 1;
	                      else  pos = 2;
      }
      int comppos = pos;
      size_t compare = part->GetFlow(comppos);
      int refpos  = 3-pos;
      size_t ref = part->GetFlow(refpos);
      bool chain(true);
      do {
	for (int j=0;j<blob->NOutP();j++) {
	  part = blob->OutParticle(j);
	  if (part->Info()=='f')      continue;
	  if (!part->Flav().Strong()) continue;
	  if (part->GetFlow(3-comppos)==compare) {
	    part->SetInfo('f');
	    compare = part->GetFlow(comppos);
            if(part->GetFlow(3-refpos)==ref) {
              chain=false;
              break;
            }
	  }
	}
      } while (chain);
    }
  }
  p_shower->CleanUp();
  DEBUG_INFO("succeeded.");
  return Return_Value::Success;
}

bool HadronDecays_Apacic_Interface::FillBlobs(ATOOLS::Blob_List *blobs)
{
  for (int i=0;i<p_blob->NInP();i++) {
    p_blob->InParticle(i)->SetStatus(part_status::decayed);
  }
  blobs->push_back(p_blob);
  return true;
}

int HadronDecays_Apacic_Interface::PerformShowers()
{
  THROW(fatal_error, "Hadron Decay Blobs must have one incoming particle only,"+
        string(" but seemingly there are two incoming particles.\n")+
        "Will abort the run.");
}

int HadronDecays_Apacic_Interface::PerformDecayShowers()
{ 
  DEBUG_FUNC("started");
  APACIC::Tree * tree = p_shower->GetFinTree();
  double scale = sqr(p_decayblob->InParticle(0)->Flav().PSMass());
  Interface_Tools::InitializeNOutGoing(p_decayblob,tree,scale);
  m_boost = Poincare(p_decayblob->InParticle(0)->Momentum());
  tree->BoRo(m_boost);
  if (p_shower->GetApacic()->FinShower()->PerformShower(tree)==-1) {
    delete p_blob;
    return -1;
  }
  p_shower->GetApacic()->FinShower()->SetAllColours(tree->GetRoot());
  m_boost.Invert();
  tree->BoRo(m_boost);
  p_shower->GetApacic()->FinShower()->ExtractFinalState(p_blob,tree->GetRoot());
  tree->Reset();

  DEBUG_INFO("succeeded.");
  return 1;
}
