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
{
}

Return_Value::code HadronDecays_Apacic_Interface::DefineInitialConditions(ATOOLS::Blob *blob) 
{
  p_blob = new Blob();
  p_blob->SetType(btp::FS_Shower);
  p_blob->SetTypeSpec("APACIC++ 2.0");
  p_blob->SetStatus(blob_status::needs_hadronization);
  p_blob->SetId();
  unsigned int pos,refpos,compare;
  int          ref,comppos;
  bool         chain(false);
  list<Particle *> * singlet;
  Particle * part;
  for (int i=0;i<blob->NOutP();i++) {
    part = blob->OutParticle(i);
    p_blob->AddToInParticles(part);
    if (part->Info()=='f')      continue;
    if (!part->Flav().Strong()) continue;
    if (!chain && part->GetFlow(1)!=0) {
      pos   = i;
      chain = true;
      part->SetInfo('f');
      singlet = new list<Particle *>;
      singlet->push_back(part);
      m_singlets.push_back(singlet);
      if (part->Flav().IsGluon())   pos = 1;
      else if (part->Flav().IsQuark()) {
	if (!part->Flav().IsAnti()) pos = 1;
	                       else pos = 2;
      }
      else if (part->Flav().IsDiQuark()) {
	if (part->Flav().IsAnti())  pos = 1;
	                      else  pos = 2;
      }
      comppos = pos;
      compare = part->GetFlow(comppos);
      refpos  = 3-pos;
      ref     = part->GetFlow(refpos);
      do {
	for (int j=0;j<blob->NOutP();j++) {
	  part = blob->OutParticle(j);
	  if (part->Info()=='f')      continue;
	  if (!part->Flav().Strong()) continue;
	  if (part->GetFlow(3-comppos)==compare) {
	    part->SetInfo('f');
	    singlet->push_back(part);
	    if (part->GetFlow(3-refpos)==ref) chain = false;
	    break;
	  }
	}
      } while (chain);
    }
  }
  return Return_Value::Success;
}

bool HadronDecays_Apacic_Interface::FillBlobs(ATOOLS::Blob_List *blobs)
{
  //cout<<METHOD<<endl<<(*blobs)<<endl;

  blobs->push_back(p_blob);
  return true;
}

int HadronDecays_Apacic_Interface::PerformShowers()
{
  msg.Error()<<"ERROR in "<<METHOD<<" : "<<endl
	     <<"   Hadron Decay Blobs must have one incoming particle only,"<<endl
	     <<"   but seemingly there are two incoming particles."<<endl
	     <<"   Will abort the run."<<endl;
  abort();
}

int HadronDecays_Apacic_Interface::PerformDecayShowers()
{ 
  //cout<<METHOD<<endl;
  APACIC::Tree * tree = p_shower->GetFinTree();
  for (list<list<Particle * > * >::iterator slit=m_singlets.begin();
       slit!=m_singlets.end();) {
    if (FillTree((*slit))) {
      //cout<<(*tree)<<endl;
      p_shower->GetApacic()->FinShower()->PerformShower(tree,false); 
      p_shower->GetApacic()->FinShower()->SetAllColours(tree->GetRoot());
      m_boost.Invert();
      tree->BoRo(m_boost);
      //cout<<(*p_shower->GetFinTree())<<endl;
      p_shower->GetApacic()->FinShower()->ExtractFinalState(p_blob,tree->GetRoot());
      tree->Reset();
      (*slit)->clear();
      slit = m_singlets.erase(slit);
    }
  }
  return 1;
}


bool HadronDecays_Apacic_Interface::FillTree(list<Particle *> * singlet)
{
  APACIC::Tree * tree = p_shower->GetFinTree();
  tree->Reset();
  if (singlet->size()>3) {
    return false;
  }
  Vec4D momom(Vec4D(0.,0.,0.,0.));
  for (list<Particle * >::iterator sit=singlet->begin();
       sit!=singlet->end();sit++) momom += (*sit)->Momentum();
  double scale(momom.Abs2());

  APACIC::Knot * mo = tree->NewKnot(Flavour(kf::photon),momom,scale,0.);
  mo->part->SetStatus(part_status::decayed);
  mo->part->SetInfo('M');
  mo->E2      = scale;
  mo->maxpt2  = scale;
  mo->z       = 0.5;
  mo->costh   = -1.; 
  mo->thcrit  = 200.;
  mo->stat    = 1;  

  std::cout<<METHOD<<" : singlet with mass "<<sqrt(scale)<<std::endl;
  list<Particle * >::iterator sit=singlet->begin();
  if (singlet->size()==2) {
    mo->left             = tree->NewKnot((*sit));
    mo->left->prev       = mo;
    mo->left->stat       = 3;    
    mo->left->t          = mo->t;
    mo->left->tout       = sqr((*sit)->Flav().PSMass());
    mo->left->maxpt2     = 0.;
    mo->left->didkin     = true;
    (*sit)->SetInfo('F');
    sit = singlet->erase(sit);

    mo->right            = tree->NewKnot((*sit));
    mo->right->prev      = mo;
    mo->right->stat      = 3;    
    mo->right->t         = mo->t;
    mo->right->tout      = sqr((*sit)->Flav().PSMass());
    mo->right->maxpt2    = 0.;
    mo->right->didkin    = true;
    (*sit)->SetInfo('F');
    sit = singlet->erase(sit);

    m_boost = Poincare(momom);
    tree->BoRo(m_boost);
    mo->left->E2         = sqr(mo->left->part->Momentum()[0]);
    mo->right->E2        = sqr(mo->right->part->Momentum()[0]);
    return true;
  }
  else {
    msg.Error()<<"ERROR in "<<METHOD<<" : "<<endl
	       <<"   Singlets consisting of three partons not considered yet, "<<endl
	       <<"   will continue and hope for the best."<<endl;
  }
  return false;
}
