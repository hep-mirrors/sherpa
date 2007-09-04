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
  p_blob = new Blob();
  p_blob->SetType(btp::FS_Shower);
  p_blob->SetTypeSpec("APACIC++2.0");
  p_blob->SetStatus(blob_status::needs_hadronization);
  p_blob->SetId();
  unsigned int pos(0),refpos(0),compare(0);
  int          comppos(0);
  size_t       ref;
  bool         chain(false);
  Particle * part;

  m_N_of_singlets=0;
  m_particles.clear();
  for (int i=0;i<blob->NOutP();i++) {
    part = blob->OutParticle(i);
    p_blob->AddToInParticles(part);
    if (part->Info()=='f')      continue;
    if (!part->Flav().Strong()) continue;
    if (!chain && (part->GetFlow(1)!=0 || part->GetFlow(2)!=0)) {
      chain = true;
      m_N_of_singlets++;
      part->SetInfo('f');
      m_particles.push_back(part);
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
	    m_particles.push_back(part);
	    compare = part->GetFlow(comppos);
	  }
	  if (part->GetFlow(3-refpos)==ref) chain = false;
	}
      } while (chain);
    }
  }
  return Return_Value::Success;
}

bool HadronDecays_Apacic_Interface::FillBlobs(ATOOLS::Blob_List *blobs)
{
  blobs->push_back(p_blob);
  return true;
}

int HadronDecays_Apacic_Interface::PerformShowers()
{
  msg_Error()<<"ERROR in "<<METHOD<<" : "<<endl
	     <<"   Hadron Decay Blobs must have one incoming particle only,"<<endl
	     <<"   but seemingly there are two incoming particles."<<endl
	     <<"   Will abort the run."<<endl;
  abort();
}

int HadronDecays_Apacic_Interface::PerformDecayShowers()
{ 
  APACIC::Tree * tree = p_shower->GetFinTree();
  if (FillTree(tree)) {
    if (p_shower->GetApacic()->FinShower()->PerformShower(tree)==-1) {
      delete p_blob;
      return -1;
    }
    p_shower->GetApacic()->FinShower()->SetAllColours(tree->GetRoot());
    m_boost.Invert();
    tree->BoRo(m_boost);
    p_shower->GetApacic()->FinShower()->ExtractFinalState(p_blob,tree->GetRoot());
    tree->Reset();
    m_particles.clear();
  }
  return 1;
}


bool HadronDecays_Apacic_Interface::FillTree(APACIC::Tree * tree)
{
  switch (m_N_of_singlets) {
  case 1:
    if (m_particles.size()==2) return FillBinaryDecayTree(tree);
    if (m_particles.size()==3) return FillTertiaryDecayTree(tree);
    break;
  case 2:
    if (m_particles.size()==4) return FillSpectatorDecayTree(tree);
    break;
  }
  msg_Error()<<"Error in METHOD:"<<std::endl
	     <<"   Do not know how to deal with a "
	     <<m_N_of_singlets<<" singlets, "
	     <<m_particles.size()<<" particles situation."<<std::endl
	     <<"   Will return .false. and hope for the best."<<std::endl;
  return false;
}


bool HadronDecays_Apacic_Interface::FillBinaryDecayTree(APACIC::Tree * tree) {
  tree->Reset();

  Vec4D momom(Vec4D(0.,0.,0.,0.));
  for (int i=0;i<2;i++) momom += m_particles[i]->Momentum();
  double scale(momom.Abs2());

  APACIC::Knot * mo = tree->NewKnot(Flavour(kf::photon),momom,scale,0.);
  mo->part->SetStatus(part_status::decayed);
  mo->part->SetInfo('M');
  mo->E2      = scale;
  mo->maxpt2  = scale;
  mo->z       = 0.5;
  mo->costh   = -1.; 
  mo->thcrit  = 200.;
  mo->stat    = 0;  
  mo->shower  = 0;

  Particle * part      = m_particles[0];
  mo->left             = tree->NewKnot(part);
  mo->left->prev       = mo;
  mo->left->stat       = 3;    
  mo->left->t          = mo->t;
  mo->left->tout       = sqr(part->Flav().PSMass());
  mo->left->maxpt2     = 0.;
  mo->left->didkin     = true;
  mo->left->part->SetInfo('H');

  part                 = m_particles[1];
  mo->right            = tree->NewKnot(part);
  mo->right->prev      = mo;
  mo->right->stat      = 3;    
  mo->right->t         = mo->t;
  mo->right->tout      = sqr(part->Flav().PSMass());
  mo->right->maxpt2    = 0.;
  mo->right->didkin    = true;
  mo->right->part->SetInfo('H');

  m_boost = Poincare(momom);
  tree->BoRo(m_boost);
  mo->left->E2         = sqr(mo->left->part->Momentum()[0]);
  mo->right->E2        = sqr(mo->right->part->Momentum()[0]);
  mo->z                = sqrt(mo->left->E2/mo->E2);
  return true;
}


bool HadronDecays_Apacic_Interface::FillTertiaryDecayTree(APACIC::Tree * tree) {
  tree->Reset();

  Vec4D momom(Vec4D(0.,0.,0.,0.));
  for (int i=0;i<3;i++) momom += m_particles[i]->Momentum();
  double scale(momom.Abs2()), sijmin(scale), sij;

  int left(0),right1(1),right2(2);
  for (int i=0;i<2;i++) {
    for (int j=i+1;j<3;j++) {
      sij = (m_particles[i]->Momentum()+m_particles[j]->Momentum()).Abs();
      if (sij<sijmin) { 
	sijmin = sij; 
	right1 = i; 
	right2 = j;
	left   = 3-i-j;
      }
    }
  }

  APACIC::Knot * mo        = tree->NewKnot(Flavour(kf::photon),momom,scale,0.);
  mo->part->SetStatus(part_status::decayed);
  mo->part->SetInfo('M');
  mo->E2                   = scale;
  mo->maxpt2               = scale;
  mo->z                    = 0.5;
  mo->costh                = -1.; 
  mo->thcrit               = 200.;
  mo->stat                 = 0;  
  mo->shower               = 0;

  Particle * part          = m_particles[left];
  mo->left                 = tree->NewKnot(part);
  mo->left->prev           = mo;
  mo->left->stat           = 3;    
  mo->left->t              = mo->t;
  mo->left->tout           = sqr(part->Flav().PSMass());
  mo->left->maxpt2         = 0.;
  mo->left->didkin         = true;
  mo->left->part->SetInfo('H');

  Vec4D rightmom = m_particles[right1]->Momentum() + m_particles[right2]->Momentum();
  scale = rightmom.Abs2();

  mo->right                = tree->NewKnot(mo->left->part->Flav().Bar(),
					   rightmom,scale,0.);
  mo->right->part->SetStatus(part_status::decayed);
  mo->right->part->SetInfo('H');
  mo->right->part->SetFlow(1,mo->left->part->GetFlow(2));
  mo->right->part->SetFlow(2,mo->left->part->GetFlow(1));
  mo->right->E2            = scale;
  mo->right->maxpt2        = scale;
  mo->right->t             = scale;
  mo->right->tout          = scale;
  mo->right->z             = 0.5;
  mo->right->costh         = -1.; 
  mo->right->thcrit        = 200.;
  mo->right->prev          = mo;
  mo->right->stat          = 0;  
  mo->right->shower        = 0;
  mo->right->didkin        = true;

  part                     = m_particles[right1];
  mo->right->right         = tree->NewKnot(part);
  mo->right->right->prev   = mo->right;
  mo->right->right->stat   = 3;    
  mo->right->right->t      = mo->t;
  mo->right->right->tout   = sqr(part->Flav().PSMass());
  mo->right->right->maxpt2 = 0.;
  mo->right->right->didkin = true;
  mo->right->right->part->SetInfo('H');

  part                     = m_particles[right2];
  mo->right->left          = tree->NewKnot(part);
  mo->right->left->prev    = mo->right;
  mo->right->left->stat    = 3;    
  mo->right->left->t       = mo->t;
  mo->right->left->tout    = sqr(part->Flav().PSMass());
  mo->right->left->maxpt2  = 0.;
  mo->right->left->didkin  = true;
  mo->right->left->part->SetInfo('H');

  m_boost = Poincare(momom);
  tree->BoRo(m_boost);
  mo->right->left->E2      = sqr(mo->right->left->part->Momentum()[0]);
  mo->right->right->E2     = sqr(mo->right->right->part->Momentum()[0]);
  mo->right->E2            = sqr(sqrt(mo->right->left->E2)+sqrt(mo->right->right->E2));
  mo->right->z             = sqrt(mo->right->left->E2/mo->right->E2);
  mo->right->zs            = mo->right->z;
  mo->left->E2             = sqr(mo->left->part->Momentum()[0]);
  mo->right->E2            = sqr(mo->right->part->Momentum()[0]);
  mo->zs = mo->z           = sqrt(mo->left->E2/mo->E2);

  return true;
}


bool HadronDecays_Apacic_Interface::FillSpectatorDecayTree(APACIC::Tree * tree) {
  tree->Reset();

  Vec4D momom(Vec4D(0.,0.,0.,0.));
  for (int i=0;i<4;i++) momom += m_particles[i]->Momentum();
  double scale(momom.Abs2());

  APACIC::Knot * mo = tree->NewKnot(Flavour(kf::photon),momom,scale,0.);
  mo->part->SetStatus(part_status::decayed);
  mo->part->SetInfo('M');
  mo->E2      = scale;
  mo->maxpt2  = scale;
  mo->z       = 0.5;
  mo->costh   = -1.; 
  mo->thcrit  = 200.;
  mo->stat    = 0;  
  mo->didkin  = true;

  momom = Vec4D(0.,0.,0.,0.);
  for (int i=0;i<2;i++) momom += m_particles[i]->Momentum();
  scale = momom.Abs2();

  mo->left = tree->NewKnot(Flavour(kf::photon),momom,scale,0.);
  mo->left->part->SetStatus(part_status::decayed);
  mo->left->part->SetInfo('M');
  mo->left->t       = scale;
  mo->left->E2      = scale;
  mo->left->maxpt2  = scale;
  mo->left->z       = 0.5;
  mo->left->costh   = -1.; 
  mo->left->thcrit  = 200.;
  mo->left->prev    = mo;
  mo->left->stat    = 0;  
  mo->left->didkin  = true;

  Particle * part            = m_particles[0];
  mo->left->left             = tree->NewKnot(part);
  mo->left->left->prev       = mo->left;
  mo->left->left->stat       = 3;    
  mo->left->left->t          = mo->left->t;
  mo->left->left->tout       = sqr(part->Flav().PSMass());
  mo->left->left->maxpt2     = 0.;
  mo->left->left->didkin     = true;
  mo->left->left->part->SetInfo('H');

  part                       = m_particles[1];
  mo->left->right            = tree->NewKnot(part);
  mo->left->right->prev      = mo->left;
  mo->left->right->stat      = 3;    
  mo->left->right->t         = mo->left->t;
  mo->left->right->tout      = sqr(part->Flav().PSMass());
  mo->left->right->maxpt2    = 0.;
  mo->left->right->didkin    = true;
  mo->left->right->part->SetInfo('H');


  momom = Vec4D(0.,0.,0.,0.);
  for (int i=2;i<4;i++) momom += m_particles[i]->Momentum();
  scale = momom.Abs2();

  mo->right = tree->NewKnot(Flavour(kf::photon),momom,scale,0.);
  mo->right->part->SetStatus(part_status::decayed);
  mo->right->part->SetInfo('M');
  mo->right->t       = scale;
  mo->right->E2      = scale;
  mo->right->maxpt2  = scale;
  mo->right->z       = 0.5;
  mo->right->costh   = -1.; 
  mo->right->thcrit  = 200.;
  mo->right->prev    = mo;
  mo->right->stat    = 0;  
  mo->right->didkin  = true;

  part                        = m_particles[2];
  mo->right->left             = tree->NewKnot(part);
  mo->right->left->prev       = mo->right;
  mo->right->left->stat       = 3;    
  mo->right->left->t          = mo->right->t;
  mo->right->left->tout       = sqr(part->Flav().PSMass());
  mo->right->left->maxpt2     = 0.;
  mo->right->left->didkin     = true;
  mo->right->left->part->SetInfo('H');

  part                        = m_particles[3];
  mo->right->right            = tree->NewKnot(part);
  mo->right->right->prev      = mo->right;
  mo->right->right->stat      = 3;    
  mo->right->right->t         = mo->right->t;
  mo->right->right->tout      = sqr(part->Flav().PSMass());
  mo->right->right->maxpt2    = 0.;
  mo->right->right->didkin    = true;
  mo->right->right->part->SetInfo('H');


  m_boost = Poincare(mo->part->Momentum());
  tree->BoRo(m_boost);
  mo->left->left->E2           = sqr(mo->left->left->part->Momentum()[0]);
  mo->left->right->E2          = sqr(mo->left->right->part->Momentum()[0]);
  mo->left->E2                 = sqr(mo->left->part->Momentum()[0]);
  mo->left->zs = mo->left->z   = sqrt(mo->left->left->E2/mo->left->E2);
  mo->right->left->E2          = sqr(mo->right->left->part->Momentum()[0]);
  mo->right->right->E2         = sqr(mo->right->right->part->Momentum()[0]);
  mo->right->E2                = sqr(mo->right->part->Momentum()[0]);
  mo->right->zs = mo->right->z = sqrt(mo->right->left->E2/mo->right->E2);
  mo->zs = mo->z               = sqrt(mo->left->E2/mo->E2);

  return true;
}
