#include "Apacic.H"
#include "Initial_State_Shower.H"
#include "Final_State_Shower.H"
#include "Tree.H"

#include "ISR_Handler.H"

#include "Run_Parameter.H"
#include "Random.H"

using namespace APACIC;
using namespace PDF;
using namespace ATOOLS;


Apacic::Apacic(ISR_Handler * _isr,MODEL::Model_Base * _model,int _maxjetnumber,
	       bool _isron,bool _fsron,Data_Read * _dataread):
  m_isron(_isron), m_fsron(_fsron), m_showers(_isron||_fsron),
  p_inishower(NULL), p_finshower(NULL), p_initrees(NULL), p_fintree(NULL)
{
  if (m_fsron) {
    p_fintree   = new Tree();
    p_finshower = new Final_State_Shower(_model,_dataread);
  }
  if (m_isron) {
    p_initrees  = new Tree*[2];
    for (int i=0;i<2;i++) p_initrees[i] = new Tree();
    p_inishower = new Initial_State_Shower(_isr,p_finshower,_model,_dataread);
  }
}
  
Apacic::~Apacic() 
{
  if (p_fintree)     { delete p_fintree; p_fintree = 0; }
  if (p_initrees)    {
    for (int i=0;i<2;i++) { delete p_initrees[i]; p_initrees[i] = 0; }
    delete p_initrees;p_initrees = 0;
  }  
  if (p_inishower)   { delete p_inishower; p_inishower = 0; }
  if (p_finshower)   { delete p_finshower; p_finshower = 0; }
}

void Apacic::PrepareTrees() {
  if (m_fsron) p_fintree->Reset(); 
  if (m_isron) for (int i=0;i<2;i++) p_initrees[i]->Reset();
}

void Apacic::SetJetvetoPt2(const double q2i, const double q2f)
{ 
  if (m_fsron) 
    p_finshower->SetJetvetoPt2(q2f); 
  if (m_isron)
    p_inishower->SetJetvetoPt2(q2i,q2f); 
}

void Apacic::SetFactorisationScale(const double scale)
{
  if (m_isron) {
    p_inishower->SetFactorisationScale(scale);
  }
}

int Apacic::PerformShowers(bool ini,bool fin,int jetveto,double x1,double x2) {
  if (!m_showers) return 1;
  if (msg.LevelIsDebugging()) {
    if (m_fsron) 
      p_finshower->OutputTree(p_fintree);
    if (m_isron) {
      p_inishower->OutputTree(p_initrees[0]);
      p_inishower->OutputTree(p_initrees[1]);
    }
  }
  if (m_fsron) {
    Vec4D sum=p_fintree->GetRoot()->part->Momentum();
    Poincare cms(sum);
    p_fintree->BoRo(cms);

    int fsrstatus = p_finshower->PerformShower(p_fintree,jetveto);
    int number=0;
    p_finshower->GetMomentum(p_fintree->GetRoot(),number);

    if (fsrstatus==0) return fsrstatus;
    // check ME if still njet ME!
    // if isr is on, this check will be performed after the initial state shower
    if (!m_isron) {
      double ycut   = rpa.gen.Ycut();
      double s  = sqr(rpa.gen.Ecms());
      double dr2 = sqr(rpa.gen.DeltaR());
      SetJetvetoPt2(s*ycut,dr2*s*ycut);
      if (!p_finshower->ExtraJetCheck()) {
	return 3;
      }
    }
    p_finshower->SetAllColours(p_fintree->GetRoot());

    if (!m_isron) {
      //      Vec4D vl =Vec4D(x1+x2,0.,0.,x2-x1);
      Vec4D vl =Vec4D(sum[0],-1.*Vec3D(sum));
      Poincare lab(vl);
      p_fintree->BoRo(lab);
    }
  }

  if (m_isron) {
    p_inishower->InitShowerPT(p_initrees[0]->GetRoot()->maxpt2);
    if (!(p_inishower->PerformShower(p_initrees,jetveto))) return 0;


    // boost in Beam-System (using x1, x2 from initiator) done in is-shower

    // determine Boost+Rotate in cms(using moms from root1 und root2)
    Vec4D mom1=p_initrees[0]->GetRoot()->part->Momentum();
    Vec4D mom2=p_initrees[1]->GetRoot()->part->Momentum();

    Vec4D vl =Vec4D(mom1[0]+mom2[0], -1.*Vec3D(mom1+mom2));
    Poincare lab(vl);
    lab.BoostBack(mom1);
    lab.BoostBack(mom2);
    Poincare rot(Vec4D::ZVEC,mom1);
    rot.RotateBack(mom1);
    rot.RotateBack(mom2);

    // rotate and boost fs in beam-system "if (fin)"
    if (fin) {
      p_fintree->BoRo(rot);
      p_fintree->BoRo(lab);

      // check ME if still njet ME!
      double ycut   = rpa.gen.Ycut();
      double s  = sqr(rpa.gen.Ecms());
      double dr2 = sqr(rpa.gen.DeltaR());
      SetJetvetoPt2(s*ycut,dr2*s*ycut);
      if (!p_finshower->ExtraJetCheck()) {
	return 3;
      }
    }
  }

  if (msg.LevelIsDebugging()) {
    if (m_fsron) 
      p_finshower->OutputTree(p_fintree);
    if (m_isron) {
      p_inishower->OutputTree(p_initrees[0]);
      p_inishower->OutputTree(p_initrees[1]);
    }
  }
  return 1;
}

bool Apacic::ExtractPartons(bool ini,bool fin,Blob_List * bl,Particle_List * pl) {
  if (fin) p_finshower->ExtractPartons(p_fintree->GetRoot(),0,bl,pl);

  if (ini) {
    for (int i=0;i<2;i++) {
      p_inishower->ExtractPartons(p_initrees[i]->GetInitiator(),i,0,bl,pl);
    }
  }
  return 1;
}

