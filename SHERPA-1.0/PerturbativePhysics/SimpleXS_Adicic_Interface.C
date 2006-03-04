#include "SimpleXS_Adicic_Interface.H"

#include "Adicic.H"
#include "Run_Parameter.H"
#include "XS_Base.H"
#include "Exception.H"

using namespace SHERPA;
using namespace ADICIC;
using namespace EXTRAXS;
using namespace ATOOLS;


SimpleXS_Adicic_Interface::
SimpleXS_Adicic_Interface(Matrix_Element_Handler* _p_mehandler,
			  Shower_Handler* _p_shower) :
  Perturbative_Interface(_p_mehandler,_p_shower),
  m_2to2id(0), m_scale(0.0), p_hard(NULL),
  p_isme(NULL), p_mefs(NULL), p_is(NULL), p_fs(NULL) {}



SimpleXS_Adicic_Interface::~SimpleXS_Adicic_Interface() {}



int SimpleXS_Adicic_Interface::DefineInitialConditions(Blob * blob) {
  if (blob==NULL) return 0;
  if ((blob->NInP()!=2) || (blob->NOutP()!=2)) {
    THROW(fatal_error,"Cannot handle blobs with more than 4 legs.");
  }
  p_hard=blob;
  //p_shower->CleanUp();
  return InitColours(blob);
}



bool SimpleXS_Adicic_Interface::InitColours(Blob * blob) {
  XS_Base* xs=p_mehandler->GetXS();
  if(!(xs->SetColours(p_mehandler->Momenta()))) return false;
  if(blob->InParticle(0)->Momentum()[3]<blob->InParticle(1)->Momentum()[3]) {
    blob->SwapInParticles(0,1);
    blob->SwapOutParticles(0,1);
  }
  for(int j=0; j<2; ++j) {
    for (int i=0; i<blob->NInP(); ++i)
      blob->InParticle(i)->SetFlow(j+1,xs->Colours()[i][j]);
    for (int i=0; i<blob->NOutP(); ++i)
      blob->OutParticle(i)->SetFlow(j+1,xs->Colours()[i+blob->NInP()][j]);
  }
  m_2to2id=Valid(blob);
  //std::cout<<"VALID="<<m_2to2id<<std::endl;//////////////////////////////////
  switch(m_2to2id) {
  case 1: {
    Dipole::Branch     q((*blob->OutParticle(0)));
    Dipole::Antibranch qbar((*blob->OutParticle(1)));
    assert(p_shower->GetCascade().AddChain(q,qbar));
    return true;
  }
  case -1: {
    m_scale=xs->Scale(PHASIC::stp::as);
    Dipole::Antibranch q(*blob->InParticle(0));    //Incoming quark!
    Dipole::Branch     qbar(*blob->InParticle(1));    //Incoming antiquark!
    dpv.evo.SetChainParticleLimit(5);
    dpv.evo.SetChainCorrelationLimit(5);
    dpv.sud.SetMaxIIScale(m_scale);
    //dpv.sud.SetMaxIIScale();/////////////////////////////////////////////////
    //std::cout<<sqrt(dpa.sud.MaxIIK2t())<<std::endl;//////////////////////////
    assert(p_shower->GetCascade().AddChain(qbar,q,
					   dpa.HighestIIInvScale(m_scale)));
    //assert(p_shower->GetCascade().AddChain(q,qbar,
    //                                       dpa.MaxIIInvScale(m_scale)));
    return true;
  }
  case -2: {
    m_scale=xs->Scale(PHASIC::stp::as);
    Dipole::Antibranch q(*blob->InParticle(1));    //Incoming quark!
    Dipole::Branch     qbar(*blob->InParticle(0));    //Incoming antiquark!
    dpv.evo.SetChainParticleLimit(5);
    dpv.evo.SetChainCorrelationLimit(5);
    dpv.sud.SetMaxIIScale(m_scale);
    //dpv.sud.SetMaxIIScale();/////////////////////////////////////////////////
    //std::cout<<sqrt(dpa.sud.MaxIIK2t())<<std::endl;//////////////////////////
    assert(p_shower->GetCascade().AddChain(qbar,q,
					   dpa.HighestIIInvScale(m_scale)));
    //assert(p_shower->GetCascade().AddChain(q,qbar,
    //                                       dpa.MaxIIInvScale(m_scale)));
    return true;
  }
  default: return false;
  }
}



int SimpleXS_Adicic_Interface::Valid(Blob* blob) {
  bool incol=false;
  bool oucol=false;
  for(int i=0; i<blob->NInP(); ++i) {
    for(int j=1; j<3; j++) {
      if(blob->InParticle(i)->GetFlow(j)!=0) { incol=true; break;}
    }
  }
  for(int i=0; i<blob->NOutP(); ++i) {
    for(int j=1; j<3; j++) {
      if(blob->OutParticle(i)->GetFlow(j)!=0) { oucol=true; break;}
    }
  }
  if(incol==false && oucol==true) {
    if(blob->OutParticle(0)->Flav().IsAnti()) blob->SwapOutParticles(0,1);
    if(blob->OutParticle(0)->GetFlow(2)==0 &&
       blob->OutParticle(1)->GetFlow(1)==0 &&
       blob->OutParticle(0)->GetFlow(1)!=0 &&
       blob->OutParticle(0)->GetFlow(1)==blob->OutParticle(1)->GetFlow(2))
      return 1;
    return 0;
  }
  if(incol==true && oucol==false) {
    if(blob->InParticle(0)->Flav().IsAnti()) {
      if(blob->InParticle(0)->GetFlow(1)==0 &&
	 blob->InParticle(1)->GetFlow(2)==0 &&
	 blob->InParticle(0)->GetFlow(2)!=0 &&
	 blob->InParticle(0)->GetFlow(2)==blob->InParticle(1)->GetFlow(1)) {
	return -2;
      }
      return 0;
    } else {
      if(blob->InParticle(0)->GetFlow(2)==0 &&
	 blob->InParticle(1)->GetFlow(1)==0 &&
	 blob->InParticle(0)->GetFlow(1)!=0 &&
	 blob->InParticle(0)->GetFlow(1)==blob->InParticle(1)->GetFlow(2)) {
	return -1;
      }
      return 0;
    }
  }
  return 0;
}



bool SimpleXS_Adicic_Interface::FillBlobs(Blob_List * blobs) {
  switch(m_2to2id) {
  case 1: {
    p_fs=new Blob(); assert(p_fs);
    p_fs->AddToInParticles(p_hard->OutParticle(0));
    p_fs->AddToInParticles(p_hard->OutParticle(1));
    p_fs->SetStatus(1);
    p_fs->SetType(btp::FS_Shower);
    p_fs->SetTypeSpec("ADICIC++0.0");
    p_fs->SetId();
    blobs->push_back(p_fs);
    return true;
  }
  case -1: {}
  case -2: {
    p_fs = new Blob(); assert(p_fs);
    for(int i=0; i<p_hard->NOutP(); ++i) {
      p_fs->AddToInParticles(p_hard->OutParticle(i));
    }
    p_fs->SetStatus(1);
    p_fs->SetType(btp::FS_Shower);
    p_fs->SetTypeSpec("Sherpa");
    p_fs->SetId();
    blobs->push_back(p_fs);
    p_isme=new Blob(); assert(p_isme);
    for(int i=0; i<p_hard->NInP(); ++i) {
      p_isme->AddToOutParticles(p_hard->InParticle(i));
    }
    p_isme->SetStatus(1);
    p_isme->SetType(btp::ME_PS_Interface_IS);
    p_isme->SetTypeSpec("Sherpa");
    p_isme->SetId();
    blobs->push_front(p_isme);
    p_is=new Blob(); assert(p_is);
    p_is->SetStatus(1);
    p_is->SetType(btp::IS_Shower);
    p_is->SetTypeSpec("ADICIC++0.0");
    p_is->SetId();
    blobs->push_front(p_is);
    //std::cout<<" => after interface:\n";
    //for(ATOOLS::Blob_List::const_iterator blit=blobs->begin();
    //	blit!=blobs->end(); ++blit)
    //    std::cout<<**blit<<"\n";
    //std::cout<<" => now leaving.\n";
    return true;
  }
  default: return false;
  }
}



int SimpleXS_Adicic_Interface::PerformShowers() {
  return p_shower->PerformShowers(false,1,1.,1.,rpa.gen.Ycut());
}
