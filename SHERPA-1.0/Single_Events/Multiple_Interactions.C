#include "Multiple_Interactions.H"
#include "Matrix_Element_Handler.H"

#ifdef PROFILE__Multiple_Interactions
#include "prof.hh"
#endif

using namespace SHERPA;

Multiple_Interactions::Multiple_Interactions(MI_Handler *_p_mihandler):
  m_one(true),
  p_mihandler(_p_mihandler)
{
  m_name=std::string("Multiple_Interactions : ")+p_mihandler->Name();
  m_type=std::string("Perturbative");
}

Multiple_Interactions::~Multiple_Interactions() {}

bool Multiple_Interactions::Treat(Blob_List *bloblist,double &weight)
{
#ifdef PROFILE__Multiple_Interactions
  PROFILE_HERE;
#endif
  if (p_mihandler->Type()==MI_Handler::None) return false;
  if (bloblist->empty()) {
    ATOOLS::msg.Error()<<"Potential error in Multiple_Interactions::Treat."<<std::endl
		       <<"   Incoming blob list is empty!"<<std::endl;
    return false;
  }
  Blob * myblob;
  bool hit=false;
  m_one=false;
  for (Blob_Iterator bit=bloblist->begin();bit!=bloblist->end();++bit) {
    if ((*bit)->Type()==std::string("Hard Subprocess : ")) {
      m_one=true;
      if ((*bit)->Status()==2) {
	myblob=*bit;
	if (p_mihandler->GenerateHardProcess(myblob)) {
	  CompleteBlob(myblob);
	  myblob->SetStatus(1);
	  // weight=p_mihandler->Weight();
	  weight=1.0;
	  ATOOLS::Blob *blob = new ATOOLS::Blob();
	  blob->SetType(std::string("Hard Subprocess : "));
	  blob->SetStatus(2);
	  blob->SetId(bloblist->size());
	  bloblist->push_back(blob);
	  hit=true;
	  break;
	}
	else {
	  delete myblob;
	  bloblist->erase(bit--);
	  break;
	}
      }
      else if ((*bit)->Status()==-1) {
	myblob=(*bit);
	p_mihandler->SameHardProcess(myblob); 
	CompleteBlob(myblob);
	myblob->SetStatus(1);
	// weight=p_mihandler->Weight();
	weight=1.0;
	hit=true;
	break;
      }
    }
  }
  if (!m_one) {
    bool found=false;
    for (Blob_Iterator bit=bloblist->begin();bit!=bloblist->end();++bit) {
      if ((*bit)->Type().find("ME PS Interface (Sherpa, FS)")!=std::string::npos) {
	myblob=*bit;
	found=true;
	break;
      }
    }
    if (!found) {
      for (Blob_Iterator bit=bloblist->begin();bit!=bloblist->end();++bit) {
	if ((*bit)->Type().find("Signal Process")!=std::string::npos) {
	  myblob=*bit;
	  break;
	}
      }
    }
    if (myblob->Status()!=0) return false;
    m_pperpmax=1.0e37; 
    m_ecmsmax=ATOOLS::rpa.gen.Ecms();
    m_xmax[1]=m_xmax[0]=1.0;
    double val=0.0;
    ATOOLS::Vec4D cur, ptot;
    for (int i=0;i<myblob->NInP();++i) {
      m_xmax[i]-=2.0*myblob->InParticle(i)->Momentum()[0]/ATOOLS::rpa.gen.Ecms();
    }
    for (int i=0;i<myblob->NOutP();++i) {
      ptot+=cur=myblob->OutParticle(i)->Momentum();
      val=ATOOLS::Max(val,sqrt(cur[1]*cur[1]+cur[2]*cur[2]));
    }
    m_ecmsmax-=sqrt(ptot.Abs2());
    m_pperpmax=ATOOLS::Min(m_pperpmax,val);
    myblob = new ATOOLS::Blob();
    myblob->SetType(std::string("Hard Subprocess : "));
    myblob->SetStatus(2);
    myblob->SetId(bloblist->size());
    bloblist->push_back(myblob);
    p_mihandler->SetScaleMax(m_pperpmax,0);
    p_mihandler->SetScaleMax(m_ecmsmax,1);
    p_mihandler->SetScaleMax(m_xmax[0],2);
    p_mihandler->SetScaleMax(m_xmax[1],3);
    p_mihandler->Reset();
    hit=true;
  }
  return hit;
}

void Multiple_Interactions::CleanUp() 
{ 
  return; 
}

void Multiple_Interactions::CompleteBlob(Blob *blob)
{
#ifdef PROFILE__Multiple_Interactions
  PROFILE_HERE;
#endif
  EXTRAXS::XS_Base *xs=NULL;
  if (blob->Type()==std::string("Hard Subprocess : ")) {
    xs=p_mihandler->HardMEHandler()->GetXS();
  }
  if (blob->Type()==std::string("Soft Subprocess : ")) {
    xs=p_mihandler->SoftMEHandler()->GetXS();
  }
  blob->SetPosition(Vec4D(0.,0.,0.,0.));
//   blob->SetType(blob->Type()+xs->ProcessName());
  blob->SetStatus(1);

  Vec4D cms = Vec4D(0.,0.,0.,0.);
  for (int i=0;i<xs->Nin();i++) cms += xs->Momenta()[i];
  blob->SetCMS(cms);
  blob->SetBeam(-1);

  Particle *particle;
  for (int i=0;i<(int)blob->NInP();i++) {
    particle=blob->InParticle(i);
    particle->SetNumber((long int)particle);
    particle->SetStatus(2);
    particle->SetInfo('G');
  }
  for (int i=0;i<(int)blob->NOutP();i++) {
    particle=blob->OutParticle(i);
    particle->SetNumber((long int)particle);
    particle->SetStatus(1);
    particle->SetInfo('H');
  }
}
