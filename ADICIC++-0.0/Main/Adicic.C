#include "Adicic.H"
#include "Chain_Handler.H"
#include "Particle_List.H"
#include "Message.H"


using namespace ADICIC;
using namespace ATOOLS;
using namespace std;

Adicic::Adicic() :
  p_chain(new Chain())
{
}

Adicic::~Adicic()
{
  if (p_chain) { delete p_chain; p_chain = NULL; }
}

int  Adicic::PerformShowers() 
{
  Chain_Handler chainhandler((*p_chain));
  while(chainhandler.EvolveChainByOneStep()) {
    Vec4D test;
    p_chain->CheckMomentumConservation(test);
  }
  return 1;
}

bool Adicic::ExtractPartons(ATOOLS::Blob_List * blobs,ATOOLS::Particle_List *)
{
  Blob * blob;
  Particle_List plist;
  int upflow,length,count;
  std::cout<<"Blob list : "<<blobs->size()<<std::endl;
  for (Blob_Iterator blit=blobs->begin();blit!=blobs->end();++blit) {
    if((*blit)->Type()==btp::FS_Shower) {
      blob = * blit;
      for (int i=0;i<blob->NInP();++i) blob->InParticle(i)->SetStatus(2);
      upflow = blob->InParticle(0)->GetFlow(1);
      if (!p_chain->ExtractPartons(plist)) {
	msg.Error()<<"Error in Adicic::ExtractPartons."<<std::endl;
	return false;
      }
      count=1;
      length=plist.size();
      if(length>2) {
	for(Particle_Iterator piter=plist.begin(); piter!=plist.end(); ++piter) {
	  blob->AddToOutParticles(*piter);
	  (*piter)->SetInfo('F');
	  if(count==1) upflow=(*piter)->GetFlow(1);
	  else if(count<length) {
	    (*piter)->SetFlow(2,upflow);
	    (*piter)->SetFlow(1,-1);
	    upflow=(*piter)->GetFlow(1);
	  }
	  else (*piter)->SetFlow(2,upflow);
	  ++count;
	}
      }
      else {
	for(Particle_Iterator piter=plist.begin(); piter!=plist.end(); ++piter) {
	  blob->AddToOutParticles(*piter);
	  (*piter)->SetInfo('F');
	}
      }
      std::cout<<"Final Blob : "<<std::endl<<blob<<std::endl;
      return true;
    }
  }
  return false;
}







