#include "Cluster.H"
#include "Clusters_2_Hadrons.H"
#include "Cluster_Decay_Handler.H"
#include "Cluster_Formation_Handler.H"


#include "Particle_List.H"
#include "Particle.H"
#include "Blob_List.H"
#include "Blob.H"
#include "Rambo.H"
#include "Vector.H"
#include "Random.H"

using namespace AHADIC;
using namespace ATOOLS;

int main(int argc,char* argv[]) 
{  
  ATOOLS::ParticleInit(std::string("./"));

  Particle_List * pl  = new Particle_List;
  Blob_List     * bl  = new Blob_List;
  Blob          * blob;
  Cluster_List  * clusters;

  AHADIC::hadpars.Init(std::string("./"),std::string("Cluster.dat"));
  Clusters_2_Hadrons        * c2hadrons    = new Clusters_2_Hadrons();
  Cluster_Formation_Handler * cformer      = new Cluster_Formation_Handler(c2hadrons);
  Cluster_Decay_Handler     * cdecayer     = new Cluster_Decay_Handler(c2hadrons);

  int nin = 2, nout = 2+8; 

  Flavour * flavs   = new Flavour[nin+nout];
  flavs[0]          = Flavour(kf::e);     
  flavs[1]          = Flavour(kf::e).Bar();
  flavs[nin]        = Flavour(kf::c);
  flavs[nin+nout-1] = Flavour(kf::c).Bar();
  for (int i=nin+1;i<nin+nout-1;i++) {
    flavs[i]        = Flavour(kf::gluon);
  }
  PHASIC::Rambo rambo(nin,nout,flavs);
  Vec4D * vectors     = new Vec4D[nin+nout];
  vectors[0]          = Vec4D(45.6,0.,0.,45.6);
  vectors[1]          = Vec4D(45.6,0.,0.,-45.6);

  for (int nevents=0;nevents<2;nevents++) {
    Blob     * hard = new Blob;
    Particle * part;
    for (int i=0;i<nin;i++) {
      part = new Particle(0,flavs[i],vectors[i]);
      part->SetNumber(0);
      hard->AddToInParticles(part);
    }
    rambo.GeneratePoint(vectors,NULL);
    int col1 = 500;
    for (int i=nin;i<nin+nout;i++) {
      part = new Particle(0,flavs[i],vectors[i]);
      part->SetNumber(0);
      if (flavs[i].IsQuark() && !flavs[i].IsAnti()) part->SetFlow(1,col1);  
      if (flavs[i].IsQuark() &&  flavs[i].IsAnti()) part->SetFlow(2,col1++);  
      if (flavs[i].IsGluon()) {
	part->SetFlow(2,col1); 
	part->SetFlow(1,++col1); 
      } 
      hard->AddToOutParticles(part);
      pl->push_back(part);
    }
    hard->SetStatus(1);
    hard->SetId();
    hard->SetType(btp::FS_Shower);
    bl->push_back(hard);
    
    blob     = cformer->FormClusters(bl);
    clusters = cformer->GetClusters();
    c2hadrons->Transition(clusters,blob);
    cdecayer->DecayThem(clusters,blob);
    std::cout<<(*bl)<<std::endl;

    Vec4D check = Vec4D(0.,0.,0.,0.);
    for (int i=0;i<blob->NInP();i++) check = check+blob->InParticle(i)->Momentum();
    std::cout<<check<<std::endl;
    check = Vec4D(0.,0.,0.,0.);
    for (int i=0;i<blob->NOutP();i++) check = check+blob->OutParticle(i)->Momentum();
    std::cout<<check<<std::endl;

    if (!(*bl).empty()) {
      for (Blob_Iterator blit=(*bl).begin();blit!=(*bl).end();++blit) delete (*blit);
      (*bl).clear();
    }

    if (Particle::Counter()>2 || Blob::Counter()!=0) 
      msg.Error()<<"Error in Main while cleaning up the event : "<<std::endl
		 <<"   After event : "<<Particle::Counter()<<" / "<<Blob::Counter()
		 <<" particles / blobs undeleted !"<<std::endl
		 <<"   Continue and hope for the best."<<std::endl;
    Blob::Reset();
    Particle::Reset();
    Flow::ResetCounter();
  }


  delete cformer;

  std::cout<<"Finished."<<std::endl;
}
