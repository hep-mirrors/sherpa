#include "Ahadic.H"

using namespace AHADIC;
using namespace ATOOLS;

Ahadic::Ahadic(std::string _m_path,std::string _m_file) : 
  m_path(_m_path),m_file(_m_file), 
  p_clusters(NULL), p_constituents(NULL), p_hadrons(NULL),
  p_c2hadrons(NULL), p_cformer(NULL), p_cdecayer(NULL)
{
  hadpars.Init(m_path,m_file);
  p_c2hadrons    = new Clusters_2_Hadrons();
  p_cformer      = new Cluster_Formation_Handler(p_c2hadrons);
  p_cdecayer     = new Cluster_Decay_Handler(p_c2hadrons);
};


bool Ahadic::Hadronize(ATOOLS::Blob_List *bloblist,
		       ATOOLS::Particle_List *pl)
{ 
  //std::cout<<"In Ahadic::Hadronize("<<bloblist->size()<<")"<<std::endl;
  Blob * blob = p_cformer->FormClusters(bloblist);
  //std::cout<<(*blob)<<std::endl;
  p_clusters  = p_cformer->GetClusters();
  p_c2hadrons->Transition(p_clusters,blob);
  p_cdecayer->DecayThem(p_clusters,blob);
  //std::cout<<(*bloblist)<<std::endl;
  /*
    Vec4D check = Vec4D(0.,0.,0.,0.);
    for (int i=0;i<blob->NInP();i++) check = check+blob->InParticle(i)->Momentum();
    std::cout<<check<<std::endl;
    check = Vec4D(0.,0.,0.,0.);
    for (int i=0;i<blob->NOutP();i++) check = check+blob->OutParticle(i)->Momentum();
    std::cout<<check<<std::endl;
  */
  return true; 
}

