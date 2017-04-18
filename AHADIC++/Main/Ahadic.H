#ifndef AHADIC_Main_Ahadic_H
#define AHADIC_Main_Ahadic_H

#include "AHADIC++/Formation/Beam_Particles_Shifter.H"
#include "AHADIC++/Formation/Singlet_Former.H"
#include "AHADIC++/Formation/Singlet_Checker.H"
#include "AHADIC++/Formation/Gluon_Decayer.H"
#include "AHADIC++/Tools/Soft_Cluster_Handler.H"
#include "AHADIC++/Decays/Cluster_Decayer.H"
#include "AHADIC++/Tools/Cluster.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Org/Return_Value.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include <list>

namespace AHADIC {
  class Ahadic {
  private:
    std::list<Singlet *>        m_singlet_list; 
    std::list<Cluster *>        m_cluster_list; 
    std::list<Proto_Particle *> m_hadron_list; 

    Beam_Particles_Shifter m_beamparticles;
    Singlet_Former         m_sformer;
    Soft_Cluster_Handler   m_softclusters;
    Singlet_Checker        m_singletchecker;
    Gluon_Decayer          m_gluondecayer;
    Cluster_Decayer        m_clusterdecayer;
    

    bool ExtractSinglets(ATOOLS::Blob * blob);
    bool ShiftBeamParticles();
    bool CheckSinglets();
    bool DecayGluons();
    bool DecayClusters();
    void Reset();

    ATOOLS::Return_Value::code Hadronize(ATOOLS::Blob * blob, int retry=0);
    void FillOutgoingParticles(ATOOLS::Blob * blob);

    bool SanityCheck(ATOOLS::Blob * blob,double norm2);
    void CleanUp(ATOOLS::Blob * blob);
  public:
    Ahadic(std::string path=std::string("./"),
	   std::string file=std::string("Cluster.dat"));
    ~Ahadic();

    ATOOLS::Return_Value::code Hadronize(ATOOLS::Blob_List *);
  };
}
#endif
