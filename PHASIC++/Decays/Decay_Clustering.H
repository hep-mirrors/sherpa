#ifndef Decay_Clustering_h
#define Decay_Clustering_h

#include "ATOOLS/Phys/Particle.H"

namespace ATOOLS {
  class Blob_List;
  class Blob;
  template<typename Scalar>
  class Lorentz_Ten2;
  typedef Lorentz_Ten2<double>  Lorentz_Ten2D;
  class Cluster_Amplitude;
}

namespace PHASIC {
  class Multi_Channel;
  class Single_Channel;

  typedef std::pair<ATOOLS::Particle *,ATOOLS::Particle *> ParticlePair;
  typedef std::vector<ParticlePair> ParticlePair_Vector;

  class Decay_Clustering {
    void AddDecayClustering(ATOOLS::Cluster_Amplitude*& ampl,
                            ATOOLS::Blob* blob,
                            size_t& imax,
                            size_t idmother);
    void AddPhotonsClustering(ATOOLS::Cluster_Amplitude*& ampl,
                              const ATOOLS::Particle_Vector daughters,
                              ParticlePair_Vector& photons,
                              size_t& imax,
                              const std::vector<size_t>& ids);
    void AssignPhotons(const ATOOLS::Particle_Vector& daughters,
                       ParticlePair_Vector& photons);
    ATOOLS::Vec4D RecombinedMomentum(const ATOOLS::Particle * daughter,
                                     const ParticlePair_Vector& photons,
                                     size_t& stat);
  public :
    Decay_Clustering() {}
    void DefineInitialShowerConditions(ATOOLS::Blob* initialblob, ATOOLS::Blob* showerblob);
  };
}

#endif
