#ifndef SHERPA_SoftPhysics_Hadron_Rescattering_Handler_H
#define SHERPA_SoftPhysics_Hadron_Rescattering_Handler_H

#include "HADRON_RESCATTERING/Main/Hadron_Rescatterings.H"

namespace SHERPA {
  struct Collision {
  public:
    double m_time, m_dist2;
    ATOOLS::Particle * p_A, * p_B;
    Collision(const double & time,const double & dist2,
	      ATOOLS::Particle * A,ATOOLS::Particle * B) :
      m_time(time), m_dist2(dist2), p_A(A), p_B(B) {}
  };
  
  class Hadron_Rescattering_Handler {
  private:
    HADRON_RESCATTERING::Hadron_Rescatterings  m_rescattering;
    std::string m_name;
    bool        m_on;

    std::set<ATOOLS::Blob * >     m_treatedblobs;
    std::set<ATOOLS::Particle * > m_particles;
    std::map<double,Collision * > m_collisions;

    void Schedule(ATOOLS::Particle * part1,
		  ATOOLS::Particle * part2);
  public :
    Hadron_Rescattering_Handler();
    ~Hadron_Rescattering_Handler();

    bool operator()();
    void HarvestParticles(ATOOLS::Blob * blob);
    void CleanUp(const size_t & mode=0);

    void Output();

    inline const std::string & Name() const { return m_name; }
  };
}

#endif

