#ifndef Blob_h
#define Blob_h

#include <string>
#include <vector>
#include <iterator> 
#include <iostream> 

#include "Vector.H"
#include "Poincare.H"


namespace ATOOLS {
  class Particle; 

  typedef std::vector<Particle *> Particle_Queue;

  class Blob {
    friend std::ostream& operator<<( std::ostream&, const Blob &);
    friend std::ostream& operator<<( std::ostream&, const Blob *);
  private:
    Vec4D     m_position;
    int                 m_id;
    double              m_weight;
    int                 m_status,m_beam;
    std::string         m_type;
    std::vector<double> m_weightcontainer;
    Particle_Queue        m_inparticles, m_outparticles;
    Vec4D     m_cms_vec;
    Poincare  m_cms_boost;
  public:
    Blob(const Vec4D _pos = Vec4D(0.,0.,0.,0.), const int _id=-1) : 
      m_position(_pos), m_id(_id), m_weight(1.) { m_beam = -1; m_status = 0; };
    ~Blob();
    void     AddToInParticles(Particle *);
    void     AddToOutParticles(Particle *);
    Particle_Queue GetOutParticles() { return m_outparticles; }
    Particle_Queue GetInParticles()  { return m_inparticles;  }
    Particle * OutParticle(int);
    Particle * InParticle(int);
    Particle * RemoveInParticle(int,bool = true);
    Particle * RemoveInParticle(Particle *,bool = true);
    Particle * RemoveOutParticle(int,bool = true);
    Particle * RemoveOutParticle(Particle *,bool = true);
    void     DeleteInParticle(Particle *);
    void     DeleteOutParticle(Particle *);
    void     RemoveOwnedParticles();
    void     DeleteOwnedParticles();

    void     BoostInCMS();
    void     BoostInLab();

    void     SetVecs();
    void     SetPosition(Vec4D pos)  { m_position = pos; }
    void     SetCMS(Vec4D _cms)      { m_cms_vec  = _cms; }
    void     SetCMS();
    void     SetId(int _id)                    { m_id       = _id; }
    void     SetStatus(int _status)            { m_status   = _status; }
    void     SetType(std::string _type)        { m_type     = _type; }
    void     SetBeam(int _beam)                { m_beam     = _beam; }
    void     SetWeight(double _weight)         { 
      m_weightcontainer.clear();
      m_weightcontainer.push_back(_weight);
      m_weight   = _weight; 
    }
    void     AddPartialWeight(double _weight)  { 
      m_weightcontainer.push_back(_weight);
      m_weight *= _weight;
    }

    int      Id()                        const { return m_id; }
    int      Status()                    const { return m_status; }
    int      Beam()                      const { return m_beam; }
    int      NInP()                      const { return m_inparticles.size(); }
    int      NOutP()                     const { return m_outparticles.size(); }
 
    double              const Weight()   const { return m_weight; }
    std::vector<double> WeightContainer()      { return m_weightcontainer; } 
    Vec4D       CheckMomentumConservation();
    Vec4D       const Position()         const { return m_position; }
    Vec4D       const CMS()              const { return m_cms_vec; }
    std::string const Type()             const { return m_type; }
  };


  /*!
    \file 
    \brief  contains the class Blob
  */
  
  /*!
    \class Blob 
    \brief This class contains a point where a given number of incomming and outgoing Particle interact
    
    A typical Blob is the hard process at a center of a collision. Another Blob
    might be the transition between a hard particles and soft particles, i.e. the particle
    shower.
  */
  
  /*!
    \var int Blob::m_id 
    \brief contains an unique number for each Blob in an event.
  */
  
  /*!
    \var   char Blob::m_type;
    \brief Classifies the type of blob.
    
    Info about the blob type:
    H = primary hard interaction - i.e. signal
    h = secondary hard interaction, for instance underlying event vertex
    D = hard decay, like for instance top decay
    d = soft decay, like in fragmentation
    F = final state shower for H
    I = initial state shower for H
    i = hadron to particle transition
  */
  
  /*!
    \fn void Blob::BoostInLab()
    \brief boost blob back in lab system
    
    \warning this is not jet implemented!
  */
}

#endif
  




