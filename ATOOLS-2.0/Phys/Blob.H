#ifndef Blob_h
#define Blob_h

#include <string>
#include <vector>
#include <iterator> 
#include <iostream> 

#include "Vector.H"
#include "Poincare.H"


namespace APHYTOOLS {
  class Parton; 

  typedef std::vector<Parton *> Parton_Queue;

  class Blob {
    friend std::ostream& operator<<( std::ostream&, const Blob &);
    friend std::ostream& operator<<( std::ostream&, const Blob *);
  private:
    int   m_id;
    int   m_status;
    int   m_beam;
    std::string         m_type;
    Parton_Queue        m_inpartons;
    Parton_Queue        m_outpartons;
    AMATOOLS::Vec4D     m_cms_vec;
    AMATOOLS::Poincare  m_cms_boost;
    AMATOOLS::Vec4D     m_position;
  public:
    Blob(const AMATOOLS::Vec4D _pos = AMATOOLS::Vec4D(0.,0.,0.,0.), const int _id=0) : 
      m_position(_pos), m_id(_id) { m_beam = -1; m_status = 0; };
    ~Blob();
    void     AddToInPartons(Parton *);
    void     AddToOutPartons(Parton *);
    Parton * OutParton(int);
    Parton * InParton(int);
    Parton * RemoveInParton(int);
    Parton * RemoveInParton(Parton *);
    Parton * RemoveOutParton(int);
    Parton * RemoveOutParton(Parton *);
    void     DeleteInParton(Parton *);
    void     DeleteOutParton(Parton *);
    void     RemoveOwnedPartons();
    void     DeleteOwnedPartons();

    void     BoostInCMS();
    void     BoostInLab();

    void     SetVecs();
    void     SetPosition(AMATOOLS::Vec4D pos) { m_position = pos; }
    void     SetCMS(AMATOOLS::Vec4D _cms)     { m_cms_vec  = _cms; }
    void     SetId(int _id)                   { m_id       = _id; }
    void     SetStatus(int _status)           { m_status   = _status; }
    void     SetType(std::string _type)       { m_type     = _type; }
    void     SetBeam(int _beam)               { m_beam     = _beam; }

    int      Id()                       const { return m_id; }
    int      Status()                   const { return m_status; }
    int      Beam()                     const { return m_beam; }
    int      NInP()                     const { return m_inpartons.size(); }
    int      NOutP()                    const { return m_outpartons.size(); }

    AMATOOLS::Vec4D    CheckMomentumConservation();
    AMATOOLS::Vec4D    const Position() const { return m_position; }
    AMATOOLS::Vec4D    const CMS()      const { return m_cms_vec; }
    std::string        const Type()     const { return m_type; }
  };


  /*!
    \file 
    \brief  contains the class APHYTOOLS::Blob
  */
  
  /*!
    \class Blob 
    \brief This class contains a point where a given number of incomming and outgoing Parton interact
    
    A typical Blob is the hard process at a center of a collision. Another Blob
    might be the transition between a hard partons and soft partons, i.e. the parton
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
    i = hadron to parton transition
  */
  
  /*!
    \fn void Blob::BoostInLab()
    \brief boost blob back in lab system
    
    \warning this is not jet implemented!
  */
}

#endif
  




