#ifndef Blob_H
#define Blob_H

#include <string>
#include <vector>
#include <map>

#include "Poincare.H"
#include "Particle.H"

namespace ATOOLS {

  class Blob_Data_Base {
  public:
    template <class Type> Type Get();
    virtual ~Blob_Data_Base();
    virtual std::ostream & operator>>(std::ostream &) const =0;
  };

  template <class Type> class Blob_Data : public Blob_Data_Base {
    Type m_data;
  public:
    Blob_Data(const Type & d);
    Type Get() { return m_data; }
    void Set(const Type & d) { m_data=d; }
    std::ostream & operator>>(std::ostream &) const;
    ~Blob_Data();
  };

  std::ostream& operator<<(std::ostream&,const Blob_Data_Base &);

  typedef std::map<std::string,Blob_Data_Base *> String_BlobDataBase_Map;

  struct btp {
    enum code {
      Unspecified         =    0,
      Signal_Process      =    1,
      Hard_Decay          =    2,
      Hard_Collision      =    3,
      Soft_Collision      =    4,
      ME_PS_Interface_IS  =    9,
      ME_PS_Interface_FS  =   10,
      FS_Shower           =   11,
      IS_Shower           =   12,
      Beam                =  101,
      Bunch               =  102,
      Fragmentation       =  201,
      Cluster_Formation   =  202,
      Cluster_Decay       =  250,
      Hadron_Decay        = 1001
    };
  };// end of struct btp

  std::ostream& operator<<(std::ostream&,const btp::code);

  class Blob {
    friend std::ostream& operator<<( std::ostream&, const Blob &);
    friend std::ostream& operator<<( std::ostream&, const Blob *);
  private:
    Vec4D                   m_position;
    int                     m_id;
    double                  m_weight;
    int                     m_status,m_beam;
    bool                    m_hasboost;
    btp::code               m_type;
    std::string             m_typespec;
    std::vector<double>     m_weightcontainer;
    String_BlobDataBase_Map m_datacontainer;
    Particle_Vector         m_inparticles, m_outparticles;
    Vec4D                   m_cms_vec;
    Poincare                m_cms_boost;
  public:
    Blob(const Vec4D _pos = Vec4D(0.,0.,0.,0.), const int _id=-1);
    ~Blob();
    void     AddToInParticles(Particle *);
    void     AddToOutParticles(Particle *);
    Particle_Vector GetOutParticles() { return m_outparticles; }
    Particle_Vector GetInParticles()  { return m_inparticles;  }
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
    void     SetPosition(Vec4D pos)            { m_position = pos; }
    void     SetCMS(Vec4D _cms)                { m_cms_vec  = _cms; }
    void     SetCMS();
    void     SetId(int _id)                    { m_id       = _id; }
    void     SetStatus(int _status)            { m_status   = _status; }
    void     SetType(btp::code _type)          { m_type     = _type; }
    void     SetTypeSpec(std::string _type)    { m_typespec = _type; }
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
    Blob_Data_Base * operator[](const std::string name) 
    {
      String_BlobDataBase_Map::const_iterator cit=m_datacontainer.find(name);
      if (cit==m_datacontainer.end()) return 0;
      return cit->second;
    } 
    void     AddData(const std::string name, Blob_Data_Base * data); 
    void     ClearAllData();

    int      Id()                        const { return m_id; }
    int      Status()                    const { return m_status; }
    int      Beam()                      const { return m_beam; }
    int      NInP()                      const { return m_inparticles.size(); }
    int      NOutP()                     const { return m_outparticles.size(); }
 
    double              const Weight()   const { return m_weight; }
    std::vector<double> WeightContainer()      { return m_weightcontainer; } 
    Vec4D        CheckMomentumConservation();
    const Vec4D& Position()              const { return m_position; }
    const Vec4D& CMS()                   const { return m_cms_vec; }
    btp::code   const Type()             const { return m_type; }
    std::string const TypeSpec()         const { return m_typespec; }
    void SwapInParticles(const size_t i,const size_t j);
    void SwapOutParticles(const size_t i,const size_t j);
  };// end of class Blob


  typedef std::map<int,ATOOLS::Blob *> Int_Blob_Map;


  template <class Type>
  Blob_Data<Type>::Blob_Data(const Type & d) : m_data(d) {}


  template <class Type>
  std::ostream & Blob_Data<Type>::operator>>(std::ostream & s) const 
  {
    s<<m_data;
    return s;
  }

  template <class Type>
  Type Blob_Data_Base::Get() 
  {
    return ((Blob_Data<Type>*)this)->Get();
  }

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
  



