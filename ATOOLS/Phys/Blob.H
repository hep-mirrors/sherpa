#ifndef ATOOLS_Phys_Blob_H
#define ATOOLS_Phys_Blob_H

#include <string>
#include <vector>
#include <map>

#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Org/STL_Tools.H"

namespace ATOOLS {

  class Blob_Data_Base {
  private:
    static long int s_number;
  public:
    template <class Type> Type &Get();
    template <class Type> void Set(const Type &data);
    Blob_Data_Base();
    Blob_Data_Base(const Blob_Data_Base &base);
    virtual ~Blob_Data_Base();
    virtual std::ostream & operator>>(std::ostream &) const =0;
  };

  template <class Type> class Blob_Data : public Blob_Data_Base {
    Type m_data;
  public:
    Blob_Data(const Type & d);
    Type &Get() { return m_data; }
    void Set(const Type & d) { m_data=d; }
    std::ostream & operator>>(std::ostream &) const;
    ~Blob_Data();
  };

  std::ostream& operator<<(std::ostream&,const Blob_Data_Base &);

  typedef std::map<std::string,Blob_Data_Base *> String_BlobDataBase_Map;

  struct blob_status {
    enum code {
      inactive            = 0,
      needs_signal        = 1,
      needs_showers       = 2,
      needs_harddecays    = 4,
      needs_beams         = 8,
      needs_softUE        = 16,
      needs_hadronization = 32,
      needs_hadrondecays  = 64,
      needs_extraQED      = 128,
      needs_smearing      = 16384,
      internal_flag       = 32768,
      fully_active        = 65535
    };
  };

  inline const blob_status::code operator|(const blob_status::code c1,const blob_status::code c2)
  { return (blob_status::code)((int)c1|(int)c2); }
  inline const blob_status::code operator&(const blob_status::code c1,const blob_status::code c2)
  { return (blob_status::code)((int)c1&(int)c2); }

  struct btp {
    enum code {
      Signal_Process      =     1,
      Hard_Decay          =     2,
      Hard_Collision      =     4,
      Soft_Collision      =     8,
      Shower              =    16,
      QED_Radiation       =    32,
      Beam                =   256,
      Bunch               =   512,
      Fragmentation       =  1024,
      Cluster_Formation   =  2048,
      Cluster_Decay       =  4096,
      Hadron_Decay        =  8192,
      Hadron_Mixing       = 16384,
      Hadron_To_Parton    = 32768,
      Unspecified         = 65535
    };
  };// end of struct btp

  inline const btp::code operator|(const btp::code c1,const btp::code c2)
  { return (btp::code)((int)c1|(int)c2); }
  inline const btp::code operator&(const btp::code c1,const btp::code c2)
  { return (btp::code)((int)c1&(int)c2); }

  std::ostream& operator<<(std::ostream&,const btp::code);

  class Blob {
    friend std::ostream& operator<<( std::ostream&, const Blob &);
  private:
    static long unsigned int s_currentnumber;
    Vec4D                   m_position;
    int                     m_id;
    double                  m_weight;
    blob_status::code       m_status;
    int                     m_beam;
    bool                    m_hasboost;
    btp::code               m_type;
    std::string             m_typespec;
    std::vector<double>     m_weightcontainer;
    String_BlobDataBase_Map m_datacontainer;
    Particle_Vector         m_inparticles, m_outparticles;
    Vec4D                   m_cms_vec;
    Poincare                m_cms_boost;
    static int              s_totalnumber;
    bool IsConnectedTo(const btp::code &type,
		       std::set<const Blob*> &checked) const;
  public:
    Blob(const Vec4D _pos = Vec4D(0.,0.,0.,0.), const int _id=-1);
    Blob(const Blob *);
    ~Blob();
    void     AddToInParticles(Particle *);
    void     AddToOutParticles(Particle *);
    Particle_Vector GetOutParticles() { return m_outparticles; }
    Particle_Vector GetInParticles()  { return m_inparticles;  }
    Particle * OutParticle(int);
    Particle * InParticle(int);
    const Particle *ConstOutParticle(const size_t i) const;
    const Particle *ConstInParticle(const size_t i) const;
    Particle * RemoveInParticle(int,bool = true);
    Particle * RemoveInParticle(Particle *,bool = true);
    Particle * RemoveOutParticle(int,bool = true);
    Particle * RemoveOutParticle(Particle *,bool = true);
    void     RemoveInParticles(const int = 1);
    void     RemoveOutParticles(const int = 1);
    void     DeleteInParticles(const int = -1);
    void     DeleteOutParticles(const int = -1);
    void     DeleteInParticle(Particle *);
    void     DeleteOutParticle(Particle *);
    void     RemoveOwnedParticles(const bool = true);
    void     DeleteOwnedParticles();

    void     Boost(Poincare boost);
    void     BoostInCMS();
    void     BoostInLab();

    void        SetVecs();
    void        SetPosition(Vec4D pos)            { m_position = pos; }
    void        SetCMS(Vec4D _cms)                { m_cms_vec  = _cms; }
    void        SetCMS();
    void        SetId(const int _id=0);
    static void ResetCounter()                    { s_totalnumber = 0; }
    static int  Counter()                         { return s_totalnumber; }
    inline static void Reset(const int number=0)  { s_currentnumber = number; }
    inline bool Has(blob_status::code status)     { return (int(status&m_status)>0); }
    bool IsConnectedTo(const btp::code &type) const;
    Blob *      UpstreamBlob() const;
    Blob *      DownstreamBlob() const;
    void        SetStatus(blob_status::code status=blob_status::inactive) {
      m_status = status;
    }
    void        UnsetStatus(blob_status::code status=blob_status::inactive) {
      m_status = blob_status::code(int(m_status) & int(~status));
    }
    void        AddStatus(blob_status::code status=blob_status::inactive) {
      m_status = m_status | status;
    }
    void        SetType(btp::code _type)          { m_type     = _type; }
    void        SetTypeSpec(std::string _type)    { m_typespec = _type; }
    void        SetBeam(int _beam)                { m_beam     = _beam; }
    void        SetWeight(double _weight)         { 
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
    const String_BlobDataBase_Map & GetData() const {return m_datacontainer;}

    void     AddData(const std::string name, Blob_Data_Base * data); 
    void     ClearAllData();

    int      Id()                        const { return m_id; }
    int      Status()                    const { return m_status; }
    int      Beam()                      const { return m_beam; }
    int      NInP()                      const { return m_inparticles.size(); }
    int      NOutP()                     const { return m_outparticles.size(); }
 
    double              const Weight()   const { return m_weight; }
    std::vector<double> WeightContainer()      { return m_weightcontainer; } 
    Vec4D        CheckMomentumConservation() const;
    bool  MomentumConserved();
    const Vec4D& Position()              const { return m_position; }
    const Vec4D& CMS()                   const { return m_cms_vec; }
    btp::code   const Type()             const { return m_type; }
    std::string const TypeSpec()         const { return m_typespec; }
    void SwapInParticles(const size_t i,const size_t j);
    void SwapOutParticles(const size_t i,const size_t j);
  };// end of class Blob


  typedef std::map<int,ATOOLS::Blob *> Int_Blob_Map;
  typedef std::map<ATOOLS::Particle *,ATOOLS::Blob *> Particle_Blob_Map;



  template <class Type>
  Blob_Data<Type>::Blob_Data(const Type & d) : m_data(d) {}


  template <class Type>
  std::ostream & Blob_Data<Type>::operator>>(std::ostream & s) const 
  {
    s<<m_data;
    return s;
  }

  template <class Type>
  Type &Blob_Data_Base::Get() 
  {
    return ((Blob_Data<Type>*)this)->Get();
  }

  template <class Type>
  void Blob_Data_Base::Set(const Type &data) 
  {
    return ((Blob_Data<Type>*)this)->Set(data);
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
  
  /*!
    \fn Blob* Blob::UpstreamBlob()
    \brief Return unique upstream blob, if it exists
    
    Will return the production blob of all inparticles, if it is equal among them.
    Otherwise it will return NULL.
  */
  
  /*!
    \fn Blob* Blob::DownstreamBlob()
    \brief Return unique downstream blob, if it exists
    
    Will return the decay blob of all outparticles, if it is equal among them.
    Otherwise it will return NULL.
  */
}

#endif
  



