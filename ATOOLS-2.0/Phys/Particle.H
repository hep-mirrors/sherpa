#ifndef Particle_H
#define Particle_H

#include "Vector.H"
#include "Flavour.H"
//#include "Blob.H"
#include "Flow.H"

#include <iostream>


// compare with       HepMC_Version1.01/HepMC/*.h !!

namespace ATOOLS {
  class Blob;
   
  class Particle {
  private : 
    int       m_number;       
    int       m_status;         
    char      m_info;   
    int       m_jetnumber; 
    Flavour   m_fl;         
    Flow    * p_flow;
    Vec4D     m_momentum; 
    Blob    * p_startblob;
    Blob    * p_endblob;
    double    m_dec_time;      
    double    m_finalmass;
    static int s_totalnumber;
  public:
    Particle();
    Particle(const Particle * );
    Particle(const Particle & );
    Particle& operator=(const Particle &);
    Particle(int,Flavour=Flavour(kf::none),Vec4D=Vec4D(0.,0.,0.,0.));
    ~Particle();
    friend std::ostream& operator<<(std::ostream &,Particle*);

    void      Copy(Particle * );
    double    ProperTime();
    double    LifeTime();
    Vec3D     Distance(double = -1.);
    Vec4D     Momentum() const;
    double    E();
    double    FinalMass();
    int       Status() const;
    char      Info() const;
    Vec4D     XProd();
    Blob    * ProductionBlob();
    Vec4D     XDec();
    Blob    * DecayBlob();
    double    Time() const ;
    Flavour   Flav() const ;
    Flow    * GetFlow() const ;
    int       GetFlow( int ) const ;
    int       Number() const;
    int       JetNumber() const;
    void   SetMomentum(const Vec4D &);
    void   SetFinalMass(const double _lower,const double _upper=-1.);
    void   SetStatus( int status = 0 );
    void   SetInfo(char);
    void   SetProductionBlob(Blob *);
    void   SetDecayBlob(Blob *);
    void   SetTime(const int);
    void   SetTime();
    void   SetFlav(const Flavour &);
    void   SetFlow(Flow *);
    void   SetFlow(int code_index,int code=0);
    void   SetNumber(const int);
    void   SetJetNumber(const int);

    bool operator==(Particle part);

    static void   ResetCounter() { s_totalnumber = 0; }
    static int    Counter()      { return s_totalnumber; }
  };

  // Format of a particle to transport via MPI:
  struct MPI_Particle {
    int     id;     
    int     m_status;         
    char    m_info;   
    int     m_fl;      
    int     m_flow[2];
    double  m_mom[4];  
    int     blob_owned_by;
    int     is_initial;           // 1 - is initial particle of this blob, 0 - is final particle of this
  };

  void Particle2MPI(const Particle * , MPI_Particle &);  
  Particle * MPI2Particle(const MPI_Particle &);

  class Is_Photon {
  public:
    bool operator() ( Particle* p ) {
      if ( p && p->Flav().IsPhoton() ) return 1;
      return 0;
    }
  };

  class Is_Gluon {
  public:
    bool operator() ( Particle* p ) {
      if ( p && p->Flav().IsGluon() ) return 1;
      return 0;
    }
  };

  class Is_Final_State {
  public:
    bool operator() ( const Particle* p ) {
      if ( p && (p->Status() == 1) ) return 1;
      return 0;
    }
  };

  class Is_Charged {
  public:
    bool operator() ( const Particle* p ) {
      if ( p && (p->Flav().IntCharge() !=0) ) return 1;
      return 0;
    }
  };

  class Is_Charged_Hadron {
  public:
    bool operator() ( const Particle* p ) {
      if ( p && (p->Flav().IntCharge() !=0) 
	   &&  p->Flav().IsHadron() &&  !p->Flav().IsDiQuark()) return 1;
      return 0;
    }
  };

  /*! 
    \class Particle
    \brief Characterising individual particles occuring in event generation.
    
    This class contains all information needed for an Event Record based 
    on particles and methods to manipulate this information. 
    The particles defined here are close to the particles of the HepMC format
    by Matt Dobbs. ( http://cern.ch/HepMC/ )
  */
  
  /*! 
    \var   int Particle::m_number     
    \brief The particle's number (integer >0)
  */
  
  /*! 
    \var   char Particle::info
    \brief What the particle really does and where it stems from.
    
    info flag, marks production step: 
    \verbatim
    G - IS for hard interaction, H - FS from hard interaction, 
    M - internal in ME, 
    I - initial particle, F - final particle, 
    i - particle from ISR, f - particle from FSR, decayed
    D - particle from hadron decays,
    P - primary hadrons
    \endverbatim
  */ 
  
  /*!
    \var   int Particle::m_status;
    \brief The particle status (HepMC 1.01).
    
    The status flag can take the following values:
    \verbatim
    0       null entry
    1       existing entry (not decayed or fragmented), i.e. final state
    2       decayed or fragmented particle
    3       documentation line
    4-10    reserved
    11-200  at disposal of each model builder
    201-... user status
    \endverbatim
  */
  
  /*!
    \var   int Particle::m_info
    \brief The particle info, just for jun!
    
    info flag, marks production step: 
    G - IS for hard interaction, H - FS from hard interaction, 
    
    I - initial particle, F - final particle, 
    i - particle from ISR, f - particle from FSR, decayed
  */
  
  
  /*!
    \var   int Particle::m_jetnumber
    \brief number of jet the particle belongs to (Maybe obsolete)
  */
  
  /*!
    \var     Flavour m_fl
    \brief The particle's flavour. 
    
    The particle's flavour. 
    
    Additional information can be obtained using methods 
    related to the class Flavour.
  */
  
  /*!
    \var   Vec4D   m_momentum
    \brief The particle's momentum.
  */
  
  /*! 
    \var  Blob * Particle::m_startblob, * Particle::m_endblob
    \brief Pointers to the particle's production and decay vertex.
  */
  
  /*!
    \var   double Particle::m_dec_time
    \brief The proper decay time of the particle.
  */
  
  /*!  
    \fn    double Particle::ProperTime()
    \brief Calculate m_dec_time.
    
    Calculate proper decay time \f$\tau^*\f$ of particle in its rest frame via
    \f[
    \tau^* = hbar \sqrt{q^2/((q^2-M^2)^2 + (\Gamma/M q^2)^2)}
    \f]
    Result is given in seconds.
  */
  
  /*!
    \fn    double Particle::LifeTime()
    \brief Calculate lifetime.
    
    Calculate lifetime of the particle according to an exponential
    decay law (like radioactive decays). Lifetime is given in cm-system of
    the hard process (where the particle list is filled), boost taken care of
    by factor gamma.
  */
  
  /*!
    \fn    Vec3D Particle::Distance()
    \brief Calculate decay length.
    
    Calculate the distance travelled by the particle during its lifetime.
    Result is given as a 3-vector in mm
  */ 
}
#endif






