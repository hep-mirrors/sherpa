#ifndef Parton_H
#define Parton_H

#include "Vector.H"
#include "Flavour.H"
//#include "Blob.H"
#include "Flow.H"

#include <iostream>


// compare with       HepMC_Version1.01/HepMC/*.h !!

namespace APHYTOOLS {
  class Blob;
   
  class Parton {
  private : 
    int               m_number;       
    int               m_status;         
    char              m_info;   
    int               m_jetnumber; 
    Flavour           m_fl;         
    Flow            * p_flow;
    AMATOOLS::Vec4D   m_momentum; 
    Blob            * p_startblob;
    Blob            * p_endblob;
    double            m_dec_time;      
  public:
    Parton();
    Parton(Parton * );
    Parton(int,Flavour=Flavour(kf::none),AMATOOLS::Vec4D=AMATOOLS::Vec4D(0.,0.,0.,0.));
    ~Parton();
    friend std::ostream& operator<<(std::ostream &,Parton*);

    void Copy(Parton * );
    double ProperTime();
    double LifeTime();
    AMATOOLS::Vec3D Distance();
    AMATOOLS::Vec4D    Momentum() const;
    double             E();
    int                Status() const;
    char               Info() const;
    AMATOOLS::Vec4D    XProd();
    Blob             * ProductionBlob();
    AMATOOLS::Vec4D    XDec();
    Blob             * DecayBlob();
    double             Time() const ;
    APHYTOOLS::Flavour Flav() const ;
    APHYTOOLS::Flow  * GetFlow() const ;
    int                GetFlow( int ) const ;
    int                Number() const;
    int                JetNumber() const;
    void   SetMomentum(const AMATOOLS::Vec4D &);
    void   SetStatus( int status = 0 );
    void   SetInfo(char);
    void   SetProductionBlob(Blob *);
    void   SetDecayBlob(Blob *);
    void   SetTime(const int);
    void   SetTime();
    void   SetFlav(APHYTOOLS::Flavour &);
    void   SetFlow(APHYTOOLS::Flow *);
    void   SetFlow(int code_index,int code=0);
    void   SetNumber(const int);
    void   SetJetNumber(const int);
  };

  // Format of a parton to transport via MPI:
  struct MPI_Parton {
    int     id;     
    int     m_status;         
    char    m_info;   
    int     m_fl;      
    int     m_flow[2];
    double  m_mom[4];  
    int     blob_owned_by;
    int     is_initial;           // 1 - is initial parton of this blob, 0 - is final parton of this
  };

  void Parton2MPI(const Parton * , MPI_Parton &);  
  Parton * MPI2Parton(const MPI_Parton &);

  class Is_Photon {
  public:
    bool operator() ( Parton* p ) {
      if ( p && p->Flav().IsPhoton() ) return 1;
      return 0;
    }
  };

  class Is_Gluon {
  public:
    bool operator() ( Parton* p ) {
      if ( p && p->Flav().IsGluon() ) return 1;
      return 0;
    }
  };

  class Is_Final_State {
  public:
    bool operator() ( const Parton* p ) {
      if ( p && (p->Status() == 1) ) return 1;
      return 0;
    }
  };

  /*! 
    \class Parton
    \brief Characterising individual partons occuring in event generation.
    
    This class contains all information needed for an Event Record based 
    on partons and methods to manipulate this information. 
    The partons defined here are close to the particles of the HepMC format
    by Matt Dobbs. ( http://cern.ch/HepMC/ )
  */
  
  /*! 
    \var   int Parton::m_number     
    \brief The parton's number (integer >0)
  */
  
  /*! 
    \var   char Parton::info
    \brief What the parton really does and where it stems from.
    
    info flag, marks production step: 
    \verbatim
    G - IS for hard interaction, H - FS from hard interaction, 
    M - internal in ME, 
    I - initial particle, F - final particle, 
    i - particle from ISR, f - particle from FSR, decayed
    \endverbatim
  */ 
  
  /*!
    \var   int Parton::m_status;
    \brief The parton status (HepMC 1.01).
    
    The status flag can take the following values:
    \verbatim
    0       null entry
    1       existing entry (not decayed or fragmented), i.e. final state
    2       decayed or fragmented parton
    3       documentation line
    4-10    reserved
    11-200  at disposal of each model builder
    201-... user status
    \endverbatim
  */
  
  /*!
    \var   int Parton::m_info
    \brief The parton info, just for jun!
    
    info flag, marks production step: 
    G - IS for hard interaction, H - FS from hard interaction, 
    
    I - initial particle, F - final particle, 
    i - particle from ISR, f - particle from FSR, decayed
  */
  
  
  /*!
    \var   int Parton::m_jetnumber
    \brief number of jet the parton belongs to (Maybe obsolete)
  */
  
  /*!
    \var     Flavour m_fl
    \brief The parton's flavour. 
    
    The parton's flavour. 
    
    Additional information can be obtained using methods 
    related to the class Flavour.
  */
  
  /*!
    \var   AMATOOLS::Vec4D   m_momentum
    \brief The parton's momentum.
  */
  
  /*! 
    \var  Blob * Parton::m_startblob, * Parton::m_endblob
    \brief Pointers to the parton's production and decay vertex.
  */
  
  /*!
    \var   double Parton::m_dec_time
    \brief The proper decay time of the parton.
  */
  
  /*!  
    \fn    double Parton::ProperTime()
    \brief Calculate m_dec_time.
    
    Calculate proper decay time \f$\tau^*\f$ of particle in its rest frame via
    \f[
    \tau^* = hbar \sqrt{q^2/((q^2-M^2)^2 + (\Gamma/M q^2)^2)}
    \f]
    Result is given in seconds.
  */
  
  /*!
    \fn    double Parton::LifeTime()
    \brief Calculate lifetime.
    
    Calculate lifetime of the particle according to an exponential
    decay law (like radioactive decays). Lifetime is given in cm-system of
    the hard process (where the parton list is filled), boost taken care of
    by factor gamma.
  */
  
  /*!
    \fn    AMATOOLS::Vec3D Parton::Distance()
    \brief Calculate decay length.
    
    Calculate the distance travelled by the particle during its lifetime.
    Result is given as a 3-vector in mm
  */ 
}
#endif






