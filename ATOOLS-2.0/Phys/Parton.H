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
    //! The parton's number (integer >0)
    int     Number;     
    /*!
      The parton status (HepMC 1.01):
       0       null entry
       1       existing entry (not decayed or fragmented), i.e. final state
       2       decayed or fragmented parton
       3       documentation line
       4-10    reserved
       11-200  at disposal of each model builder
       201-... user status
    */
    int     m_status;         
    /*! 
      info flag, marks production step: 
      G - IS for hard interaction, H - FS from hard interaction, 

      I - initial particle, F - final particle, 
      i - particle from ISR, f - particle from FSR, decayed
    */ 
    char    m_info;     // just for fun!
    int     JetNumber;  //! Maybe obsolete
    /*! 
      The parton's flavour. Additional information can 
      be obtained using methods related to the class Flavour.
    */
    Flavour m_fl;         
    Flow    m_flow;
    AMATOOLS::vec4d   m_momentum; //! The parton's momentum.
    //! The parton's production vertex
    Blob * startblob;
    //! The parton's decay vertex as a fourvector in the system of the hard process
    Blob * endblob;
    //! The parton's lifetime in the system of the hard process
    double            T_Dec;      
  public:
    //! Standard constructor
    Parton();
    //! Copy constructor
    Parton(Parton * );
    //! Explicit constructor to set number, flavour and momentum of the parton.
    Parton(int,Flavour,AMATOOLS::vec4d);
    //! Standard Destructor
    ~Parton();
    //! Copy 
    void Copy(Parton * );
    /*
      Calculate proper time tau^* of particle in its rest frame via
      tau^* = hbar sqrt[q^2/((q^2-M^2)^2 + (Gamma/M q^2)^2)]
      Result is given in seconds
    */
    double ProperTime();
    /*
      Calculate lifetime of the particle according to an exponential
      decay law (like radioactive decays). Lifetime is given in cm-system of
      the hard process (where the parton list is filled), boost taken care of
      by factor gamma.
    */
    double Lifetime();
    /*
      Calculate the distance travelled by the particle during its lifetime.
      Result is given as a 3-vector in mm
    */
    AMATOOLS::vec3d Distance();

    friend std::ostream& operator<<(std::ostream &,Parton*);

    //! (almost) HepMC syntax:
    /*! 
      Some self-explanatory methods to extract or manipulate information above 
     */
    AMATOOLS::vec4d    momentum() const;
    double             E();
    int                status() const;
    char               info() const;
    AMATOOLS::vec4d    xprod();
    Blob             * prod();
    AMATOOLS::vec4d    xdec();
    Blob             * dec();
    double             time() const ;
    APHYTOOLS::Flavour flav() const ;
    APHYTOOLS::Flow    flow() const ;
    int                flow( int ) const ;
    int                Get_Numb() const;
    int                Get_JetN() const;

    void   set_momentum(const AMATOOLS::vec4d &);
    void   set_status( int status = 0 );
    void   set_info(char);
    void   set_prod(Blob *);
    void   set_dec(Blob *);
    void   set_time(const int);
    void   set_time();
    void   set_flav(APHYTOOLS::Flavour &);
    void   set_flow( const APHYTOOLS::Flow & );
    void   set_flow( int code_index, int code=0 );
    void   Set_Numb(const int);
    void   Set_JetN(const int);
  };

  // Format of a parton to transport via MPI:
  struct MPI_Parton {
    //
    int     id;     
    int     m_status;         
    char    m_info;   
    int     m_fl;      // kfcode ( -kfcode for anti)
    int     m_flow[2];
    double  m_mom[4]; //! The parton's momentum.

    // Blob information
    int     blob_owned_by;
    int     is_initial;           // 1 - is initial parton of this blob, 0 - is final parton of this
  };


  void Parton2MPI(const Parton * , MPI_Parton &);
  
  Parton * MPI2Parton(const MPI_Parton &);


  // some simple predecate operators:
  // where is_photon() is a predicate like:
  class Is_Photon {
  public:
    bool operator() ( Parton* p ) {
      if ( p && p->flav().isphoton() ) return 1;
      return 0;
    }
  };

  class Is_Gluon {
  public:
    bool operator() ( Parton* p ) {
      if ( p && p->flav().isgluon() ) return 1;
      return 0;
    }
  };

  class Is_Final_State {
  public:
    bool operator() ( const Parton* p ) {
      if ( p && (p->status() == 1) ) return 1;
      return 0;
    }
  };

}

  

/*! 
  \class Parton
  \brief Characterising individual partons occuring in event generation.

  This class contains all information needed for an Event Record based 
  on partons and methods to manipulate this information. 
  The partons defined here are close to the particles of the HepMC format
  by Matt Dobbs. ( http://cern.ch/HepMC/ )
*/

/*! 
  \var   int Number     
  \brief The parton's number (integer >0)
*/

/*! 
  \var   char info
  \brief What the parton really does and where it stems from.
  info flag, marks production step: 
  G - IS for hard interaction, H - FS from hard interaction, 
  M - internal in ME, 
  I - initial particle, F - final particle, 
  i - particle from ISR, f - particle from FSR, decayed
*/ 

/*!
  \var   int m_status;
  \brief The parton status (HepMC 1.01).
  
  The status flag can take the following values:
  0       null entry
  1       existing entry (not decayed or fragmented), i.e. final state
  2       decayed or fragmented parton
  3       documentation line
  4-10    reserved
  11-200  at disposal of each model builder
  201-... user status
*/

/*! 
  Blob * startblob, *endblob
  \brief Pointers to the parton's production and decay vertex.
*/

/*
  \var   double T_Dec
  \brief The proper decay time of the parton.
*/

/*  
  \fn    double ProperTime()
  \brief Calculate T_Dec.
  Calculate proper decay time tau^* of particle in its rest frame via
  tau^* = hbar sqrt[q^2/((q^2-M^2)^2 + (Gamma/M q^2)^2)]
  Result is given in seconds.
*/

/*
  \fn    double Lifetime()
  \brief Calculate lifetime.
  Calculate lifetime of the particle according to an exponential
  decay law (like radioactive decays). Lifetime is given in cm-system of
  the hard process (where the parton list is filled), boost taken care of
  by factor gamma.
*/

/*
  \fn    AMATOOLS::vec3d Distance()
  \brief Calculate decay length.
  Calculate the distance travelled by the particle during its lifetime.
  Result is given as a 3-vector in mm
*/

#endif






