#ifndef Lund_Fortran_Interface_h
#define Lund_Fortran_Interface_h

#include "Blob_List.H"
#include "Particle_List.H"
#include "Particle.H"

namespace SHERPA {
  class Lund_Fortran_Interface {
    static const int maxentries = 2000;
    double  * phep, * vhep;
    int     * jmohep, * jdahep, * isthep, * idhep;

    double  * pjet, * xjet;
    int       nk, * kfjet, * mothers, * daughters;

    void Init();
    void AddPartonToString(ATOOLS::Particle *,int &);
    void FillPrimaryHadronsInBlob(ATOOLS::Blob *,ATOOLS::Blob_List *,
				  ATOOLS::Particle_List * = NULL);
    void FillSecondaryHadronsInBlob(ATOOLS::Blob *,ATOOLS::Blob_List *,int,int,
				    ATOOLS::Particle_List * = NULL);
  public:
    Lund_Fortran_Interface(double,double,double);
    ~Lund_Fortran_Interface();
    bool Hadronize(ATOOLS::Blob *,ATOOLS::Blob_List *,
		   ATOOLS::Particle_List * = NULL);
  };

  extern Lund_Fortran_Interface LFInterface;

  /*!
    \class Lund_Fortran_Interface
    \brief Contains an interface to the Fortran routines of Pythia responsible for
           hadronization and hadron decays.

    This is the interface to the Lund string model of Pythia 6.124 and its hadron decays.	   
  */
  /*!
    \fn static const int Lund_Fortran_Interface::maxentries = 2000
    The maximal number of parton entries for the interface. This is the size for
    all vectors entering the Fortran routines.
  */
  /*!
    \fn double * Lund_Fortran_Interface::phep
    Four momentum vectors for partons to be hadronized. They are of the form 
    \f$(p_x,\, p_y,\, p_z,\, E)^T\f$ and enter the fortran code through the external
    routine finterf_.
  */
  /*!
    \fn double * Lund_Fortran_Interface::vhep
    Position four vectors for partons to be hadronized. They are of the form 
    \f$x,\, y,\, z,\, t)^T\f$ and enter the fortran code through the external
    routine finterf_.
  */
  /*!
    \fn int * Lund_Fortran_Interface::jmohep
    Pointers to the mother particles for the partons to be hadronized. They are all set 
    to zero.
  */
  /*!
    \fn int * Lund_Fortran_Interface::jdahep
    Pointers to the daughter particles for the partons to be hadronized. They are all set 
    to zero.
  */
  /*!
    \fn int * Lund_Fortran_Interface::isthep
    Status flags for the partons to be hadronized, set to 1.
  */
  /*!
    \fn int * Lund_Fortran_Interface::idhep
    Particle identity flag for the partons to be hadronized, set to the corresponding PDG code
    through the method ATOOLS::Flavour::HepEvt().
  */
  /*!
    \fn double * Lund_Fortran_Interface::pjet
    Four momenta of the particles that are read out through the subroutine fhawface_.
    Again, the form is \f$(p_x,\, p_y,\, p_z,\, E)^T\f$.
  */
  /*!
    \fn double * Lund_Fortran_Interface::xjet
    Production positions of the particles that are read out through the subroutine fhawface_.
    Again, the form is \f$(x,\, y,\, z,\, t)^T\f$.
  */
  /*!
    \fn int Lund_Fortran_Interface::nk
    The number of particles in the HepEvt common block that is read out.
  */
  /*!
    \fn int * Lund_Fortran_Interface::kfjet
    The flavour label of the particles.
  */
  /*!
    \fn int * Lund_Fortran_Interface::mothers
    Pointers to the position(s) of the mother particle(s) in the list that is read out. 
    Usually, this is either the string entry signifying a primary hadron, or the hadron that 
    produced the actual one in its decay. 
  */
  /*!
    \fn int * Lund_Fortran_Interface::daughters
    Pointers to the position of the first and the last daughter in the list.
  */
  /*!
    \var void Lund_Fortran_Interface::Init()
    Instantiates some of the variables.
  */
  /*!
    \var void Lund_Fortran_Interface::AddPartonToString(ATOOLS::Parton *,int &)
    Translates an additional single parton into the HepEvt standard and sets its
    phep, vhep, jmohep, jdahep, idhep, isthep variables.
  */
  /*!
    \var void Lund_Fortran_Interface::FillPrimaryHadronsInBlob(ATOOLS::Blob *,
                                                               ATOOLS::Blob_List *,
				                               ATOOLS::Parton_List * = NULL)
    Reads out the result given by Pythia, i.e. selects only hadrons to deal with and:
    - defines the flavour;
    - reads out four momentum and production position;
    - instamntiates a new corresponding parton and generates and sets its number;
    - sets the link to the production blob and its position and adds the new parton to the
      list of outgonig partons of this blob;
    - if daughter flags are different from zero (i.e. if the particle decays),
      a decay blob is instantiated and suitable links are established, i.e. the link to the
      decay blob of the parton and the parton itself is added to the list of incoming partons 
      of its decay blob;
    - the full decay chain is spanned by a call to FillSecondaryHadronsInBlob with suitable
      arguments.
    Finally the fragmentation blob is set inactive.  
  */
  /*!
    \var void Lund_Fortran_Interface::FillSecondaryHadronsInBlob(ATOOLS::Blob *,
                                                                 ATOOLS::Blob_List *,int,int,
				                                 ATOOLS::Parton_List * = NULL)
    This method is to fill in daughters of decaying hadrons in the decay blob and continue in the
    decay chain by calling itself with modified arguments. The procedure of filling in
    partons etc. is very much in parallel to what happens in FillPrimaryHadronsInBlob.
  */
  /*!
    \var Lund_Fortran_Interface::Lund_Fortran_Interface(double,double,double)
    Initializes the interface with three arguments for the Lund parameters \f$a\f$, \f$b\f$,
    and \f$\sigma\f$.
  */
  /*!
    \var Lund_Fortran_Interface::~Lund_Fortran_Interface()
    Deletes all pointers to the variables instantiated in it().
  */
  /*!
    \var bool Lund_Fortran_Interface::Hadronize(ATOOLS::Blob *,ATOOLS::Blob_List *,
		                                ATOOLS::Parton_List * = NULL)
    For a specific blob the Lund hadronization is performed. First, the tuype of the blob is
    set to "Fragmentation (Lund : Pythia 6.163)", then a photon pseudoparticle is created
    with appropriate production position, momentum, etc.. By iterating over all incoming
    partons of the blob, and translating them into the HepoEvt standard through AddPartonToString,
    the parameters of the common block are prepared for hadronization. This is performed
    by a call of the subroutine finterf_ to be foound in the Fortran interface 
    file finterf.f.
  */
}

#endif
