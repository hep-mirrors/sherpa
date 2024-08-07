#ifndef ATOOLS_Phys_Particle_H
#define ATOOLS_Phys_Particle_H

#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/Flow.H"

#include <iostream>
#include <list>

namespace ATOOLS {
  class Blob;

  struct part_status {
    enum code {
      undefined     = 0,
      active        = 1,
      decayed       = 2,
      documentation = 3,
      fragmented    = 4,
      internal      = 5
    };
  };

  class Particle {
  private :
    static long unsigned int s_currentnumber;
    int               m_number, m_beam;
    size_t            m_meid;
    part_status::code m_status;
    char              m_info;
    Flavour           m_fl;
    Vec4D             m_momentum, m_position;
    Flow              m_flow;
    Blob            * p_startblob, * p_endblob;
    Particle        * p_originalpart;
    double            m_dec_time, m_finalmass;
    bool              m_ownpos;
    bool              m_fromdec;
    static int        s_totalnumber;
  public:
    Particle();
    Particle(const Particle & );
    Particle& operator=(const Particle &);
    Particle(int, Flavour=Flavour(kf_none), Vec4D=Vec4D(0.,0.,0.,0.),
	     char='a');
    ~Particle();
    friend std::ostream& operator<<(std::ostream &,const Particle &);

    void              Copy(Particle * );
    double            ProperTime();
    double            LifeTime();
    Vec3D             Distance(double = -1.);
    const Vec4D&      Momentum() const;
    double            E() const;
    double            FinalMass() const;
    part_status::code Status() const;
    char              Info() const;
    Vec4D             XProd() const;
    Blob            * ProductionBlob() const;
    Vec4D             XDec() const;
    Blob            * DecayBlob() const;
    double            Time() const ;
    Flavour           Flav() const ;
    const Flavour   & RefFlav() const;
    unsigned int      GetFlow( unsigned int ) const ;
    int               Number() const;
    int               Beam() const;
    Particle        * OriginalPart() const;
    void   SetMomentum(const Vec4D &);
    void   SetPosition(const Vec4D &);
    void   SetFinalMass(const double=-1,const double=-1);
    void   SetStatus(part_status::code=part_status::undefined);
    void   SetInfo(char);
    void   SetProductionBlob(Blob *);
    void   SetDecayBlob(Blob *);
    void   SetTime(const double);
    void   SetTime();
    void   SetFlav(const Flavour &);
    void   SetFlow(int code_index,int code=0);
    void   SetNumber(const int number=0);
    void   SetBeam(const int);
    void   SetOriginalPart(Particle*);
    bool operator==(Particle part) const;
    inline void   SetFromDec(const bool &val){m_fromdec=val;}
    inline bool   FromDec() const {return m_fromdec;}
    static void   ResetCounter(const int number=0)  { s_totalnumber = number; }
    static int    Counter()                         { return s_totalnumber; }
    inline static void Reset(const int number=0)    { s_currentnumber = number; }
    inline static long unsigned int CurrentNumber() { return s_currentnumber; }
    inline size_t MEId() const { return m_meid; }
    inline void   SetMEId(const size_t &id) { m_meid=id; }
  };

  typedef std::vector<Particle *>   Particle_Vector;
  typedef std::map<Particle *,int>  Particle_Int_Map;
  typedef std::map<int,Particle *>  Int_Particle_Map;

  typedef std::list<Particle*>      Part_List;
  typedef Part_List::iterator       Part_Iterator;
  typedef Part_List::const_iterator Part_Const_Iterator;
  typedef Particle_Vector::iterator PVIt;

  std::ostream & operator<<(std::ostream &,const Particle &);
  std::ostream & operator<<(std::ostream & s, const Part_List & pl);

  typedef std::map<int,std::pair<Particle *,bool> > Translation_Map;

  class PMMSort {
  public:
    bool operator()(const Particle * a,const Particle * b) const {
      return (a->Number()>b->Number());
    }

  };// end of class String_Sort

  typedef std::map<ATOOLS::Particle *,ATOOLS::Vec4D,PMMSort> ParticleMomMap;
  typedef std::set<Particle *,PMMSort>       ParticleSet;



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
    \var   bool Particle::m_fromdec
    \brief True, if this particle is emerged from a hard decay. This property is set in the Hard Decay Handler and propagated through the CSS Shower.
  */

  /*!
    \var   char Particle::m_info
    \brief What the particle really does and where it stems from.

    info flag, marks production step:
    \verbatim
    G - IS for hard interaction, H - FS from hard interaction,
    M - internal in ME,
    I - initial particle, F - final particle,
    i - particle from ISR, f - particle from FSR, decayed
    D - particle from hadron decays,
    P - primary hadrons
    S - soft photons, s - soft photons splitting products
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
    \var     Flavour Particle::m_fl
    \brief The particle's flavour.

    The particle's flavour.

    Additional information can be obtained using methods
    related to the class Flavour.
  */

  /*!
    \var   Vec4D   Particle::m_momentum
    \brief The particle's momentum.
  */

  /*!
    \var  Blob * Particle::p_startblob
    \brief Pointers to the particle's production vertex.
    \see Particle::m_endblob
  */

  /*!
    \var  Blob * Particle::p_endblob
    \brief Pointers to the particle's decay vertex.
    \see Particle::m_startblob
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
