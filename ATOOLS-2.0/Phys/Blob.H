#ifndef Blob_h
#define Blob_h

#include <string>
#include <deque>
#include <iterator> 
#include <iostream> 

#include "Vector.H"
#include "Poincare.H"
//#include "Parton.H"

namespace APHYTOOLS {
  class Parton; 
  class Blob {
    friend std::ostream& operator<<( std::ostream&, const Blob &);
    friend std::ostream& operator<<( std::ostream&, const Blob *);
  private:
    int                   id;
    int                   status;
    int                   beam;
    AMATOOLS::vec4d       cms;
    AMATOOLS::Poincare    cms_boost;
    AMATOOLS::vec4d       position;
    int                   Nin, Nout;
    std::deque<Parton *>  InPartons;
    std::deque<Parton *>  OutPartons;
    std::string           type;
  public:
    Blob(const AMATOOLS::vec4d _pos = AMATOOLS::vec4d(0.,0.,0.,0.), const int _id=0) : 
      position(_pos), id(_id) { Nin = 0; Nout=0; beam = -1; status = 0; };
    ~Blob();
    Parton * OutParton(int);
    Parton * InParton(int);
    void     AddToInPartons(Parton *);
    void     AddToOutPartons(Parton *);
    void     DeleteOwnedPartons();
    void     BoostInCMS();
    void     BoostInLab();

    void     SetVecs();
    void     SetPosition(AMATOOLS::vec4d pos) { position = pos; }
    void     SetCMS(AMATOOLS::vec4d _cms)     { cms      = _cms; }
    void     SetId(int _id)                   { id       = _id; }
    void     SetStatus(int _status)           { status   = _status; }
    void     SetType(std::string _type)       { type     = _type; }
    void     SetBeam(int _beam)               { beam     = _beam; }

    AMATOOLS::vec4d    CheckMomentumConservation();
    AMATOOLS::vec4d    const Position() const { return position; }
    AMATOOLS::vec4d    const CMS()      const { return cms; }
    int Id()                            const { return id; }
    int Status()                        const { return status; }
    int Beam()                          const { return beam; }
    int NInP()                          const { return Nin; }
    int NOutP()                         const { return Nout; }
    std::string const Type()            const { return type; }
    /*
  public:
    class Parton_Iterator {
    public:
      Parton_Iterator();
      Parton_Iterator(const Parton_Iterator& );
    };
    */
  };
}


/*!
  \var   char type;
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

#endif
  




