#ifndef Reconstructed_Object_h
#define Reconstructed_Object_h

#include "Cell.H"
#include "Flavour.H"
#include "Particle.H"

namespace ANALYSIS {
  class Primitive_Analysis;

  class Reconstructed_Object {
  private:
    bool            m_includetracks;
    ATOOLS::Flavour m_flav;
    double          m_E,m_eta,m_phi;
    ATOOLS::Vec4D   m_mom;

    std::vector<Cell *>  m_cells;
    std::vector<Track *> m_tracks;
  public:
    Reconstructed_Object(ATOOLS::Flavour flav,const double E,const double eta,const double phi);
    Reconstructed_Object(ATOOLS::Flavour flav,std::vector<Cell *> &);
    Reconstructed_Object(Track * track);
    ~Reconstructed_Object();

    double E()              const { return m_E; }
    double Eta()            const { return m_eta; }
    double Phi()            const { return m_phi; }
    ATOOLS::Flavour Flav()  const { return m_flav; }
    ATOOLS::Vec4D Mom()     const { return m_mom; }
    ATOOLS::Vec4D TrueMom() const;
    
    void   SetFlav(const ATOOLS::Flavour flav) { m_flav = flav; }

    void   CorrectE(const double inv) { m_E/=inv; }
    void   Update();
    void   SetUsed(const bool);

    void   SetIncludeTracks(const bool it) { m_includetracks = it; }

    void SetCells(const std::vector<Cell *> & cells);
    void AddCell(Cell * cell);
    std::vector<Cell *> & GetCells()     { return m_cells; }
    bool IsIncluded(const Cell *) const;

    void SetTracks(const std::vector<Track *> & tracks);
    void AddTrack(Track * track);
    std::vector<Track *> & GetTracks()   { return m_tracks; }
    bool IsIncluded(const Track *) const;

    int NumberOfParticles() const;

    ATOOLS::Particle * CreateParticle();
  };

  typedef std::list<Reconstructed_Object *>           ObjectList;
  typedef std::list<Reconstructed_Object *>::iterator ObjectListIterator;
}

#endif
