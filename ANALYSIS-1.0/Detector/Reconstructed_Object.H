#ifndef Reconstructed_Object_h
#define Reconstructed_Object_h

#include "Cell.H"
#include "Flavour.H"
#include "Particle.H"

namespace ANALYSIS {
  class Primitive_Analysis;

  class Reconstructed_Object {
  private:
    ATOOLS::Flavour m_flav;
    double m_E,m_eta,m_phi;
    ATOOLS::Vec4D m_mom;

    std::vector<Cell *>  * p_cells;
    std::vector<Track *> * p_tracks;
  public:
    Reconstructed_Object(ATOOLS::Flavour flav,const double E,const double eta,const double phi);
    Reconstructed_Object(ATOOLS::Flavour flav,std::vector<Cell *> *);
    ~Reconstructed_Object();

    double E()             const { return m_E; }
    double Eta()           const { return m_eta; }
    double Phi()           const { return m_phi; }
    ATOOLS::Flavour Flav() const { return m_flav; }
    ATOOLS::Vec4D Mom()    const { return m_mom; }
    
    void   CorrectE(const double inv) { m_E/=inv; }

    void SetCells(std::vector<Cell *> * cells) { p_cells=cells; }
    void AddCell(Cell * cell)                  { p_cells->push_back(cell); }
    std::vector<Cell *> * GetCells()     const { return p_cells; }
    bool IsIncluded(const Cell *) const;

    void SetTracks(std::vector<Track *> * tracks) { p_tracks=tracks; }
    void AddTrack(Track * track)                  { p_tracks->push_back(track); }
    std::vector<Track *> * GetTracks()     const  { return p_tracks; }
    bool IsIncluded(const Track *) const;

    ATOOLS::Particle * CreateParticle();
  };

  typedef std::list<Reconstructed_Object *>           ObjectList;
  typedef std::list<Reconstructed_Object *>::iterator ObjectListIterator;
}

#endif
