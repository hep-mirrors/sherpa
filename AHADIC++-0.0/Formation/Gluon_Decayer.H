#ifndef Gluon_Decayer_H
#define Gluon_Decayer_H

#include "Cluster.H"
#include "Blob.H"

namespace AHADIC {
  struct DecaySpecs { double weight,zmin,zmax; };
  typedef std::map<ATOOLS::Flavour,DecaySpecs *> FlavDecayMap;
  typedef FlavDecayMap::iterator                 FDIter;

  class Gluon_Decayer {
  private:
    FlavDecayMap        m_options;
    double              m_z,m_phi,m_mass;
    ATOOLS::Flavour     m_flav;
    ATOOLS::Vec4D       m_p1vec, m_p2vec;

    void                Shift(Part_List * pl);
    bool                DecayIt(ATOOLS::Particle * part,
				ATOOLS::Particle *& part1,ATOOLS::Particle *& part2);
    void                SelectDecay(const ATOOLS::Vec4D &);
    double              ZWeight(const double, const double);
    void                SelectZ(const DecaySpecs*);
    void                SelectPhi();
    void                BuildKinematics(const ATOOLS::Vec4D &);
  public:
    Gluon_Decayer();
    ~Gluon_Decayer();
    bool DecayList(Part_List * plin=NULL);
  };
}
#endif
