#ifndef Analysis_Detector_Cell_H
#define Analysis_Detector_Cell_H

#include "ATOOLS/Phys/Particle.H"

namespace ANALYSIS {
  class Track {
  public:
    Track() : 
      used(false), eta(0.), phi(0.), mom(ATOOLS::Vec4D(0.,0.,0.,0.)), 
      flav(ATOOLS::Flavour(kf_none)) {}

    bool   used;
    double eta,phi;
    ATOOLS::Vec4D mom;
    ATOOLS::Flavour flav;
  };


  class Etastrip;
  class Detector_Segment;

  class Cell {
    double        m_etamin,m_etamax,m_phimin,m_phimax,m_eta,m_phi,m_costheta,m_sintheta;
    double        m_summedE;
    bool          m_used;
    Etastrip    * p_etastrip;
    Cell        * p_up, * p_down;
    ATOOLS::Vec4D m_direction;
    std::map<ATOOLS::Particle *,double> m_energydeposits;
  public:
    Cell(const double*,Etastrip *);
    ~Cell();

    void          Reset();
    void          Dimensions(double *) const;
    void          Centroid(double &,double &) const;
    bool          InPhiRange(const double phi) const { return (phi>=m_phimin && phi<=m_phimax); }

    void          AddDeposit(const double dep)       { m_summedE += dep; }
    void          AddParticle(ATOOLS::Particle *,double=-1.);
    Cell *        AddDeposit(const double,const double);
    Cell *        AddParticle(const double,ATOOLS::Particle *,double=-1.);
    Cell *        LocateCell(double);
    void          SetUsed(const bool used)   { m_used = used; }

    void          MultiplyDeposit(const double fac) { m_summedE *= fac; }

    bool          Used()               const { return m_used; }
    double        TotalDeposit()       const { return m_summedE; }
    double        EPerp()              const { return m_summedE*m_sintheta; }
    ATOOLS::Vec4D Direction()          const { return m_direction; }
    ATOOLS::Vec4D TrueMom()            const;

    Etastrip *    GetEtastrip()        const { return p_etastrip; }
    Cell *        GetUp()              const { return p_up; }
    Cell *        GetDown()            const { return p_down; }
    void          SetUp(Cell * up)           { p_up = up; }
    void          SetDown(Cell * down)       { p_down = down; }

    double        R2(const double eta,const double phi);

    std::map<ATOOLS::Particle *,double> * ParticleEntries() {
      return (&m_energydeposits);
    }

    void Print() const;
  };

  class Etastrip {
    double             m_etamin,m_etamax;
    long int           m_nphi;
    Detector_Segment * p_segment;
    Etastrip         * p_plus, * p_minus;
    Cell             * p_zero;
  public:
    Etastrip(const double,const double,const long int,Detector_Segment *);
    ~Etastrip();

    void       Reset();
    void       Dimensions(double &,double &) const;
    bool       InEtaRange(const double eta) const { return (eta>=m_etamin && eta<=m_etamax); }

    Cell     * AddDeposit(const double,const double,const double);
    Cell     * AddParticle(const double,const double,ATOOLS::Particle *,double=-1.);
    Cell     * LocateCell(const double,const double);
    Etastrip * LocateEtastrip(const double);

    Detector_Segment * GetSegment()         const { return p_segment; }
    Cell             * GetZero()            const { return p_zero; }
    Etastrip         * GetPlus()            const { return p_plus; }
    Etastrip         * GetMinus()           const { return p_minus; }
    void               SetPlus(Etastrip * es)     { p_plus = es; }
    void               SetMinus(Etastrip * es)    { p_minus = es; }

  };
}

#endif
