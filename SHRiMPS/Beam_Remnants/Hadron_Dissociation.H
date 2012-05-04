#ifndef SHRIMPS_Beam_Remnants_Hadron_Dissociation_H
#define SHRIMPS_Beam_Remnants_Hadron_Dissociation_H

#include "SHRiMPS/Beam_Remnants/Continued_PDF.H"
#include "SHRiMPS/Eikonals/Form_Factors.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Math/Histogram_2D.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include <vector>

namespace SHRIMPS {
  class Hadron_Dissociation {
  private:
    bool                             m_elastic;
    std::vector<ATOOLS::Particle * > m_particles;
    std::vector<double>              m_xs;
    std::vector<ATOOLS::Vec4D>       m_qtvecs;

    Continued_PDF * p_pdf;
    ATOOLS::Flavour m_bunch, m_quark, m_remnant;
    ATOOLS::Vec4D   m_bunchmom;
    ATOOLS::Blob *  p_blob;
    double          m_ycut;

    bool                                           m_analyse;
    std::map<std::string, ATOOLS::Histogram * >    m_histomap;
    std::map<std::string, ATOOLS::Histogram_2D * > m_histomap2D;

    void FixFlavourConstituents();
    void FillParticleList(const int & N);
    void DefineTransverseMomenta(Form_Factor * ff);
  public:
    Hadron_Dissociation(Continued_PDF *const pdf);
    ~Hadron_Dissociation();
    bool DefineDissociation(const int & Nladders,const double & xcut,
			    const double & eta,Form_Factor * ff);
    void Reshuffle(const size_t & N);
    void FillBeamBlob();
    void AddParticlesToBlob(ATOOLS::Blob * blob,int beam);
    bool MustReplaceColour(const size_t & pos,
			   const size_t & c1,const size_t & c2);
    inline void Reset(const ATOOLS::Vec4D & mom) {
      m_quark    = m_remnant = ATOOLS::Flavour(kf_none);
      m_bunchmom = mom;
      m_particles.clear();
      m_xs.clear();
      m_qtvecs.clear();
      m_elastic  = true;
    }
    inline const bool & Elastic() const { return m_elastic; } 
    inline void SetBeamBlob(ATOOLS::Blob * blob) { p_blob = blob; }
    inline ATOOLS::Blob * GetBeamBlob() const { return p_blob; }
    inline const std::vector<ATOOLS::Particle * > & GetParticles() const {
      return m_particles;
    }
    inline ATOOLS::Particle * GetParticle(const size_t & n) const {
      if (n<m_particles.size()) return m_particles[n];
      return NULL;
    }
    inline void GetXs(const size_t & n,double & xp,double & xm,double & xt2,
		      const double & shat) {
      if (n<m_particles.size()) {
	if (m_bunchmom[3]>0.) { xp = m_xs[n]; xm = 0.; }
	                 else { xm = m_xs[n]; xp = 0.; }
	xt2 = ATOOLS::dabs(m_qtvecs[n].Abs2())/shat;
	return;
      }
      msg_Error()<<"Error in "<<METHOD<<"("<<n<<"): out of bounds.\n";
    }
    inline const ATOOLS::Vec4D Kperp(const size_t & n) const {
      if (n<m_particles.size()) return m_qtvecs[n];
      msg_Error()<<"Error in "<<METHOD<<"("<<n<<"): out of bounds.\n";
      return ATOOLS::Vec4D(0.,0.,0.,0.);
    }
    inline size_t Size() const {
      return m_particles.size();
    }
    inline void DeleteParticles() {
      while (!m_particles.empty()) {
	delete m_particles.back();
	m_particles.pop_back();
      }
      m_elastic = true;
    }

    void PrintParticles() const;
  };
}

#endif
