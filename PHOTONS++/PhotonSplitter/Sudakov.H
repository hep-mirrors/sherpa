#ifndef PHOTONS_Photon_Splitter_Sudakov_H
#define PHOTONS_Photon_Splitter_Sudakov_H

#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Blob.H"
#include "PHOTONS++/PhotonSplitter/Splitting_Functions.H"
#include "PHOTONS++/PhotonSplitter/Kinematics_Base.H"
#include <cstddef>
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Math/Histogram_2D.H"

namespace PHOTONS {

  class Sudakov {
  private:
    SF_Vector           m_splitters;
    YFS_Particle_Vector m_splittingphotons;
    YFS_Particle_Vector m_spectators;

    std::vector<size_t> m_splitterIds;
    int                 m_NInP;
    bool                m_decay;

    int m_mode;
    int m_virtualityOrdering;
    int m_spectatorScheme;
    bool m_debug_initProbabilistic;
    double m_masscutoff;

    double m_t0;
    double m_t;

    // flavour-dependent enhancements
    std::map<kf_code,double> m_enhancefac;

    // for storing splittings
    YFS_Particle_List m_remainingphotons;
    YFS_Particle_Vector m_addedparticles;

    bool m_addedanything;

    static std::string          s_histo_base_name;
    static ATOOLS::Histogram    s_histo_dipole;
    static ATOOLS::Histogram_2D s_histo_tdR;

    void RegisterDefaults();

  public:
    Sudakov();
    Sudakov(int mode);
    ~Sudakov();

    void SetNInParticles(int n);
    void AddSplitter(ATOOLS::Particle *softphoton, const size_t&);
    void AddChargedParticle(ATOOLS::Particle* p, const size_t& id);
    inline YFS_Particle_Vector GetSpectators() { return m_spectators; }
    bool ClearAll();

    void SetCutoff();
    YFS_Particle* DefineInitialConditions(double &t,ATOOLS::Vec4D pphoton);

    bool Run();
    bool Generate();

    std::string Info() const { return std::string(""); }
    
    inline YFS_Particle_List GetRemainingSoftPhotons() { return m_remainingphotons; }
    inline YFS_Particle_Vector GetAddedParticles() { return m_addedparticles; }
    inline bool AddedAnything() { return m_addedanything; }
  };

  /*!
    \file Sudakov.H
    \brief Contains the class PHOTONS::Sudakov
  */

  /*! 
    \class Sudakov
    \brief Implements the veto algorithm for a 1-step parton shower to split soft photons.
  */
  /////////////////////////////////////////////////////////////////////////////////////////
  // Description of member variables for Photon_Splitter
  /////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var int m_mode
    \brief encodes which splittings to include. Settings are additive. 
    1 = electrons, 2 = muons, 4 = taus, 8 = hadrons up to mass cutoff.
  */
  /*!
    \var int m_virtualityOrdering
    \brief stores the ordering scheme to be used. 0 is transverse momentum ordering, 
    1 is virtuality ordering and 2 is mixed (default).
  */
  /*!
    \var int m_spectatorScheme
    \brief stores the spectator scheme to be used: 
    0 = all final-state charged particles that exist prior to this module being called (default), 
    1 = only the final-state charged particle that the soft photon is calculated to be emitted off.
  */
  /*!
    \var bool m_debug_initProbabilistic
    \brief if off, chooses starting scale by winner-takes all. Default is on.
  */
  /*!
    \var double m_masscutoff
    \brief the mass of the heaviest hadron which can be produced by photon splittings. Note only 
    fermion and scalar splitting functions are currently implemented.
  */
  /*!
    \var double m_t0
    \brief the IR cutoff of the algorithm. Defaults to 4m^2 where m is the mass of the lightest 
    particle which can be produced in photon splittings. 
  */
  /*!
    \var YFS_Particle_List m_remainingphotons
    \brief stores the photons which have not (yet) split.
  */
  /*!
    \var YFS_Particle_Vector m_addedparticles
    \brief stores the particles which are to be added to the blob. 
  */
  /*!
    \var bool m_addedanything 
    \brief flags whether anything has changed so the blob can be updated. 
  */
  /*!
    \var std::vector<size_t> m_splitterIds
    \brief keeps track of the splitter photon associated with each splitting function.
  */
  /*!
    \var int m_NInP
    \brief the number of incoming particles in the blob. Used to classify particles as 
    incoming or outgoing by their ID without referring to the blob. 
  */

  /////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for Photon_Splitter
  /////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn void AddSplitter(Particle* photon, const size_t& id)
    \brief Creates the splitting functions for this photon.
  */
  /*!
    \fn void AddChargedParticle(Particle* p, const size_t& id)
    \brief Assesses whether the charged particle can be a spectator and adds it if so.
  */
  /*!
    \fn bool ClearAll()
    \brief resets the Sudakov for the next event.
  */
  /*!
    \fn void SetCutoff()
    \brief sets the IR cutoff for the algorithm as described above. 
  */
  /*!
    \fn YFS_Particle* DefineInitialConditions(double &t, Vec4D pphoton)
    \brief Reconstructs the emission history of the given photon and returns 
    the recontructed spectator and starting scale. 
  */
  /*!
    \fn bool Run(Blob *)
    \brief Runs the algorithm until the scale is below the cutoff. 
  */
  /*!
    \fn bool Generate (Blob *)
    \brief Generates trial emissions, lets the splitting functions compete and corrects 
    the splitting probability with a veto. 
  */
}

#endif
