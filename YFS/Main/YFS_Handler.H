#ifndef YFS_Main_YFS_Handler_H
#define YFS_Main_YFS_Handler_H

#include "YFS/Main/YFS_Base.H"
#include "YFS/Main/YFS_FB_Analysis.H"
#include "YFS/Main/ISR.H"
#include "YFS/Main/FSR.H"
#include "YFS/Main/Coulomb.H"
#include "YFS/Main/Debug.H"
#include "YFS/Main/Photon.H"
#include "YFS/Main/Emission.H"
#include "YFS/NLO/NLO_Base.H"
// #include "BEAM/Main/Beam_Base.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"
// #include "MODEL/Main/Model_Base.H"
#include "YFS/Main/Define_Dipoles.H"
#include "YFS/CEEX/Ceex_Base.H"
#include "YFS/NLO/NLO_Base.H"
#include "ATOOLS/Phys/Momenta_Stretcher.H"
#include "ATOOLS/Phys/Weights.H"
#include "PHOTONS++/PhotonSplitter/Photon_Splitter.H"


// #include <iostream>
// #include <fstream>
namespace ATOOLS { class Blob_Data_Base; }
namespace YFS    {class Real_ff;         }
namespace YFS    {class Coulomb;         }


namespace YFS {

    class YFS_Handler : public YFS_Base {

    public:
         // constructor
        YFS_Handler();

        // destructor
        ~YFS_Handler();
        double m_splimits[3];
        double m_yfsweight, m_form, m_fsrW;
        bool   m_setparticles, m_asymbeams;
        double m_real, m_virt;
        double m_nlo_real, m_nlo_virtual, m_nlo_rv, m_nlo_rr;
        double m_nlo_real_hardest, m_nlo_rv_hardest, m_nlo_rr_2hardest;
        double m_nlo_real_2hardest, m_nlo_rv_2hardest;
        ATOOLS::Weights_Map m_nlo_weightsmap;
        double m_nbar, m_formfactor, m_betatWW, m_coulSub, m_ww_formfact,m_eex;
        int m_negskip;
        YFS::Coulomb *p_coulomb;
        PHASIC::Process_Base *proc;
        ATOOLS::Particle_Vector m_particles;
        ATOOLS::ParticleMomMap m_inparticles, m_outparticles;

        // std::string m_name, m_type, m_debugDIR_ISR, m_debugDIR_FSR;
        ATOOLS::Poincare m_cmsboost;
        ATOOLS::Vec4D   m_Wm, m_Wp;
        Vec4D_Vector m_plab, m_reallab, m_FSRPhotons, m_ISRPhotons, m_bornMomenta;
        // Every photon of the event, each tagged with the dipole it came from,
        // so its ISR/FSR origin never has to be recovered positionally.
        YFS::Photon_Vector m_photons;
        // The photons the NLO layer sees: ISR plus the FSR *ME* photons
        // (captured before HidePhotons), which is what InitNLO hands it and is
        // NOT the same list as m_photons.
        YFS::Photon_Vector m_me_photons;
        ATOOLS::Flavour_Vector m_flavs;
        std::vector<double> m_mass;
        ATOOLS::Momenta_Stretcher m_stretcher;
        BEAM::Beam_Spectra_Handler *p_beams;
        PHOTONS::Photon_Splitter *p_splitter;

    public:

        YFS::ISR *p_isr;
        YFS::Debug *p_debug;
        ATOOLS::nlo_type::code m_nlotype;
        // True only while m_nlo_real/m_nlo_virtual/m_nlo_rv/m_nlo_rr belong to
        // the event currently being weighted. CalculateBeta() bails out early
        // during integration when Integrate_NLO is off (m_rmode is set to 1
        // only in Phase_Space_Handler::OneEvent), leaving m_real at its Reset()
        // value of 1 and the m_nlo_* accumulators holding the PREVIOUS event's
        // numbers. The named YFS.* weights must not be built from those.
        bool m_nlo_current = false;
        YFS::Define_Dipoles *p_dipole;
        Vec4D m_photonSumISR, m_photonSumFSR;
        Vec4D_Vector m_fsrphotonsforME, m_isrphotonsforME;
        Flavour_Vector m_inflav;
        Flavour_Vector m_outflav;
        
        YFS::NLO_Base *p_nlo;
        YFS::YFS_FB_Analysis *p_fb;
        YFS::FSR *p_fsr;
        YFS::YFS_Form_Factor *p_yfsFormFact;
        YFS::Define_Dipoles *p_dipoles;
        YFS::Ceex_Base *p_ceex;

        void Init(std::vector<double> &splimits, const Flavour_Vector &fl);

        // bool On();
        bool MakeYFS(ATOOLS::Vec4D_Vector &p);
        bool MakeYFS();
        void SetLimits(const double &smin);
        void SetFlavours(const Flavour_Vector &fl);
        void SetMomenta(const ATOOLS::Vec4D_Vector &plab);
        void SetBornMomenta(const ATOOLS::Vec4D_Vector &plab);
        void InitializeCEEX(const ATOOLS::Flavour_Vector &fl);
        void CalculateBeta();
        void GenerateWeight();
        //! Fills m_nlo_weightsmap with the YFS.* ratios. Defined next to
        //! GenerateWeight, which supplies the two weights it compares.
        void BuildNamedWeights(double w_lo, double w_full);
        void CreatMomentumMap();
        void MakeCEEX();
        void InitNLO();
        Vec4D_Vector GetPhotons();
        inline double GetWeight() const      { return m_yfsweight;  }
        inline const ATOOLS::Weights_Map& GetNLOWeightsMap() const { return m_nlo_weightsmap; }
        inline ATOOLS::Vec4D GetPhotonSum() const      { return m_photonSumISR;  }
        inline ATOOLS::Vec4D GetPhotonSumFSR() const      { return m_photonSumFSR;  }
        // The event momenta/photons are stored in the boosted incoming-pair
        // rest frame (see SetBornMomenta). These getters feed the event blob,
        // which must be in the original lab frame, so boost the returned copy
        // back. Boosting the copy (not the member) leaves the internal COM
        // state intact for the weight/ME machinery and is independent of when
        // the blob is filled relative to the weight calculation.
        inline ATOOLS::Vec4D_Vector GetISRPhotons() const {
          Vec4D_Vector k(m_ISRPhotons); ToLab(k); return k;
        }
        inline ATOOLS::Vec4D_Vector GetFSRPhotons() const {
          Vec4D_Vector k(m_FSRPhotons); ToLab(k); return k;
        }
        inline ATOOLS::Flavour GetFlavour(const int i) const { return m_flavs[i];  }
        inline void SetBorn(double b) { m_born = b; }
        inline double GetBorn() { return m_born; }

        inline yfsmode::code  Mode()   const {return m_mode; }

        inline bool UsingLoopTool()   const {return m_looptool; }
        inline bool FillBlob()   const {return m_fillblob; }
        // Undo the incoming-pair-rest-frame boost that SetBornMomenta applied
        // (m_cmsboost: lab->COM), returning the blob-facing momenta to the
        // original lab frame. m_cmsboost is the identity when no boost was
        // applied (already-at-rest incoming, or pure-FSR mode), so this is a
        // safe no-op in those cases - a standard at-rest COM collider included.
        inline void ToLab(Vec4D &q) const { m_cmsboost.BoostBack(q); }
        inline void ToLab(Vec4D_Vector &v) const { for (auto &q : v) ToLab(q); }
        inline Vec4D_Vector GetMomenta() { Vec4D_Vector p(m_plab); ToLab(p); return p; }
        inline Vec4D GetMomenta(int i) { Vec4D q(m_plab[i]); ToLab(q); return q; }
        inline Vec4D_Vector BornMomenta() { Vec4D_Vector p(m_bornMomenta); ToLab(p); return p; }
        // out-particle map for FSR event blobs, boosted back to the lab frame
        inline ATOOLS::ParticleMomMap GetOutParticles() const {
          ATOOLS::ParticleMomMap m(m_outparticles);
          for (auto &kv : m) m_cmsboost.BoostBack(kv.second);
          return m;
        }

        inline void SetInFlav(Flavour &fl) {m_inflav.push_back(fl);}
        inline Flavour GetInFlav(int i) {return m_inflav[i];}
        inline YFS::Photon_Vector       &Photons()       { return m_photons; }
        inline const YFS::Photon_Vector &Photons() const { return m_photons; }
        inline void SetNLOType(ATOOLS::nlo_type::code nlo_type) {m_nlotype = nlo_type;}
        // Per-process request to drop the Born from the YFS weight, used when a
        // process is declared without 'B' in NLO_Part (e.g. an extra explicit
        // one-photon process supplying RR/RV, whose Born is already the real
        // correction of the Born-level process). The card setting No_Born stays
        // a global override. The handler is shared, so YFS_Process re-applies
        // this in MakeActive().
        inline void SetNoBorn(bool nb) {m_no_born = nb || m_no_born_setting;}
        inline NLO_Base *NLO() {return p_nlo; }
        NLO_Base *EnsureNLO();
        inline void resetparticles() {m_setparticles=false;}

        void SetBeam(BEAM::Beam_Spectra_Handler *beam);
        inline void SetSprime(double sp) { m_sp = sp; }
        inline void SetRunMode(const int &rmode) { m_rmode = rmode; }

        double Eikonal(Vec4D k);
        double CalculateNLO();
        bool CalculateFSR();
        bool CalculateISR();
        bool CalculateFSR(Vec4D_Vector &k);
        bool CheckMomentumConservation();
        void CheckMasses();
        void Reset();
        void AddFormFactor();
        void MakeWWVecs(ATOOLS::Vec4D_Vector p);
        void CalculateCoulomb();
        void CalculateWWForm();
        void YFSDebug(double W);
        void SplitPhotons(ATOOLS::Blob * blob);
        void CheckResonance();
    };// end of class YFS_Handler



}// end of namespace YFS

#endif
