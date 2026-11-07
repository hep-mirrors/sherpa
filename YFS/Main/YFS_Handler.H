#ifndef YFS_Main_YFS_Handler_H
#define YFS_Main_YFS_Handler_H

#include "YFS/Main/YFS_Base.H"
#include "YFS/Main/ISR.H"
#include "YFS/Main/FSR.H"
#include "YFS/Main/Coulomb.H"
#include "YFS/Main/Debug.H"
#include "YFS/Main/Photon.H"
#include "YFS/Main/Emission.H"
#include <map>
#include "YFS/NLO/NLO_Base.H"
// #include "BEAM/Main/Beam_Base.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"
// #include "MODEL/Main/Model_Base.H"
#include "YFS/Main/Define_Dipoles.H"
#include "YFS/ACRAIC/Acraic_Base.H"
#include "ATOOLS/Phys/Weights.H"
#include <set>
#include <memory>


// #include <iostream>
// #include <fstream>
namespace ATOOLS { class Blob_Data_Base; }


namespace YFS {

    class YFS_Weight_Builder;

    class YFS_Handler : public YFS_Base {

        friend class YFS_Weight_Builder;

    public:
         // constructor
        YFS_Handler();

        // destructor
        ~YFS_Handler();

    private:

        struct YFS_Event {
          // -- kinematics -----------------------------------------------
          Vec4D_Vector m_plab, m_reallab, m_bornMomenta;
          Vec4D_Vector m_acraiclegs;
          Vec4D_Vector m_sprimeBeams;
          ATOOLS::Poincare m_cmsboost;
          ATOOLS::Vec4D m_Wm, m_Wp;
          ATOOLS::ParticleMomMap m_outparticles;

          // -- photons --------------------------------------------------
          YFS::Photon_Vector m_photons;
          YFS::Photon_Vector m_me_photons;
          Vec4D_Vector m_ISRPhotons, m_FSRPhotons;
          Vec4D_Vector m_fsrphotonsforME;
          ATOOLS::Vec4D m_photonSumISR, m_photonSumFSR;

          // -- weights and form factors ---------------------------------
          double m_yfsweight = 0.;
          double m_formfactor = 1.;
          double m_formfactor_sum = 0.;
          double m_coulSub = 0.;
          double m_coulho = 1.;
          double m_ww_formfact = 1.;
          double m_real = 1.;
          double m_eex = 0.;

          // -- fixed-order accumulators ---------------------------------
          double m_nlo_real = 0., m_nlo_virtual = 0., m_nlo_rv = 0., m_nlo_rr = 0.;
          //! Sum of the fixed-order real corrections above the double real
          double m_nlo_rn = 0.;
          //! YFS: VIRTUAL_COMBINE 1: v times the real correction, (V/B)(R+RR+RN)
          double m_nlo_vxr = 0.;
          //! YFS: BETA0_2: the squared-level beta_0^(2) of YFS.NLO, B v2 times
          //! (1 + real/B) where v multiplies the real (VIRTUAL_COMBINE product)
          double m_nlo_beta02 = 0.;
          //! YFS: DYSON_VP: the squared-level VP resummation of YFS.NLO, placed as v2
          double m_nlo_vp = 0.;
          double m_nlo_real_hardest = 0., m_nlo_rv_hardest = 0.;
          //! YFS: IFI_Real: B (prod_j r_IF(k_j) - 1), the per-photon IF
          //! restoration inside m_nlo_real (0 without IFI_Real)
          double m_nlo_real_ifisoft = 0.;
          double m_nlo_rr_2hardest = 0., m_nlo_real_2hardest = 0.;
          double m_nlo_rv_2hardest = 0.;
          ATOOLS::Weights_Map m_nlo_weightsmap;
          bool m_nlo_current = false;
          std::vector<std::pair<ATOOLS::Vec4D,double> > m_extrealphot;

          // -- ACRAIC, per event ------------------------------------------
          double m_acraicfactor = 0.;
          //! named-weight ratio. Set in GenerateWeight once IFI_Real is known.
          double m_corr_nominal = 1., m_corr_acraic = 0., m_corr_eex = 1.;
          //! ACRAIC: ORDER 2 only - the O(alpha^1) ACRAIC correction of the
          //! same event (named weight ACRAIC_O1), composed exactly as ORDER 1
          double m_corr_acraic_o1 = 0.;
          //! ACRAIC: STRICT_O1_STREAM, (rho_1 + v Re<C_0,C_1>)/rho_crude
          double m_corr_acraic_o1s = 0.;
          double m_corr_acraic_o0 = 0.;   //!< ACRAIC: WEIGHT_O0, rho_0/rho_crude
          //! ACRAIC: FIXED_ORDER_STREAM, the O(alpha) truncation (see GenerateWeight)
          double m_corr_acraic_fo = 0.;
          //! ... and with (1+v/2) on the whole of A_1 (ACRAIC_O1V1)
          double m_corr_acraic_o1v1 = 0.;
          //! ORDER 2, both external-virtual compositions: (v/2) on A_0
          //! (ACRAIC_O2V0), the closest to KKMC CEEX2
          double m_corr_acraic_v0 = 0., m_corr_acraic_v1 = 0.;
          //! ... and with (1+v/2) on the whole of A_2 (ACRAIC_O2V2): the
          //! nominal plus the (v/2) beta_2 terms, an O(alpha^3) estimate
          double m_corr_acraic_v2 = 0.;
          //! ACRAIC: BETA0_2: the ACRAIC correction without beta_0^(2) (ACRAIC_NO_B02)
          double m_corr_acraic_nob02 = 0.;
          //! YFS: DYSON_VP: the ACRAIC correction without the VP resummation (ACRAIC_NO_VP)
          double m_corr_acraic_novp = 0.;
          //! ACRAIC: TNP_STREAMS: the class templates T_c/rho_crude (with the IF
          //! cut) at O(alpha) (ACRAIC_TNP1_*) and at ORDER 2 (ACRAIC_TNP2_*)
          std::vector<double> m_tnp1, m_tnp2;
        };
        YFS_Event m_ev;

        void StartEvent() { m_ev = YFS_Event(); }

        bool   m_setparticles, m_asymbeams;
        int m_negskip;
        long m_nonfinite = 0;
        long m_regcheck_n[2] = {0, 0}, m_regcheck_bad[2] = {0, 0};   //!< [ISR, FSR]
        double m_regcheck_min[2] = {1e300, 1e300};
        long m_regcheck_zero = 0;
        std::unique_ptr<YFS::Coulomb> p_coulomb;
        ATOOLS::Particle_Vector m_particles;
        std::vector<std::unique_ptr<ATOOLS::Particle> > m_particle_store;
        ATOOLS::Flavour_Vector m_flavs;
        std::vector<double> m_mass;
        BEAM::Beam_Spectra_Handler *p_beams = nullptr;
        std::unique_ptr<YFS::ISR> p_isr;
        std::unique_ptr<YFS::Debug> p_debug;
        ATOOLS::nlo_type::code m_nlotype{ATOOLS::nlo_type::born};
        Flavour_Vector m_inflav;
        
        std::unique_ptr<YFS::NLO_Base> p_nlo;
        std::unique_ptr<YFS::FSR> p_fsr;
        std::unique_ptr<YFS::YFS_Form_Factor> p_yfsFormFact;
        std::unique_ptr<YFS::Define_Dipoles> p_dipoles;
        std::unique_ptr<YFS::Acraic_Base> p_acraic;
        std::map<std::string, std::unique_ptr<YFS::Acraic_Base> > m_acraicstore;
        std::string m_acraickey;
        //! CalculateBeta through BuildNamedWeights and their probes (YFS_Weight_Builder.C)
        std::unique_ptr<YFS::YFS_Weight_Builder> p_weights;

        struct Acraic_Stats {
          //! ACRAIC vs the external providers, per point: running mean and worst.
          long   m_cmp_n = 0;
          double m_vsum = 0., m_vworst = 0.;   //!< virtual
          double m_rsum = 0., m_rworst = 0.;   //!< real
          double m_tsum = 0., m_tworst = 0.;   //!< total
          //! (ACRAIC factor)/(EEX factor) over the events where ACRAIC supplied one.
          double m_oesum = 0.;
          long   m_oen = 0;
          //! Events where the ACRAIC weight was unusable and EEX was used instead.
          long   m_bad = 0;
          void AddCompare(double dv, double dr, double dt) {
            m_vsum += dv; m_vworst = std::max(m_vworst, dv);
            m_rsum += dr; m_rworst = std::max(m_rworst, dr);
            m_tsum += dt; m_tworst = std::max(m_tworst, dt);
            ++m_cmp_n;
          }
          void AddOverEex(double ratio) { m_oesum += ratio; ++m_oen; }
        };
        Acraic_Stats m_acraicstats;
        std::map<std::string, long> m_fsrfail;
        long m_fsrtrials = 0;
        void CountFSRFailure(const std::string &what);
        bool ExchangeLinesPassTMin(const ATOOLS::Vec4D_Vector &p) const
        { return ExchangeLinesPassTMin(p, m_exchange_tmin); }
        bool ExchangeLinesPassTMin(const ATOOLS::Vec4D_Vector &p, double tmin) const;
        //! This event's s, handed to every YFS object the handler owns.
        void SetEventS(const double s);
        //! Whether the kept initial-initial dipole was built on this event's
        //! beams (it is built once, and a beam spectrum moves the beams).
        bool IIDipoleIsCurrent();
        PHASIC::Process_Base *m_acraicborn = nullptr, *m_acraicreal = nullptr;
        //! the Born and real process groups the channels above are drawn from
        PHASIC::Process_Base *m_borngroup = nullptr, *m_realgroup = nullptr;
        //! the process group each n-photon real provider was built with, the
        //! channels of a container being drawn from it
        std::map<YFS::Real_Correction*, PHASIC::Process_Base*> m_realprovgroups;
        //! YFS: REAL_BORN_PHOTON_MULTICHANNEL: this event's G (ACRAIC crude)
        double m_bpmc_evtG = 1.;  int m_bpmc_evtn = 0;

    public:

        void Init(std::vector<double> &splimits, const Flavour_Vector &fl);

        // bool On();
        bool MakeYFS(ATOOLS::Vec4D_Vector &p);
        bool MakeYFS();
        void SetLimits(const double &smin, const double &s);
        void SetFlavours(const Flavour_Vector &fl);
        /*!
         * The flavours of the channel of a process group whose weight is
         * being evaluated (e.g. one of the 93 93 children). The group's phase
         * point was made with the container flavours, so only the final-state
         * flavours, masses and the ACRAIC object change: no particle store
         * rebuild, which would discard the momenta. ISR mode only, as the
         * final-state masses enter the photon kinematics otherwise.
         */
        void SetChannelFlavours(const Flavour_Vector &fl);
        void SetMomenta(const ATOOLS::Vec4D_Vector &plab);
        void SetBornMomenta(const ATOOLS::Vec4D_Vector &plab);
        void InitializeACRAIC(const ATOOLS::Flavour_Vector &fl);
        //! Point ACRAIC (and the NLO base) at the Born and real channels of fl
        void SetAcraicChannelProcs(const ATOOLS::Flavour_Vector &fl);
        //! Hand ACRAIC the Born and one-photon real processes, so it can take
        //! its spin amplitudes from Comix rather than hand-coded spinors.
        void SetAcraicProcs(PHASIC::Process_Base *born, PHASIC::Process_Base *real);
        void GenerateWeight();
        //! The YFS.* column names this configuration emits. Fixed once per
        //! event from configuration alone, and the only names emit() will
        //! write -- a per-event code path must not be able to conjure a column
        //! that other events lack, since HepMC binds the set from event one.
        std::set<std::string> m_wnames;

        void CreatMomentumMap();
        void MakeACRAIC();
        Vec4D_Vector GetPhotons();
        inline double GetWeight() const      { return m_ev.m_yfsweight;  }
        inline const ATOOLS::Weights_Map& GetNLOWeightsMap() const { return m_ev.m_nlo_weightsmap; }
        inline ATOOLS::Vec4D GetPhotonSum() const      { return m_ev.m_photonSumISR;  }
        // The event momenta/photons are stored in the boosted incoming-pair
        // rest frame (see SetBornMomenta). These getters feed the event blob,
        // which must be in the original lab frame, so boost the returned copy
        // back. Boosting the copy (not the member) leaves the internal COM
        // state intact for the weight/ME machinery and is independent of when
        // the blob is filled relative to the weight calculation.
        inline ATOOLS::Vec4D_Vector GetISRPhotons() const {
          Vec4D_Vector k(m_ev.m_ISRPhotons); ToLab(k); return k;
        }
        inline ATOOLS::Vec4D_Vector GetFSRPhotons() const {
          Vec4D_Vector k(m_ev.m_FSRPhotons); ToLab(k); return k;
        }
        inline ATOOLS::Flavour GetFlavour(const int i) const { return m_flavs[i];  }
        inline void SetBorn(double b) { m_born = b; }

        inline yfsmode::code  Mode()   const {return m_mode; }
        //! YFS: Integrate_NLO - whether CalculateBeta() also runs while integrating
        inline bool IntegrateNLO() const { return m_int_nlo != 0; }

        // Undo the incoming-pair-rest-frame boost that SetBornMomenta applied
        // (m_ev.m_cmsboost: lab->COM), returning the blob-facing momenta to the
        // original lab frame. m_ev.m_cmsboost is the identity when no boost was
        // applied (already-at-rest incoming, or pure-FSR mode), so this is a
        // safe no-op in those cases - a standard at-rest COM collider included.
        inline void ToLab(Vec4D &q) const { m_ev.m_cmsboost.BoostBack(q); }
        inline void ToLab(Vec4D_Vector &v) const { for (auto &q : v) ToLab(q); }
        inline Vec4D_Vector GetMomenta() { Vec4D_Vector p(m_ev.m_plab); ToLab(p); return p; }
        inline Vec4D GetMomenta(int i) { Vec4D q(m_ev.m_plab[i]); ToLab(q); return q; }
        inline Vec4D_Vector BornMomenta() { Vec4D_Vector p(m_ev.m_bornMomenta); ToLab(p); return p; }
        //! the event before final-state radiation (after ISR), in the lab:
        //! the momenta CalculateFSR started from
        inline Vec4D_Vector PreFSRMomenta() const { Vec4D_Vector p(m_ev.m_reallab); ToLab(p); return p; }
        // out-particle map for FSR event blobs, boosted back to the lab frame
        inline ATOOLS::ParticleMomMap GetOutParticles() const {
          ATOOLS::ParticleMomMap m(m_ev.m_outparticles);
          for (auto &kv : m) m_ev.m_cmsboost.BoostBack(kv.second);
          return m;
        }

        inline void SetInFlav(Flavour &fl) {m_inflav.push_back(fl);}
        inline Flavour GetInFlav(int i) {return m_inflav[i];}
        inline YFS::Photon_Vector       &Photons()       { return m_ev.m_photons; }
        inline const YFS::Photon_Vector &Photons() const { return m_ev.m_photons; }
        inline void SetNLOType(ATOOLS::nlo_type::code nlo_type) {m_nlotype = nlo_type;}
        // Per-process request to drop the Born from the YFS weight, used when a
        // process is declared without 'B' in NLO_Part (e.g. an extra explicit
        // one-photon process supplying RR/RV, whose Born is already the real
        // correction of the Born-level process). The card setting No_Born stays
        // a global override. The handler is shared, so YFS_Process re-applies
        // this in MakeActive().
        inline void SetNoBorn(bool nb) {m_no_born = nb || m_no_born_setting;}
        //! Raw observer: the handler keeps ownership.
        inline NLO_Base *NLO() {return p_nlo.get(); }
        //! The three members reached from outside YFS. Accessors rather than
        //! public data so the storage can change without touching SHERPA and
        //! PHASIC++; Dipoles() is non-const because YFS_Process configures the
        //! resonance map through it.
        inline Define_Dipoles *Dipoles() { return p_dipoles.get(); }
        inline const ATOOLS::Particle_Vector &Particles() const { return m_particles; }
        NLO_Base *EnsureNLO();
        inline void resetparticles() {m_setparticles=false;}

        void SetBeam(BEAM::Beam_Spectra_Handler *beam);
        inline void SetSprime(double sp) { m_sp = sp; }
        inline void SetRunMode(const int &rmode) { m_rmode = rmode; }

        double Eikonal(Vec4D k);
        long m_acrv_n = 0, m_acrv_nomatch = 0, m_acrv_nom1 = 0;
        bool CalculateFSR();
        bool CalculateISR();
        bool CalculateFSR(Vec4D_Vector &k);
        void Reset();
        void AddFormFactor();
        /*!
          The W+ and W- momenta of the final state: an on-shell W+ W- pair, or
          the daughters of a four-fermion final state paired by flavour (CKM
          diagonal): W+ -> nu_l l+ or u-type quark + d-type antiquark of the
          same generation, W- the charge conjugates. False if neither exists.
          wp is BVV_WW's p1.
        */
        bool WPairMomenta(const Vec4D_Vector &p, Vec4D &wp, Vec4D &wm) const;
        //! hand ACRAIC its W decay groups (ACRAIC: W_STAGES); false = flat groups
        bool SetACRAICWStageGroups();
        void CalculateCoulomb();
        void CalculateWWForm();
        void YFSDebug(double W);
        void CheckResonance();
    };// end of class YFS_Handler



}// end of namespace YFS

#endif
