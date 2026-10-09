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
#include <map>
#include "YFS/NLO/NLO_Base.H"
// #include "BEAM/Main/Beam_Base.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"
// #include "MODEL/Main/Model_Base.H"
#include "YFS/Main/Define_Dipoles.H"
#include "YFS/ACRAIC/Acraic_Base.H"
#include "ATOOLS/Phys/Momenta_Stretcher.H"
#include "ATOOLS/Phys/Weights.H"
#include <set>
#include <memory>


// #include <iostream>
// #include <fstream>
namespace ATOOLS { class Blob_Data_Base; }


namespace YFS {

    class YFS_Handler : public YFS_Base {

    public:
         // constructor
        YFS_Handler();

        // destructor
        ~YFS_Handler();

    private:

        /*!
          Everything that belongs to ONE event, in one object so that the
          event boundary is a construction rather than a list of assignments.

          The handler used to hold these flat alongside its run- and
          process-scope state, and Reset() cleared eight of them; the other
          thirty were simply overwritten before the next read, an invariant
          nothing enforced. m_ev.m_nlo_current is the surviving trace of where that
          went wrong - the named YFS.* weights were once built from the
          PREVIOUS event's accumulators. Pure-FSR mode was worse: it enters at
          CalculateFSR() and never reached MakeYFS(), so Reset() never ran at
          all.

          StartEvent() now assigns a default-constructed YFS_Event at the one
          point both modes pass through (SetBornMomenta), so a field added
          later starts fresh without anyone having to remember to clear it.
          The defaults below are the values the constructor and Reset() used
          to establish - in particular m_ev.m_real, m_ev.m_formfactor and m_ev.m_ww_formfact
          are 1, not 0.
        */
        struct YFS_Event {
          // -- kinematics -----------------------------------------------
          Vec4D_Vector m_plab, m_reallab, m_bornMomenta;
          /*!
            m_plab with the final pairs that hid photons replaced by their
            record momenta (AbsorbHiddenPhotons): legs that balance against
            the VISIBLE photons. Empty when nothing was hidden. ACRAIC takes its
            final legs from here, since it gets only the visible photons.
          */
          Vec4D_Vector m_acraiclegs;
          /*!
            The ISR-reduced beams BEFORE ToLab, so they are in the same frame
            as m_plab[2..] (which stay in the incoming-pair rest frame). With
            those final legs they are a phase-space point that balances and
            whose invariant is s'. m_plab[0..1] are the same momenta boosted to
            the lab, so {m_plab[0],m_plab[1],m_plab[2..]} is frame-mixed and
            balances only for symmetric beams.
          */
          Vec4D_Vector m_sprimeBeams;
          ATOOLS::Poincare m_cmsboost;
          ATOOLS::Vec4D m_Wm, m_Wp;
          ATOOLS::ParticleMomMap m_inparticles, m_outparticles;

          // -- photons --------------------------------------------------
          //! Every photon of the event, each tagged with the dipole it came
          //! from, so its ISR/FSR origin never has to be recovered positionally.
          YFS::Photon_Vector m_photons;
          //! The photons the NLO layer sees: ISR plus the FSR *ME* photons
          //! (captured before HidePhotons). NOT the same list as m_photons.
          YFS::Photon_Vector m_me_photons;
          Vec4D_Vector m_ISRPhotons, m_FSRPhotons;
          Vec4D_Vector m_isrphotonsforME, m_fsrphotonsforME;
          ATOOLS::Vec4D m_photonSumISR, m_photonSumFSR;

          // -- weights and form factors ---------------------------------
          double m_yfsweight = 0.;
          double m_formfactor = 1.;
          //! FormFactorSum() as evaluated by AddFormFactor, i.e. the exponent
          //! m_formfactor was actually built from. Cached rather than
          //! recomputed so YFS.NLO_FixedOrder's (1+Y)/exp(Y) is provably the
          //! same Y in numerator and denominator.
          double m_formfactor_sum = 0.;
          double m_coulSub = 0.;
          double m_ww_formfact = 1.;
          double m_real = 1.;
          //! The EEX beta-expansion estimate, assigned in CalculateBeta only
          //! when BETA > 0; read by the YFS.EEX named weight, where an
          //! unassigned read produced a NaN column rather than a missing one.
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
          /*!
            True only once the accumulators above belong to THIS event.
            CalculateBeta() bails out early during integration when
            Integrate_NLO is off (m_rmode is set to 1 only in
            Phase_Space_Handler::OneEvent), and the named YFS.* weights must
            not be built from accumulators it never filled. Structurally false
            at the start of every event now, rather than by assignment.
          */
          bool m_nlo_current = false;
          /*!
            Per-photon real from the EXTERNAL provider, captured the moment
            CalculateReal() returns. It cannot be read later: CalculateNLO()
            calls InitNLO() again before the virtual, and that REBUILDS
            p_nlo->m_photons from scratch, so every m_beta10 is back to zero
            by the time anything downstream looks.
          */
          std::vector<std::pair<ATOOLS::Vec4D,double> > m_extrealphot;

          // -- ACRAIC, per event ------------------------------------------
          //! ACRAIC's O(alpha) factor for this event; 0 = ACRAIC produced none.
          double m_acraicfactor = 0.;
          //! The O(alpha) correction actually applied, and ACRAIC's, for the
          //! named-weight ratio. Set in GenerateWeight once IFI_Real is known.
          double m_corr_nominal = 1., m_corr_acraic = 0., m_corr_eex = 1.;
          //! ACRAIC: ORDER 2 only - the O(alpha^1) ACRAIC correction of the
          //! same event (named weight ACRAIC_O1), composed exactly as ORDER 1
          double m_corr_acraic_o1 = 0.;
          double m_corr_acraic_o0 = 0.;   //!< ACRAIC: WEIGHT_O0, rho_0/rho_crude
          //! ... and with (1+v/2) on the whole of A_1 (ACRAIC_O1V1)
          double m_corr_acraic_o1v1 = 0.;
          //! ORDER 2, both external-virtual compositions: (v/2) on A_0
          //! (ACRAIC_O2V0) and on A_1, KKMC's (1+v/2) on beta_1 (ACRAIC_O2V1)
          double m_corr_acraic_v0 = 0., m_corr_acraic_v1 = 0.;
          //! ... and with (1+v/2) on the whole of A_2 (ACRAIC_O2V2): the
          //! nominal plus the (v/2) beta_2 terms, an O(alpha^3) estimate
          double m_corr_acraic_v2 = 0.;
          //! ACRAIC: BETA0_2: the ACRAIC correction without beta_0^(2) (ACRAIC_NO_B02)
          double m_corr_acraic_nob02 = 0.;
          //! YFS: DYSON_VP: the ACRAIC correction without the VP resummation (ACRAIC_NO_VP)
          double m_corr_acraic_novp = 0.;
        };
        YFS_Event m_ev;

        //! Begin a new event: discard everything the previous one left behind.
        //! Called from SetBornMomenta, the one entry point both the ISR/ISRFSR
        //! and the pure-FSR paths pass through (Phase_Space_Point.C:229,241).
        void StartEvent() { m_ev = YFS_Event(); }

        bool   m_setparticles, m_asymbeams;
        // m_nbar and m_betatWW removed: both shadowed YFS_Base's. m_betatWW
        // is set from the card (WW_BETAT) on the base, but the test in
        // CalculateBeta resolved to this copy, which nothing ever wrote.
        //! The EEX beta-expansion estimate of the correction, assigned in
        //! CalculateBeta only when BETA > 0. Initialised and reset per event
        //! because it is read by the YFS.EEX named weight, and an unassigned
        //! read there produced a NaN column rather than a missing one.
        //! FormFactorSum() as evaluated by AddFormFactor, i.e. the exponent
        //! that m_ev.m_formfactor was actually built from. Cached rather than
        //! recomputed so YFS.NLO_FixedOrder's (1+Y)/exp(Y) is provably the same
        //! Y in numerator and denominator, whatever happens to the dipoles in
        //! between.
        int m_negskip;
        std::unique_ptr<YFS::Coulomb> p_coulomb;
        /*!
          The externally visible view (SHERPA's Signal_Processes copies it as a
          Particle_Vector), holding observing pointers into m_particle_store.
          The store is what owns them; the two are filled and cleared together
          in SetFlavours.
        */
        ATOOLS::Particle_Vector m_particles;
        std::vector<std::unique_ptr<ATOOLS::Particle> > m_particle_store;

        // std::string m_name, m_type, m_debugDIR_ISR, m_debugDIR_FSR;
        // Every photon of the event, each tagged with the dipole it came from,
        // so its ISR/FSR origin never has to be recovered positionally.
        // The photons the NLO layer sees: ISR plus the FSR *ME* photons
        // (captured before HidePhotons), which is what InitNLO hands it and is
        // NOT the same list as m_ev.m_photons.
        ATOOLS::Flavour_Vector m_flavs;
        std::vector<double> m_mass;
        ATOOLS::Momenta_Stretcher m_stretcher;
        //! Not owned: handed in by SetBeam, lives in Beam_Spectra_Handler.
        BEAM::Beam_Spectra_Handler *p_beams = nullptr;
        std::unique_ptr<YFS::ISR> p_isr;
        std::unique_ptr<YFS::Debug> p_debug;
        // born unless YFS_Process (NLO_Mode YFS) says otherwise: a plain LO
        // process never calls SetNLOType, and an uninitialised value sent it
        // into CalculateNLO with no NLO_Base (null p_nlo, crash in InitNLO)
        ATOOLS::nlo_type::code m_nlotype{ATOOLS::nlo_type::born};
        // True only while m_ev.m_nlo_real/m_ev.m_nlo_virtual/m_ev.m_nlo_rv/m_ev.m_nlo_rr belong to
        // the event currently being weighted. CalculateBeta() bails out early
        // during integration when Integrate_NLO is off (m_rmode is set to 1
        // only in Phase_Space_Handler::OneEvent), leaving m_ev.m_real at its Reset()
        // value of 1 and the m_nlo_* accumulators holding the PREVIOUS event's
        // numbers. The named YFS.* weights must not be built from those.
        //! Beam flavours, set once at initialisation; REMNANTS reads them
        //! back through GetInFlav. Not per-event.
        Flavour_Vector m_inflav;
        
        std::unique_ptr<YFS::NLO_Base> p_nlo;
        std::unique_ptr<YFS::YFS_FB_Analysis> p_fb;
        std::unique_ptr<YFS::FSR> p_fsr;
        std::unique_ptr<YFS::YFS_Form_Factor> p_yfsFormFact;
        std::unique_ptr<YFS::Define_Dipoles> p_dipoles;
        std::unique_ptr<YFS::Acraic_Base> p_acraic;
        //! One Acraic_Base per process (flavour set): the handler is shared by
        //! every process of the card, and a ACRAIC object holds process-specific
        //! state (stages, Comix/loop helicity maps, calibration). p_acraic is
        //! the current process's; the others wait here (InitializeACRAIC).
        std::map<std::string, std::unique_ptr<YFS::Acraic_Base> > m_acraicstore;
        std::string m_acraickey;
        //! Bookkeeping for the ACRAIC-vs-EEX comparison at end of run.
        //! ACRAIC's O(alpha) factor for the current event; 0 = ACRAIC produced none.
        //! The O(alpha) correction actually applied, and ACRAIC's, for the
        //! named-weight ratio. Set in GenerateWeight once IFI_Real is known.
        /*!
          YFS: ACRAIC_Compare - point-by-point ACRAIC against the EXTERNAL
          providers. Prints the first N events and accumulates the spread, so
          a disagreement can be attributed to the virtual or the real rather
          than just seen in the total.
        */
        void   AcraicCompare();
        //! YFS: FSR_DENSITY_GATE - the FSR generator density against the eikonals.
        void ProbeFSRDensity(YFS::Dipole &D);
        //! ACRAIC: GOLDEN - bit-exact reference stream of the ACRAIC factor (@@@ ACRAICGOLD).
        void ProbeGolden(const double r0, const double r1);
        //! ACRAIC: WEIGHT_PROBE - what makes a heavy ACRAIC event (@@@ ACRAICWT).
        void ProbeAcraicWeight();
        //! ACRAIC: WEIGHT_PROBE - the ACRAIC factor next to the YFS.NLO one (@@@ ACRAICCMP).
        void ProbeAcraicCompare(const bool haveacraic);
        //! ACRAIC: REF_PROBE - every momentum of a one-photon event (@@@ REFP).
        void ProbeReference(const bool haveacraic);
        //! YFS: REAL_STAB - the pre/post-emission momentum residual of each resonant pair (@@@ PREPOST).
        void ProbeRealStab();
        //! YFS: RR_PROBE - the O(alpha^2) weight pieces and exactness identities (@@@ NNLOEV/RRID/RVID).
        void ProbeNNLOEvent();
        //! YFS: LO_PROBE - how CalculateBeta moved the LO pieces (@@@ LOPRE).
        void ProbeLOPre(const double pre_yfsw, const double pre_ff,
                        const double pre_v, const double pre_born);
        //! YFS: LO_PROBE - the Born-level weight pieces (@@@ LOPROBE).
        void ProbeLO();
        //! ACRAIC: WEIGHT_PROBE - the final per-column factors (@@@ WPROBE).
        void ProbeFinalWeights(const double corr_eex, const double corr_acraic,
                               const double wif, const bool acraico2,
                               const double corr_acraic_o1,
                               const double corr_acraic_v0,
                               const double corr_acraic_v1);
        //! YFS: LO_PROBE - the final LO and full weights (@@@ LOPROBE2).
        void ProbeLOFinal(const double w_lo, const double w_full,
                          const bool acraic_nom, const double corr_acraic,
                          const double corr_eex);
        //! YFS: SIGN_PROBE - the hardest photons, m_ll and both weights, per event.
        void ProbeSign(const double w_full);
        /*!
          Per-photon real from the EXTERNAL provider, captured the moment
          CalculateReal() returns. It cannot be read later: CalculateNLO()
          calls InitNLO() again before the virtual, and that REBUILDS
          p_nlo->m_photons from scratch, so every m_beta10 is back to zero by
          the time anything downstream looks.
        */
        /*!
          End-of-run bookkeeping for the two ACRAIC cross-checks. Grouped rather
          than left as ten loose members because none of it takes part in the
          event: it accumulates and is printed once by the destructor, so it
          has no business sitting beside the per-event momenta.
        */
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
        /*!
          YFS: FSR_FAILURE_STATS - why CalculateFSR returned false or set a
          zero weight, per dipole and failure kind, reported at the end of
          the run. A trial that fails here is a zero-weight event, so a high
          failure rate lowers every column; the pole scheme's factor-5
          deficit on cc_em_mup (NOTES-w-stages-2026-09-27.md 1.4) is what
          this was added to find.
        */
        std::map<std::string, long> m_fsrfail;
        long m_fsrtrials = 0;
        void CountFSRFailure(const std::string &what);
        bool ExchangeLinesPassTMin(const ATOOLS::Vec4D_Vector &p) const;
        //! This event's s, handed to every YFS object the handler owns.
        void SetEventS(const double s);
        //! Whether the kept initial-initial dipole was built on this event's
        //! beams (it is built once, and a beam spectrum moves the beams).
        bool IIDipoleIsCurrent();
        PHASIC::Process_Base *m_acraicborn = nullptr, *m_acraicreal = nullptr;
        //! YFS: REAL_BORN_PHOTON_MULTICHANNEL: this event's G (ACRAIC crude)
        double m_bpmc_evtG = 1.;  int m_bpmc_evtn = 0;

    public:

        void Init(std::vector<double> &splimits, const Flavour_Vector &fl);

        // bool On();
        bool MakeYFS(ATOOLS::Vec4D_Vector &p);
        bool MakeYFS();
        void SetLimits(const double &smin, const double &s);
        void SetFlavours(const Flavour_Vector &fl);
        void SetMomenta(const ATOOLS::Vec4D_Vector &plab);
        void SetBornMomenta(const ATOOLS::Vec4D_Vector &plab);
        void InitializeACRAIC(const ATOOLS::Flavour_Vector &fl);
        //! Hand ACRAIC the Born and one-photon real processes, so it can take
        //! its spin amplitudes from Comix rather than hand-coded spinors.
        void SetAcraicProcs(PHASIC::Process_Base *born, PHASIC::Process_Base *real);
        void CalculateBeta();
        void GenerateWeight();
        //! Fills m_ev.m_nlo_weightsmap with the YFS.* ratios. Defined next to
        //! GenerateWeight, which supplies the two weights it compares.
        void BuildNamedWeights(double w_lo, double w_full);
        bool BuildLadderWeights(ATOOLS::Weights &w);
        //! The YFS.* column names this configuration emits. Fixed once per
        //! event from configuration alone, and the only names emit() will
        //! write -- a per-event code path must not be able to conjure a column
        //! that other events lack, since HepMC binds the set from event one.
        std::set<std::string> m_wnames;

        //! Per-event assertions for the invariants that could not be made
        //! structural. Each one corresponds to a bug that reached a physics
        //! result before being noticed, so the point is to fail where the
        //! invariant breaks rather than several frames downstream. Off unless
        //! YFS: CHECK_INVARIANTS is set.
        void CheckInvariants() const;
        void CreatMomentumMap();
        void MakeACRAIC();
        void InitNLO();
        Vec4D_Vector GetPhotons();
        inline double GetWeight() const      { return m_ev.m_yfsweight;  }
        inline const ATOOLS::Weights_Map& GetNLOWeightsMap() const { return m_ev.m_nlo_weightsmap; }
        inline ATOOLS::Vec4D GetPhotonSum() const      { return m_ev.m_photonSumISR;  }
        inline ATOOLS::Vec4D GetPhotonSumFSR() const      { return m_ev.m_photonSumFSR;  }
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
        inline double GetBorn() { return m_born; }

        inline yfsmode::code  Mode()   const {return m_mode; }

        inline bool UsingLoopTool()   const {return m_looptool; }
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
        inline const ATOOLS::ParticleMomMap &OutParticles() const { return m_ev.m_outparticles; }
        NLO_Base *EnsureNLO();
        inline void resetparticles() {m_setparticles=false;}

        void SetBeam(BEAM::Beam_Spectra_Handler *beam);
        inline void SetSprime(double sp) { m_sp = sp; }
        inline void SetRunMode(const int &rmode) { m_rmode = rmode; }

        double Eikonal(Vec4D k);
        double CalculateNLO();
        //! ACRAIC: REAL_VIRTUAL 2 - sum_j dv_j Re<A_1, M_1j>/4 with dv_j from
        //! YFS.NLO's RV (NLO_Base::m_rvdv), matched to ACRAIC photons by lab
        //! momentum; 0 without YFS: RV_MODE 1
        double AcraicRealVirtualRemainder();
        //! ACRAIC: REAL_VIRTUAL helicity: the same, resolved by helicity
        //! (Acraic_RealVirtualHelicity.C); per-photon fallback to dv/2
        double AcraicRealVirtualHelicityRemainder();
        /*!
          ACRAIC: BETA0_2 with an external virtual: d2 on A_0 alone,
            |A_X + (v/2) A_v + d2 A_0|^2 - |A_X + (v/2) A_v|^2
              = 2 d2 Re<A_0, A_X + (v/2) A_v> + d2^2 |A_0|^2
          (per 4 initial helicities), from the overlaps r0x = Re<A_0,A_X>/4
          and r0v = Re<A_0,A_v>/4.
        */
        double AcraicBeta02Rho(double d2, double v, double r0x, double r0v) const;
        //! BETA0_2: the per-event pieces, on the event ACRAIC dumped (CHECK_XS)
        void Beta02DumpProbe(double v) const;
        long m_acrv_n = 0, m_acrv_nomatch = 0, m_acrv_nom1 = 0;
        bool CalculateFSR();
        bool CalculateISR();
        bool CalculateFSR(Vec4D_Vector &k);
        bool CheckMomentumConservation();
        void CheckMasses();
        void Reset();
        void AddFormFactor();
        void MakeWWVecs(ATOOLS::Vec4D_Vector p);
        //! hand ACRAIC its W decay groups (ACRAIC: W_STAGES); false = flat groups
        bool SetACRAICWStageGroups();
        void CalculateCoulomb();
        void CalculateWWForm();
        void YFSDebug(double W);
        void CheckResonance();
    };// end of class YFS_Handler



}// end of namespace YFS

#endif
