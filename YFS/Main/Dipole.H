#ifndef YFS_Main_Dipole_H
#define YFS_Main_Dipole_H

#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Math/Poincare.H"
#include <optional>
#include "METOOLS/Loops/Divergence_Array.H"
#include "YFS/NLO/Virtual.H"
#include "YFS/Main/Soft_Region.H"
#include <vector>
#include <cmath>

using namespace ATOOLS;

namespace YFS{

    class ISR;
    class FSR;
    struct EmissionResult;


    // Kaellen function lambda(x,y,z) = x^2+y^2+z^2-2(xy+xz+yz), returned as
    // |lambda| so sqrt(Lambda(...)/s) stays real below threshold.
    inline double Lambda(double x, double y, double z) {
      return std::abs(x * x + y * y + z * z - 2. * x * y - 2. * x * z - 2. * y * z);
    }

    struct dipoletype {
    enum code {
      initial  = 0,
      final    = 1,
      ifi      = 2
    };
  };




    class Dipole
    {
    public:
        Dipole(Flavour_Vector const &fl, Vec4D_Vector const &mom,  
                Vec4D_Vector const &born, dipoletype::code d, const double alpha);
        ~Dipole();

        std::vector<double> m_masses, m_charges, m_theta, m_Q;
        Vec4D_Vector m_dipolePhotons, m_dipolePhotonsEEX;
        std::vector<std::string> m_names;
        double m_QiQj, m_Qi, m_Qj, m_thetai, m_thetaj, m_thetaij, m_sp, m_b1, m_b2, m_gamma, m_gammap, m_alp, m_alpi,m_rescale_alpha;
        SoftRegion m_genregion, m_expregion;
        bool m_hasgenregion = false, m_hasexpregion = false;
        bool m_irfinite;

        // m_real_only and m_resonance are initialised as the constructor reads them.
        int m_Nphotons = 0, m_resonance = 0, m_leftfl, m_rightfl, m_real_only = 0;
        double m_ranPhi, m_ranTheta, m_born, m_photonscale, m_u,m_irscale;
        dipoletype::code m_type;
        // m_oldmomenta is the momenta as constructed, never rewritten.
        Vec4D_Vector m_momenta, m_bornmomenta, m_oldmomenta, m_newmomenta, m_ghost, m_eikmomentum,m_mephotons;
        Flavour_Vector m_flavs;
        Vec4D m_photonSum, m_QFrame;
        //! The boost into the dipole rest frame, rebuilt per emission and
        //! cleared at its end (optional so Dipole stays copyable).
        std::optional<Poincare> p_Pboost;
        Poincare p_rotate, p_boost;
        // Radiate, dispatching on m_type: initial-initial samples ISR, final-final
        // samples FSR, initial-final carries the interference and does not radiate.
        // me_acc is the caller's running list of ME photons across dipoles.
        EmissionResult GenerateEmissions(ISR *isr, FSR *fsr,
                                         double born, double v,
                                         Vec4D_Vector &me_acc);

        void PrintInfo();
        void AddPhotonsToDipole(Vec4D_Vector &Photons);
        void AddPhotonToDipole(Vec4D &Photon);
        inline Vec4D Momentum(int i) { return m_momenta[i]; }
        inline void SetMomentum(int i, Vec4D mom) { m_momenta[i] = mom; }
        inline void SetEikMomentum(int i, Vec4D mom) { m_eikmomentum[i] = mom; }
        void AddToGhosts(Vec4D &p);
        double Mass();
        double Eikonal(const Vec4D &k,const Vec4D &p1,const Vec4D &p2);
        double EikonalMassless(const Vec4D &k,const Vec4D &p1,const Vec4D &p2);
        double Eikonal(const Vec4D &k);
        double EikonalInterferance(const Vec4D &k);
        // The EEX beta expansion, truncated at `betaorder` (passed explicitly as
        // the recursion evaluates subtractions at lower orders).
        double EEX(const int betaorder);
        double VirtualEEX(const int betaorder);
        double Hard(const Vec4D & k, int order);
        double Hard(const Vec4D & k1, const Vec4D & k2);
        double Hard(const Vec4D & k1, const Vec4D & k2, const Vec4D & k3);
        double Beta1(const Vec4D & k1, int order);
        double Beta2(const Vec4D & k1, const Vec4D & k2, int order);
        double Beta3(const Vec4D & k1, const Vec4D & k2,const Vec4D & k3, int order);
        double xi(const double &alp, const double &beta, const double &gamma);
        double xi(const double &alp, const double &a1, const double &b1, const double &a2, const double &b2);
        void Boost();
        bool BoostNLO();
        void MakePair(double cms, Vec4D &p1, Vec4D &p2);
        void Boost(Vec4D &p);
        void BoostNLO(Vec4D &p);
        void RandomRotate(Vec4D &p);
        void BoostToQFM(bool boostback);
        void CalculateGamma();
        void Clean();
        bool IsDecayAllowed();
        void SetFlavLab(int i, int j);
        Vec4D Sum();
        inline int size() {return m_momenta.size();}
        //! The legs as they were BEFORE this dipole radiated, which is what
        //! enters the form factor. DipoleSet maintains that invariant: II is
        //! built with the born momenta, FF after ISR and before FSR, IF with
        //! born initial legs and post-ISR final legs.
        inline void LegsBeforeRadiation(Vec4D &p1, Vec4D &p2) const {
          p1 = m_bornmomenta[0];
          p2 = m_bornmomenta[1];
        }
        inline Vec4D GetBornMomenta(int i) { return m_bornmomenta[i];}
        inline Vec4D_Vector GetBornMomenta() { return m_bornmomenta;}
        inline Vec4D GetMomenta(int i) { return m_momenta[i];}
        inline Vec4D_Vector GetMomenta() { return m_momenta;}
        inline Vec4D GetNewMomenta(int i) { return m_newmomenta[i];}
        inline Vec4D_Vector GetNewMomenta() { return m_newmomenta;}
        inline Vec4D GetOldMomenta(int i) { return m_oldmomenta[i];}
        const inline double GetMass(int i) { return m_masses[i];}
        inline void SetMass(int i, double m) { m_masses[i] = m; }
        inline Flavour GetFlav(int i) const { return m_flavs[i];}
        inline Flavour_Vector GetFlavors() const { return m_flavs;}
        inline Vec4D_Vector GetPhotons() const { return m_dipolePhotons;}
        inline Vec4D GetPhotonSum() const { return m_photonSum;}
        inline Vec4D GetGhost(int i) const { return m_ghost[i];}
        inline void SetRotate(Poincare &rot) { p_rotate = rot;}
        inline void SetBoost(Poincare &boost) { p_boost = boost;}
        inline dipoletype::code Type() const {return m_type;}
        inline double Sprime() const {return m_sp;}
        /*!
          The region ledger (Soft_Region.H): the soft region the dipole's
          photons were generated outside of (ISR: beam-CMS sphere of
          IR_CUTOFF/2; FSR: Q-frame sphere of EminQ) and that of its exponent.
        */
        inline void SetGenRegion(const SoftRegion &r) { m_genregion = r; m_hasgenregion = true; }
        inline void SetExpRegion(const SoftRegion &r) { m_expregion = r; m_hasexpregion = true; }
        inline bool HasGenRegion() const { return m_hasgenregion; }
        inline bool HasExpRegion() const { return m_hasexpregion; }
        inline const SoftRegion &GenRegion() const { return m_genregion; }
        inline const SoftRegion &ExpRegion() const { return m_expregion; }
        inline void SetAlpha(double alp) { m_alp=alp; m_alpi=alp/M_PI;}
        inline void SetNPhoton(int N) {m_Nphotons = N; }
        inline void SetBorn(double B) {m_born = B; }
      inline bool IsResonance() const {return m_resonance;}
        inline void SetResonance(bool res) { m_resonance=res;}
        inline double ChargeNorm() const {return -m_QiQj*m_thetaij;}
        inline void SetPhotonScale(double &scale) {m_photonscale = scale;}
        inline void SetMEPhotons(Vec4D_Vector &k) {m_mephotons = k;}
        inline Vec4D_Vector GetMEPhotons() {return m_mephotons;}
        inline bool IsFinite() {return m_irfinite;}
        inline int Left() const {return m_leftfl;}
        inline int Right() const {return m_rightfl;}
        inline bool RealOnly() const {return m_real_only;}
        inline void SetRealOnly(int real) {m_real_only=real;}

    };

    std::ostream& operator<<(std::ostream &out, const Dipole &Dip);
    std::ostream &operator<<(std::ostream &out,const dipoletype::code &c);
}
#endif
