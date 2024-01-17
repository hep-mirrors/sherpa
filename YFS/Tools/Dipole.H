#ifndef YFS_Main_Dipole_H
#define YFS_Main_Dipole_H

#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Math/Poincare.H"
// #include "YFS/Main/YFS_Handler.H"

using namespace ATOOLS;

namespace YFS{

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
        Dipole(Flavour_Vector const &fl, Vec4D_Vector const &mom,  Vec4D_Vector const &born, dipoletype::code d);
        ~Dipole();

        std::vector<double> m_masses, m_charges;
        Vec4D_Vector m_dipolePhotons;
        std::vector<std::string> m_names;
        double m_QiQj, m_Qi, m_Qj, m_thetai, m_thetaj, m_thetaij, m_sp, m_b1, m_b2, m_gamma, m_gammap, m_alp, m_alpi,m_rescale_alpha;
        bool m_setBorn, m_setbeams=1,m_sameflav;
        int m_Nphotons, m_resonance;
        dipoletype::code m_type;
        Vec4D_Vector m_momenta, m_bornmomenta, m_oldmomenta, m_newmomenta, m_ghost, m_beams, m_eikmomentum;
        Flavour_Vector m_flavs;
        Vec4D m_photonSum, m_QFrame;
        Poincare *p_boost, *p_rotate, *p_Pboost, *p_Qboost;
        void PrintInfo();
        void AddPhotonsToDipole(Vec4D_Vector &Photons);
        void AddPhotonsToDipole(Vec4D &Photon);
        inline void ClearPhotons() { m_dipolePhotons.clear(); };
        inline Vec4D Momentum(int i) { return m_momenta[i]; }
        inline void SetMomentum(int i, Vec4D mom) { m_momenta[i] = mom; }
        inline void SetNewMomentum(int i, Vec4D mom) { m_newmomenta[i] = mom; }
        inline void SetEikMomentum(int i, Vec4D mom) { m_eikmomentum[i] = mom; }
        inline void SetBeams(int i, Vec4D mom) { m_beams[i] = mom; }
        void AddToGhosts(Vec4D &p);
        void SetBornMomentum(Vec4D_Vector k);
        double Mass();
        double Eikonal(Vec4D k, Vec4D p1, Vec4D p2);
        double EikonalMassless(Vec4D k, Vec4D p1, Vec4D p2);
        double Eikonal(Vec4D k);
        double EEX(const Vec4D & k);
        void Boost();
        void Boost(Vec4D &p);
        void RandomRotate(Vec4D &p);
        void BoostToQFM(bool boostback);
        void CalculateBeta();
        void CalculateGamma();
        void Clean();
        bool IsDecayAllowed();
        Vec4D Sum();
        inline int size() {return m_momenta.size();}
        inline Vec4D GetBornMomenta(int i) { return m_bornmomenta[i];}
        inline Vec4D GetMomenta(int i) { return m_momenta[i];}
        inline Vec4D GetNewMomenta(int i) { return m_newmomenta[i];}
        inline Vec4D GetOldMomenta(int i) { return m_oldmomenta[i];}
        inline Vec4D GetEikMomenta(int i) { return m_eikmomentum[i];}
        inline Flavour GetFlav(int i) { return m_flavs[i];}
        inline Vec4D_Vector GetPhotons() { return m_dipolePhotons;}
        inline Vec4D GetPhotonSum() { return m_photonSum;}
        inline Vec4D GetGhost(int i) { return m_ghost[i];}
        inline dipoletype::code Type() {return m_type;}
        inline void AddPhotonsSumToDipole(Vec4D k) {m_photonSum = k;}
        inline double Sprime() {return m_sp;}
        inline void SetAlpha(double alp) { m_alp=alp; m_alpi=alp/M_PI;}
        inline void SetNPhoton(int N) {m_Nphotons = N; }
        inline bool IsResonance() {return m_resonance;}
        inline void SetResonance(bool res) { m_resonance=res;}
        inline bool IsSame() {return m_sameflav;}

    };

    std::ostream& operator<<(std::ostream &out, const Dipole &Dip);
    std::ostream &operator<<(std::ostream &out,const dipoletype::code &c);
}
#endif