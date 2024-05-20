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
        Dipole(Flavour_Vector const &fl, Vec4D_Vector const &mom,  
                Vec4D_Vector const &born, dipoletype::code d, const double alpha);
        ~Dipole();

        std::vector<double> m_masses, m_charges, m_theta, m_Q;
        Vec4D_Vector m_dipolePhotons, m_dipolePhotonsEEX;
        std::vector<std::string> m_names;
        double m_QiQj, m_Qi, m_Qj, m_thetai, m_thetaj, m_thetaij, m_sp, m_b1, m_b2, m_gamma, m_gammap, m_alp, m_alpi,m_rescale_alpha;
        bool m_setBorn, m_setbeams=1,m_sameflav;
        int m_Nphotons, m_resonance, m_betaorder;
        double m_ranPhi, m_ranTheta, m_born;
        dipoletype::code m_type;
        Vec4D_Vector m_momenta, m_bornmomenta, m_oldmomenta, m_newmomenta, m_ghost, m_beams, m_eikmomentum;
        Flavour_Vector m_flavs;
        Vec4D m_photonSum, m_QFrame;
        Poincare *p_boost, *p_rotate, *p_rotatex, *p_rotatey, *p_Pboost, *p_Qboost;
        void PrintInfo();
        void AddPhotonsToDipole(Vec4D_Vector &Photons);
        void AddPhotonToDipole(Vec4D &Photon);
        inline void ClearPhotons() { m_dipolePhotons.clear(); };
        inline Vec4D Momentum(int i) { return m_momenta[i]; }
        inline void SetMomentum(int i, Vec4D mom) { m_momenta[i] = mom; }
        inline void SetNewMomentum(int i, Vec4D mom) { m_newmomenta[i] = mom; }
        inline void SetEikMomentum(int i, Vec4D mom) { m_eikmomentum[i] = mom; }
        inline void SetBornMomentum(int i, Vec4D mom) {m_bornmomenta[i] = mom; }
        inline void SetBeams(int i, Vec4D mom) { m_beams[i] = mom; }
        void AddToGhosts(Vec4D &p);
        double Mass();
        double Eikonal(const Vec4D &k,const Vec4D &p1,const Vec4D &p2);
        double EikonalMassless(const Vec4D &k,const Vec4D &p1,const Vec4D &p2);
        double Eikonal(const Vec4D &k);
        double EikonalInterferance(const Vec4D &k);
        double EEX(const int betaorder);
        double VirtualEEX(const int betaorder);
        double Hard(const Vec4D & k);
        double Hard(const Vec4D & k1, const Vec4D & k2);
        double Hard(const Vec4D & k1, const Vec4D & k2, const Vec4D & k3);
        double Beta1(const Vec4D & k1);
        double Beta2(const Vec4D & k1, const Vec4D & k2);
        double Beta3(const Vec4D & k1, const Vec4D & k2,const Vec4D & k3);
        double xi(const double &alp, const double &beta, const double &gamma);
        double xi(const double &alp, const double &a1, const double &b1, const double &a2, const double &b2);
        void Boost();
        void BoostLab();
        void Boost(Vec4D &p);
        void RandomRotate(Vec4D &p);
        void BoostToQFM(bool boostback);
        void BoostToCMS(Vec4D_Vector &k, bool boostback=false);
        void CalculateBeta();
        void CalculateGamma();
        void Clean();
        bool IsDecayAllowed();
        Vec4D Sum();
        inline int size() {return m_momenta.size();}
        inline Vec4D GetBornMomenta(int i) { return m_bornmomenta[i];}
        inline Vec4D_Vector GetBornMomenta() { return m_bornmomenta;}
        inline Vec4D GetMomenta(int i) { return m_momenta[i];}
        inline Vec4D_Vector GetMomenta() { return m_momenta;}
        inline Vec4D GetNewMomenta(int i) { return m_newmomenta[i];}
        inline Vec4D_Vector GetNewMomenta() { return m_newmomenta;}
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
        inline void SetBorn(double B) {m_born = B; }
        inline bool IsResonance() {return m_resonance;}
        inline void SetResonance(bool res) { m_resonance=res;}
        inline bool IsSame() {return m_sameflav;}
        inline double ChargeNorm() {return -m_QiQj*m_thetaij;}

    };

    std::ostream& operator<<(std::ostream &out, const Dipole &Dip);
    std::ostream &operator<<(std::ostream &out,const dipoletype::code &c);
}
#endif
