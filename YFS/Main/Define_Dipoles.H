#ifndef YFS_Main_Define_Dipoles_H
#define YFS_Main_Define_Dipoles_H

#include "ATOOLS/Phys/Particle.H"
#include "YFS/Tools/Dipole.H"
#include "YFS/Main/YFS_Base.H"
#include "YFS/Main/YFS_Form_Factor.H"


namespace MODEL {
  class Single_Vertex;
  typedef std::vector<Single_Vertex *> Vertex_List;
}



namespace YFS {
  typedef std::vector<YFS::Dipole> Dipole_Vector;

 


  class Define_Dipoles : YFS_Base {
    public:
      bool                    m_success;
      bool                    m_photonsadded;
      int m_N;
      std::map<ATOOLS::Flavour, ATOOLS::Vec4D> m_test_dip, m_born_dip;
      std::map<ATOOLS::Flavour, int> m_flav_label;
      // std::map<ATOOLS::Vec4D, int> m_mom_label;
      ATOOLS::Vec4D_Vector m_chargedinparticles;
      ATOOLS::Vec4D_Vector m_neutralinparticles;
      ATOOLS::Vec4D_Vector m_chargedoutparticles;
      ATOOLS::Vec4D_Vector m_neutraloutparticles;
      ATOOLS::Vec4D_Vector m_softphotons, m_bornmomenta;
      ATOOLS::Vec4D m_softphotonSum;
      Dipole_Vector m_dipoles, m_olddipoles, m_newdipoles, m_dipolesII, m_dipolesFF, m_dipolesIF;
      std::vector<std::vector<int>> m_dip;
      std::map<std::string,MODEL::Vertex_List> m_proc_restab_map;
      YFS::YFS_Form_Factor *p_yfsFormFact;


      YFS::dipoletype     m_type;
      Define_Dipoles();
      ~Define_Dipoles();

      int m_in, m_out;
      std::vector<double> m_massInC, m_massInN;
      std::vector<double> m_massOutC, m_massOutN;

      ATOOLS::Flavour_Vector m_dipoleFlav;
      ATOOLS::Vec4D_Vector m_dipoleMom;
      void AddRadiation();
      void MakeDipolesII(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const &mom, ATOOLS::Vec4D_Vector const &born);
      void MakeDipolesIF(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const &mom, ATOOLS::Vec4D_Vector const &born);
      void MakeDipoles(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const &mom,  ATOOLS::Vec4D_Vector const &born);
      void MakeDipolesFF(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const &mom,ATOOLS::Vec4D_Vector const &born);
      void Dipole_II(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const &mom);
      void Dipole_FF(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const &mom);
      void Dipole_IF(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const &mom, ATOOLS::Vec4D_Vector const &born);
      void CleanInParticles();
      void CleanOutParticles();
      void Get4Mom(ATOOLS::Flavour_Vector const &fl,ATOOLS::Vec4D_Vector mom, std::vector<int> v);
      void Get4Mom(ATOOLS::Flavour_Vector const &fl,ATOOLS::Vec4D_Vector mom);
      bool CheckMasses();
      void IsResonant(YFS::Dipole &D);
      double CalculateRealSub(const Vec4D &k);
      double CalculateRealSubEEX(const Vec4D &k);
      double CalculateVirtualSub();
      double CalculateRealVirtualSub(const Vec4D & k);
      double CalculateEEXReal(const Vec4D & k);
      double CalculateEEXVirtual();
      double CalculateEEX(const Vec4D & k1, const Vec4D & k2);
      void generate_pairings(std::vector<std::vector<int>>& pairings, std::vector<int>& curr_pairing, std::vector<int>& available_nums);
      // inline void SetOldMomentum(int i, ATOOLS::Vec4D mom) { m_olddipoles.m_momenta[i] = mom;}
      // inline void SetNewMomentum(int i, ATOOLS::Vec4D mom) { m_newdipoles.m_momenta[i] = mom;}
      // inline ATOOLS::Vec4D OldMomentum(int i) { return m_olddipoles.m_momenta[i]; } 
      // inline ATOOLS::Vec4D NewMomentum(int i) { return m_newdipoles.m_momenta[i]; } 
      void CleanUp();
      std::ostream& operator<<(std::ostream &out);
      inline bool DoneSuccessfully()     { return m_success; }
      inline bool AddedAnything()        { return m_photonsadded; }
      inline void AddSoftPhotons(ATOOLS::Vec4D k) { m_softphotonSum +=k; }
      inline Dipole_Vector Old() { return m_olddipoles; }
      inline Dipole_Vector New() { return m_newdipoles; }
      inline Dipole * GetDipoleII() { return &m_dipolesII[0]; }
      inline Dipole_Vector * GetDipoleFF() { return &m_dipolesFF; }
      inline Dipole_Vector * GetDipoleIF() { return &m_dipolesIF; }
      inline void SetProcResMap(std::map<std::string,MODEL::Vertex_List> map) {m_proc_restab_map=map;}
  };


}
#endif
