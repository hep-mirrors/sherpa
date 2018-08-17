#include "EXTAMP/CS_Dipole.H"

namespace EXTAMP {

  ///////////////////////////////////////////////////////////////
  ////////// FINAL FINAL ////////////////////////////////////////
  ///////////////////////////////////////////////////////////////

  
  /* Kinematic variables of hep-ph/9605323v3 eq. (5.3) - (5.6) */
  struct FF_Dipole_Kinematics : public Dipole_Kinematics {
    ATOOLS::Vec4D m_pi;
    ATOOLS::Vec4D m_pj;
    ATOOLS::Vec4D m_pk;
    ATOOLS::Vec4D m_pl;
    ATOOLS::Vec4D m_n;         // for id case
    ATOOLS::Vec4D m_pa;         // for id case
    ATOOLS::Vec4D m_pij_tilde;
    ATOOLS::Vec4D m_pk_tilde;
    ATOOLS::Vec4D m_pl_tilde;  // for res case
    ATOOLS::Vec4D m_n_tilde;   // for id case   //TODO: discard if not needed
    double m_zi;
    double m_zj;
    double m_zain;   // for id case
    double m_y;
    double m_viab;   // for id case
    double m_xain;   // for id case
    double m_vitilde;// for id case
    double m_alphamin; 

    double Alpha() const { return m_alphamin; }

    double ShowerX()  const { return m_zi; }
    double ShowerY()  const { return m_y ; }
    double ShowerQ2() const { return (m_pi+m_pj+m_pk).Abs2(); }
  };


  class FF_Dipole : virtual public CS_Dipole {

  public:

    FF_Dipole(const Dipole_Info& di) : CS_Dipole(di) {};
    virtual ~FF_Dipole() {};

    void  CalcKinematics(const ATOOLS::Vec4D_Vector& p);
    const ATOOLS::Vec4D_Vector& Momenta() const { return m_kin.m_born_mom; }

  protected:

    double CalcB() const;
    double CalcA() const;
    double CalcKinDependentPrefac() const;
    ATOOLS::Vec4D CalcPtilde() const;

    const Dipole_Kinematics* const LastKinematics() const { return &m_kin; }

    /* Legacy structure: NLO_subevts only maintain pointers to const
       momenta and flavours. These pointers need to point to
       persistent copies of momenta/flavours that do not run out of
       scope. Therefore have to maintain a copy of the last calculated
       momentum configuration in this class and let the momentum
       pointer of one NLO_subevt point to that instance. */
    FF_Dipole_Kinematics m_kin;

  };


  ///////////////////////////////////////////////////////////////
  ////////// FINAL INITIAL //////////////////////////////////////
  ///////////////////////////////////////////////////////////////


  /* Kinematic variables of hep-ph/9605323v3 eq. (5.3) - (5.6) */
  struct FI_Dipole_Kinematics : public Dipole_Kinematics {
    ATOOLS::Vec4D m_pi;
    ATOOLS::Vec4D m_pj;
    ATOOLS::Vec4D m_pa;
    ATOOLS::Vec4D m_pij_tilde;
    ATOOLS::Vec4D m_pa_tilde;
    ATOOLS::Vec4D m_pb;        // for id case
    ATOOLS::Vec4D m_n;         // for id case
    double m_x;
    double m_zi;
    double m_zj;
    double m_zain;   // for id case
    double m_viab;   // for id case
    double m_xain;   // for id case
    double m_vitilde;// for id case
    double m_alphamin; 

    double Alpha() const { return m_alphamin; }

    double ShowerX()  const { return m_zi; }
    double ShowerY()  const { return m_x ; }
    double ShowerQ2() const { return (m_pi+m_pj-m_pa).Abs2(); }
  };


  class FI_Dipole : virtual public CS_Dipole {

  public:

    FI_Dipole(const Dipole_Info& di) : CS_Dipole(di) {};
    virtual ~FI_Dipole() {};

    void  CalcKinematics(const ATOOLS::Vec4D_Vector& p);
    const ATOOLS::Vec4D_Vector& Momenta() const { return m_kin.m_born_mom; }

  protected:

    double CalcB() const;
    double CalcA() const;
    double CalcKinDependentPrefac() const;
    ATOOLS::Vec4D CalcPtilde() const;

    const Dipole_Kinematics* const LastKinematics() const { return &m_kin; }

    FI_Dipole_Kinematics m_kin;

  };


  ///////////////////////////////////////////////////////////////
  ////////// INITIAL FINAL //////////////////////////////////////
  ///////////////////////////////////////////////////////////////
  

  /* Kinematic variables of hep-ph/9605323v3 eq. (5.3) - (5.6) */
  struct IF_Dipole_Kinematics : public Dipole_Kinematics {
    ATOOLS::Vec4D m_pa;
    ATOOLS::Vec4D m_pi;
    ATOOLS::Vec4D m_pk;
    ATOOLS::Vec4D m_pai_tilde;
    ATOOLS::Vec4D m_pk_tilde;
    ATOOLS::Vec4D m_n;         // for id case
    double m_ui;
    double m_x;
    double m_xain;             // for id case
    double m_viab;             // for id case
    double m_vitilde;          // for id case
    double m_alphamin;

    double Alpha() const { return m_alphamin; }

    double ShowerX()  const { return m_x; }
    double ShowerY()  const { return m_ui; }
    double ShowerQ2() const { return (-m_pa+m_pi+m_pk).Abs2(); }
  };


  class IF_Dipole : virtual public CS_Dipole {

  public:

    IF_Dipole(const Dipole_Info& di) : CS_Dipole(di) {};
    virtual ~IF_Dipole() {};

    void  CalcKinematics(const ATOOLS::Vec4D_Vector& p);
    const ATOOLS::Vec4D_Vector& Momenta() const { return m_kin.m_born_mom; }

  protected:

    double CalcB() const;
    double CalcA() const;
    double CalcKinDependentPrefac() const;
    ATOOLS::Vec4D CalcPtilde() const;

    const Dipole_Kinematics* const LastKinematics() const { return &m_kin; }

    IF_Dipole_Kinematics m_kin;

  };


  ///////////////////////////////////////////////////////////////
  ////////// INITIAL INITIAL ////////////////////////////////////
  ///////////////////////////////////////////////////////////////
  

  /* Kinematic variables of hep-ph/9605323v3 eq. (5.3) - (5.6) */
  struct II_Dipole_Kinematics : public Dipole_Kinematics {
    ATOOLS::Vec4D m_pa;
    ATOOLS::Vec4D m_pi;
    ATOOLS::Vec4D m_pb;
    ATOOLS::Vec4D m_pai_tilde;
    ATOOLS::Vec4D m_pb_tilde;
    double m_ui;
    double m_x;
    double m_v;

    double Alpha() const { return m_v; }

    double ShowerX()  const { return m_x; }
    double ShowerY()  const { return m_v; }
    double ShowerQ2() const { return (-m_pa+m_pi-m_pb).Abs2(); }
  };


  class II_Dipole : virtual public CS_Dipole {

  public:

    II_Dipole(const Dipole_Info& di) : CS_Dipole(di) {};
    virtual ~II_Dipole() {};

    void  CalcKinematics(const ATOOLS::Vec4D_Vector& p);
    const ATOOLS::Vec4D_Vector& Momenta() const { return m_kin.m_born_mom; }

  protected:

    double CalcB() const;
    double CalcA() const;
    double CalcKinDependentPrefac() const;
    ATOOLS::Vec4D CalcPtilde() const;

    const Dipole_Kinematics* const LastKinematics() const { return &m_kin; }
    
    II_Dipole_Kinematics m_kin;

  };

}
