//bof
//Version: 3 ADICIC++-0.0/2005/09/30    (surely preliminarily)

//The Sudakov classes are provided, which manage the Sudakov calculation needed
//for emissions of II dipoles.
//                 ==========



#ifndef _IISudakov_Group_H_
#define _IISudakov_Group_H_ _IISudakov_Group_H_


#include <string>
#include "Sudakov_Calculator.H"





namespace ADICIC {



  template<Dipole::Type _Diptype>
  class IISudakov_Group : public Sudakov_Calculator {

  private:

    //Blocked methods:

    IISudakov_Group(const IISudakov_Group<_Diptype>&);
    IISudakov_Group<_Diptype>& operator=(const IISudakov_Group<_Diptype>&);
    const bool operator==(const IISudakov_Group<_Diptype>&) const;

  private:

    const Radiation::Type m_radtype;

    double m_sdip;
    double m_a;
    double m_s;
    double m_x2tmin;
    double m_x2tmax;

    double m_x2t;
    double m_ymax;
    double m_rap;
    double m_corr;

    std::list<Sudakov_Base*> l_sud;    //Put the Sudakov's together.

    //Private methods:

    inline const bool TestEfracs(const double, const double) const;
    inline       void Reset();

    const bool InitRadiation() const;
    const bool InitWithCurrentDipole();

  public:

    IISudakov_Group(const Radiation::Type=Radiation::g);
    ~IISudakov_Group();

    //Public (inline) methods:

    inline void Which() const;
    inline void ShowSpecification() const;

    inline const double A() const;
    inline const double Scale() const;
    inline const double X2tmin() const;
    inline const double X2t() const;
    inline const double Ymax() const;

    const bool GenerateVariablesFor(const Dipole&, Sudakov_Result&);

  };



  //###########################################################################



  //================
  //Blocked classes:
  //================

  template<Dipole::Type _Diptype, Radiation::Group _Radgroup>
  class IISudakov : public Sudakov_Base {
    IISudakov();
    void ShowSpecification() const;
    const double CalculateRapLimit();
    const double GenerateRap();
    const double GenerateCorr(const Multidouble&);
    const double GenerateX2t();
  };
  template<>
  class IISudakov<Dipole::iigg,Radiation::igluon> : public Sudakov_Base {
    IISudakov();
    void ShowSpecification() const;
    const double CalculateRapLimit();
    const double GenerateRap();
    const double GenerateCorr(const Multidouble&);
    const double GenerateX2t();
  };
  template<>
  class IISudakov<Dipole::iiqbarq,Radiation::quark> : public Sudakov_Base {
    IISudakov();
    void ShowSpecification() const;
    const double CalculateRapLimit();
    const double GenerateRap();
    const double GenerateCorr(const Multidouble&);
    const double GenerateX2t();
  };





  //===============================================
  //Specified class: gluon emission (Annihilation.)
  //===============================================

  template<Dipole::Type _Diptype>
  class IISudakov<_Diptype,Radiation::gluon> : public Sudakov_Base {

  private:

    IISudakov();
    IISudakov(const IISudakov<_Diptype,Radiation::gluon>&);
    IISudakov<_Diptype,Radiation::gluon>&
    operator=(const IISudakov<_Diptype,Radiation::gluon>&);
    const bool operator==(const IISudakov<_Diptype,Radiation::gluon>&) const;

  private:

    //Compile time information.
    static const short  s_x1pow;
    static const short  s_x3pow;
    static const double s_colfac;
    static const double s_iieffbas;
    static const double s_pdfapprox;

    double                           m_genx2tfac;
    Multiflavour                     m_mfl;
    const IISudakov_Group<_Diptype>& m_sgroup;

  public:

    IISudakov(const IISudakov_Group<_Diptype>&, const Sudakov_Flavour&);
    ~IISudakov();

    inline void Which() const;
    inline void ShowSpecification() const;

    inline void   SetPDFFlavs(bool);
    inline void   SetGenX2tFac();
    inline double GeneratePDFCorr(const Multidouble&);

    inline const double CalculateRapLimit();
    inline const double GenerateRap();
    inline const double GenerateCorr(const Multidouble&);
    const double GenerateX2t();

  };





  //===========================================
  //Specified class: igluon emission (Compton.)
  //===========================================

  template<Dipole::Type _Diptype>
  class IISudakov<_Diptype,Radiation::igluon> : public Sudakov_Base {

  private:

    IISudakov();
    IISudakov(const IISudakov<_Diptype,Radiation::igluon>&);
    IISudakov<_Diptype,Radiation::igluon>&
    operator=(const IISudakov<_Diptype,Radiation::igluon>&);
    const bool operator==(const IISudakov<_Diptype,Radiation::igluon>&) const;

  private:

    //Compile time information.
    static const double s_colfac;
    static const double s_iieffbas;
    static const double s_pdfapprox;

    double                           m_genx2tfac;
    Multiflavour                     m_mfl;
    const IISudakov_Group<_Diptype>& m_sgroup;

  public:

    IISudakov(const IISudakov_Group<_Diptype>&, const Sudakov_Flavour&);
    ~IISudakov();

    inline void Which() const;
    inline void ShowSpecification() const;

    inline bool   NfVeto(const double) const;
    inline void   SetPDFFlavs(bool);
    inline void   SetGenX2tFac();
    inline double GeneratePDFCorr(const Multidouble&);

    inline const double CalculateRapLimit();
    inline const double GenerateRap();
    inline const double GenerateCorr(const Multidouble&);
    const double GenerateX2t();

  };





  //===========================================
  //Specified class: quark-antiquark splitting.
  //===========================================

  template<Dipole::Type _Diptype>
  class IISudakov<_Diptype,Radiation::quark> : public Sudakov_Base {

  private:

    IISudakov();
    IISudakov(const IISudakov<_Diptype,Radiation::quark>&);
    IISudakov<_Diptype,Radiation::quark>&
    operator=(const IISudakov<_Diptype,Radiation::quark>&);
    const bool operator==(const IISudakov<_Diptype,Radiation::quark>&) const;

  private:

    //Compile time information.
    static const double s_colfac;
    static const double s_iieffbas;
    static const double s_pdfapprox;

    double                           m_genx2tfac;
    const Sudakov_Flavour            m_inicode;
    std::vector<Multiflavour>        v_mfl;
    const IISudakov_Group<_Diptype>& m_sgroup;

  public:

    IISudakov(const IISudakov_Group<_Diptype>&, const Sudakov_Flavour&);
    ~IISudakov();

    inline void Which() const;
    inline void ShowSpecification() const;

    inline bool   NfVeto(const double) const;
    inline void   SetPDFFlavs(bool);
    inline void   SetGenX2tFac();
    inline double GeneratePDFCorr(const Multidouble&);

    inline const double CalculateRapLimit();
    inline const double GenerateRap();
    inline const double GenerateCorr(const Multidouble&);
    const double GenerateX2t();

  };



  //###########################################################################



  class IISudakov_Stats {
  private:
    bool         out;
    unsigned     incs;
    unsigned     bviols;
    unsigned     tviols;
    double       blim;
    double       tlim;
    double       min;
    double       max;
    Dipole::Type type;
    std::string  name;
  public:
    IISudakov_Stats(Dipole::Type, const double&,
		    const std::string& n="NoName");
    IISudakov_Stats(Dipole::Type, const std::string& n="NoName",
		    bool o=true,
		    const double& bot=0.0, const double& top=1.0);
    ~IISudakov_Stats();
    const bool Include(const double& weight);
  };



}    //eo namespace ADICIC





#include "IISudakov_Group.inl.hh"


#endif    //eo _IISudakov_Group_H_



//eof
