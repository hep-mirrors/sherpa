//bof
//Version: 1 ADICIC++-0.0/2004/06/07    (surely preliminarily)

//Possibility of having information at compile time.
//Globally defined parameter sets belonging to Sudakov_Calculator.H.





namespace ADICIC {



  template<Dipole::Type _DipType>
  struct Sudakov_Info {
    static const Dipole::Type  Dipoletype=Dipole::qqbar;
    static const short         X1power=2;
    static const short         X3power=2;
    static const double        Colourfactor=/*1.5*/0.75*M_PI;
  };

  template<>
  struct Sudakov_Info<Dipole::qg> {
    static const Dipole::Type  Dipoletype=Dipole::qg;
    static const short         X1power=2;
    static const short         X3power=3;
    static const double        Colourfactor=2.0*M_PI/3.0;
  };

  template<>
  struct Sudakov_Info<Dipole::gqbar> {
    static const Dipole::Type  Dipoletype=Dipole::gqbar;
    static const short         X1power=3;
    static const short         X3power=2;
    static const double        Colourfactor=2.0*M_PI/3.0;
  };

  template<>
  struct Sudakov_Info<Dipole::gg> {
    static const Dipole::Type  Dipoletype=Dipole::gg;
    static const short         X1power=3;
    static const short         X3power=3;
    static const double        Colourfactor=2.0*M_PI/3.0;
  };



}





//eof
