//bof
//Version: 1 ADICIC++-0.0/2004/07/12    (surely preliminarily)

//Possibility of having information at compile time.
//Globally defined parameter sets belonging to Sudakov_Calculator.H.





namespace ADICIC {



  template<Dipole::Type _DipType>
  struct Sudakov_Info {
    static const Dipole::Type  Dipoletype;
    static const short         X1power;
    static const short         X3power;
    static const double        Colourfactor;
  };



}





//eof
