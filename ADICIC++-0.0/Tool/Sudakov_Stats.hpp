//bof
//Version: 4 ADICIC++-0.0/2006/07/30    (surely preliminarily)

//Helpful tool to crosscheck the veto algorithm in the Sudakov calculations.



#ifndef _Sudakov_Stats_hpp_
#define _Sudakov_Stats_hpp_ _Sudakov_Stats_hpp_


#include <string>
#include "Sudakov_Calculator.H"





namespace ADICIC {



  class Sudakov_Stats {
    static const bool        s_chk;
    static const std::size_t s_mode;
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
    std::vector<Multidouble> entries;
    std::vector<Multidouble> wgtsum;
  private:
    template<std::size_t n>
    double CalcEstimate(const std::size_t dx, const double& resort);
  public:
    Sudakov_Stats(Dipole::Type, const double&, const std::string& n="NoName");
    Sudakov_Stats(Dipole::Type, const std::string& n="NoName", bool out=true,
		  const double& bot=0.0, const double& top=1.0,
		  double exa=-1.0);
    ~Sudakov_Stats();
    const bool Include(const double& weight);
    const bool Include(const double& weight, const double& estim,
		       const double& x);
    double GiveEstimate(const double& x, const double& resort);
  };



}    //eo namespace ADICIC





//#include "....inl.hh"


#endif    //eo _Sudakov_Stats_hpp_



//eof
