//bof
//Version: 2 ADICIC++-0.0/2004/08/09

//Possibility of having information at compile time.
//Globally defined parameters influencing Chain_Handler.H.



#ifndef _Evolution_Strategy_hpp_
#define _Evolution_Strategy_hpp_ _Evolution_Strategy_hpp_





namespace ADICIC {



  namespace Chain_Evolution_Strategy {



    //Set the corresponding value:
    static const int Label=1;



    struct Unknown {};       //All other Labels.
    struct Production {};    //Label==1.
    struct Emission {};      //Label==2.
    struct Mass {};          //Label==3.



    template<int label> struct Map {
      typedef Unknown Ret;
    };
    template<> struct Map<1> {
      typedef Production Ret;
    };
    template<> struct Map<2> {
      typedef Emission Ret;
    };
    template<> struct Map<3> {
      typedef Mass Ret;
    };



    //typedef Map<Label>::Ret Ret;



  };    //eo namespace Chain_Evolution_Strategy



}    //eo namespace ADICIC





#endif    //eo _Evolution_Strategy_hpp_



//eof
