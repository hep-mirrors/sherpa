//bof
//Version: 2 ADICIC++-0.0/2004/08/27

//Possibility of having information at compile time.
//Globally defined parameters influencing Chain_Handler.H.



#ifndef _Sudakov_Strategy_hpp_
#define _Sudakov_Strategy_hpp_ _Sudakov_Strategy_hpp_





namespace ADICIC {



  namespace Radiation {
    enum Group {
      quark    = 50,
      qtop     = 51,
      qbot     = 52,
      gluon    = 100,
      incorrect=-9999
    };
    enum Type {
      d      = 1,
      du     = 2,
      dus    = 3,
      dusc   = 4,
      duscb  = 5,
      g      = 10,
      gd     = 11,
      gdu    = 12,
      gdus   = 13,
      gdusc  = 14,
      gduscb = 15
    };
  }





  namespace Sudakov_Strategy {



    //Set the corresponding value:
    static const int Label=1;



    struct Unknown {};          //All other Labels.
    struct Factorization {};    //Label==1.    Factorization of Sudakov's.
    struct Distribution {};     //Label==2.    Distribution of Sudakov's.



  }    //eo namespace Sudakov_Strategy



}    //eo namespace ADICIC





#endif    //eo _Sudakov_Strategy_hpp_



//eof
