//bof
//Version: 3 ADICIC++-0.0/2005/07/22

//Possibility of having information at compile time.
//Globally defined parameters influencing the evolution of chains and cascades.



#ifndef _Evolution_Strategy_hpp_
#define _Evolution_Strategy_hpp_ _Evolution_Strategy_hpp_


#include <cstdlib>





namespace ADICIC {



  namespace Chain_Evolution_Strategy {

    enum Type {
      Unknown    = 0,
      Production = 1,
      Emission   = 11,
      Mass       = 21,
      stop       = -77
    };

    static const size_t NumberOfTypes = 5;

    static const Type List[NumberOfTypes]={
      Unknown,
      Production,
      Emission,
      Mass,
      stop    //Keep this at the bottom.
    };

  }    //eo namespace Chain_Evolution_Strategy





  //Enhance readability -> chain evolution labels. Do not change numbers.
  namespace cel {
    enum code {
      def  = 0,
      stop = 1
    };
  }



}    //eo namespace ADICIC





#endif    //eo _Evolution_Strategy_hpp_



//eof
