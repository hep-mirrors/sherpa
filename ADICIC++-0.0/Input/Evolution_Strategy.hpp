//bof
//Version: 4 ADICIC++-0.0/2006/06/01

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





  //Enhance readability! Do not change numbers.
  //Dipole shower mode dsm (j i.e. just).
  namespace dsm {
    enum code {
      off  = 0,
      jff  = 1,
      jfi  = 2,
      jif  = 4,
      jii  = 8,
      iiff = 9,
      all  = 15
    };
  }


  //Enhance readability! Do not change numbers.
  //Factorization scale type fascat.
  namespace fascat {
    enum code {
      p2t  = 0,    //the evolution variable p2t.
      k2t  = 1,    //the lab squared transverse momentum.
      m2t  = 2,    //the lab squared transverse mass.
      shat = 3,    //for testing purposes, e.g. single-emission test.
      stop = 4
    };
  }



}    //eo namespace ADICIC





#endif    //eo _Evolution_Strategy_hpp_



//eof
