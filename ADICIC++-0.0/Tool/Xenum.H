//bof
//Version: 1 ADICIC++-0.0/2004/03/11

//A collection of specific enum structures.



#ifndef _Xenum_H_
#define _Xenum_H_ _Xenum_H_


//#include





//namespace ... {



//===================================================================
//  Enumerators for an extended bool, a Signed bool and a Bool type
//===================================================================


enum xbool {
  negative=-1,    nil=0,       positive=1,
  minus=-1,                    plus=1,
  onn=-1,         off=0,       on=1,
  end=-1,         both=0,      front=1,
                  between=0
};



enum Sbool {
  Minus=-1,    Plus=1,
  Left=-1,     Right=1
};



enum Bool {
  Blocked=-1,    Off=0,    On=1,
  Dead=-1,       Old=0,    New=1
};










//======================================================================
//  The enumerator types Duo and Trio - two and three stage bool types
//======================================================================


enum Duo {
  Down=-1,
  Up=1
};



enum Trio {
  Negative=-1,
  Nil=0,
  Positive=1
};










//===================================================
//  The booln enumerator types - n stage bool types
//===================================================


namespace bool3 {
  enum level {
    zero=0,
    one=1,
    two=2
  };
}



namespace bool4 {
  enum level {
    zero=0,
    one=1,
    two=2,
    three=3
  };
}



namespace bool5 {
  enum level {
    zero=0,
    one=1,
    two=2,
    three=3,
    four=4
  };
}





#endif    //eo _Xenum_H_



//eof










