// defines Random-Class

#ifndef random_h
#define random_h

#include <fstream>
#include "MathTools.H"

namespace AMATOOLS {
  class Random {
    long id;
    int inext,inextp;
    long ma[56];
    
    double ran2(long *idum);
    double ran1(long *idum);
    
    double Ran3();
    void   Init_Ran3(long*);
    // for read and write status
    long int written;
    std::fstream * myoutstream;
    char outname[150];
  public:
    Random(long nid){Set_Seed(nid);}
    
    inline double get() {return Ran3();}
    inline double getNZ() {      // not zero random number
      double ran1;
      do {ran1=get();
      } while (ran1==0.0); 
      return ran1;
    }
    inline long Get_Seed() {return id;}
    inline void Set_Seed(long nid) {
      id = nid<0 ? nid : -nid;
      Init_Ran3(&id);
      // for status recovery
      written=0;    
      myoutstream=0;
      strcpy(outname,"");
    }
    
    inline double theta() {
      return acos(2.*get()-1.);
    }  
    int WriteOutStatus(char*);
    void ReadInStatus(char* filename, long int index=0 );
    ~Random(){ 
      if (myoutstream!=0) {
	myoutstream->close();
	delete myoutstream;
      }
    }; 
  };

  extern Random Ran;
}


#endif
