#ifndef Color_Function_h
#define Color_Function_h

namespace AMEGIC {

  class cf {
  public:
    enum code {T=0, F=1, D=2, None=3,G=4, Unknown=99}; 
  };
  
  class Color_Function {
  public:
    Color_Function()  {
      type = cf::Unknown;
      Next = 0;
      for (short int i=0;i<3;++i) partarg[i] = -1;
      for (short int i=0;i<3;++i) strarg[i] = '?';
    }
    Color_Function(cf::code _type, 
		   int _partarg0 = -1, int _partarg1 = -1, int _partarg2 = -1, 
		   char _strarg0 = '?', char _strarg1 = '?', char _strarg2 = '?') 
      :type(_type) { 
      
      Next = 0;
      partarg[0] = _partarg0;
      partarg[1] = _partarg1;
      partarg[2] = _partarg2;
      strarg[0]  = _strarg0;
      strarg[1]  = _strarg1;
      strarg[2]  = _strarg2;
      }
    
    cf::code  type; // Number of args, multi gluons etc.
    int  partarg[3];
    char strarg[3];
    Color_Function* Next;

    void SetParticleArg(int a, int b,int c=-1) {
      partarg[0] = a;partarg[1] = b;
      if (type!=cf::D && type!=cf::G) partarg[2] = c;
    }

    void SetStringArg(char a, char b, char c='?') {
      strarg[0] = a;strarg[1] = b;
      if (type!=cf::D && type!=cf::G) strarg[2] = c;
    }

    std::string String() {
      
      if (type==cf::None) return std::string("1");
      
      std::string help;
      switch (type) {
      case 0: help = std::string("T[ , , ]");break;
      case 1: help = std::string("F[ , , ]");break;
      case 2: help = std::string("D[ , ]");break;
      case 4: help = std::string("G[ , ]");break;
      }
      for (short int i=0;i<3;i++) {
	if ((type==cf::D || type==cf::G) && i==2) break;
	help[2+i*2] = strarg[i];
      }
      return help;
    }
   
  };  

  // Format of a Color_Function to transport via MPI:
  struct MPI_Color_Function {
    
    //    cf::code m_type;
    int    m_type;
    int    m_partarg[3];         
    char   m_strarg[4];         
  };
  
  //for MPI usage
  
  void Color_Function2MPI(const Color_Function * , MPI_Color_Function &);
  
  Color_Function * MPI2Color_Function(const MPI_Color_Function &);

  std::ostream & operator<<(std::ostream &, const MPI_Color_Function &);

}
#endif










