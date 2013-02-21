//%module Cluster_Amplitude
%{
#include <ATOOLS/Phys/Cluster_Leg.H>
#include <ATOOLS/Phys/Decay_Info.H>
#include <ATOOLS/Phys/Cluster_Amplitude.H>
#include <vector>
#include <map> 
%}


// Tell SWIG about size_t
typedef unsigned int size_t;

namespace ATOOLS {

  // ColorID is not part of the Cluster_Amplitude header but is included in this file for convenience
  struct ColorID {
    int m_i, m_j;
    inline ColorID(const int &i=-1,const int &j=-1): m_i(i), m_j(j) {}
    inline ColorID Conj() { return ColorID(m_j,m_i); }
    inline bool operator==(const ColorID &c) const 
    { return c.m_i==m_i && c.m_j==m_j; }
    inline bool Singlet() { return m_i==m_j; }
  };// end of class ColorID

  
  typedef std::map<int,int> CI_Map;

  class Cluster_Amplitude {
  protected:

    Cluster_Amplitude(Cluster_Amplitude *const prev=NULL);

    ~Cluster_Amplitude();

  public:

    static Cluster_Amplitude* New(Cluster_Amplitude *const prev=NULL);

    void CreateLeg(const ATOOLS::Vec4D &p,const ATOOLS::Flavour &fl,
		   const ATOOLS::ColorID &col=ColorID(),
		   const size_t &id=std::string::npos);
    
    inline void SetNIn(const size_t &nin)  { m_nin=nin; }

    // The default CreateLeg method takes a Vec4D as the first argument. Vec4D is a Vec4<double> typedef
    // that is not available in python and therefore, the class is extended by a function that simply
    // casts a Vec4<double> into a Vec4D before calling CreateLeg
    %extend{
      void CreateLegFromPyVec4D(const ATOOLS::Vec4<double> &p,const ATOOLS::Flavour &fl,
    			  const ATOOLS::ColorID &col=ATOOLS::ColorID(),
    			  const size_t &id=std::string::npos){
    	$self->CreateLeg((ATOOLS::Vec4D&) p, fl,col,id);
      }
    }    

  };// end of class Cluster_Amplitude

  // SWIG needs to rename the following operator to succesfully wrap the functionality
  // In Python, the operator is used like Exception.Input_Exception( STR, EXC),
  // where STR is of std::ostream type and EXC is of ATOOLS::Exception type.
  // The usaage of this method requires std::ostream to be at least minimally wrapped (see iostream.i)
  %rename(Stream_Cluster_Amplitude) &operator<<(std::ostream &ostr,const ATOOLS::Cluster_Amplitude &ampl); 
  std::ostream &operator<< 
    (std::ostream &ostr,const ATOOLS::Cluster_Amplitude &ampl); 

}// end of namespace ATOOLS
