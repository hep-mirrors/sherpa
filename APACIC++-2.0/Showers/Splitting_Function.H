#ifndef Splitting_Function_H
#define Splitting_Function_H

#include "Flavour.H"
#include "Sudakov_Tools.H"
#include <iostream>

namespace APACIC {
  /*! 
    Purely virtual class -
    All specific splitting functions are derived from this class.
  */
  class Splitting_Function {
  protected:
    const static double Nc;
    const static double CF;
    const static double CA;
    const static double TR;
    // Flavours of this particular splitting function
    APHYTOOLS::Flavour flavs[3]; // a -> b c
    // Domain of last intergral
    double zmin,zmax;
    // Last integral
    double lastint;
  public:
    //! returns the exact value of a splitting function at some given z
    virtual double operator() (double) = 0;

    /*! 
      returns some z chosen according to logarithmic accuracy 
      dependent on the lower cut off zmin                         
    */
    virtual double GetZ()=0;             

    //! returns the maximal coupling constant
    virtual double GetCoupling()=0;             
    //! returns the actual coupling constant at scale t
    virtual double GetCoupling(double t)=0;             

    /*! 
      returns the correction weight of the exact splitting function w.r.t. the
      crude (logarithmic) estimate.
    */
    virtual double GetWeight(double,double,bool)=0;

    /*! 
      returns integral of the splitting function to logarithmic
      accuracy with symmetric boundaries zmin=z_0 and zmax=1-z_0  
    */
    virtual double CrudeInt(double z0)  {return CrudeInt(z0,1.-z0);};
    /*! 
      returns integral of the splitting function to logarithmic
      accuracy with asymmetric boundaries zmin and zmax               
    */
    virtual double CrudeInt(double,double) = 0;        

    //! returns the last value of the integral.
    virtual double GetLastInt() {return lastint;};

    //! adds a splitting function (used in groups only) 
    virtual void Add(Splitting_Function *) {
      std::cerr<<" Error in Spliting_Function: something nasty is going on"<<std::endl;
      std::cout<<" Error in Spliting_Function: something nasty is going on"<<std::endl;
    };
    //! select one splitting function (used in Groups only)
    virtual void SelectOne() {};
  
    // determines APHYTOOLS::Flavours
    virtual APHYTOOLS::Flavour & GetFlA() { return flavs[0];};
    virtual APHYTOOLS::Flavour & GetFlB() { return flavs[1];};
    virtual APHYTOOLS::Flavour & GetFlC() { return flavs[2];};

    virtual void PrintStat(int mode=0) {
      if (mode>0) for(int i=0;i<mode;++i) std::cout<<' ';
      std::cout<<"Splitting Function: "
	  <<GetFlA()<<" -> "<<GetFlB()<<" + "<<GetFlC()<<std::endl;
    };
  };
} 

#endif 
