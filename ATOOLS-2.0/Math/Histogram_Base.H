#ifndef Histogram_Base_H
#define Histogram_Base_H

#include <string>


namespace ATOOLS {
  class Root_Histogram;
  class Histogram;

  class Histogram_Base {
  public :
    Histogram_Base()                                          {}
    Histogram_Base(const std::string,const std::string)       {}
    Histogram_Base(const int,const double,const double,const int,
		   const std::string=std::string(""))         {}
    Histogram_Base(const int,const int,const double *,
		   const std::string=std::string(""))         {}
    Histogram_Base(const Histogram *)                         {}
    Histogram_Base(const Root_Histogram *)                    {}
    virtual ~Histogram_Base()                                 {}

    virtual void Reset()                                               =0;
    virtual void Insert(double x)                                      =0;
    virtual void Insert(double x, double weight, int ntimes)           =0;
    virtual double Mean()     const                                    =0;
    virtual double Integral() const                                    =0;
    virtual double Ymax()     const                                    =0; 
    virtual double Ymin()     const                                    =0; 
    virtual std::string Title() const                                  =0;

    virtual void   InsertRange(double start, double end, double value) =0;
    virtual double LogCoeff() const                                    =0;

    virtual void   Scale(double factor)                                =0;
    virtual void   Output()                                            =0;
    virtual void   Output(const std::string)                           =0;
    virtual void   Finalize()                                          =0;
    virtual void   Restore()                                           =0;
    virtual void   SetBin(int i, double x)                             =0;
    virtual void   Extrapolate(double x, double * ys,int mode)         =0;
    virtual void   SetTitle(const std::string)                         =0;
    virtual void   SetLogY()                                           =0;
    virtual void   SetLegend(const bool,
			     const double,const double,
			     const double,const double)                =0;

    // basic access methods
    virtual int    Type()           const=0;
    virtual int    Depth()          const=0;
    virtual int    Nbins()          const=0;
    virtual double Xmin()           const=0;
    virtual double Xmax()           const=0;
    virtual double Value(int i)     const=0;
    virtual double Fills()          const=0;
    virtual double BinSize(int i=1) const=0;

    virtual Histogram_Base & operator+=(const Root_Histogram & histo)=0;
    virtual Histogram_Base & operator+=(const Histogram & histo)=0;
  };
}



#endif 
