#ifndef Histogram_H
#define Histogram_H

#include <string>


namespace AMATOOLS {

  class Histogram {
  private :
    double      lower,upper,binsize,logbase;
    int         nbins,depth;
    long int    fills;
    double   ** bins;
    int         type;
    bool        active;
    int         logarithmic;
  public :
    /* 
       type = logarithmic*10 + (depth-1)  
       (depth==1) means normal histo 
       (depth==2) means hist with maximum store
       left, right limits, number of bins
    */
    Histogram(int type,double,double,int);
    Histogram(Histogram *);
    Histogram(std::string);
    ~Histogram();
    void Reset();
    // add one to corresponding bin
    void Insert(double);
    // add value (and possibly the maximum) to corresponding bin
    void Insert(double,double);
    void InsertRange(double,double,double);
    void Scale(double);
    void Output();
    void Output(std::string);
    void Finalize();
    double * Bin(int);
    double * Bin(double);
    /*!
      This method extrapolates an histogram in the first entry and
      takes the maximum of the two adjacent entries for all other entries,
      therefore it is suited for some kind of weighted (by the first entry)
      hit or miss method.
      The integer in the method Extrapolate labels how the histogram is treated.
      mode usage
      0    Take only the bin in question
      +-1  Add all bins to the left (-) or right (+)
      +-2  Add all bins apart from the overflow ones to the left (-) or right (+)

      Tested is only mode "+1"!!!
    */
    void Extrapolate(double,double *,int);
    int Depth() { return depth; }


    // basic access methods
    int    Nbin() const { return nbins-2; }
    double Xmin() const { return lower; }
    double Xmax() const { return upper; }
    double Value(int i) const { return bins[i][0]; }
  };
}


#endif 
