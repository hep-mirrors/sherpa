#ifndef Histogram_H
#define Histogram_H

#include <string>


namespace AMATOOLS {

  class Histogram {
  private :
    int        m_type, m_nbin;
    double     m_lower, m_upper;
    double  ** m_bins;
    long int   m_fills;
    double     m_binsize, m_logbase;
    int        m_depth;
    bool       m_active;
    int        m_logarithmic;
  public :
    Histogram(int type,double xmin, double xmax, int nbins);
    Histogram(Histogram *);
    Histogram(std::string);
    ~Histogram();
    void Reset();
    void Insert(double x);
    void Insert(double x, double weight);
    void InsertRange(double start, double end, double value);
    void Scale(double factor);
    void Output();
    void Output(std::string);
    void Finalize();
    double * Bin(int i);
    double * Bin(double x);
    void Extrapolate(double x, double * ys,int mode);

    // basic access methods
    int Depth() { return m_depth; }
    int    Nbin() const { return m_nbin-2; }
    double Xmin() const { return m_lower; }
    double Xmax() const { return m_upper; }
    double Value(int i) const { return m_bins[i][0]; }
  };


  /*! 
    \file 
    \brief contains class AMATOOLS::Histogram
  */

  /*!
    \class Histogram
    \brief simple histogramm class
    
    This class can be used to create, file, and write out simple
    histograms. It is used in the APHYTOOLS::Primitive_Analysis.
    It can store values, and maxima. Equidistant and
    logarithmic bining is supported.
  */

  /*!
    \fn   Histogram::Histogram(int type,double xmin,double xmax, int nbins)
    \brief Constructor

    \param type   type of histogram (see below)
    \param xmin   start of histogram
    \param xmax   end of histogram
    \param nbins  number of bins

    Initialise a new Histogram acording to the given parameters
    \verbatim
       type = logarithmic*10 + (depth-1)  
       (depth==1) means normal histo 
       (depth==2) means hist with maximum store
    \endverbatim
  */

  /*!
    \fn Histogram::Histogram(Histogram *)
    \brief Copy constructor
  */

  /*!
    \fn Histogram::Histogram(std::string)
    \brief Constructor reading in values from a file
  */    

  /*!
    \fn Histogram::~Histogram()
    \brief Destructor
  */

  /*!
    \fn void Histogram::Reset()
    \brief set all bins of the histogram to zero
  */
  
  /*! 
    \fn void Histogram::Insert(double x)
    \brief add one to corresponding bin
  */
    
  /*!
    \fn void  Histogram::Insert(double x, double weight)
    \brief add value (and possibly the maximum) to corresponding bin
  */
    
  /*!
    \fn  void  Histogram::InsertRange(double start, double end, double value)
    \brief fill all bin in range (bins only partly covered by range are filled coresspondingly less)
  */
    
  /*!
    \fn void  Histogram::Scale(double factor)
    \brief multiply all bins with a given scale factor
  */
    
  /*!
    \fn void  Histogram::Output()
    \brief print content of histogram on screen
  */
    
  /*!
    \fn  void  Histogram::Output(std::string)
    \brief write out content of histogram in a file of given name
  */
    
  /*!
    \fn void  Histogram::Finalize()
    \brief normalize the historam to one
  */
    
  /*!
    \fn  double *  Histogram::Bin(int)
    \brief return the value array of a specified bin
  */
    
  /*!
    \fn  double *  Histogram::Bin(double)
    \brief return the value array of the bin where the given x falls in
  */
    
  /*!
    \fn void Histogram::Extrapolate(double x ,double * ys ,int mode)
    \brief extrapolate, returns array with value and maximum

      This method extrapolates an histogram in the first entry and
      takes the maximum of the two adjacent entries for all other entries,
      therefore it is suited for some kind of weighted (by the first entry)
      hit or miss method.
      The integer in the method Extrapolate labels how the histogram is treated.
      mode usage
      \verbatim
      0    Take only the bin in question
      +-1  Add all bins to the left (-) or right (+)
      +-2  Add all bins apart from the overflow ones to the left (-) or right (+)

      Tested is only mode "+1"!!!
      \endverbatim
  */
    
  /*!
    \fn int     Histogram::Depth() 
    \brief returns number of values per bin

    Depth is 1 if only the value is stored, and 2 if also a maximum is obtainable.
  */
      
  /*!
    \fn int     Histogram::Nbin() const
    \brief  returns number of bins (not counting the overflow bins)
  */
      
  /*!
    \fn double  Histogram::Xmin() const 
    \brief  returns start of histogram
  */
      
  /*!
    \fn double  Histogram::Xmax() const 
    \brief  return end of histogram
  */
      
  /*!
    \fn double  Histogram::Value(int i) const
    \brief  returns the value of bin with number i.
  */
}


#endif 
