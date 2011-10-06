#ifndef ATOOLS_Math_Histogram_H
#define ATOOLS_Math_Histogram_H

#include <string>


namespace ATOOLS {

  class Histogram {
  private :
    int        m_type, m_nbin;
    double     m_lower, m_upper;
    double   * m_yvalues, * m_y2values;
    double   * m_psvalues, * m_ps2values;
    double  ** m_mvalues;
    double   * m_tmp;
    double     m_fills,m_psfills;
    double     m_mfills,m_mpsfills;
    double     m_binsize, m_logbase,m_mcb;
    int        m_depth;
    bool       m_active, m_finished, m_initialized;
    int        m_logarithmic,m_fuzzyexp;
    std::string m_name;
    void MPIInit();
  public :
    Histogram(int type,double xmin, double xmax, int nbins,
	      const std::string &name="");
    Histogram(const Histogram *);
    Histogram(const std::string &);
    ~Histogram();
    void MPISync();
    void Reset();
    void CopyFrom(const Histogram *histo);
    void Insert(double x);
    void Insert(int i, double weight, double n=1.0);
    void Insert(double x, double weight, double n=1.0);
    void InsertMCB(double x, double weight, double =1.);
    void InsertMCBIM(double coordinate,double value);
    void FinishMCB();
    void InsertRange(double start, double end, double value);
    void Scale(double factor);
    void Output();
    void Output(const std::string);
    void Finalize();
    void Restore();
    double Mean() const;
    void SetBin(int i, double y)  { m_yvalues[i]=y; }
    void SetBin2(int i, double y2)  { m_y2values[i]=y2; }
    void SetBinPS(int i, double ps)  { m_psvalues[i]=ps; }
    void SetBinPS2(int i, double ps)  { m_ps2values[i]=ps; }
    double Bin(int i) const;
    double Bin2(int i) const;
    double Bin(double x) const;
    double BinOrInterpolate(int i) const;
    double BinOrInterpolate(double x) const;
    void Extrapolate(double x, double * ys,int mode);
    inline void SetFills(double fills) {  m_fills=fills; }

    // basic access methods
    int Type() { return m_type; }
    int Depth() { return m_depth; }
    int    Nbin() const { return m_nbin-2; }
    double Xmin() const { return m_lower; }
    double Xmax() const { return m_upper; }
    double Value(int i) const { return m_yvalues[i]; }
    double Value2(int i) const { return m_y2values[i]; }
    double Fills() const { return m_fills; }
    double BinSize() const { return m_binsize; }

    double Integral() const;
    double Ymax() const;
    double Ymin() const;
    double LogCoeff() const;


    Histogram & operator+=(const Histogram & histo);
    Histogram & operator=(const Histogram & histo);
    void Addopt(const Histogram & histo);
    void AddGeometric(const Histogram & histo);
    void BinMin(const Histogram & histo);
    void BinMax(const Histogram & histo);
    int CheckStatistics(const Histogram & histo,double& avgs,double& maxs);

    inline const std::string &Name() const { return m_name; }
  };


  /*! 
    \file 
    \brief contains class Histogram
  */

  /*!
    \class Histogram
    \brief simple histogramm class
    
    This class can be used to create, file, and write out simple
    histograms. It is used in the Primitive_Analysis.
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
