#ifndef ATOOLS_Math_Histogram_2D_H
#define ATOOLS_Math_Histogram_2D_H

#include <string>


namespace ATOOLS {

  class Histogram_2D {
  private :
    int        m_type, m_nbin, m_nbinx, m_nbiny;
    double     m_lowerx, m_upperx, m_lowery, m_uppery;
    double   * m_zvalues;
    double   * m_z2values;
    double   * m_psvalues;
    double   * m_tmp;
    double     m_fills,m_psfills;
    double     m_binsizex, m_binsizey, m_logbasex, m_logbasey, m_mcb;
    int        m_depth;
    bool       m_active, m_finished, m_initialized;
    int        m_logarithmicx, m_logarithmicy,m_fuzzyexp;
  public :
    Histogram_2D(int type,double xmin, double xmax, int nbinsx,
                          double ymin, double ymax, int nbinsy);
    Histogram_2D(const Histogram_2D *);
    Histogram_2D(const std::string &);
    ~Histogram_2D();
    void Reset();
    void Insert(double x, double y);
    void Insert(int ix, int iy, double weight, double n=1.);
    void Insert(double x, double y, double weight, double n=1.);
    void InsertMCB(double x, double y, double weight, double n=1.);
    void InsertMCBIM(double coordinatex, double coordinatey, double value);
    void FinishMCB();
    void InsertRange(double startx, double endx,
                     double starty, double endy, double value);
    void Scale(double factor);
    void Output();
    void Output(const std::string);
    void Finalize();
    void Restore();
    double Mean() const;
    inline void SetBin(int i, double entry)  { m_zvalues[i]=entry; }
    inline void SetBin(int ix, int iy, double entry)  { m_zvalues[iy+ix*m_nbiny+1]=entry; }
    inline double Bin(int i) { return m_zvalues[i]; }
    inline double Bin(int ix, int iy) { return m_zvalues[iy+ix*m_nbiny+1]; }
    double Bin(double x, double y);
//     void Extrapolate(double x, double * ys,int mode); //??
    inline void SetFills(double fills) {  m_fills=fills; }

    // basic access methods
    int Type() { return m_type; }
    int Depth() { return m_depth; }
    int    Nbin() const { return m_nbin-2; }
    int    NbinX() const { return m_nbinx; }
    int    NbinY() const { return m_nbiny; }
    double Xmin() const { return m_lowerx; }
    double Xmax() const { return m_upperx; }
    double Ymin() const { return m_lowery; }
    double Ymax() const { return m_uppery; }
    double Value(int i) const { return m_zvalues[i]; }
    double Value(int ix, int iy) const { return m_zvalues[iy+ix*m_nbiny+1]; }
    double Value2(int i) const { return m_z2values[i]; }
    double Value2(int ix, int iy) const { return m_z2values[iy+ix*m_nbiny+1]; }
    double Fills() const { return m_fills; }
    double BinSize() const { return m_binsizex*m_binsizey; }
    double BinSizeX() const { return m_binsizex; }
    double BinSizeY() const { return m_binsizey; }

    double Integral() const;
    double Integral(int xminbin, int xmaxbin, int yminbin, int ymaxbin) const;
    double Integral(double xmin, double xmax, double ymin, double ymax) const;
    double Zmax() const;
    double Zmin() const;
    double LogCoeff() const;


    Histogram_2D & operator+=(const Histogram_2D & histo);
    void Addopt(const Histogram_2D & histo);
  };


  /*!
    \file
    \brief contains class Histogram_2D
  */

  /*!
    \class Histogram_2D
    \brief simple histogramm class

    This class can be used to create, file, and write out simple
    2D histograms. It is not yet used in the Primitive_Analysis.
    It can store values, and maxima. Equidistant and
    logarithmic bining is supported in both x and y.
  */

  /*!
    \fn   Histogram_2D::Histogram_2D(int type,double xmin,double xmax, int nbinsx, double ymin, double ymax, int nbinsy)
    \brief Constructor

    \param type   type of histogram (see below)
    \param xmin   start of histogram in x
    \param xmax   end of histogram in x
    \param nbinsy number of bins in x
    \param ymin   start of histogram in y
    \param ymax   end of histogram in y
    \param nbinsx number of bins in y

    Initialise a new Histogram_2D acording to the given parameters
    \verbatim
       type = logarithmicx*100 + logarithmicy*10 + (depth-1)
       (depth==1) means normal histo
       (depth==2) means hist with maximum store
    \endverbatim
  */

  /*!
    \fn Histogram_2D::Histogram_2D(Histogram_2D *)
    \brief Copy constructor
  */

  /*!
    \fn Histogram_2D::Histogram_2D(std::string)
    \brief Constructor reading in values from a file
  */

  /*!
    \fn Histogram_2D::~Histogram_2D()
    \brief Destructor
  */

  /*!
    \fn void Histogram_2D::Reset()
    \brief set all bins of the histogram to zero
  */

  /*!
    \fn void Histogram_2D::Insert(double x)
    \brief add one to corresponding bin
  */

  /*!
    \fn void  Histogram_2D::Insert(double x, double weight)
    \brief add value (and possibly the maximum) to corresponding bin
  */

  /*!
    \fn  void  Histogram_2D::InsertRange(double start, double end, double value)
    \brief fill all bin in range (bins only partly covered by range are filled coresspondingly less)
  */

  /*!
    \fn void  Histogram_2D::Scale(double factor)
    \brief multiply all bins with a given scale factor
  */

  /*!
    \fn void  Histogram_2D::Output()
    \brief print content of histogram on screen
  */

  /*!
    \fn  void  Histogram_2D::Output(std::string)
    \brief write out content of histogram in a file of given name
  */

  /*!
    \fn void  Histogram_2D::Finalize()
    \brief normalize the historam to one
  */

  /*!
    \fn  double *  Histogram_2D::Bin(int)
    \brief return the value array of a specified bin
  */

  /*!
    \fn  double *  Histogram_2D::Bin(double)
    \brief return the value array of the bin where the given x falls in
  */

  /*!
    \fn void Histogram_2D::Extrapolate(double x ,double * ys ,int mode)
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
    \fn int     Histogram_2D::Depth()
    \brief returns number of values per bin

    Depth is 1 if only the value is stored, and 2 if also a maximum is obtainable.
  */

  /*!
    \fn int     Histogram_2D::Nbin() const
    \brief  returns number of bins (not counting the overflow bins)
  */

  /*!
    \fn double  Histogram_2D::Xmin() const
    \brief  returns start of histogram
  */

  /*!
    \fn double  Histogram_2D::Xmax() const
    \brief  return end of histogram
  */

  /*!
    \fn double  Histogram_2D::Value(int i) const
    \brief  returns the value of bin with number i.
  */
}


#endif
