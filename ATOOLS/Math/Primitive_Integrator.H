#ifndef ATOOLS_Math_Primitive_Integrator_H
#define ATOOLS_Math_Primitive_Integrator_H
/*!
  \file Primitive_Integrator.H
  \brief Declares the classes Primitive_Integrand,
  Primitive_Channel and Primitive_Integrator
*/

#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/My_Limits.H"
#include <fstream>
#include <string>
#include <vector>
#include <deque>
#include <map>

namespace ATOOLS {

  class Primitive_Integrand {
  public:

    // destructor
    virtual ~Primitive_Integrand();

    // member functions
    virtual double operator()(const std::vector<double> &point) = 0;
    /*!
      \fn virtual double operator()(const std::vector<double> &point) const
      \brief Evaluates the integrand at the given point
    */
    virtual void AddPoint(const double value,const double weight,
			  const int mode=0);
    /*!
      \fn virtual void AddPoint(const double value,const double weight,
                                const int mode=0)
      \brief Updates the integrand with the approved phasespace points.
    */
    virtual void FinishConstruction(const double apweight);
    /*!
      \fn virtual void FinishConstruction(const double apweight)
      \brief Informs the integrand that the cell buildup is finished
      and updates it with the new a priori weight.
    */

  };// end of class Primitive_Integrand
  /*!
    \class Primitive_Integrand
    \brief Base class for any integrand to be used with
    Primitive_Integrator.
    
    This class is the base class for all integrands which are to be used
    with the Primitive_Integrator. The sole member function operator()
    acts as the function call and evaluates the integrand 
    at the given point.
  */

  struct rmc {
    
    enum code {
      none      = 0,
      construct = 1,
      shuffle   = 2,
      run       = 4
    };

  };// end of struct rmc

  class Primitive_Integrator;

  class Primitive_Channel {
  private:
    
    Primitive_Integrator *p_integrator;

    double m_alpha, m_oldalpha, m_weight, m_splitpos;
    double m_sum, m_sum2, m_max, m_np, m_ssum, m_ssum2, m_snp;

    std::vector<double>             m_this;
    std::vector<Primitive_Channel*> m_next;

    std::vector<std::pair<std::vector<double>,double> > m_points;

    size_t m_pos, m_split;

    static long unsigned int s_npoints;

    friend std::ostream &
    operator<<(std::ostream &str,const Primitive_Channel &channel);

    Primitive_Channel(Primitive_Integrator *const integrator);
    
  public:
    
    // constructor
    Primitive_Channel(Primitive_Integrator *const integrator,
		      Primitive_Channel *const prev,const size_t i,
		      const double &pos=std::numeric_limits<double>::max());
    /*!
      \fn Primitive_Channel(Primitive_Channel *const prev,const size_t i,
      const double &pos=std::numeric_limits<double>::max())
      \brief The standard constructor. It creates a new Primitive_Channel
      by splitting the given channel in the given dimesion into two 
      commensurate parts.
    */
    
    // destructor
    ~Primitive_Channel();

    // member functions
    double Point(Primitive_Integrand *const function,
		 std::vector<double> &point);
    /*!
      \fn double Point(Primitive_Integrand *const function,
      std::vector<double> &point)
      \brief Generates one point inside the integration domain
      employing a uniform weight distrbution.
    */

    bool Find(const std::vector<double> &point) const;
    /*!
      \fn bool Find(const std::vector<double> &point) const
      \brief Determines whether the given point lies inside the 
      integration domain.
    */
    
    void Reset();
    /*!
      \fn void Reset()
      \brief Resets the weight, the squared weight and the number of points.
    */
    void DeletePoints(const int mode=0);
    /*!
      \fn void DeletePoints(const int mode=0)
      \brief Deletes the stored phase space points.
    */
    void SetWeight();
    /*!
      \fn void SetWeight()
      \brief Sets the weight of the channel, i.e. determines the 
      volume of the hypercube the channel represents.
    */
    void SetAlpha(const double &alpha);
    /*!
      \fn void SetAlpha(const double &alpha)
      \brief Sets the external weight of the channel.
    */
    void Store();
    /*!
      \fn void Store()
      \brief Adds the the squared weight, and the number of points 
      to the stored data employing the weight of the channel.
    */
    void SelectSplitDimension();
    /*!
      \fn void SelectSplitDimension(Primitive_Channel *const current)
      \brief Selectes a dimesion for the splitting of the active
      integration channel. 
      
      This method determines the dimesion for the splitting 
      of the active integration channel. During the previous 
      integration step 
      the corresponding channel has stored information about the weight 
      distribution inside itself into a double container. This information
      is evaluated and the dimesion where the splitting results 
      in the largest 
      reduction of the variance or the maximum (according to the given 
      integration mode) is selected.
    */

    bool WriteOut(std::fstream *const file,
		  std::map<Primitive_Channel*,size_t> &pmap) const;
    /*!
      \fn bool WriteOut(std::fstream *const file,
      std::map<Primitive_Channel*,size_t> &pmap) const
      \brief Writes out the channel to the given file.
    */
    bool ReadIn(std::fstream *const file,
		std::map<size_t,Primitive_Channel*> &pmap);
    /*!
      \fn bool ReadIn(std::fstream *const file,
      std::map<size_t,Primitive_Channel*> &pmap)
      \brief Reads in the channel from the given file.
    */

    static void CreateRoot(Primitive_Integrator *const integrator,
			   const std::vector<double> &min,
			   const std::vector<double> &max,
			   std::vector<Primitive_Channel*> &channels);
    /*!
      \fn static void CreateRoot(const std::vector<double> &min,
      const std::vector<double> &max,std::vector<Primitive_Channel*> &channels)
      \brief Creates a root node, i.e. a channel which covers 
      the whole integration domain specified by the parameters min and max.
      Stores the root node and the nodes determining the boundaries 
      to the parameter channels.
    */

    // inline functions
    inline void SetPosition(const size_t pos) { m_pos=pos; }
    /*!
      \fn inline void SetPosition(const size_t pos)
      \brief Sets the position of the channel in the 
      information container which is used to determine the dimension 
      for the splitting of the integration domain.
    */
    inline void SetSplitDimension(const size_t split) { m_split=split; }
    /*!
      \fn inline void SetSplitDimension(const size_t split)
      \brief Sets the splitting dimension for the next splitting 
      of the channel.
    */

    inline double Alpha() const { return m_alpha; }
    /*!
      \fn inline double Alpha() const
      \brief Returns the weight of the channel.
    */
    inline double OldAlpha() const { return m_oldalpha; }
    /*!
      \fn inline double OldAlpha() const
      \brief Returns the previous weight of the channel.
    */
    inline double Weight() const { return m_weight; }
    /*!
      \fn inline double Weight() const
      \brief Returns the volume of the integration domain.
    */

    inline size_t Position() const { return m_pos; }
    /*!
      \fn inline size_t Position() const
      \brief Returns the position of the channel in the 
      information container which is used to determine the dimension 
      for the splitting of the integration domain.
    */
    inline size_t SplitDimension() const { return m_split; }
    /*!
      \fn inline size_t SplitDimension() const
      \brief Returns the dimension for the next splitting 
      of the channel.
    */

    inline double Sum() const { return m_sum; }
    /*!
      \fn inline double Sum() const
      \brief Returns \f$\sum w\f$.
    */
    inline double SSum() const { return m_ssum; }
    /*!
      \fn inline double SSum() const
      \brief Returns \f$\sum_{stored} \frac{w}{\alpha}\f$.
    */
    inline double Sum2() const { return m_sum2; }
    /*!
      \fn inline double Sum2() const
      \brief Returns \f$\sum w^2\f$.
    */
    inline double SSum2() const { return m_ssum2; }
    /*!
      \fn inline double SSum2() const
      \brief Returns \f$\sum_{stored} \frac{w^2}{\alpha^2}\f$.
    */
    inline double Max() const { return m_max; }
    /*!
      \fn inline double Max() const
      \brief Returns \f$max\{w\}\f$.
    */
    inline double Points() const { return m_np; }
    /*!
      \fn inline double Points() const
      \brief Returns \f$n\f$.
    */
    inline double SPoints() const { return m_snp; }
    /*!
      \fn inline double SPoints() const
      \brief Returns \f$n_{stored}\f$.
    */

    inline double Mean() const { return m_sum/m_np; }
    /*!
      \fn inline double Mean() const
      \brief Returns \f$\frac{\sum w}{n}\f$.
    */
    inline double SMean() const { return m_ssum/m_snp; }
    /*!
      \fn inline double SMean() const
      \brief Returns \f$\frac{\sum_{stored} \frac{w}{\alpha}}{n_{stored}}\f$.
    */
    inline double Variance() const    
    { return (m_sum2-m_sum*m_sum/m_np)/(m_np-1.0); }
    /*!
      \fn inline double Variance() const
      \brief Returns \f$\frac{\sum (w^2-\left<w\right>^2)}{n-1}\f$.
    */
    inline double SVariance() const    
    { return (m_ssum2-m_ssum*m_ssum/m_snp)/(m_snp-1.0); }
    /*!
      \fn inline double SVariance() const
      \brief Returns \f$\frac{\sum_{stored} 
      (\frac{w^2}{\alpha^2}-\left<\frac{w}{\alpha}\right>^2)}
      {n_{stored}-1}\f$.
    */
    inline double Sigma() const  
    { return sqrt(Variance()/m_np); }
    /*!
      \fn inline double Sigma() const
      \brief Returns \f$\sqrt{\frac{Variance()}{n}}\f$
    */
    inline double SSigma() const  
    { return sqrt(SVariance()/m_snp); }
    /*!
      \fn inline double SSigma() const
      \brief Returns \f$\sqrt{\frac{SVariance()}{n_{stored}}}\f$
    */

    inline bool Boundary() const { return m_next.front()==NULL; }
    /*!
      \fn inline bool Boundary() const
      \brief Determines whether the channel is a boundary of the 
      integration domain.
    */
    inline void SaveAlpha() { m_oldalpha=m_alpha; }
    /*!
      \fn inline void SaveAlpha()
      \brief Stores the weight.
    */

  };// end of class Primitive_Channel
  /*!
    \class Primitive_Channel
    \brief Represents one integration channel of a Primitive_Integrator.

    This class represents one integration channel of a Primitive_Integrator.
    Its integration domain is a hypercube in the n-dimensional parameter space 
    of the integrator. Inside this hypercube the weight is distributed uniformly.
    The channel stores the sum of the weights \f$\sum w\f$, the sum of the 
    squared weights \f$\sum w^2\f$, the maximum \f$max\{w\}\f$ 
    and the number of points \f$n\f$ which have been diced inside its 
    integration domain. Additionally these data may be added
    to an array of saved data from previous runs or previous optimization 
    steps employing the current weight, \f$\alpha\f$ of the channel.
  */

  std::ostream &
  operator<<(std::ostream &str,const Primitive_Channel &channel);

  struct imc {
    
    enum code {
      varopt = 0,
      maxopt = 1
    };

  };// end of struct imc

  class Primitive_Integrator {
  public:

    typedef std::pair<size_t,std::pair<size_t,size_t> > Position_Pair;
    typedef std::map<std::string,Position_Pair>         Position_Map;

  private:

    long unsigned int m_nopt, m_nmax, m_ndiced;

    double m_error, m_scale, m_apweight;
    double m_sum, m_sum2, m_max, m_np, m_nrealp;

    std::deque<double> m_smax;

    std::vector<double>   m_rmin, m_rmax, m_point, m_asum;
    std::map<size_t,bool> m_nosplit;

    std::vector<Primitive_Channel*> m_channels;

    Primitive_Integrand *p_function;

    size_t m_ncells, m_split, m_shuffle, m_last;

    imc::code m_mode;
    rmc::code m_rmode;

    std::string  m_vname, m_uname;
    Position_Map m_reserved;

    void Split();
    /*!
      \fn void Split()
      \brief Splits the active integration channel.
    */
    bool Shuffle();
    /*!
      \fn bool Shuffle()
      \brief Shuffles the weights of the integration channels.

      This method re-distributes the weights of the integration channels
      according to their variance.
    */

    double Update(const int mode);
    /*!
      \fn double Update(const int mode)
      \brief Updates the integrator.

      This method collects the integrals and the errors of the single 
      channels and determines the total integral and the total error via
      \f[
      I\,=\;\left<w\right>\,=\;\frac{\sum\limits_{i\in channels} \sum_i w}
      {\sum\limits_{i\in channels} n_i}
      \f]
      and
      \f[
      \sigma\,=\;\sqrt{\frac{\sum\limits_{i\in channels} 
      \sum_i (w^2-\left<w\right>^2)}{\sum\limits_{i\in channels} n_i-1}}
      \f]
    */

    void CheckTime() const;
    /*!
      \fn void CheckTime() const
      \brief Checks whether the user defined time limit has beeen reached.
    */

  public:

    // constructor
    Primitive_Integrator();
    /*!
      \fn Primitive_Integrator()
      \brief The standard constructor.
    */

    // destructor
    ~Primitive_Integrator();

    // member functions
    void SetDimension(const size_t dim);
    /*!
      \fn void SetDimension(const size_t dim)
      \brief Sets the dimension of the integration domain.
    */
    void SetMin(const std::vector<double> &min);
    /*!
      \fn void SetMin(const std::vector<double> &min)
      \brief Sets the lower edge of the integration domain.
    */ 
    void SetMax(const std::vector<double> &max);
    /*!
      \fn void SetMax(const std::vector<double> &max)
      \brief Sets the upper edge of the integration domain.
    */ 

    void   Reset();
    /*!
      \fn void   Reset()
      \brief Resets the integrator by deleting all cells.
    */
    void   Initialize();
    /*!
      \fn void   Initialize()
      \brief Initializes the integrator by creating a root cell 
      of type Primitive_Channel inside the integration domain.
    */
    double Integrate(Primitive_Integrand *const function);
    /*!
      \fn double Integrate(Primitive_Integrand *const function)
      \brief Integrates the integrand.
      
      This method iterates the splitting of the root cell until 
      either a sufficiently small error or the maximum number of cells 
      or the maximum number of integration points has been reached. 
    */

    void Point();
    /*!
      \fn void Point()
      \brief Generates one point inside the integration domain according 
      to the weights of the single channels.
    */
    void   Point(std::vector<double> &x);
    /*!
      \fn void   Point(std::vector<double> &x)
      \brief Generates one point inside the integration domain.
    */
    double Weight(const std::vector<double> &x) const;
    /*!
      \fn double Weight(const std::vector<double> &x) const
      \brief Determines the weight of the given point.
    */

    bool WriteOut(const std::string &filename) const;
    /*!
      \fn bool WriteOut(const std::string &filename) const
      \brief Writes out the status of the integrator to the given file.
    */
    bool ReadIn(const std::string &filename);
    /*!
      \fn bool ReadIn(const std::string &filename)
      \brief Reads in the status of the integrator from the given file.
    */

    void Reserve(const std::string &key,
		 const size_t n,const size_t nprev=0);
    /*!
      \fn void Reserve(const std::string &key,
      const size_t n,const size_t nprev=0)
      \brief Reserves the points generated in the given dimension 
      for obejcts with the given key.

      This method is used to tell the integrator that points which are 
      generated inside the n'th - (n+nprev)'th dimesion shall be reserved 
      for the purpose of objects which ask for it using the given key. 
      This is particularly useful if the integration domain shall contain
      the weight distribution of an external multi-channel. The functionality
      of this multi channel may then be taken over by the Primitive_Integrator 
      but the multi-channel may get its discriminators from the integrator
      anyway by simply employing a unique key.
    */
    void Split(const std::string &key,const size_t nprev,
	       const std::vector<double> &pos,const bool nosplit=true);
    /*!
      \fn void Split(const std::string &key,const size_t nprev,
      const std::vector<double> &pos,const bool nosplit=true)
      \brief Splits an indicated dimension according to the given positions.
      
      This method is to be used after the Reserve(const std::string &key,
      const size_t n,const size_t nprev=0) method. After nprev dimesions have
      been reserved for external objects using the given key, the integration 
      domain may be splitted in the nprev'th of these dimensions according to
      the positions given by the parameters in pos. The split flag nosplit 
      determines, whether the integration domain is allowed to be split in 
      this dimension again.
    */
    const double *const Reserved(const std::string &key) const;
    /*!
      \fn const double *const Reserved(const std::string &key) const
      \brief Returns the reserved points for the given key.
    */

    // inline functions
    inline void SetNOpt(const long unsigned int &nopt) { m_nopt=nopt; }
    /*!
      \fn inline void SetNOpt(const long unsigned int &nopt)
      \brief Sets the number of points between two optimization steps.
    */
    inline void SetNMax(const long unsigned int &nmax) { m_nmax=nmax; }
    /*!
      \fn inline void SetNMax(const long unsigned int &nmax)
      \brief Sets the maximum number of points.
    */

    inline void SetNCells(const size_t &ncells) { m_ncells=ncells; }
    /*!
      \fn inline void SetNCells(const size_t &ncells)
      \brief Sets the maximum number of cells.
    */

    inline void SetError(const double &error) { m_error=error; }
    /*!
      \fn inline void SetError(const double &error)
      \brief Sets the maximum error.
    */
    inline void SetScale(const double &scale) { m_scale=scale; }
    /*!
      \fn inline void SetScale(const double &scale)
      \brief Sets a possible scale the result is multiplied with.
    */
    inline void SetMode(const imc::code mode) { m_mode=mode; }
    /*!
      \fn inline void SetMode(const imc::code mode)
      \brief Sets the integration mode.

      This method sets the integration mode. Mode 0 corresponds to variance 
      minimization, mode 1 corresponds to maximum minimization.
    */

    inline void SetVariableName(const std::string &vname) { m_vname=vname; }
    /*!
      \fn inline void SetVariableName(const std::string &vname)
      \brief Sets the name of the variable the the integral is identified with.
    */
    inline void SetUnitName(const std::string &uname) { m_uname=uname; } 
    /*!
      \fn inline void SetUnitName(const std::string &uname)
      \brief Sets the unit of the variable the the integral is identified with.
    */

    inline void SetSplitMode(const size_t split) { m_split=split; }
    /*!
      \fn inline void SetSplitMode(const size_t split)
      \brief Sets the split mode. If set to 0 no channel will be split.
    */
    inline void SetShuffleMode(const size_t shuffle) { m_shuffle=shuffle; }
    /*!
      \fn inline void SetShuffleMode(const size_t shuffle)
      \brief Sets the shuffle mode. If set to 0 weights will not be shuffled.
    */

    inline void SetFunction(Primitive_Integrand *const function)
    { p_function=function; }
    /*!
      \fn inline void SetFunction(Primitive_Integrand *const function)
      \brief Sets the integrand.
    */
    inline Primitive_Integrand *const Function() const
    { return p_function; }
    /*!
      \fn inline Primitive_Integrand *const Function() const
      \brief Returns the integrand.
    */

    inline imc::code Mode() const { return m_mode; }
    /*!
      \fn inline double Mode() const
      \brief Returns the current mode.
    */
    inline rmc::code RunMode() const { return m_rmode; }
    /*!
      \fn inline double RunMode() const
      \brief Returns the current integration mode.
    */
    inline double Mean() const { return m_sum/m_np; }
    /*!
      \fn inline double Mean() const
      \brief Returns \f$\frac{\sum w}{n}\f$.
    */
    inline double Max() const { return m_max; }
    /*!
      \fn inline double Max() const
      \brief Returns \f$max\{w\}\f$.
    */
    inline double Points() const { return m_np; }
    /*!
      \fn inline double Points() const
      \brief Returns \f$n\f$.
    */
    inline double Variance() const    
    { return (m_sum2-m_sum*m_sum/m_np)/(m_np-1.0); }
    /*!
      \fn inline double Variance() const
      \brief Returns \f$\frac{\sum (w^2-\left<w\right>^2)}{n-1}\f$.
    */
    inline double Sigma() const  
    { return sqrt(Variance()/m_np); }
    /*!
      \fn inline double Sigma() const
      \brief Returns \f$\sqrt{\frac{Variance()}{n}}\f$
    */

    inline double APrioriWeight() const { return m_apweight; }
    /*!
      \fn inline double APrioriWeight() const
      \brief Returns the a priori weight of a single channel.
    */

  };// end of class Primitive_Integrator
  /*!
    \class Primitive_Integrator
    \brief Simple integrator class.
    
    This class implements a unitary probability decomposition algorithm 
    for an almost arbitrary integration domain in n-dimensional space.
    The integration domain must be given as a n-dimensional hypercube.
    The integrand is repesented by an instance of the abstract class
    Primitive_Integrand.
  */

}// end of namespace ATOOLS

#endif
