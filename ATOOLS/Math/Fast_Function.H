#ifndef ATOOLS_Math_Fast_Function_H
#define ATOOLS_Math_Fast_Function_H

#include "ATOOLS/Math/Function_Base.H"
#include <vector>
#include <iostream>

namespace ATOOLS {

  class Intervall {
    double m_minval;
    double m_maxval;
  public:
    Intervall(double a=0, double e=0) :
      m_minval(a), m_maxval(e) {}
    inline double   Min() const { return m_minval; }
    inline double   Max() const { return m_maxval; }
    inline double&  Min() { return m_minval; }
    inline double&  Max() { return m_maxval; }
    inline double   Length() const { return m_maxval-m_minval; }
    inline bool     IsIncluded(double x) {
      return m_minval<=x && m_maxval >=x; 
    }

    friend  std::ostream & operator<<(std::ostream &, const Intervall & );
  };

  class  Fast_Function : public Function_Base {
  public:
    class Pair {
    public:
      double x, y;
      Pair(double _x=0, double _y =0) : x(_x), y(_y) {};
    };

    typedef std::vector<Pair> Data_List;
    typedef Data_List::iterator  Data_Iter;
  private:
    Data_List m_data;
    double m_ymin, m_ymax;

    double LinInter(Data_Iter & it, double x);
    double LinInterInv(Data_Iter & it, double y);
    double LinInter(int, double x);
    double LinInterInv(int, double y);

    int Nesting(int,int,double);
  public:
    Fast_Function();
    Fast_Function(int size);
    
    Data_Iter Insert(double x, double y);

    double operator() (double x);

    // return a x value to a given y value, using a linear interpolation method
    // assuming a monotonous function
    double Invers(double y);

    void   Init(Function_Base &, double xmin, double xmax, int min_points);

    const Intervall XRange() const {return Intervall(Xmin(),Xmax());}
    const Intervall YRange() const {return Intervall(m_ymin,m_ymax);}

    double Xmin() const {return m_data.front().x;}
    double Xmax() const {return m_data.back().x;}

    double Ymin() const {return m_ymin;}
    double Ymax() const {return m_ymax;}

    void   Clear() { 
      m_data.clear();
      m_ymin =  1.e99;
      m_ymax = -1.e99;
    }

    void   WriteOut(const char * name);
    bool   ReadIn(const char * name);

    //
    int    Length()  const {return m_data.size();}
    double GetX(int i) const {return m_data[i].x;}
    double GetY(int i) const {return m_data[i].y;}
    void   SetX(int i,double value) {m_data[i].x = value;}
    void   SetY(int i,double value) {m_data[i].y = value;}

    friend  std::ostream & operator<<(std::ostream & s, const Fast_Function & ff);
  };

  std::ostream & operator<<(std::ostream &, const Fast_Function & );

  std::ostream & operator<<(std::ostream &, const Intervall & );

  std::ostream & operator<<(std::ostream &, const Fast_Function::Pair &);


  /*!
    \file 
    \brief Contains the class Fast_Function
  */
  
  /*!
    \class Fast_Function
    \brief provides an interpolated function on the basis of number of data points

    This class can be use to obtain values for a function \f$f(x)\f$ which is only known 
    on certain points \f${x_1,\dots,x_n}\f$.
  */

  /*!
   \fn Fast_Function::Fast_Function()
   \brief Default constructor
  */

  /*!
   \fn Fast_Function::Fast_Function(int size)
   \brief Special constructor initializing a given number of data pairs.
  */

  /*!
   \fn Data_Iter Fast_Function::Insert(double x, double y)
   \brief adds a pair \f$(x,y)\f$ of data points, and return an iterator to the new element.
  */

  /*!
    \fn double Fast_Function::operator() (double x)
    \brief returns \f$y=f(x)\f$ using a linear interpolation method
  */

  /*!
   \fn double Fast_Function::Invers(double y)
   \brief returns \f$x=f^{-1}(x)\f$ using a linear interplation method 
  */

  /*!
    \fn void Fast_Function::Init(Function_Base &, double xmin, double xmax, int min_points)
    \brief initailize a Fast_Function, with a given (slow) function
  */

  /*!
    \fn const Intervall Fast_Function::XRange() const
    \brief returns domain of data points.
  */

  /*!
    \fn  const Intervall Fast_Function::YRange() const
    \brief returns size of y domain.
  */

  /*!
    \fn const double Fast_Function::Xmin() const
    \brief    return minimum x (of data points)
  */

  /*!
    \fn const double Fast_Function::Xmax() const
    \brief    return maximum x (of data points)
  */

  /*!
    \fn const double Fast_Function::Ymin() const
    \brief    return minimum y (of data points)
  */

  /*!
    \fn const double Fast_Function::Ymax() const
    \brief    return maximum y (of data points)
  */

  /*!
    \fn    void Fast_Function::Clear() 
    \brief delete all stored data pairs.
  */

  /*!
    \fn    void Fast_Function::WriteOut(char * name)
    \brief write out all data points to a file of given name
  */

  /*!
    \fn    bool Fast_Function::ReadIn(char * name)
    \brief read in data points from a file of given name
  */

  /*!
    \fn   int    Fast_Function::Length()  const 
    \brief return number of data points
  */

  /*!
    \fn    double Fast_Function::GetX(int i) const 
    \brief access x value of \em i th data point
  */

  /*!
    \fn    double Fast_Function::GetY(int i) const 
    \brief access y value of \em i th data point
  */

  /*!
    \fn    double Fast_Function::SetX(int i, double value) 
    \brief define x value of \em i th data point
  */

  /*!
    \fn    double Fast_Function::SetY(int i, double value) 
    \brief define y value of \em i th data point
  */
}



#endif
