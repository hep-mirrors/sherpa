#ifndef Fast_Function_H
#define Fast_Function_H

#include "Function_Base.H"
#include <vector>
#include <iostream>

namespace AMATOOLS {

  // =======================================
  // adapted from fastfunc by R. Kuhn
  // =======================================

  class Intervall {
    double minval;
    double maxval;
  public:
    Intervall(double a=0, double e=0) :
      minval(a), maxval(e) {}
    inline ~Intervall() {}
    inline double Min() const { return minval; }
    inline double Max() const { return maxval; }
    inline double& Min() { return minval; }
    inline double& Max() { return maxval; }
    inline double Length() const { return maxval-minval; }
    inline bool IsIncluded(double x) {
      return minval<=x && maxval >=x; }

     friend  std::ostream & operator<<(std::ostream &, const Intervall & );
  };

  class  Fast_Function : public Function_Base {
  public:
    class Pair {
    public:
      double x, y;
      Pair(double _x=0, double _y =0) : x(_x), y(_y) {};
    };

    
    typedef std::vector<Pair> DList;
    typedef DList::iterator  DIter;
  private:
    DList data;
    double ymin, ymax;

    double LinInter(DIter & it, double x);
    double LinInterInv(DIter & it, double y);
  public:
    Fast_Function();
    Fast_Function(int size);
    
    // insert a x-y-pair int the function base (returns iterator to the new element)
    DIter Insert(double x, double y);
    // return a y value to a given x value, using a linear interpolation method
    double operator() (double x);
    // return a x value to a given y value, using a linear interpolation method
    // assuming a monotonous function
    double Invers(double y);

    void Init(Function_Base &, double xmin, double xmax, int min_points);

    // return domain and y range
    const Intervall XRange() const {return Intervall(Xmin(),Xmax());}
    const Intervall YRange() const {return Intervall(ymin,ymax);}
    // return minimum and maximum x value
    const double Xmin() const {return data.front().x;}
    const double Xmax() const {return data.back().x;}
    // return minimum and maximum x value
    const double Ymin() const {return ymin;}
    const double Ymax() const {return ymax;}
    // delete all stored data pairs
    void Clear() { 
      data.clear();
      ymin=1.e99;
      ymax=-1.e99;
    }

    //    void output(double xmin,double xmax,char* name);
    //    int  input(char* name);
    //    void   as_output(double xmin,double xmax,char* name);
    //    double as_input(char* name,int);
    void WriteOut(char * name);
    bool ReadIn(char * name);


    //
    int    Length()  const {return data.size();}
    double GetX(int i) const {return data[i].x;}
    double GetY(int i) const {return data[i].y;}
    void   SetX(int i,double value) {data[i].x = value;}
    void   SetY(int i,double value) {data[i].y = value;}

    friend  std::ostream & operator<<(std::ostream & s, const Fast_Function & ff);
  };

  std::ostream & operator<<(std::ostream &, const Fast_Function & );

  std::ostream & operator<<(std::ostream &, const Intervall & );

  std::ostream & operator<<(std::ostream &, const Fast_Function::Pair &);



}



#endif
