#ifndef Grid_H
#define Grid_H

#include "Node.H"

namespace ATOOLS {

  class Grid {// adopted from numerical recipes
  private:

    size_t m_dim, m_points;

    Node<double> *p_root;

    void Allocate(Node<double> *current,
		  const size_t node);

    void SetX(const size_t position,const std::vector<double> &x,
	      Node<double> *current,const size_t node);
    double Y(const std::vector<size_t> &position,
	     Node<double> *current,const size_t node) const;
    void SetY(const std::vector<size_t> &position,const double value,
	      Node<double> *current,const size_t node);
    
    double Interpolate(const std::vector<double> &xgrid,
		       const std::vector<double> &ygrid,
		       const double &x) const;
    void Interpolate(const std::vector<double> &x,const size_t i,
		     Node<double> *const node,std::vector<double> &y,
		     const size_t j);

  public:

    // constructor
    Grid(const size_t dim,const size_t points);

    // destructor
    ~Grid();

    // member functions
    double Interpolate(const std::vector<double> &x) const;

    // inline functions
    inline void SetX(const size_t position,
		     const std::vector<double> &x) 
    { SetX(position,x,p_root,1); }
    inline void SetY(const std::vector<size_t> &position,
		     const double value) 
    { SetY(position,value,p_root,1); }

    inline size_t Points() const    { return m_points; }
    inline size_t Dimension() const { return m_dim;    }

    inline double Y(const std::vector<size_t> &position) const
    { return Y(position,p_root,1); }

    inline double operator()(const std::vector<double> &x) const
    { return Interpolate(x); }

  };// end of class Interpolate

};// end of namespace ATOOLS

#endif
