#ifndef Grid_H
#define Grid_H

#include "Node.H"
#include <iostream>
#include <fstream>
#include <string>

namespace ATOOLS {

  class Grid {// adopted from numerical recipes
  private:

    size_t m_dim, m_points, m_order;

    Node<double> *p_root;

    std::string   m_file;
    std::fstream *p_file;

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
		     const size_t j) const;

    size_t Find(const double &x,ATOOLS::Node<double> *const node,
		size_t &l,size_t &r) const;
    double Data(const std::vector<double> &x,
		ATOOLS::Node<double> *const node,const size_t i) const;

    bool WriteOut(Node<double> *current,const size_t node);
    bool ReadIn(Node<double> *current,const size_t node);

  public:

    // constructor
    Grid(const size_t dim,const size_t points);

    // destructor
    ~Grid();

    // member functions
    double Interpolate(const std::vector<double> &x) const;
    double Data(const std::vector<double> &x) const;

    bool WriteOut(const std::string &filename);
    bool ReadIn(const std::string &filename);

    // inline functions
    inline void SetX(const size_t position,
		     const std::vector<double> &x) 
    { SetX(position,x,p_root,1); }
    inline void SetY(const std::vector<size_t> &position,
		     const double value) 
    { SetY(position,value,p_root,1); }

    inline void SetOrder(const size_t order) { m_order=order; }

    inline size_t Points() const    { return m_points; }
    inline size_t Dimension() const { return m_dim;    }
    inline size_t Order() const     { return m_order;  }

    inline double Y(const std::vector<size_t> &position) const
    { return Y(position,p_root,1); }

    inline double operator()(const std::vector<double> &x) const
    { return Interpolate(x); }
    inline double operator[](const std::vector<double> &x) const
    { return Data(x); }
    
  };// end of class Interpolate

}// end of namespace ATOOLS

#endif
