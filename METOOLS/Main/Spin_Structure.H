#ifndef METOOLS_Main_Spin_Structure_H
#define METOOLS_Main_Spin_Structure_H

#include "METOOLS/Main/Polarization_Index.H"
#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Math/Matrix.H"
#include "ATOOLS/Org/Message.H"
#include <algorithm>
#include <iomanip>
#include "ATOOLS/Org/Exception.H"

namespace METOOLS
{
  bool SortByFirst(const std::pair<int,int> p1, const std::pair<int,int> p2);

  template<class Value>
  class Spin_Structure : public std::vector<Value>, public Polarization_Index {
  protected:
    size_t GetNumber(std::vector<std::pair<int,int> >& spins) const
    {
      sort(spins.begin(),spins.end(),SortByFirst);
      
      if(spins.size()!=m_spins.size()) {
	msg_Error()<<METHOD<<" Error: wrong size of spin std::vector."<<std::endl;
        abort();
      }
      int mult(1);
      size_t num(0);
      for(size_t i=0; i<spins.size(); i++) {
	num += mult * spins[i].second;
	mult *= m_spins[i];
      }
      if(num>this->size()) {
	msg_Error()<<METHOD<<" Error: tried to access value out of bounce. "
	  <<"num="<<num<<" > "<<this->size()<<std::endl;
        abort();
      }
      return num;
    }

  public:
    Spin_Structure() {}
    Spin_Structure(const std::vector<int>& spins, const Value& value):
      Polarization_Index(spins)
    {
      this->resize(m_n, value);
    }
    
    Spin_Structure(const ATOOLS::Flavour_Vector& flavs, const Value& value)
    {
      m_spins = std::vector<int>(flavs.size());
      m_n=1;
      for(size_t i=0;i<flavs.size();i++) {
        if(flavs[i].IsVector() && flavs[i].IsMassive()==0) m_spins[i] = 2;
        else m_spins[i] = flavs[i].IntSpin()+1;
	m_n*=m_spins[i];
      }
      this->resize(m_n,value);
    }
    
    Spin_Structure(const ATOOLS::Flavour_Vector& flavs,
                   const std::vector<int>& indices)
    {
      m_spins = std::vector<int>(indices.size());

      m_n=1;
      for(size_t i=0;i<indices.size();i++) {
        if(flavs[indices[i]].IsVector() && flavs[indices[i]].IsMassive()==0)
          m_spins[i] = 2;
        else m_spins[i] = flavs[indices[i]].IntSpin()+1;
	m_n*=m_spins[i];
      }
      this->resize(m_n);
    }
    
    Spin_Structure(const ATOOLS::Particle_Vector& particles)
    {
      m_spins = std::vector<int>(particles.size());

      m_n=1;
      for(size_t i=0;i<particles.size();i++) {
        if(particles[i]->Flav().IsVector() && particles[i]->Flav().IsMassive()==0) m_spins[i] = 2;
        else m_spins[i] = particles[i]->Flav().IntSpin()+1;
	m_n*=m_spins[i];
      }
      this->resize(m_n);
    }
    
    ~Spin_Structure() {}

    inline size_t GetNumber(const std::vector<int> &spins) const { return (*this)(spins); }

    inline std::vector<int> GetSpinCombination(size_t number) const { return (*this)(number); }

    void    Insert(Value value, std::vector<std::pair<int,int> >& spins)
    {
      (*this)[GetNumber(spins)]=value;
    }

    void    Insert(Value value, const std::vector<int>& spins)
    {
      (*this)[GetNumber(spins)]=value;
    }

    void    Add(Value value, std::vector<std::pair<int,int> >& spins)
    {
      (*this)[GetNumber(spins)]=value;
    }
    
    void    Insert(Value value, size_t index)
    {
      (*this)[index]=value;
    }

    Value   Get(const std::vector<int>& spins) const
    {
      return (*this)[GetNumber(spins)];
    }
    
    Value   Get(size_t index) const
    {
      return (*this)[index];
    }

    void    CreateTrivial(Value value)
    {
      size_t n = this->size();
      this->clear();
      this->resize(n,value);
    }

    Spin_Structure<Value>& operator+= (const Spin_Structure<Value>& addend)
    {
      for(size_t i=0;i<this->size();i++) {
	(*this)[i]+=addend[i];
      }
      return *this;
    }
  };

  template<class Value>
  std::ostream& operator<<(std::ostream& ostr, const Spin_Structure<Value>& s) {
    ostr<<"   Spin_Structure with "<<s.size()<<" spin combinations:"<<std::endl;
    for(size_t i=0;i<s.size();i++) {
      ostr<<std::setw(3)<<i;
      std::vector<int> spins = s.GetSpinCombination(i);
      for(size_t j=0;j<spins.size();j++) {
        ostr<<std::setw(8)<<spins[j]<<" | ";
      }
      ostr<<s[i]<<std::endl;
    }
    return ostr;
  }

  class Spin_Amplitudes : public Spin_Structure<Complex> {
  public:
    Spin_Amplitudes(const std::vector<int>& spins, const Complex& value);
    Spin_Amplitudes(const ATOOLS::Flavour_Vector& flavs, const Complex& value);
    Spin_Amplitudes(const ATOOLS::Flavour_Vector& flavs,
                    const std::vector<int>& indices);
    Spin_Amplitudes(const ATOOLS::Particle_Vector& particles);
    virtual ~Spin_Amplitudes();
    double  SumSquare() const;
    virtual void Calculate(const ATOOLS::Vec4D_Vector& momenta,bool anti=false);
  };

  /*! 
    \class Spin_Structure
    \brief Storing objects by spin combination.
    
    This class provides methods to access and manipulate arbitrary objects, by
    specifying the corresponding spin combination. It inherits from STL vector
    which contains one 'object' for each spin combination.
  */

  /*!
    \var std::vector<int> Spin_Structure::m_spins
    \brief vector containing the number of spin combinations for each particle "node".
  */

  /*!
    \var std::vector<Value>
    \brief storage of objects
    
    Objects are ordered as in this example (where the node part1 has
    2 spin combinations, part2 has 1, part3 has 3, part4 has 2):
        |  part1 |    part2 |    part3 |    part4 |
      0 |      0 |        0 |        0 |        0 |
      1 |      1 |        0 |        0 |        0 |
      2 |      0 |        0 |        1 |        0 |
      3 |      1 |        0 |        1 |        0 |
      4 |      0 |        0 |        2 |        0 |
      5 |      1 |        0 |        2 |        0 |
      6 |      0 |        0 |        0 |        1 |
      7 |      1 |        0 |        0 |        1 |
      8 |      0 |        0 |        1 |        1 |
      9 |      1 |        0 |        1 |        1 |
     10 |      0 |        0 |        2 |        1 |
     11 |      1 |        0 |        2 |        1 |

     The polarisation index for each particle corresponds to:
      - <var>lambda</var>
        - \f$0 \to \epsilon^+\f$
        - \f$1 \to \epsilon^-\f$
        - \f$2 \to \epsilon^0\f$
        .
  */

  /*!
    \fn size_t Spin_Structure::GetNumber(const vector<int>& spins) const
    \brief Determine number of given combination in the STL vector.

    Here, the spins have to be specified in the same order as in m_spins (i.e.
    as the flavours or particles in the constructor where specified).
  */

  /*!
    \fn size_t Spin_Structure::GetNumber(vector<pair<int,int> >& spins) const
    \brief Determine number of given combination in the STL vector.

    Here, the spins (second int in pair) can be in arbitrary order, as long as
    the first int in the pair specifies the number in m_spins it belongs to.
  */

  /*!
    \fn std::vector<int> Spin_Structure::GetSpinCombination(size_t number) const
    \brief Determine spin combination from number of combination in the vector.

    The other way around from the GetNumber methods. Useful to
    traverse through all spin combinations.
  */

  /*!
    \fn Spin_Structure::Spin_Structure(ATOOLS::Flavour* flavs, size_t size)
    \brief Constructor
  */

  /*!
    \fn Spin_Structure::Spin_Structure(ATOOLS::Particle_Vector particles)
    \brief Constructor

    ATOOLS::Particles' flavours determine the number of spin combinations of each node.
  */

  /*!
    \fn Spin_Structure::~Spin_Structure()
    \brief Destructor

    Doesn't do anything, because there are no pointer members.
  */

  /*!
    \fn void Spin_Structure::Insert(Value value, vector<pair<int,int> >& spins)
    \brief Inserts value at the right position determined by GetNumber.
  */

  /*!
    \fn void Spin_Structure::Insert(Complex amp, size_t index)
    \brief Inserts value at the given position.
  */

  /*!
    \fn Value Spin_Structure::Get(const vector<int>& spins) const
    \brief Retrieves the value from the right position determined by GetNumber.
  */

  /*!
    \fn vector<Complex> Spin_Structure::Get(size_t index) const;
    \brief Retrieves the value from the given position.
  */

  /*!
    \fn void Spin_Structure::CreateTrivial(Value value)
    \brief Inserts the given object for all spin combinations.
  */

}


#endif
