#ifndef Color_Integrator_H
#define Color_Integrator_H

#include <vector>
#include <iostream>

namespace PHASIC {

  class Representation {
  protected:

    size_t m_id, m_i, m_j;
    int    m_type;

  public:

    // constructor
    inline Representation(const size_t &id,const int &type): 
      m_id(id), m_i(0), m_j(0), m_type(type) {}

    // inline functions
    inline void SetI(const size_t &i) { m_i=i; }
    inline void SetJ(const size_t &j) { m_j=j; }

    inline size_t I() const { return m_i; }
    inline size_t J() const { return m_j; }

    inline size_t Id() const { return m_id; }

    inline int Type() const { return m_type; }

  };// end of class Representation

  std::ostream &operator<<(std::ostream &ostr,const Representation &v);
  
  class Color_Integrator {
  private:

    std::vector<Representation*> m_ids;

    std::vector<std::vector<size_t> > m_orders, m_confs;

    std::vector<double> m_weights, m_cweights, m_asum;

    size_t m_pairs;
    bool   m_check;

    double m_weight, m_n, m_nv;

    double Factorial(const double &n) const;
    double NOverM(const double &n,const double &m) const;

    size_t GenerateIndex();

    bool DiceColours();
    
    int ConstructConfigurations(std::vector<size_t> ids,
				std::vector<size_t> perm,
				bool sing,double weight);
    int ConstructConfigurations();

    bool TrivialCheck();
    bool CheckPermutation(const std::vector<size_t> &perm);

    bool AddConfiguration(const size_t &i);

  public:

    // constructor
    Color_Integrator();

    // destructor
    ~Color_Integrator();

    // member functions
    bool ConstructRepresentations(const std::vector<size_t> &ids,
				  const std::vector<int> &types);

    bool GeneratePoint();
    bool Initialize();

    std::vector<int> I() const;
    std::vector<int> J() const;

    // inline functions
    inline void SetCheck(const bool check) { m_check=check; }

    inline bool Check() const { return m_check; }

    inline double GlobalWeight() const { return m_weight; }

    inline const std::vector<std::vector<size_t> > &Orders() const
    { return m_orders; }
    inline const std::vector<double> &Weights() const 
    { return m_weights; }

  };// end of class Color_Integrator

  std::ostream &operator<<(std::ostream &ostr,const std::vector<size_t> &v);
  std::ostream &operator<<(std::ostream &ostr,
			   const std::vector<std::vector<size_t> > &v);

}// end of namespace PHASIC

#endif
