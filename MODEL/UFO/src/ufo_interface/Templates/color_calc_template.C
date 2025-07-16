#include "MODEL/Main/Single_Vertex.H"
#include "METOOLS/Explicit/Vertex.H"

namespace METOOLS {

  class ${color_name}_Calculator: public Color_Calculator {
  private:

    static std::complex<double> m_cfacs${array_declaration};

    // Outgoing index
    size_t m_out;

    // Mapping of indices, taking care of the 'rotation', in which the
    // vertex occurs such that j[m_inds[i]] = ji with i \in
    // {0,1,2,...,n_external-1}
    std::vector<size_t> m_inds;

  public:

    inline ${color_name}_Calculator(const Vertex_Key &key): 
      Color_Calculator(key) 
    {
      m_out          = p_v->V()->id.back();
      size_t n_ext   = p_v->V()->id.size();
      size_t max_ind = n_ext-1;

      // Set up index mapping
      m_inds.resize(n_ext);
      for(size_t i(0); i<n_ext; i++)
	m_inds[i] = (i + (max_ind - m_out) ) % n_ext;
    }

    std::string Label() const
    {
      return "${color_name}";
    }

    bool Evaluate(const CObject_Vector &j)
    {
      m_c.clear();
      switch(m_out){
${get_cfs}
      }
      return !m_c.empty();
    }

    static void init_cfacs()
    {

${array_vals}
    }

  };

  std::complex<double> ${color_name}_Calculator::m_cfacs${array_declaration} = ${array_init};

}

using namespace METOOLS;
using namespace ATOOLS;

DECLARE_GETTER(${color_name}_Calculator,"${color_name}",
	       Color_Calculator,Vertex_Key);

Color_Calculator *ATOOLS::Getter
<Color_Calculator,Vertex_Key,${color_name}_Calculator>::
operator()(const Vertex_Key &key) const
{
  static int init(0);
  if (init==0) ${color_name}_Calculator::init_cfacs();
  init=1;
  return new ${color_name}_Calculator(key);
}

void ATOOLS::Getter<Color_Calculator,Vertex_Key,${color_name}_Calculator>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"${color_name}";
}
