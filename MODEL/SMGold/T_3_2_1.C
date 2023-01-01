#include "MODEL/Main/Single_Vertex.H"
#include "METOOLS/Explicit/Vertex.H"

namespace METOOLS {

  class T_3_2_1_Calculator: public Color_Calculator {
  private:

    static std::complex<double> m_cfacs[3][3][3][3];

    // Outgoing index
    size_t m_out;

    // Mapping of indices, taking care of the 'rotation', in which the
    // vertex occurs such that j[m_inds[i]] = ji with i \in
    // {0,1,2,...,n_external-1}
    std::vector<size_t> m_inds;

  public:

    inline T_3_2_1_Calculator(const Vertex_Key &key): 
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
      return "T_3_2_1";
    }

    bool Evaluate(const CObject_Vector &j)
    {
      m_c.clear();
      switch(m_out){

case 2:
for(size_t i(0); i<3; i++)
for(size_t m(0); m<3; m++)
if(m_cfacs[i][m][(*j[m_inds[1]])(0)-1][(*j[m_inds[0]])(1)-1]!=0.0)
m_c.push_back(CInfo(i+1,m+1,m_cfacs[i][m][(*j[m_inds[1]])(0)-1][(*j[m_inds[0]])(1)-1]));
break;

case 1:
for(size_t i(0); i<3; i++)
if(m_cfacs[(*j[m_inds[2]])(1)-1][(*j[m_inds[2]])(0)-1][i][(*j[m_inds[0]])(1)-1]!=0.0)
m_c.push_back(CInfo(0,i+1,m_cfacs[(*j[m_inds[2]])(1)-1][(*j[m_inds[2]])(0)-1][i][(*j[m_inds[0]])(1)-1]));
break;

case 0:
for(size_t i(0); i<3; i++)
if(m_cfacs[(*j[m_inds[2]])(1)-1][(*j[m_inds[2]])(0)-1][(*j[m_inds[1]])(0)-1][i]!=0.0)
m_c.push_back(CInfo(i+1,0,m_cfacs[(*j[m_inds[2]])(1)-1][(*j[m_inds[2]])(0)-1][(*j[m_inds[1]])(0)-1][i]));
break;

default:
if(m_cfacs[(*j[m_inds[2]])(1)-1][(*j[m_inds[2]])(0)-1][(*j[m_inds[1]])(0)-1][(*j[m_inds[0]])(1)-1]!=0.0)
m_c.push_back(CInfo(0,0,m_cfacs[(*j[m_inds[2]])(1)-1][(*j[m_inds[2]])(0)-1][(*j[m_inds[1]])(0)-1][(*j[m_inds[0]])(1)-1]));

      }
      return !m_c.empty();
    }

    static void init_cfacs()
    {

m_cfacs[0][0][0][0]=std::complex<double>(4.71404520791031733662e-01,0.00000000000000000000e+00);
m_cfacs[0][0][1][1]=std::complex<double>(-2.35702260395515839075e-01,0.00000000000000000000e+00);
m_cfacs[0][0][2][2]=std::complex<double>(-2.35702260395515894587e-01,0.00000000000000000000e+00);
m_cfacs[0][1][0][1]=std::complex<double>(7.07106781186547572737e-01,0.00000000000000000000e+00);
m_cfacs[0][2][0][2]=std::complex<double>(7.07106781186547572737e-01,0.00000000000000000000e+00);
m_cfacs[1][0][1][0]=std::complex<double>(7.07106781186547572737e-01,0.00000000000000000000e+00);
m_cfacs[1][1][0][0]=std::complex<double>(-2.35702260395515839075e-01,0.00000000000000000000e+00);
m_cfacs[1][1][1][1]=std::complex<double>(4.71404520791031733662e-01,0.00000000000000000000e+00);
m_cfacs[1][1][2][2]=std::complex<double>(-2.35702260395515894587e-01,0.00000000000000000000e+00);
m_cfacs[1][2][1][2]=std::complex<double>(7.07106781186547572737e-01,0.00000000000000000000e+00);
m_cfacs[2][0][2][0]=std::complex<double>(7.07106781186547572737e-01,0.00000000000000000000e+00);
m_cfacs[2][1][2][1]=std::complex<double>(7.07106781186547572737e-01,0.00000000000000000000e+00);
m_cfacs[2][2][0][0]=std::complex<double>(-2.35702260395515894587e-01,0.00000000000000000000e+00);
m_cfacs[2][2][1][1]=std::complex<double>(-2.35702260395515894587e-01,0.00000000000000000000e+00);
m_cfacs[2][2][2][2]=std::complex<double>(4.71404520791031789173e-01,0.00000000000000000000e+00);

    }

  };

  std::complex<double> T_3_2_1_Calculator::m_cfacs[3][3][3][3] = {{{{ std::complex<double>(0.0,0.0) }}}};

}

using namespace METOOLS;
using namespace ATOOLS;

DECLARE_GETTER(T_3_2_1_Calculator,"T_3_2_1",
	       Color_Calculator,Vertex_Key);

Color_Calculator *ATOOLS::Getter
<Color_Calculator,Vertex_Key,T_3_2_1_Calculator>::
operator()(const Vertex_Key &key) const
{
  static int init(0);
  if (init==0) T_3_2_1_Calculator::init_cfacs();
  init=1;
  return new T_3_2_1_Calculator(key);
}

void ATOOLS::Getter<Color_Calculator,Vertex_Key,T_3_2_1_Calculator>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"T_3_2_1";
}
