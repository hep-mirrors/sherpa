#include "METOOLS/HadronCurrents/Line_Shapes.H"
#include "METOOLS/HadronCurrents/F0_Decays.H"
#include "METOOLS/HadronCurrents/Rho_Decays.H"
#include "METOOLS/HadronCurrents/Omega_Decays.H"
#include "METOOLS/HadronCurrents/A1_Decays.H"
#include "ATOOLS/Org/Message.H"


using namespace METOOLS;
using namespace ATOOLS;
using namespace std;


/////////////////////////////////////////////////////////////
//
// Wrapper class around all lineshapes (i.e. total widths)
//
/////////////////////////////////////////////////////////////

namespace METOOLS {
  Line_Shapes * LineShapes(NULL);
}

Line_Shapes::Line_Shapes() {}

void Line_Shapes::Init() {
  msg_Out()<<"=================================================\n";
  Add(Flavour(kf_f_0_600),       new F_0_500_Lineshape());
  Add(Flavour(kf_rho_770),       new Rho_770_0_Lineshape());
  Add(Flavour(kf_rho_770_plus),  new Rho_770_plus_Lineshape());
  Add(Flavour(kf_omega_782),     new Omega_782_Lineshape());
  Add(Flavour(kf_rho_1450),      new Rho_1450_0_Lineshape());
  Add(Flavour(kf_rho_1450_plus), new Rho_1450_plus_Lineshape());
  Add(Flavour(kf_rho_1700),      new Rho_1700_0_Lineshape());
  Add(Flavour(kf_rho_1700_plus), new Rho_1700_plus_Lineshape());
  Add(Flavour(kf_omega_1420),    new Omega_1420_Lineshape());
  Add(Flavour(kf_omega_1600),    new Omega_1600_Lineshape());
  Add(Flavour(kf_a_1_1260_plus), new A1_1260_plus_Lineshape());
  Add(Flavour(kf_a_1_1260),      new A1_1260_0_Lineshape());
  msg_Out()<<"Initialised "<<m_lineshapes.size()<<" lineshapes.\n"
	   <<"=================================================\n";
}

Line_Shapes::~Line_Shapes() {
  while (!m_lineshapes.empty()) {
    m_lineshapes.begin()->second->OutputLineshape(0.,3.,300);
    delete m_lineshapes.begin()->second;
    m_lineshapes.erase(m_lineshapes.begin());
  }
}

Total_Width_Base * Line_Shapes::Get(const ATOOLS::Flavour & inflav) {
  long int kfcode = inflav.Kfcode();
  map<long int,Total_Width_Base *>::iterator lsit=m_lineshapes.find(kfcode);
  if (lsit==m_lineshapes.end()) return NULL;
  return lsit->second;
}

void Line_Shapes::Add(const ATOOLS::Flavour & inflav,
		      Total_Width_Base * lineshape) {
  long int kfcode = inflav.Kfcode();
  map<long int,Total_Width_Base *>::iterator lsit=m_lineshapes.find(kfcode);
  if (lsit==m_lineshapes.end()) m_lineshapes[kfcode] = lineshape;
  else {
    delete lsit->second;
    lsit->second = lineshape;
  }
}
