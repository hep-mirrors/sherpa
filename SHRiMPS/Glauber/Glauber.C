#include "SHRiMPS/Glauber/Glauber.H"
#include "SHRiMPS/Event_Generation/Soft_Diffractive_Event_Generator.H"
#include "SHRiMPS/Tools/Special_Functions.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"
#include "SHRiMPS/Cross_Sections/Sigma_Inelastic.H"

using namespace SHRIMPS;
using namespace ATOOLS;
//using namespace ATOOLS;
//using namespace Glauber;
using namespace std;

Glauber::Glauber(int N1, int N2) :
  m_numNucleons1(N1), m_numNucleons2(N2),
  m_radius1(1.12*pow(N1,1./3.) - 0.86*pow(N1,-1./3.)),
  m_radius2(1.12*pow(N2,1./3.) - 0.86*pow(N2,-1./3.)),
  m_rho0(0.17), m_w(0.), m_a(0.54), 
  m_nucleus_1_position(ATOOLS::Vec4D(0.,0.,0.,0.)),
  m_nucleus_2_position(ATOOLS::Vec4D(0.,0.,0.,0.)),
  m_type_abs(1), m_type_SD1(2), m_type_SD2(3), m_type_DD(4), m_type_elas(5)
{
  msg_Out() << "creating glauber object" << endl;
  //set nuclei size and position
  double radius1 = 1.12*pow(m_numNucleons1,1./3.) - 0.86*pow(m_numNucleons1,-1./3.);
  double radius2 = 1.12*pow(m_numNucleons2,1./3.) - 0.86*pow(m_numNucleons2,-1./3.);
  double maxB = ( radius1 > radius2 )? 2*radius1 : 2*radius2;
  double impact_parameter = sqrt(ran->Get())*maxB;
  m_nucleus_1_position[1] = -impact_parameter/2.;
  m_nucleus_2_position[1] = impact_parameter/2.;
  //std::ofstream nuclei_impact_par;
  //nuclei_impact_par.open("./nuclei_impact_par.txt", std::ios::app);
  //nuclei_impact_par << impact_parameter << endl;
  //nuclei_impact_par.close();
  
  //distribute nucleons positions inside the nuclei
  //std::ofstream positions1;
  //positions1.open("./position_list_1.txt", std::ios::app);
  for(int i = 0; i < m_numNucleons1; i++) {
    m_nucleus1.push_back(Glauber::distributeNucleon(radius1,m_nucleus_1_position));
    //positions1 << m_nucleus1.at(i)[1] << "\t" << m_nucleus1.at(i)[2] << "\t" << m_nucleus1.at(i)[3] << endl;
  }
  //positions1.close();
  //std::ofstream positions2;
  //positions2.open("./position_list_2.txt", std::ios::app);
  for(int i = 0; i < m_numNucleons2; i++) {
    m_nucleus2.push_back(Glauber::distributeNucleon(radius2,m_nucleus_2_position));
    //positions2 << m_nucleus2.at(i)[1] << "\t" << m_nucleus2.at(i)[2] << "\t" << m_nucleus2.at(i)[3] << endl;
  }
  //double test = XS_inelastic.GetValue(0.1);
  //msg_Out() << test << "<<<<<<<<<<<<<<<<<<<" << endl;
  //positions2.close();
}

Glauber::~Glauber() {}

void Glauber::DoCollision(std::vector<ATOOLS::Vec4D> pos_N1,std::vector<ATOOLS::Vec4D> pos_N2) {
  for(int i = 0; i < m_numNucleons1; i++) {
    for(int j = 0; j < m_numNucleons2; j++) {
      double distance = sqrt(pow(pos_N1.at(i)[1] - pos_N2.at(j)[1],2) + pow(pos_N1.at(i)[2] - pos_N2.at(j)[2],2) + pow(pos_N1.at(i)[3] - pos_N2.at(j)[3],2));
      //ouble prob_in = 
    }
  }
}

ATOOLS::Vec4D Glauber::distributeNucleon(double radius, ATOOLS::Vec4D nuc_pos) {
  double accept = false;
  double R,theta,phi, x, y, z;
  while (!accept) {
    R = 2*radius*ran->Get();
    phi = 2*M_PI*ran->Get();
    theta = acos(1-2*ran->Get());
    x = R*sin(theta)*cos(phi);
    y = R*sin(theta)*sin(phi);
    z = R*cos(theta);
    if (WoodsSaxon(sqrt(pow(x,2) + pow(y,2) + pow(z,2)),radius) > ran->Get()){
      accept = true;
    }
  }
  return Vec4D(0.,x+nuc_pos[1],y+nuc_pos[2],z+nuc_pos[3]);
}

double Glauber::WoodsSaxon(double r, double Radius) {
  return m_rho0*(1 + pow(m_w*r,2)/pow(Radius,2))/(1 + exp((r - Radius)/m_a));
}
