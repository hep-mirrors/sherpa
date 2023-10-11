#include "SHRiMPS/Glauber/Glauber.H"
#include "SHRiMPS/Event_Generation/Soft_Diffractive_Event_Generator.H"
#include "SHRiMPS/Tools/Special_Functions.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"
#include "SHRiMPS/Cross_Sections/Sigma_D.H"
#include "SHRiMPS/Cross_Sections/Sigma_DD.H"
#include "SHRiMPS/Cross_Sections/Sigma_Inelastic.H"
#include "SHRiMPS/Cross_Sections/Sigma_Elastic.H"
#include "SHRiMPS/Cross_Sections/Sigma_Total.H"
#include "SHRiMPS/Cross_Sections/Cross_Sections.H"
#include "SHRiMPS/Cross_Sections/Sigma_Base.H"
#include "SHRiMPS/Tools/MinBias_Parameters.H"

using namespace SHRIMPS;
using namespace ATOOLS;
using namespace std;

Glauber::Glauber(Cross_Sections * xsecs, int N1, int N2) :
  m_numNucleons1(N1), m_numNucleons2(N2),
  m_radius1( (1.12*pow(N1,1./3.) - 0.86*pow(N1,-1./3.))),
  m_radius2((1.12*pow(N2,1./3.) - 0.86*pow(N2,-1./3.))),
  m_rho0(0.17), m_w(0.), m_a(0.54), 
  m_nucleus_1_position(ATOOLS::Vec4D(0.,0.,0.,0.)),
  m_nucleus_2_position(ATOOLS::Vec4D(0.,0.,0.,0.)),
  m_type_abs(3), m_type_SD1(0), m_type_SD2(1), m_type_DD(2), m_type_elas(4),
  m_pp(true), m_pA(false), m_AA(false), m_projState(0)
{
  double maxB = ( m_radius1 > m_radius2 )? 2*m_radius1 : 2*m_radius2;
  double impact_parameter = sqrt(ran->Get())*maxB;
  m_nucleus_1_position[1] = -impact_parameter/2.;
  m_nucleus_2_position[1] = impact_parameter/2.;
  if((N1 == 1 && N2 > 1) || (N1 > 1 && N2 == 1)) {m_pA = true; m_pp = false;}
  else if(N1 > 1 && N2 > 1) {m_AA = true; m_pp = false;}

  for(int i = 0; i < m_numNucleons1; i++) {
    m_nucleus1.push_back(Glauber::distributeNucleon(m_radius1,m_nucleus_1_position));
  }
  for(int i = 0; i < m_numNucleons2; i++) {
    m_nucleus2.push_back(Glauber::distributeNucleon(m_radius2,m_nucleus_2_position));
  }
  Glauber::DoCollision(xsecs,m_nucleus1,m_nucleus2);
  //Glauber::SaveXSs(xsecs,0.,maxB,500);
  //Glauber::SaveNucleonPositions();
}

Glauber::~Glauber() {}

void Glauber::SaveNucleonPositions() {
  std::ofstream pos_nuc_1;
  pos_nuc_1.open("./nucleons_1.txt", std::ios::app);
  for(int i = 0; i < m_numNucleons1; i++) {
    double x = m_nucleus1[i][1] - m_nucleus_1_position[1];
    double y = m_nucleus1[i][2] - m_nucleus_1_position[2];
    double z = m_nucleus1[i][3] - m_nucleus_1_position[3];
    pos_nuc_1 << x << "\t" << y << "\t" << z << "\t" << endl;
  }
  pos_nuc_1.close();

  std::ofstream pos_nuc_2;
  pos_nuc_2.open("./nucleons_2.txt", std::ios::app);
  for(int i = 0; i < m_numNucleons2; i++) {
    double x = m_nucleus2[i][1] - m_nucleus_2_position[1];
    double y = m_nucleus2[i][2] - m_nucleus_2_position[2];
    double z = m_nucleus2[i][3] - m_nucleus_2_position[3];
    pos_nuc_2 << x << "\t" << y << "\t" << z << "\t" << endl;
  }
  pos_nuc_2.close();
}

void Glauber::SaveXSs(Cross_Sections * xsecs,double bmin, double bmax,int num) {
  double b_increment = (bmax - bmin)/num;
  double current_b = bmin;
  std::ofstream cross_sections, cross_sections_ch, abs_xs_ch;
  cross_sections.open("./xsecs.dat");//, std::ios::app);
  cross_sections_ch.open("./xsecs_channels_QE.dat");//, std::ios::app);
  abs_xs_ch.open("./xsecs_channels_abs.dat");//, std::ios::app);
  std::vector<std::vector<Omega_ik *> > * p_eikonals(MBpars.GetEikonals());
  while (current_b < bmax) {
    double b_GeVmin1 = current_b/0.197;
    double inel(xsecs->GetSigmaInelastic()->GetCombinedValue(b_GeVmin1));
    double el(xsecs->GetSigmaElastic()->GetCombinedValue(b_GeVmin1));
    double QE(xsecs->GetSigmaD()->GetCombinedValue(b_GeVmin1));
    double SD0(xsecs->GetSigmaD()->GetCombinedValueSD0(b_GeVmin1));
    double SD1(xsecs->GetSigmaD()->GetCombinedValueSD1(b_GeVmin1));
    double DD(xsecs->GetSigmaD()->GetCombinedValueDD(b_GeVmin1));
    cross_sections << current_b << "\t" << inel << "\t" << SD0 << "\t" << SD1 << "\t" << DD << "\t" << el << "\t" << QE << endl;
    
    std::vector< std::vector<double> > QEchannels;
    std::vector< std::vector<double> > ABSchannels;
    for (size_t i=0;i<p_eikonals->size();i++) {
      std::vector<double> QErow;
      std::vector<double> ABSrow;
      for (size_t j=0;j<(*p_eikonals)[i].size();j++) {
        QErow.push_back(xsecs->GetSigmaD()->GetValuePerChannel(i,j,b_GeVmin1));
        ABSrow.push_back(xsecs->GetSigmaInelastic()->GetValuePerChannel(i,j,b_GeVmin1));
      }
      QEchannels.push_back(QErow);
      ABSchannels.push_back(ABSrow);
    }
    //msg_Out() << QEchannels.size() << " SIZE QE CHANNESL! " << std::endl;
    //msg_Out() << QEchannels[0].size() << " SIZE QE CHANNESL! 0" << std::endl;
    //msg_Out() << QEchannels[1].size() << " SIZE QE CHANNESL! 1" << std::endl;
    cross_sections_ch << current_b << "\t";
    abs_xs_ch << current_b << "\t";
    //msg_Out() << "B!!!!!!!!!!!!!!!" << endl;
    for(size_t i = 0; i < QEchannels.size(); i++) {
      for(size_t j = 0; j < QEchannels[i].size(); j++) {
        cross_sections_ch << QEchannels[i][j] << "\t";
        abs_xs_ch << ABSchannels[i][j] << "\t";
      }
    }
    for(size_t i = 0; i < QEchannels.size(); i++) {
      double QEfixedCol(xsecs->GetSigmaD()->GetValuePerChannel(-1,i,b_GeVmin1));
      double QEfixedRow(xsecs->GetSigmaD()->GetValuePerChannel(i,-1,b_GeVmin1));
      double ABSfixedCol(xsecs->GetSigmaInelastic()->GetValuePerChannel(-1,i,b_GeVmin1));
      double ABSfixedRow(xsecs->GetSigmaInelastic()->GetValuePerChannel(i,-1,b_GeVmin1));
      cross_sections_ch << QEfixedCol << "\t" << QEfixedRow << "\t";
      abs_xs_ch << ABSfixedCol << "\t" << ABSfixedRow << "\t";
    }
    cross_sections_ch << std::endl;
    abs_xs_ch << std::endl;
    current_b = current_b + b_increment;
  }
  cross_sections.close();
  cross_sections_ch.close();
  abs_xs_ch.close();
}

int Glauber::FixProjectileState() {
  std::vector<std::vector<Omega_ik *> > * p_eikonals(MBpars.GetEikonals());
  double disc = ran->Get();
  for (size_t j=0;j<(*p_eikonals)[0].size();j++) {
    double a0j_sqr = ( (*p_eikonals)[0][j]->Prefactor() ) / sqrt((*p_eikonals)[0][0]->Prefactor());
    disc -= a0j_sqr;
    if (disc <= 1.e-12) return j;
  }
  return -1;
}

void Glauber::DoCollision(Cross_Sections * xsecs, std::vector<ATOOLS::Vec4D> pos_N1,std::vector<ATOOLS::Vec4D> pos_N2) {
  //std::vector<std::pair<std::pair<int,int>,int>> interactions;
  for(int i = 0; i < m_numNucleons1; i++) {
    for(int j = 0; j < m_numNucleons2; j++) {
      double distance = sqrt(sqr(pos_N1.at(i)[1] - pos_N2.at(j)[1]) + sqr(pos_N1.at(i)[2] - pos_N2.at(j)[2])); // fm
      double inel(xsecs->GetSigmaInelastic()->GetCombinedValue(distance/0.197));
      double el(xsecs->GetSigmaElastic()->GetCombinedValue(distance/0.197));
      double QE(xsecs->GetSigmaD()->GetCombinedValue(distance/0.197));
      double SD0(xsecs->GetSigmaD()->GetCombinedValueSD0(distance/0.197));
      double SD1(xsecs->GetSigmaD()->GetCombinedValueSD1(distance/0.197));
      double DD(xsecs->GetSigmaD()->GetCombinedValueDD(distance/0.197));

      double prob_SD0(0.);  
      double prob_SD1(0.);  
      double prob_DD(0.);  
      double prob_el(0.);
      if (QE>0.) {
        prob_SD0 = SD0/QE;  
        prob_SD1 = SD1/QE;  
        prob_DD = DD/QE;  
        prob_el = el/QE;
        //msg_Out() << prob_SD0 + prob_SD1 + prob_DD + prob_el << endl;
      }
      if(m_pA) {
        if (i == 0 && j ==0) {m_projState=FixProjectileState();}
        else {
          if(m_numNucleons1==1) {
            QE = xsecs->GetSigmaD()->GetValuePerChannel(m_projState,-1,distance/0.197);
            inel = xsecs->GetSigmaInelastic()->GetValuePerChannel(m_projState,-1,distance/0.197);
          }
          else {
            QE = xsecs->GetSigmaD()->GetValuePerChannel(-1,m_projState,distance/0.197);
            inel = xsecs->GetSigmaInelastic()->GetValuePerChannel(-1,m_projState,distance/0.197);
          }
        }
      }
      double prob_inel(inel), prob_QE(QE);
      double rand1(ran->Get()), rand2(ran->Get());
      std::pair<int,int> current_pair;
      current_pair.first = i;
      current_pair.second = j;
      //std::ofstream interactions_file;
      //interactions_file.open("./interactions_postfix.dat", std::ios::app);
      
      bool inelastic_happened = false;
      bool QE_happened = false;
      if (rand1 < prob_inel) {
        inelastic_happened = true;
        std::pair<std::pair<int,int>,int> interaction1;
        interaction1.first = current_pair;
        interaction1.second = m_type_abs; //interaction 1 = inelastic
        m_list_of_interactions.push_back(interaction1);
        //interactions_file << interaction1.second << "\t" << distance << endl;
      }
      if (rand2 < prob_QE) {
        QE_happened = true;
        double rand3(ran->Get());
        std::pair<std::pair<int,int>,int> interaction2;
        interaction2.first = current_pair;
        if (rand3 < prob_SD0) interaction2.second = m_type_SD1; //interaction 2 = single diff in 0
        else if (rand3 < prob_SD0 + prob_SD1) interaction2.second = m_type_SD2; //interaction 3 = single diff in 1
        else if (rand3 < prob_SD0 + prob_SD1 + prob_DD) interaction2.second = m_type_DD; //interaction 4 = double diff
        else interaction2.second = m_type_elas; //interaction 5 = elastic
        m_list_of_interactions.push_back(interaction2);
        //interactions_file << interaction2.second << "\t" << distance << endl;
      }
      if (!inelastic_happened && !QE_happened) {
        std::pair<std::pair<int,int>,int> interaction3;
        interaction3.first = current_pair;
        interaction3.second = -1;
      //interactions_file << -1 << "\t" << distance << endl;
      }
      //interactions_file.close();
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
