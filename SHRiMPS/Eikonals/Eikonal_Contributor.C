#include "SHRiMPS/Eikonals/Eikonal_Contributor.H"
#include "SHRiMPS/Tools/MinBias_Parameters.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/MathTools.H"

using namespace SHRIMPS;
using namespace ATOOLS;

Eikonal_Contributor::
Eikonal_Contributor(Form_Factor * ff1,Form_Factor * ff2) :
  p_ff1(ff1), p_ff2(ff2),
  m_originalY(MBpars("originalY")), m_cutoffY(MBpars("deltaY")), 
  m_Y(m_originalY-m_cutoffY),
  m_b1min(0.), m_b2min(0.), m_b1max(p_ff1->Bmax()), m_b2max(p_ff2->Bmax()),
  m_ff1max(p_ff1->FourierTransform(m_b1min)), 
  m_ff2max(p_ff2->FourierTransform(m_b2min)), 
  m_ff1steps(-1), m_ff2steps(-1), m_ysteps(-1),
  m_b1(-1.), m_b2(-1.)
{ 
  msg_Tracking()<<"In "<<METHOD<<"("<<m_Y<<"):"<<std::endl
		<<"   FF1("<<m_b1min<<") = "<<m_ff1max<<", "
		<<"FF2("<<m_b2min<<") = "<<m_ff2max<<"."<<std::endl;
}

Eikonal_Contributor::~Eikonal_Contributor() {
  for (int i=0;i<m_ff1steps+1;i++) {
    for (int j=0;j<m_ff2steps+1;j++) {
      m_grid[i][j].clear();
    }
    m_grid[i].clear();
  }
  m_grid.clear();
}

void Eikonal_Contributor::
PrepareGrid(const int & ff1steps, const int & ff2steps)
{
  m_ff1steps = ff1steps-1;
  m_ff2steps = ff2steps-1;
  m_grid.resize(m_ff1steps+1);
  for (int i=0;i<m_ff1steps+1;i++) m_grid[i].resize(m_ff2steps+1);
  m_deltaff1 = p_ff1->FourierTransform(m_b1min)/double(m_ff1steps);  
  m_deltaff2 = p_ff2->FourierTransform(m_b2min)/double(m_ff2steps);  
  msg_Tracking()<<"In "<<METHOD<<"("<<m_ff1steps<<", "<<m_ff2steps<<") "
		<<"--> "<<m_grid.size()<<" "<<m_grid[0].size()<<std::endl
		<<" --> "<<p_ff1->FourierTransform(m_b1min)
		<<" in steps of "<<m_deltaff1
		<<" and "<<p_ff2->FourierTransform(m_b2min)
		<<" in steps of "<<m_deltaff2
		<<"."<<std::endl;
}

void Eikonal_Contributor::
InsertValues(const int & i,const int & j,const std::vector<double> & values) {
  copy(values.begin(),values.end(),back_inserter(m_grid[i][j]));
  if (m_ysteps<0) {
    m_ysteps  = m_grid[i][j].size();
    m_deltay  = 2.*m_Y/double(m_ysteps-1);
  }
}

double Eikonal_Contributor::operator()(const double & y) const 
{ 
  if (IsNan(y)) return 0.;
  double value(0.);
  if (m_b1>=m_b1max || m_b2>=m_b2max || m_b1<m_b1min || m_b2<m_b2min) {
    return 0.;
  }
  double ff1(p_ff1->FourierTransform(m_b1));
  double ff2(p_ff2->FourierTransform(m_b2));
  if (ff1<=0. || ff2<=0.) {
    //msg_Error()<<METHOD<<" ff's < 0."<<std::endl;
    return 0.;
  }
  size_t ff1bin((m_ff1max-ff1)/m_deltaff1),  ff2bin((m_ff2max-ff2)/m_deltaff2);
  double ff1low(m_ff1max-ff1bin*m_deltaff1), ff1up(m_ff1max-(ff1bin+1)*m_deltaff1);
  double ff2low(m_ff2max-ff2bin*m_deltaff2), ff2up(m_ff2max-(ff2bin+1)*m_deltaff2);
  double d1up(ff1up-ff1), d1low(ff1-ff1low), d2up(ff2up-ff2), d2low(ff2-ff2low);
  if (ff1bin>=m_grid.size()-1 || ff2bin>=m_grid[0].size()-1) {
    msg_Error()<<"Error in "<<METHOD<<"[1]("<<m_b1<<", "<<m_b2<<", "<<y<<"):"
	       <<std::endl
	       <<"   ff1 = "<<ff1<<" --> "<<ff1bin
	       <<"("<<m_grid.size()<<"/"<<m_ff1steps<<")"<<std::endl
	       <<"   ff2 = "<<ff2<<" --> "<<ff2bin
	       <<"("<<m_grid[0].size()<<"/"<<m_ff2steps<<")"<<std::endl
	       <<"   y   = "<<y
	       <<" ("<<m_grid[0][0].size()<<"/"<<m_ysteps<<")"<<std::endl
	       <<" sizes: "
	       <<m_grid[ff1bin+0][ff2bin+0].size()<<" * "
	       <<m_grid[ff1bin+0][ff2bin+1].size()<<" * "
	       <<m_grid[ff1bin+1][ff2bin+0].size()<<" * "
	       <<m_grid[ff1bin+1][ff2bin+1].size()<<"."<<std::endl;
    return 0.;
  }
  if (y<-m_Y) {
    //if (y>-m_originalY) {
    value = 
      (d1low * d2low * m_grid[ff1bin+1][ff2bin+1][0]+
       d1low * d2up  * m_grid[ff1bin+1][ff2bin+0][0]+
       d1up  * d2low * m_grid[ff1bin+0][ff2bin+1][0]+
       d1up  * d2up  * m_grid[ff1bin+0][ff2bin+0][0])/
      (m_deltaff1*m_deltaff2);
  }
  //} 
  else if (y>m_Y) {
    //else if (y<m_originalY) {
    size_t ylast(m_grid[0][0].size()-1);
    value = 
      (d1low * d2low * m_grid[ff1bin+1][ff2bin+1][ylast]+
       d1low * d2up  * m_grid[ff1bin+1][ff2bin+0][ylast]+
       d1up  * d2low * m_grid[ff1bin+0][ff2bin+1][ylast]+
       d1up  * d2up  * m_grid[ff1bin+0][ff2bin+0][ylast])/
      (m_deltaff1*m_deltaff2);
  }
  //}
  else {
    size_t ybin((y+m_Y)/m_deltay);
    double ylow(-m_Y+ybin*m_deltay), yup(-m_Y+(ybin+1)*m_deltay);
    double dyup(yup-y), dylow(y-ylow);
    if (ff1bin>=m_grid.size()-1 || ff2bin>=m_grid[0].size()-1 || 
	ybin>=m_grid[0][0].size()-1) {
      msg_Error()<<"Error in "<<METHOD<<"[2]"
		 <<"("<<m_b1<<", "<<m_b2<<", "<<y<<"):"<<std::endl
		 <<"   ff1 = "<<ff1<<" --> "<<ff1bin
		 <<"("<<m_grid.size()<<"/"<<m_ff1steps<<")"<<std::endl
		 <<"   ff2 = "<<ff2<<" --> "<<ff2bin
		 <<"("<<m_grid[0].size()<<"/"<<m_ff2steps<<")"<<std::endl
		 <<"   y   = "<<y<<" --> "<<ybin
		 <<"("<<m_grid[0][0].size()<<"/"<<m_ysteps<<")"<<std::endl;
      return 0.;
    }
    value = (d1low * d2low * dylow * m_grid[ff1bin+1][ff2bin+1][ybin+1]+
	     d1low * d2low * dyup  * m_grid[ff1bin+1][ff2bin+1][ybin+0]+
	     d1low * d2up  * dylow * m_grid[ff1bin+1][ff2bin+0][ybin+1]+
	     d1low * d2up  * dyup  * m_grid[ff1bin+1][ff2bin+0][ybin+0]+
	     d1up  * d2low * dylow * m_grid[ff1bin+0][ff2bin+1][ybin+1]+
	     d1up  * d2low * dyup  * m_grid[ff1bin+0][ff2bin+1][ybin+0]+
	     d1up  * d2up  * dylow * m_grid[ff1bin+0][ff2bin+0][ybin+1]+
	     d1up  * d2up  * dyup  * m_grid[ff1bin+0][ff2bin+0][ybin+0])/
      (m_deltaff1*m_deltaff2*m_deltay);
  }
  return value;
}
