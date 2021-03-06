#include "SHRiMPS/Eikonals/Grid.H"
#include "ATOOLS/Org/Message.H"

using namespace SHRIMPS;
using namespace ATOOLS;

Grid::Grid(const double & Ymax) : m_Ymax(Ymax) {}

Grid::~Grid() {
  for (int i=0;i<m_ff1steps+1;i++) {
    for (int j=0;j<m_ff2steps+1;j++) {
      m_grid[i][j].clear();
    }
    m_grid[i].clear();
  }
  m_grid.clear();
}

void Grid::Initialise(const size_t & ff1steps,const size_t & ff2steps,
		      const double & ff1max,const double & ff2max) 
{
  m_ff1steps = ff1steps-1;
  m_ff2steps = ff2steps-1;
  m_grid.resize(m_ff1steps+1);
  for (int i=0;i<m_ff1steps+1;i++) m_grid[i].resize(m_ff2steps+1);
  m_ff1max   = ff1max;
  m_ff2max   = ff2max;
  m_deltaff1 = m_ff1max/double(m_ff1steps);  
  m_deltaff2 = m_ff2max/double(m_ff2steps);  
}

void Grid::InsertValues(const size_t & i,const size_t & j,
			const std::vector<double> & values) {
  //msg_Out()<<METHOD<<"("<<i<<", "<<j<<" -> "<<values.size()<<") : "
  //	   <<values[0]<<" "<<values[values.size()/2]<<" "
  //	   <<values[values.size()-1]<<".\n";
  m_grid[i][j].resize(values.size());
  copy(values.begin(),values.end(),m_grid[i][j].begin());
  m_deltay = 2.*m_Ymax/double(m_grid[i][j].size()-1);
}

double Grid::
operator()(const double & ff1,const double & ff2,const double & y) {
  if (ff1==0. || ff2==0.) return 0.;
  if (!FixBins(ff1,ff2))  return 0.;
  if (y<=-m_Ymax) return ValueAtLowerYEdge();
  if (y>=m_Ymax)  return ValueAtUpperYEdge();
  return Value(y);
}

bool Grid::FixBins(const double & ff1,const double & ff2) {
  m_ff1bin = (m_ff1max-ff1)/m_deltaff1;
  m_ff2bin = (m_ff2max-ff2)/m_deltaff2;
  if (m_ff1bin>=m_grid.size()-1 || m_ff2bin>=m_grid[0].size()-1) {
    msg_Error()<<"Error in "<<METHOD<<"("<<ff1<<", "<<ff2<<"):\n"
	       <<"   "<<m_ff1bin<<"/"<<m_ff2bin<<" from "
	       <<m_ff1max<<"/"<<m_ff2max<<" and "
	       <<m_deltaff1<<"/"<<m_deltaff2<<" vs. sizes "
	       <<m_grid.size()<<"/"<<m_grid[0].size()<<".\n";
    return false;
  }
  m_ff1low = m_ff1max-m_ff1bin*m_deltaff1; 
  m_ff1up  = m_ff1max-(m_ff1bin+1)*m_deltaff1;
  m_ff2low = m_ff2max-m_ff2bin*m_deltaff2;
  m_ff2up  = m_ff2max-(m_ff2bin+1)*m_deltaff2;
  m_d1up   = m_ff1up-ff1;
  m_d1low  = ff1-m_ff1low;
  m_d2up   = m_ff2up-ff2;
  m_d2low  = ff2-m_ff2low;
  return true;
}

double Grid::ValueAtLowerYEdge() {
  return
    (m_d1low * m_d2low * m_grid[m_ff1bin+1][m_ff2bin+1][0]+
     m_d1low * m_d2up  * m_grid[m_ff1bin+1][m_ff2bin+0][0]+
     m_d1up  * m_d2low * m_grid[m_ff1bin+0][m_ff2bin+1][0]+
     m_d1up  * m_d2up  * m_grid[m_ff1bin+0][m_ff2bin+0][0])/
    (m_deltaff1*m_deltaff2);
}

double Grid::ValueAtUpperYEdge() {
  size_t ylast(m_grid[0][0].size()-1);
  return
    (m_d1low * m_d2low * m_grid[m_ff1bin+1][m_ff2bin+1][ylast]+
     m_d1low * m_d2up  * m_grid[m_ff1bin+1][m_ff2bin+0][ylast]+
     m_d1up  * m_d2low * m_grid[m_ff1bin+0][m_ff2bin+1][ylast]+
     m_d1up  * m_d2up  * m_grid[m_ff1bin+0][m_ff2bin+0][ylast])/
    (m_deltaff1*m_deltaff2);
}

double Grid::Value(const double & y) {
  size_t ybin((y+m_Ymax)/m_deltay);
  double ylow(-m_Ymax+ybin*m_deltay), yup(-m_Ymax+(ybin+1)*m_deltay);
  double dyup(yup-y), dylow(y-ylow);
  if (m_ff1bin>=m_grid.size()-1 || m_ff2bin>=m_grid[0].size()-1 || 
      ybin>=m_grid[0][0].size()-1) {
    msg_Error()<<"Error in "<<METHOD<<".\n";
    return 0.;
  }
  return 
    (m_d1low * m_d2low * dylow * m_grid[m_ff1bin+1][m_ff2bin+1][ybin+1]+
     m_d1low * m_d2low * dyup  * m_grid[m_ff1bin+1][m_ff2bin+1][ybin+0]+
     m_d1low * m_d2up  * dylow * m_grid[m_ff1bin+1][m_ff2bin+0][ybin+1]+
     m_d1low * m_d2up  * dyup  * m_grid[m_ff1bin+1][m_ff2bin+0][ybin+0]+
     m_d1up  * m_d2low * dylow * m_grid[m_ff1bin+0][m_ff2bin+1][ybin+1]+
     m_d1up  * m_d2low * dyup  * m_grid[m_ff1bin+0][m_ff2bin+1][ybin+0]+
     m_d1up  * m_d2up  * dylow * m_grid[m_ff1bin+0][m_ff2bin+0][ybin+1]+
     m_d1up  * m_d2up  * dyup  * m_grid[m_ff1bin+0][m_ff2bin+0][ybin+0])/
    (m_deltaff1*m_deltaff2*m_deltay);
}
