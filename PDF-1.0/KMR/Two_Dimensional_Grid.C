#include "Two_Dimensional_Grid.H"

using namespace PDF;

Two_Dimensional_Grid::Two_Dimensional_Grid(const size_t inipoints):
  p_root(new ATOOLS::Node<double>(0.0,true)),
  p_xl(NULL), p_xr(NULL), p_yll(NULL), p_ylr(NULL), p_yrl(NULL), p_yrr(NULL),
  m_intx(std::vector<double>(inipoints,0.0)), m_inty(0.0)
{
  for (size_t i=0;i<inipoints;++i) {
    (*p_root)()[i] = new ATOOLS::Node<double>(0.0,true);
    for (size_t j=0;j<inipoints;++j) {
      (*(*p_root)()[i])()[j] = new ATOOLS::Node<double>(0.0,true);
      (*(*(*p_root)()[i])()[j])()[0] = new ATOOLS::Node<double>(0.0,false);
    }
  }
}

Two_Dimensional_Grid::~Two_Dimensional_Grid()
{
  delete p_root;
}

size_t Two_Dimensional_Grid::FindX(const double x) 
{
  size_t l=0, r=p_root->size()-1, i=(r+l)/2;
  double xi=(*(*p_root)()[i])[0];
  while (r-l>1) {
    if (x<xi) r=i;
    else l=i;
    i=(r+l)/2;
    xi=(*(*p_root)()[i])[0];
  }
  p_xl=(*p_root)()[l];
  p_xr=(*p_root)()[r];
  if (x<xi) --i;
  return i;
}

size_t Two_Dimensional_Grid::FindY(const double y) 
{
  size_t l=0, r=p_root->size()-1, i=(r+l)/2;
  double yi=(*(*p_xl)()[i])[0];
  while (r-l>1) {
    if (y<yi) r=i;
    else l=i;
    i=(r+l)/2;
    yi=(*(*p_xl)()[i])[0];
  }
  p_yll=(*p_xl)()[l];
  p_ylr=(*p_xl)()[r];
  p_yrl=(*p_xr)()[l];
  p_yrr=(*p_xr)()[r];
  if (y<yi) --i;
  return i;
}

void Two_Dimensional_Grid::InterpolateX(const double x) const
{
}

void Two_Dimensional_Grid::InterpolateY(const double x) const
{
}


double Two_Dimensional_Grid::Value(const double x,const double y) 
{
  FindX(x);
  FindY(y);
  double x1=(*p_xl)[0], x2=(*p_xr)[0];
  double y1=(*p_yll)[0], y2=(*p_ylr)[0];
  double zx1y1=(*(*p_yll)()[0])[0], zx1y2=(*(*p_ylr)()[0])[0];
  double zx2y1=(*(*p_yrl)()[0])[0], zx2y2=(*(*p_yrr)()[0])[0];
  double fx=(x-(*p_xl)[0])/((*p_xr)[0]-(*p_xl)[0]);
  double fy=(y-(*p_yll)[0])/((*p_ylr)[0]-(*p_yll)[0]);
}

double Two_Dimensional_Grid::SetValue(const double x,const double y,
				      const double z)
{
  FindX(x);
  FindY(y);
}

