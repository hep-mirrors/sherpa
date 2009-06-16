#include "ATOOLS/Math/Fast_Function.H"
#include "ATOOLS/Math/MathTools.H"
#include <list>
#include <fstream>

using namespace ATOOLS;


Fast_Function::Fast_Function() {
  m_ymin = 1.e99;
  m_ymax = -1.e99;
}

Fast_Function::Fast_Function(int size) {
  m_data = Data_List(size);
  m_ymin =  1.e99;
  m_ymax = -1.e99;
}


void Fast_Function::Init(Function_Base & fun, double xmin, double xmax, int min_points) {
  Clear();

  int mode=1;
  switch (mode) {
  case 0:
    // min_points equidistant
    for (int i=0; i<min_points; ++i) {
      double x= xmin + (xmax-xmin)*double(i)/double(min_points-1);
      m_data.push_back(Pair(x,fun(x)));
    }
    break;
  case 1:
    // adaptive method, distributes more points in region where Interpolation is worst.
    std::list<Pair> testpoints;

    // initialize data with two points
    m_data.push_back(Pair(xmin,fun(xmin)));
    m_data.push_back(Pair(xmax,fun(xmax)));

    // initialize test points with one point
    double x=(xmin+xmax)/2.;
    double y=fun(x);
    testpoints.push_back(Pair(x,y));

    for (int i=3; i<min_points; i=i+2) {

      // loop over testpoints, tries to find worst point
      double diff=0.;
      std::list<Pair>::iterator it=testpoints.begin(), win=testpoints.begin();
      for (;it!=testpoints.end();++it) {
	y=(*this)(it->x);  // interpolation
	double newdiff=dabs(1.-y/it->y);
	if (newdiff>=diff) {
	  diff=newdiff;
	  win=it;
	}
      }
      // insert winner testpoint in data field
      Data_Iter dit=Insert(win->x,win->y);

      // and generates two new test points as neighbors to the last one
      --dit;
      x=dit->x;
      ++dit;

      x=(x+dit->x)/2.;
      y=fun(x);
      testpoints.insert(win,Pair(x,y));

      x=dit->x;
      ++dit;

      x=(x+dit->x)/2.;
      y=fun(x);
      testpoints.insert(win,Pair(x,y));
      
      // delete winner testpoint in testpoints
      testpoints.erase(win);
    }

    // transfer all remaining testpoints to data 
    std::list<Pair>::iterator it=testpoints.begin();
    for (;it!=testpoints.end();) {
      Insert(it->x,it->y);
      ++it;
    }

    // delete testpoint is done automatically at end of scope

    break;
  }
}


Fast_Function::Data_Iter Fast_Function::Insert(double x, double y) { 
  if (y>m_ymax) m_ymax = y;
  if (y<m_ymin) m_ymin = y;

  if (m_data.empty()) {
    m_data.push_back(Pair(x,y));
    Data_Iter it =m_data.end();
    return --it;
  } 
  else if (m_data.back().x<x) {
    m_data.push_back(Pair(x,y));
    Data_Iter it =m_data.end();
    return --it;
  }
  else {
    Data_Iter it=m_data.begin();
    while ((*it).x < x) {++it; }

    Data_Iter win=m_data.insert(it,Pair(x,y));
    return win;
  }
}

int Fast_Function::Nesting(int a, int b, double x) {
  if (b <=a+1) return a;
  int c=(a+b)/2;
  if (m_data[c].x<=x) 
    return Nesting(c,b,x);
  else 
    return Nesting(a,c,x);
}

double Fast_Function::Invers(double y) {
  if (m_data.empty()) {
    std::cerr<<"ERROR: Fast_Function::Invers() called for empty function!!!"<<std::endl;
    return 0;
  }
  if (m_data.size()==1) {
    if (m_data.front().y==y) {
      return m_data.front().x;
    }
    else {
      std::cerr<<"ERROR: Fast_Function::Invers() called for almost empty function!!!"<<std::endl;
      return 0;
    }
  }

  // at least two elements
  Data_Iter it=m_data.begin();
  for (;;) {
    double y1=it->y;
    ++it;
    if (it==m_data.end()) break;
    double y2=it->y;
    if (((y1<y)&&(y<=y2))||((y2<y)&&(y<=y1))) break;
  }
  if (it==m_data.end()) {
    // x is bigger than or smaller than all stored values
    std::cout<<"ERROR: Fast_Function::Invers() "<<std::endl
	     <<" given y="<<y<<" is not in range "<<YRange()<<std::endl;
    return 0; 
  }
  return LinInterInv(it,y);
}

double Fast_Function::operator()(double x) {
  if (m_data.empty()) {
    std::cout<<"ERROR: Fast_Function::opertor() called for empty function!!!"<<std::endl;
    return 0;
  }
  if (m_data.size()==1) {
    if (m_data.front().x==x) {
      return m_data.front().y;
    }
    else {
      std::cout<<"ERROR: Fast_Function::opertor() called for almost empty function!!!"<<std::endl;
      return 0;
    }
  }
  // at least two elements
  // new
  if (m_data.front().x>=x) return LinInter(0,x);
  if (m_data.back().x<=x) return LinInter(m_data.size()-1,x);
  return LinInter(Nesting(0,m_data.size()-1,x),x);
  
  /* 
  //old
  Data_Iter it=m_data.begin();
  while (((*it).x < x)&&(it!=m_data.end())) {++it; }

  if (it==m_data.end()) {
    // x is bigger than all stored values
    --it;
  }
  return LinInter(it,x);
  */
}


void Fast_Function::WriteOut(const char * name) {
  std::ofstream to(name);
  to.precision(10);
  
  for (Data_Iter it=m_data.begin();it!=m_data.end();++it) 
    to<<it->x<<"    "<<it->y<<std::endl;

}


bool Fast_Function::ReadIn(const char * name) {
  std::ifstream from(name);
  if (!from) return 0; // fail

  Clear();

  double x,y;
  for(;from;) {
    from>>x>>y;
    if ((m_data.empty())||(x!=m_data.back().x))
      m_data.push_back(Pair(x,y));
  }
  from.close();

  return 1; // success
}


double Fast_Function::LinInter(int i, double x) 
{
  double x1 = m_data[i].x;
  double y1 = m_data[i].y;

  if (i<(int)m_data.size()-1) ++i; else --i;
  double x2 = m_data[i].x;
  double y2 = m_data[i].y;

  return y1+(y2-y1)*(x-x1)/(x2-x1);
}


double Fast_Function::LinInter(Data_Iter & it, double x) {
  double x1 = it->x;
  double y1 = it->y;

  if (it!=m_data.begin()) --it; else ++it;
  double x2 = it->x;
  double y2 = it->y;

  return y1+(y2-y1)*(x-x1)/(x2-x1);
}

double Fast_Function::LinInterInv(Data_Iter & it, double y) {
  double x1 = it->x;
  double y1 = it->y;

  --it; 
  double x2 = it->x;
  double y2 = it->y;

  return x1 + (x2-x1)*(y-y1)/(y2-y1);
}


std::ostream & ATOOLS::operator<<(std::ostream & s, const Fast_Function & ff) {
  s<<"----------------"<<std::endl;
  for (Fast_Function::Data_List::const_iterator it=ff.m_data.begin();it!=ff.m_data.end();++it) 
    s<<(*it);
  return s;
}


std::ostream & ATOOLS::operator<<(std::ostream & s, const  Fast_Function::Pair & p) {
  s<<'('<<p.x<<','<<p.y<<')'<<std::endl;
  return s;
}


std::ostream & ATOOLS::operator<<(std::ostream & s, const Intervall & i) {
  s<<'['<<i.m_minval<<','<<i.m_maxval<<']'<<std::endl;
  return s;
}
