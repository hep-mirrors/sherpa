#include "Fast_Function.H"
#include "MathTools.H"
#include <list>
#include <fstream>

using namespace AMATOOLS;


Fast_Function::Fast_Function() {
  ymin=1.e99;
  ymax=-1.e99;
}

Fast_Function::Fast_Function(int size) {
  data = DList(size);
  ymin=1.e99;
  ymax=-1.e99;
}


void Fast_Function::Init(Function_Base & fun, double xmin, double xmax, int min_points) {
  Clear();

  int mode=1;
  switch (mode) {
  case 0:
    // min_points equidistant
    for (int i=0; i<min_points; ++i) {
      double x= xmin + (xmax-xmin)*double(i)/double(min_points-1);
      data.push_back(Pair(x,fun(x)));
    }
    break;
  case 1:
    // adaptive method, distributes more points in region where Interpolation is worst.
    std::list<Pair> testpoints;

    // initialize data with two points
    data.push_back(Pair(xmin,fun(xmin)));
    data.push_back(Pair(xmax,fun(xmax)));

    // initialize test points with one point
    double x=(xmin+xmax)/2.;
    double y=fun(x);
    testpoints.push_back(Pair(x,y));

    /*  // --- output ---
    std::cout<<data.size()<<" Data: "<<std::endl;
    std::cout<<(*this)<<std::endl;
    std::cout<<testpoints.size()<<" Testpoints: "<<std::endl;
    for (std::list<Pair>::iterator it=testpoints.begin(); it!=testpoints.end();++it) {
      std::cout<<(*it);
    }
    */

    for (int i=3; i<min_points; i=i+2) {

      // loop over testpoints, tries to find worst point
      double diff=0;
      std::list<Pair>::iterator it=testpoints.begin(), win;
      for (;it!=testpoints.end();++it) {
	y=(*this)(it->x);  // interpolation
	double newdiff=dabs(1.-y/it->y);
	if (newdiff>=diff) {
	  diff=newdiff;
	  win=it;
	}
      }
      // insert winner testpoint in data field
      DIter dit=Insert(win->x,win->y);

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
      /* // --- output --
      std::cout<<data.size()<<" Data: "<<std::endl;
      std::cout<<(*this)<<std::endl;
      std::cout<<testpoints.size()<<" Testpoints: "<<std::endl;
      for (std::list<Pair>::iterator it=testpoints.begin(); it!=testpoints.end();++it) {
	std::cout<<(*it);
      }
      */

    }

    // transfer all remaining testpoints to data 
    std::list<Pair>::iterator it=testpoints.begin();
    DIter dit=data.begin();
    for (;it!=testpoints.end();) {
      /*  // does not work! using our own Insert() routine
      std::cout<<" insert into fastfunc "<<Pair(it->x,it->y);
      ++dit;
      std::cout<<" before  "<<(*dit)<<std::endl;
      data.insert(dit,(*it));
      ++it;
      std::cout<<" before2 "<<(*dit)<<std::endl;
      */

      Insert(it->x,it->y);
      ++it;
    }

    // delete testpoint is done automatically at end of scope

    break;
  }
}


Fast_Function::DIter Fast_Function::Insert(double x, double y) { 
  if (y>ymax) ymax=y;
  if (y<ymin) ymin=y;

  if (data.empty()) {
    //    std::cout<<" insert in empty fastfunc "<<Pair(x,y)<<std::endl;
    data.push_back(Pair(x,y));
    DIter it =data.end();
    return --it;
  } 
  else if (data.back().x<x) {
    //    std::cout<<" append to fastfunc "<<Pair(x,y)<<std::endl;
    data.push_back(Pair(x,y));
    DIter it =data.end();
    return --it;
  }
  else {
    //    std::cout<<" insert into fastfunc "<<Pair(x,y);
    DIter it=data.begin();
    while ((*it).x < x) {++it; }

    DIter win=data.insert(it,Pair(x,y));
    //    std::cout<<" before "<<(*it)<<std::endl;
    return win;
  }
}

double Fast_Function::Invers(double y) {
  if (data.empty()) {
    std::cout<<"ERROR: Fast_Function::Invers() called for empty function!!!"<<std::endl;
    return 0;
  }
  if (data.size()==1) {
    if (data.front().y==y) {
      return data.front().x;
    }
    else {
      std::cout<<"ERROR: Fast_Function::Invers() called for almost empty function!!!"<<std::endl;
      return 0;
    }
  }
  // at least two elements
  DIter it=data.begin();
  for (;;) {
    double y1=it->y;
    ++it;
    if (it==data.end()) break;
    double y2=it->y;
    if (((y1<y)&&(y<=y2))||((y2<y)&&(y<=y1))) break;
  }
  if (it==data.end()) {
    // x is bigger than or smaller than all stored values
    std::cout<<"ERROR: Fast_Function::Invers() "<<std::endl;
    std::cout<<" given y="<<y<<" is not in range "<<YRange()<<std::endl;
    return 0; 
  }
  return LinInterInv(it,y);
}

double Fast_Function::operator()(double x) {
  if (data.empty()) {
    std::cout<<"ERROR: Fast_Function::opertor() called for empty function!!!"<<std::endl;
    return 0;
  }
  if (data.size()==1) {
    if (data.front().x==x) {
      return data.front().y;
    }
    else {
      std::cout<<"ERROR: Fast_Function::opertor() called for almost empty function!!!"<<std::endl;
      return 0;
    }
  }
  // at least two elements
  DIter it=data.begin();
  while (((*it).x < x)&&(it!=data.end())) {++it; }

  if (it==data.end()) {
    // x is bigger than all stored values
    --it;
  }
  return LinInter(it,x);
}


void Fast_Function::WriteOut(char * name) {
  std::ofstream to(name);
  to.precision(10);
  
  for (DIter it=data.begin();it!=data.end();++it) 
    to<<it->x<<"    "<<it->y<<std::endl;

  std::cout<<"File "<<name<<" "<<data.size()<<" entries saved."<<std::endl;
}


bool Fast_Function::ReadIn(char * name) {
  std::ifstream from(name);
  if (!from) return 0; // fail

  Clear();

  double x,y;
  for(;from;) {
    from>>x>>y;
    if ((data.empty())||(x!=data.back().x))
      data.push_back(Pair(x,y));
  }
  from.close();

  std::cout<<"File "<<name<<" read."<<std::endl;
  return 1; // success
}



double Fast_Function::LinInter(DIter & it, double x) {
  double x1=it->x;
  double y1=it->y;

  if (it!=data.begin()) --it; else ++it;
  double x2=it->x;
  double y2=it->y;

  return y1+(y2-y1)*(x-x1)/(x2-x1);
}

double Fast_Function::LinInterInv(DIter & it, double y) {
  double x1=it->x;
  double y1=it->y;

  --it; 
  double x2=it->x;
  double y2=it->y;

  return x1 +(x2-x1)*(y-y1)/(y2-y1);
}


std::ostream & AMATOOLS::operator<<(std::ostream & s, const Fast_Function & ff) {
  s<<"----------------"<<std::endl;
  for (Fast_Function::DList::const_iterator it=ff.data.begin();it!=ff.data.end();++it) 
    s<<(*it);
  return s;
}


std::ostream & AMATOOLS::operator<<(std::ostream & s, const  Fast_Function::Pair & p) {
  s<<'('<<p.x<<','<<p.y<<')'<<std::endl;
  return s;
}


std::ostream & AMATOOLS::operator<<(std::ostream & s, const Intervall & i) {
  s<<'['<<i.minval<<','<<i.maxval<<']'<<std::endl;
  return s;
}
