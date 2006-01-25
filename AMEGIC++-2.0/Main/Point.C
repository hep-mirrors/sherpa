#include "Point.H"

using namespace AMEGIC;


Point::Point(const Point& copy) { 
  extrafl = 0;
  Color   = new Color_Function;
  Lorentz = new Lorentz_Function;
  middle  = 0;
  ncpl = 4;
  cpl = new Complex[4];
  nextra = 0;

  *this = copy;
} 

Point::Point(int extra) : nextra(extra)  { 
  extrafl = 0;
  v       = 0;
  Color   = new Color_Function;
  Lorentz = new Lorentz_Function;
  middle  = 0;
  ncpl = 4;
  cpl = new Complex[4];
  if (nextra>0) extrafl = new ATOOLS::Flavour[nextra]; 
}

Point& Point::operator=(const Point& p) {
  if (this!=&p) {
    number = p.number;
    b      = p.b;
    t      = p.t;
    m      = p.m;
    fl     = p.fl;
      
    *Color = *p.Color; 
    *Lorentz = *p.Lorentz; 
 
    if (nextra>0) delete[] extrafl;
    nextra = p.nextra;
    if (nextra>0) {
      extrafl = new ATOOLS::Flavour[nextra]; 
      for(short int i=0;i<nextra;i++) extrafl[i] = p.extrafl[i];
    }
    left   = p.left;
    right  = p.right;
    middle = p.middle;
    prev  = p.prev;
    v = p.v;
    //cpl's
    if (ncpl!=p.ncpl) {
      delete[] cpl;
      ncpl = p.ncpl;
      cpl = new Complex[ncpl];
    } 

    for(short int i=0;i<ncpl;i++) cpl[i] = p.cpl[i];
  }
  return *this;
}

void Point::ResetExternalNumbers(int os)
{
  if (number<100 && b==1) number+=os;
  if (left) {
    left->ResetExternalNumbers(os);
    right->ResetExternalNumbers(os);
    if (middle) middle->ResetExternalNumbers(os);
  }
}

void Point::ResetFlag()
{
  t = 0;
  if (left) {
    left->ResetFlag();
    right->ResetFlag();
    if (middle) middle->ResetFlag();
  }
}

void Point::ResetProps()
{
  int st = 0;
  ResetProps(st);
}

void Point::ResetProps(int &st)
{
  if (b==2) b=1;
  if (left) {
    if (number!=0){
      st++;
      number = st;
      if (fl.IsFermion()) number+=100;
      if (fl.IsBoson())   number+=200;
    }
    left->ResetProps(st);
    right->ResetProps(st);
    if (middle) middle->ResetProps(st);
  }
}

Point* Point::CopyList(Point* p)
{
  *this = p[0];
  Point* nx = this;
  if (p[0].left) {
    left = nx + 1;
    right = left->CopyList(p[0].left) + 1;
    nx = right->CopyList(p[0].right);
    if (p[0].middle) {
      middle = nx + 1;
      nx = middle->CopyList(p[0].middle);
    }
  }
  return nx;
}

void Point::Print()
{
  std::cout<<" "<<fl<<"("<<b<<","<<number<<")";
  if (left) {
    std::cout<<"[->";
    left->Print();
    right->Print();
    if (middle) middle->Print();
    std::cout<<"]"<<std::endl;
  }
}
