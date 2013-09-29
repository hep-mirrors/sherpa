#include "AMEGIC++/Main/Point.H"

#include "ATOOLS/Org/Message.H"

using namespace AMEGIC;
using namespace MODEL;


Point::Point(const Point& copy) { 
  extrafl = 0;
  Color   = new Color_Function;
  Lorentz = NULL;
  middle  = 0;
  nextra = 0;

  *this = copy;
} 

Point::Point(int extra) : nextra(extra)  { 
  zwf     = 0;
  propid  = 0;
  extrafl = 0;
  v       = 0;
  Color   = new Color_Function;
  Lorentz = NULL;
  middle  = 0;
  if (nextra>0) extrafl = new ATOOLS::Flavour[nextra]; 
}

Point& Point::operator=(const Point& p) {
  if (this!=&p) {
    number = p.number;
    b      = p.b;
    t      = p.t;
    zwf    = p.zwf;
    propid = p.propid;
    m      = p.m;
    fl     = p.fl;
      
    *Color = *p.Color; 
    if (Lorentz) delete Lorentz;
    Lorentz=NULL;
    if (p.Lorentz) Lorentz = p.Lorentz->GetCopy(); 
 
    if (nextra>0) delete[] extrafl;
    nextra = p.nextra;
    if (nextra>0) {
      extrafl = new ATOOLS::Flavour[nextra]; 
      for(int i=0;i<nextra;i++) extrafl[i] = p.extrafl[i];
    }
    left   = p.left;
    right  = p.right;
    middle = p.middle;
    prev  = p.prev;
    v = p.v;
    //cpl's
    cpl.clear();

    for(size_t i=0;i<p.Ncpl();i++) cpl.push_back(p.cpl[i]);
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

namespace AMEGIC {

  std::ostream &operator<<(std::ostream &str,const Point &p)
  {
    str<<p.fl<<"("<<p.b<<","<<p.number<<")";
    if (p.left) {
      str<<"[->"<<*p.left<<","<<*p.right;
      if (p.middle) str<<","<<*p.middle;
      str<<"]";
    }
    return str;
  }

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

int Point::CountKK()
{
  int KKnum=0;
  if (left) {
    KKnum+=left->CountKK();
    KKnum+=right->CountKK();
    if (middle) KKnum+=middle->CountKK();
  }
  if (fl.IsKK()) KKnum++;
  return KKnum;
}

bool Point::CountT(int & tchan,const long unsigned int & kfcode) {
  long unsigned int comp;
  if (left) {
    if (left->CountT(tchan,kfcode)) { 
      comp=left->fl.Kfcode();
      if ((comp==kfcode || kfcode==0) && !fl.Strong()) tchan++; 
      return true;
    }
    if (right->CountT(tchan,kfcode)) {
      comp=right->fl.Kfcode();
      if ((comp==kfcode || kfcode==0) && !fl.Strong()) tchan++; 
      return true;
    }
    if (middle && middle->CountT(tchan,kfcode)) {
      comp=middle->fl.Kfcode();
      if ((comp==kfcode || kfcode==0) && !fl.Strong()) tchan++; 
      return true;
    }
  }
  else if (b==-1) return true;
  return false;
}

void Point::GeneratePropID()
{
  propid=0;
  if (!left) {
    propid=(1<<number);
    return;
  }
  left->GeneratePropID();
  propid+=left->propid;
  right->GeneratePropID();
  propid+=right->propid;
  if (middle) {
    middle->GeneratePropID();
    propid+=middle->propid;    
  }
}

std::string Point::GetPropID() const
{
  return fl.IDName()+ATOOLS::ToString(propid);
}

std::ostream & operator<<(std::ostream & s, const Point & p)
{
//   s<<p;
//   return s;
  s<<" t="<<p.t<<" ";
  if ((p.left==0) && (p.right==0)) {
    s<<"EndPoint : "<<p.fl<<"("<<p.b<<")"<<std::endl;
    return s;
  }
  s<<" ["<<p.fl<<"("<<p.b<<")]"<<std::endl;
  s<<"left : ";
  s<<p.left;
  s<<"right : ";
  s<<p.right;
  if(p.middle){
    s<<" middle : ";
    s<<p.middle;
  }
  return s;
}

int Point::FindQCDOrder(int & oqcd) {
 if (!this) return oqcd;
  if (v) oqcd+=v->oqcd;
  left->FindQCDOrder(oqcd);
  right->FindQCDOrder(oqcd);
  if (middle) middle->FindQCDOrder(oqcd);
  return oqcd;
}

int Point::FindQEDOrder(int & oqed) {
if (!this) return oqed;
  if (v) oqed+=v->oew;
  left->FindQEDOrder(oqed);
  right->FindQEDOrder(oqed);
  if (middle) middle->FindQEDOrder(oqed);
  return oqed;
}

