#include "Primitive_Observable_Base.H"

using namespace ATOOLS;
using namespace std;

Primitive_Observable_Base::Primitive_Observable_Base() :
  type(0), nbins(0), xmin(0.), xmax(0.), histo(NULL), sel(NULL),
  nout(0), flavs(NULL), moms(NULL), name(std::string("noname")) { };


Primitive_Observable_Base::Primitive_Observable_Base(int _type,double _xmin,double _xmax,
			  int _nbins, Selector_Base * _sel) :
  type(_type), nbins(_nbins), xmin(_xmin), xmax(_xmax), sel(_sel) 
{ 
  histo = new Histogram(type,xmin,xmax,nbins);
};


Primitive_Observable_Base::Primitive_Observable_Base(Primitive_Observable_Base * old) :
  type(old->type), nbins(old->nbins), xmin(old->xmin), xmax(old->xmax), 
  sel(old->sel), name(old->name) 
{ 
  histo = new Histogram(old->histo);
}


Primitive_Observable_Base::~Primitive_Observable_Base() {
  if (histo!=0) { delete histo; histo = 0; }
}

void Primitive_Observable_Base::SetBlobType(std::string _btype) 
{ 
  m_blobtype = _btype;
  m_blobdisc = false;
  if (_btype!=std::string("")) m_blobdisc = true;
}

void Primitive_Observable_Base::Evaluate(const Blob_List & _blist,double _weight) 
{
  size_t pos;
  bool   take;
  Particle_List plist;
  plist.clear();
  for (Blob_Const_Iterator blit=_blist.begin();blit!=_blist.end();++blit) {
    take = true;
    if (m_blobdisc) {
      take = false;
      pos  = (*blit)->Type().find(m_blobtype);
      if (pos!=std::string::npos) take = true;
    }
    if (take) {
      for (int i=0;i<(*blit)->NOutP();i++) plist.push_back((*blit)->OutParticle(i));
    }
  }
  if (plist.size()>0) Evaluate(plist,_weight);
}


void Primitive_Observable_Base::EndEvaluation(double _scale) {
  histo->Finalize();
  if (_scale!=1.) histo->Scale(_scale);
  histo->Output();
}

void Primitive_Observable_Base::SetFlavInfo(int _nout,Vec4D * _moms,Flavour * _flavs) {
  nout = _nout; moms = _moms; flavs = _flavs;
}

void Primitive_Observable_Base::Output(std::string _pname) {
  int  mode_dir = 448;
  mkdir((_pname).c_str(),mode_dir); 
  histo->Output((_pname+std::string("/")+name).c_str());
}

int             Primitive_Observable_Base::Type()  { return type;  }
int             Primitive_Observable_Base::Nbins() { return nbins; }
double          Primitive_Observable_Base::Xmin()  { return xmin;  }
double          Primitive_Observable_Base::Xmax()  { return xmax;  }
std::string     Primitive_Observable_Base::Name()  { return name;  }
Histogram     * Primitive_Observable_Base::Histo() { return histo; }
Selector_Base * Primitive_Observable_Base::Sel()   { return sel;   }
