#include "Primitive_Observable_Base.H"

#include "Primitive_Analysis.H"
#include "Shell_Tools.H"

using namespace ANALYSIS;

#define COMPILE__Getter_Function
#define OBJECT_TYPE Primitive_Observable_Base
#define PARAMETER_TYPE String_Matrix
#include "Getter_Function.C"

using namespace ATOOLS;

Primitive_Observable_Base::Primitive_Observable_Base() :
  m_type(0), m_nbins(0), m_xmin(0.), m_xmax(0.),
  m_name(std::string("noname")), m_listname(std::string("Analysed")),
  p_histo(NULL), m_nout(0), p_flavs(NULL), p_moms(NULL), 
  m_splitt_flag(true) ,p_ana(NULL), p_sel(NULL) 
{ 
  m_blobdisc = false;
}

Primitive_Observable_Base::Primitive_Observable_Base(int type,double xmin,double xmax,int nbins,
						     Selector_Base * sel) :
  m_type(type), m_nbins(nbins), m_xmin(xmin), m_xmax(xmax), 
  m_listname(std::string("Analysed")), m_splitt_flag(true), p_ana(NULL), p_sel(sel)
{ 
  p_histo = new Histogram(m_type,m_xmin,m_xmax,m_nbins);
}


Primitive_Observable_Base::Primitive_Observable_Base(const Primitive_Observable_Base & old) :
  m_type(old.m_type), m_nbins(old.m_nbins), m_xmin(old.m_xmin), m_xmax(old.m_xmax), 
  m_name(old.m_name), m_listname(old.m_listname), p_sel(old.p_sel)
{ 
  msg.Out()<<"LEGACY WARNING:  copy constructor Primitive_Observable_Base::Primitive_Observable_Base called"<<std::endl
	   <<"                 use Copy() method instead!"<<std::endl;
  if (old.p_histo) {
    p_histo = new Histogram(old.p_histo);
  }
  p_histo=NULL;
}


Primitive_Observable_Base::~Primitive_Observable_Base() 
{
  if (p_histo!=0) { delete p_histo; p_histo = 0; }
}

void Primitive_Observable_Base::SetBlobType(const std::string & btype) 
{ 
  m_blobtype = btype;
  m_blobdisc = false;
  if (btype!=std::string("")) m_blobdisc = true;
}

void Primitive_Observable_Base::Evaluate(int,const Vec4D *,const Flavour *,double, int) 
{
  msg.Error()<<"ERROR virtual function Primitive_Observable_Base::Evaluate (vecs) called "<<m_name<<std::endl;
}

void Primitive_Observable_Base::Evaluate(const Particle_List & pl,double weight,int ncount) 
{
  if (ncount>1) {
    msg.Out()<<"WARNING: "<<Name()
	     <<"::Evaluate(const Particle_List & pl,const double weight,"<<ncount<<") "<<std::endl;
    Evaluate(pl,weight,ncount);
    return;
  }
  msg.Error()<<"ERROR virutal function Primitive_Observable_Base::Evaluate (pl) called "<<m_name<<std::endl;
}

void Primitive_Observable_Base::Evaluate(const Blob_List & blobs, double value, int ncount)
{
  Particle_List * pl=p_ana->GetParticleList(m_listname);
  Evaluate(*pl,value, ncount);
}


void Primitive_Observable_Base::EndEvaluation(double scale) {
  if (p_histo) {
    p_histo->Finalize();
    if (scale!=1.) p_histo->Scale(scale);
    p_histo->Output();
  }
}

/*
void Primitive_Observable_Base::SetFlavInfo(int _nout,const Vec4D * _moms,const Flavour * _flavs) {
  m_nout = _nout; p_moms = _moms; p_flavs = _flavs;
}
*/

void Primitive_Observable_Base::Output(const std::string & pname) {
  if (p_histo) {
    int  mode_dir = 448;
    ATOOLS::MakeDir((pname).c_str(),mode_dir); 
    p_histo->Output((pname+std::string("/")+m_name).c_str());
  }
}

void Primitive_Observable_Base::SetAnalysis(Primitive_Analysis * ana) 
{
  p_ana=ana;
}

void Primitive_Observable_Base::Reset()
{
  if (p_histo) p_histo->Reset();
}

Primitive_Observable_Base & Primitive_Observable_Base::operator+=(const Primitive_Observable_Base & ob)
{
  if (p_histo) {
    (*p_histo)+=(*ob.p_histo);
  }
  else {
    msg.Out()<<"Warning in Primitive_Observable_Base::operator+= :"<<std::endl<<"   "
	     <<Name()<<" has not overloaded the operator+="<<std::endl;
  }
  return *this;
}

