#include "Primitive_Observable_Base.H"

#include "Primitive_Analysis.H"
#include "Shell_Tools.H"
#include "MyStrStream.H"

using namespace ANALYSIS;

#define COMPILE__Getter_Function
#define OBJECT_TYPE Primitive_Observable_Base
#define PARAMETER_TYPE Argument_Matrix
#include "Getter_Function.C"

using namespace ATOOLS;

int ANALYSIS::HistogramType(const std::string &scale)
{
  if (scale.length()>0 && scale[0]>47 && scale[0]<58) 
    return ToType<int>(scale); 
  if (scale=="Log") return 10;
  if (scale=="LogErr") return 11;
  if (scale=="LinErr") return 1;
  if (scale=="LogPS") return 12;
  if (scale=="LinPS") return 2;
  return 0;
}

Primitive_Observable_Base::Primitive_Observable_Base() :
  m_type(0), m_nbins(0), m_xmin(0.), m_xmax(0.),
  m_listname(std::string("Analysed")),
  p_histo(NULL), m_nout(0), p_flavs(NULL), p_moms(NULL), 
  m_splitt_flag(true) ,p_sel(NULL), m_copied(false)
{ 
  m_name="noname";
  m_blobdisc = false;
  m_obs=true;
}

Primitive_Observable_Base::Primitive_Observable_Base(int type,double xmin,double xmax,int nbins,
						     Selector_Base * sel) :
  m_type(type), m_nbins(nbins), m_xmin(xmin), m_xmax(xmax),
  m_listname(std::string("Analysed")), m_splitt_flag(true), 
  p_sel(sel), m_copied(false)
{ 
  p_histo = new Histogram(m_type,m_xmin,m_xmax,m_nbins);
  m_obs=true;
}


Primitive_Observable_Base::Primitive_Observable_Base(const Primitive_Observable_Base & old) :
  m_type(old.m_type), m_nbins(old.m_nbins), m_xmin(old.m_xmin), m_xmax(old.m_xmax),
  m_listname(old.m_listname), p_sel(old.p_sel), m_copied(false)
{ 
  m_name=old.m_name;
  msg.Out()<<"LEGACY WARNING:  copy constructor Primitive_Observable_Base::Primitive_Observable_Base called"<<std::endl
	   <<"                 use Copy() method instead!"<<std::endl;
  if (old.p_histo) {
    p_histo = new Histogram(old.p_histo);
  }
  p_histo=NULL;
  m_obs=true;
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

Analysis_Object & Primitive_Observable_Base::operator+=(const Analysis_Object & ob)
{
  return (*this)+=*static_cast<const Primitive_Observable_Base*>(&ob);
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

Analysis_Object *Primitive_Observable_Base::GetCopy() const
{
  return Copy();
}
