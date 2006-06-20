#include "Primitive_Observable_Base.H"

#include "Primitive_Analysis.H"
#include "Histogram.H"
#include "CXXFLAGS.H"
#ifdef USING__ROOT
#include "Root_Histogram.H"
#endif
#include "Shell_Tools.H"

using namespace ANALYSIS;

#define COMPILE__Getter_Function
#define OBJECT_TYPE Primitive_Observable_Base
#define PARAMETER_TYPE String_Matrix
#include "Getter_Function.C"

using namespace ATOOLS;

int ANALYSIS::HistogramType(const std::string &scale)
{
  if (scale=="Log")    return 10;
  if (scale=="Lin")    return 0;
  if (scale=="LogErr") return 11;
  if (scale=="LinErr") return 1;
  if (scale=="LogPS")  return 12;
  if (scale=="LinPS")  return 2;
  msg.Error()<<"ERROR in ANALYSIS::HistogramType:"<<std::endl
	     <<"   Do not know this scale option : "<<scale<<"."<<std::endl
	     <<"   Pretend, it's linear ('Lin')."<<std::endl;
  return 0;
}

Primitive_Observable_Base::Primitive_Observable_Base() :
  m_pobtype(POBType::Unknown),
  m_name(std::string("noname")), m_listname(std::string("Analysed")),
  p_histo(NULL), m_nout(0), p_flavs(NULL), p_moms(NULL), 
  m_splitt_flag(true) ,p_ana(NULL), p_sel(NULL), m_copied(false)
{ 
  m_blobdisc = false;
}

Primitive_Observable_Base::Primitive_Observable_Base(const int type,
						     const double xmin,const double xmax,
						     const int nbins) :
  m_pobtype(POBType::ASCII),
  m_type(type), m_xmin(xmin), m_xmax(xmax),
  m_listname(std::string("Analysed")), m_splitt_flag(true), 
  p_ana(NULL), p_sel(NULL), m_copied(false)
{ 
  p_histo = new Histogram(type,xmin,xmax,nbins);
}


Primitive_Observable_Base::Primitive_Observable_Base(const std::string infile,
						     const std::string outfile,
						     std::vector<std::string> & datanames) :
  m_pobtype(POBType::ROOT),
  m_type(100), 
  m_listname(std::string("Analysed")), m_splitt_flag(true), 
  p_ana(NULL), p_sel(NULL), m_copied(false)
{
#ifdef USING__ROOT
  p_histo = new Root_Histogram(infile,outfile,datanames);
  m_xmin  = p_histo->Xmin(); 
  m_xmax  = p_histo->Xmax(); 
#else
    msg.Error()<<"ERROR in Primitive_Observable_Base::Primitive_Observable_Base:"<<std::endl
	       <<"   Asked for root-histogram without Root being enabled."<<std::endl
	       <<"   Reconfigure with '--enable-root' and run again."<<std::endl;
    abort();
#endif
}


Primitive_Observable_Base::Primitive_Observable_Base(Histogram_Base * histo) :
  m_type(histo->Type()), m_xmin(histo->Xmin()), m_xmax(histo->Xmax()),
  m_listname(std::string("Analysed")), m_splitt_flag(true), 
  p_ana(NULL), p_sel(NULL), m_copied(false)
{ 
  if (histo->Type()<100) {
    p_histo   = new Histogram(static_cast<Histogram *>(histo));
    m_pobtype = POBType::ASCII;
  }
  else {
#ifdef USING__ROOT
    p_histo   = new Root_Histogram(static_cast<Root_Histogram *>(histo));
    m_pobtype = POBType::ROOT;
#else
    msg.Error()<<"ERROR in Primitive_Observable_Base::Primitive_Observable_Base:"<<std::endl
	       <<"   Asked for root-histogram (type = "<<histo->Type()
	       <<") without Root being enabled."<<std::endl
	       <<"   Reconfigure with '--enable-root' and run again."<<std::endl;
    abort();
#endif
  }
}

Primitive_Observable_Base::Primitive_Observable_Base(const Primitive_Observable_Base & old) :
  m_pobtype(old.m_pobtype),
  m_name(old.m_name), m_listname(old.m_listname), 
  p_sel(old.p_sel), m_copied(false)
{ 
  msg.Out()<<"LEGACY WARNING:  "
	   <<"copy constructor Primitive_Observable_Base::Primitive_Observable_Base called"<<std::endl
	   <<"                 use Copy() method instead!"<<std::endl;
  if (old.p_histo) {
    if (m_pobtype==POBType::ASCII) {
      p_histo = new Histogram(static_cast<Histogram*>(old.p_histo));
      return;
    }
    else if (m_pobtype==POBType::ROOT) {
#ifdef USING__ROOT
      p_histo = new Root_Histogram(static_cast<Root_Histogram*>(old.p_histo));
      return;
#else
      msg.Error()<<"ERROR in copy constructor of Primitive_Observable_Base:"<<std::endl
		 <<"   Asked for root-histogram (type = "<<old.p_histo->Type()
		 <<") without Root being enabled."<<std::endl
		 <<"   Reconfigure with '--enable-root' and run again."<<std::endl;
      abort();
#endif
    }
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
  //std::cout<<METHOD<<" "<<m_name<<std::endl;
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
  //std::cout<<METHOD<<":"<<pname<<" + "<<m_name<<std::endl;
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
    if (ob.m_pobtype==POBType::ASCII) (*p_histo)+=(*(static_cast<Histogram *>(ob.p_histo)));
    else if (ob.m_pobtype==POBType::ROOT) {
#ifdef USING__ROOT
      (*p_histo)+=(*(static_cast<Root_Histogram *>(ob.p_histo)));
#else
      msg.Error()<<"ERROR in Primitive_Observable_Base::operator+="<<std::endl
		 <<"   Asked for root-histogram (type = "<<ob.p_histo->Type()
		 <<") without Root being enabled."<<std::endl
		 <<"   Reconfigure with '--enable-root' and run again."<<std::endl;
      abort();
#endif    
    }
  }
  else {
    msg.Out()<<"Warning in Primitive_Observable_Base::operator+= :"<<std::endl<<"   "
	     <<Name()<<" has not overloaded the operator+="<<std::endl;
  }
  return *this;
}

