#include "Primitive_Observable_Base.H"

#include "CXXFLAGS.H"
#include "Primitive_Analysis.H"
#include "Histogram.H"
#include "Run_Parameter.H"
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
using namespace std;

int ANALYSIS::HistogramType(const string &scale)
{
  if (scale=="Log")    return 10;
  if (scale=="Lin")    return 0;
  if (scale=="LogErr") return 11;
  if (scale=="LinErr") return 1;
  if (scale=="LogPS")  return 12;
  if (scale=="LinPS")  return 2;
  msg.Error()<<"ERROR in ANALYSIS::HistogramType:"<<endl
	     <<"   Do not know this scale option : "<<scale<<"."<<endl
	     <<"   Pretend, it's linear ('Lin')."<<endl;
  return 0;
}

Primitive_Observable_Base::Primitive_Observable_Base() :
  m_pobtype(POBType::Unknown),
  m_name(string("noname")), m_listname(string("Analysed")),
  p_histo(NULL), /*p_flavs(NULL), p_moms(NULL),*/
  m_blobtype(""), m_blobdisc(false),
  m_splitt_flag(true) ,p_ana(NULL), p_sel(NULL), m_copied(false)
{ 
  m_blobdisc = false;
}

Primitive_Observable_Base::Primitive_Observable_Base(const int type,
						     const double xmin,const double xmax,
						     const int nbins) :
  m_pobtype(POBType::ASCII),
  m_type(type), m_xmin(xmin), m_xmax(xmax),
  m_listname(string("Analysed")), 
  p_histo(NULL), /*p_flavs(NULL), p_moms(NULL),*/
  m_blobtype(""), m_blobdisc(false),
  m_splitt_flag(true), 
  p_ana(NULL), p_sel(NULL), m_copied(false)
{ 
  p_histo = new Histogram(type,xmin,xmax,nbins);
}

Primitive_Observable_Base::Primitive_Observable_Base(const String_Matrix & parameters) :
  m_pobtype(POBType::Unknown),
  m_listname(""),
  p_histo(NULL), /*p_flavs(NULL), p_moms(NULL),*/
  m_blobtype(""), m_blobdisc(false),
  m_splitt_flag(true), 
  p_ana(NULL), p_sel(NULL), m_copied(false)
{
  if(parameters.size()==1) {
    m_pobtype=POBType::ASCII;
    m_blobtype=""; m_blobdisc=false;
    m_splitt_flag=true;
    p_ana=NULL; p_sel=NULL; m_copied=false;
    return;
  }
  else {
    m_pobtype=POBType::ROOT;
    m_type=100;
#ifdef USING__ROOT
    string inputfile("SherpaDefault"), title("SherpaDefault");
    m_name="SherpaDefault";
    vector<string> * datanames = new vector<string>;
    bool legend(false), logy(false);
    double lxmax(-1.), lymax(-1.),
      lxmin(-1.), lymin(-1.);
  
    for (size_t i=0;i<parameters.size();++i) {
      if (parameters[i][0]=="INPUTFILE")       inputfile  = parameters[i][1];
      else if (parameters[i][0]=="OUTPUTFILE") m_name     = parameters[i][1];
      else if (parameters[i][0]=="LIST")       m_listname = parameters[i][1];
      else if (parameters[i][0]=="TITLE") {
        title="";
        for (size_t j=1;j<parameters[i].size();j++) title+=parameters[i][j]+" ";
        int findpos = title.find("\\#");
        while ( findpos != -1) {
          title.replace( findpos, 2 , "#");
          findpos = title.find("\\#",findpos+1);
        }
      }
      else if (parameters[i][0]=="DATANAME") {
        for (size_t j=1;j<parameters[i].size();j++) datanames->push_back(parameters[i][j]);
      }
      else if (parameters[i][0]=="SETLOGARITHMIC") {
        for (size_t j=1;j<parameters[i].size();j++) logy = (parameters[i][j]=="On");
      }
      else if (parameters[i][0]=="LEGENDPOSITION") {
        if (parameters[i].size()!=5) legend="Off";
        lxmin=ATOOLS::ToType<double>(parameters[i][1]);
        lymin=ATOOLS::ToType<double>(parameters[i][2]);
        lxmax=ATOOLS::ToType<double>(parameters[i][3]);
        lymax=ATOOLS::ToType<double>(parameters[i][4]);
      }
      else if (parameters[i][0]=="LEGEND")     legend     = (parameters[i][1]=="On");
    }
    if (datanames->empty()) {
      ATOOLS::msg.Error()<<"Potential ERROR in GetObservable for Multis:"<<endl 
                        <<"   No datanames specified, use default 'Data' instead"
                          <<" and hope for the best."<<endl;
      datanames->push_back("Data");
    }
    string histofile       = "Histos/"+inputfile;
    if (system(("test -f "+histofile).c_str())) // if file does not exist in $PWD
      histofile=rpa.gen.Variable("SHERPA_SHARE_PATH")+"/Histos/"+inputfile;
    p_histo = new Root_Histogram(histofile,m_name,datanames);
    p_histo->SetTitle(title);
    p_histo->SetLegend(legend,lxmin,lymin,lxmax,lymax);
    if (logy) p_histo->SetLogY();
  
    m_xmin  = p_histo->Xmin(); 
    m_xmax  = p_histo->Xmax(); 
#else
    msg.Error()<<"ERROR in Primitive_Observable_Base::Primitive_Observable_Base:"<<endl
              <<"   Asked for root-histogram without Root being enabled."<<endl
              <<"   Reconfigure with '--enable-root' and run again."<<endl;
    abort();
#endif
  }
  std::cout<<METHOD<<": m_name="<<m_name<<std::endl;
}

// Primitive_Observable_Base::Primitive_Observable_Base(Histogram_Base * histo) :
//   m_type(histo->Type()), m_xmin(histo->Xmin()), m_xmax(histo->Xmax()),
//   m_listname(string("Analysed")), m_splitt_flag(true),
//   p_ana(NULL), p_sel(NULL), m_copied(false)
// {
//   if (histo->Type()<100) {
//     p_histo   = new Histogram(static_cast<Histogram *>(histo));
//     m_pobtype = POBType::ASCII;
//   }
//   else {
// #ifdef USING__ROOT
//     p_histo   = new Root_Histogram(static_cast<Root_Histogram *>(histo));
//     m_pobtype = POBType::ROOT;
// #else
//     msg.Error()<<"ERROR in Primitive_Observable_Base::Primitive_Observable_Base:"<<endl
// 	       <<"   Asked for root-histogram (type = "<<histo->Type()
// 	       <<") without Root being enabled."<<endl
// 	       <<"   Reconfigure with '--enable-root' and run again."<<endl;
//     abort();
// #endif
//   }
// }

Primitive_Observable_Base::Primitive_Observable_Base(const Primitive_Observable_Base & old) :
  m_pobtype(old.m_pobtype),
  m_name(old.m_name), m_listname(old.m_listname), 
  p_sel(old.p_sel), m_copied(false)
{ 
//   msg.Out()<<"LEGACY WARNING:  "
// 	   <<"copy constructor Primitive_Observable_Base::Primitive_Observable_Base called"<<endl
// 	   <<"                 use Copy() method instead!"<<endl;
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
      msg.Error()<<"ERROR in copy constructor of Primitive_Observable_Base:"<<endl
		 <<"   Asked for root-histogram (type = "<<old.p_histo->Type()
		 <<") without Root being enabled."<<endl
		 <<"   Reconfigure with '--enable-root' and run again."<<endl;
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

void Primitive_Observable_Base::SetBlobType(const string & btype) 
{ 
  m_blobtype = btype;
  m_blobdisc = false;
  if (btype!=string("")) m_blobdisc = true;
}

void Primitive_Observable_Base::Evaluate(int,const Vec4D *,const Flavour *,double, int) 
{
  msg.Error()<<"ERROR virtual function Primitive_Observable_Base::Evaluate (vecs) called "
	     <<m_name<<endl;
}

void Primitive_Observable_Base::Evaluate(const Particle_List & pl,double weight,int ncount) 
{
  if (ncount>1) {
    msg.Out()<<"WARNING: "<<Name()
	     <<"::Evaluate(const Particle_List & pl,const double weight,"<<ncount<<") "<<endl;
    Evaluate(pl,weight,ncount);
    return;
  }
  msg.Error()<<"ERROR virutal function Primitive_Observable_Base::Evaluate (pl) called "
	     <<m_name<<endl;
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

void Primitive_Observable_Base::Output(const string & pname) {
  if (p_histo) {
    int  mode_dir = 448;
    ATOOLS::MakeDir((pname).c_str(),mode_dir);
    p_histo->Output((pname+string("/")+m_name).c_str());
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
      msg.Error()<<"ERROR in Primitive_Observable_Base::operator+="<<endl
		 <<"   Asked for root-histogram (type = "<<ob.p_histo->Type()
		 <<") without Root being enabled."<<endl
		 <<"   Reconfigure with '--enable-root' and run again."<<endl;
      abort();
#endif    
    }
  }
  else {
    msg.Out()<<"Warning in Primitive_Observable_Base::operator+= :"<<endl<<"   "
	     <<Name()<<" has not overloaded the operator+="<<endl;
  }
  return *this;
}
