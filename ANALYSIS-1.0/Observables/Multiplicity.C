#include "Multiplicity.H"
#include "Message.H"
#include "MyStrStream.H"


using namespace ANALYSIS;

template <class Class>
Primitive_Observable_Base *const GetObservable(const String_Matrix &parameters)
{									
  //std::cout<<"Getter for multi"<<std::endl;
  Primitive_Observable_Base * obs(NULL);
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<4) return NULL;
    std::string list=parameters[0].size()>4?parameters[0][4]:"Analysed";
    obs = new Class(HistogramType(parameters[0][3]),
		    ATOOLS::ToType<double>(parameters[0][0]),
		    ATOOLS::ToType<double>(parameters[0][1]),
		    ATOOLS::ToType<int>(parameters[0][2]),list);
    return obs;
  }
  else {
    //std::cout<<"Getter for multi: root"<<std::endl;
    std::string inputfile(""),outputfile(""),list("Analysed"), title("Multiplicity"),
      legend("On");
    double xmin(-1.), ymin(-1.), xmax(-1.), ymax(-1.);
    std::vector<std::string> datanames;
    for (size_t i=0;i<parameters.size();++i) {
      if (parameters[i][0]=="INPUTFILE")       inputfile  = parameters[i][1];
      else if (parameters[i][0]=="DATANAME") {
	for (int j=1;j<parameters[i].size();j++) datanames.push_back(parameters[i][j]);
      }
      else if (parameters[i][0]=="TITLE") {
	title="";
	for (int j=1;j<parameters[i].size();j++) title+=parameters[i][j]+" ";
      }
      else if (parameters[i][0]=="LEGENDPOSITION") {
	if (parameters[i].size()!=5) legend="Off";
	xmin=ATOOLS::ToType<double>(parameters[i][1]);
	ymin=ATOOLS::ToType<double>(parameters[i][2]);
	xmax=ATOOLS::ToType<double>(parameters[i][3]);
	ymax=ATOOLS::ToType<double>(parameters[i][4]);
      }
      else if (parameters[i][0]=="LEGEND")     legend     = parameters[i][1];
      else if (parameters[i][0]=="LIST")       list       = parameters[i][1];
      else if (parameters[i][0]=="OUTPUTFILE") outputfile = parameters[i][1];
    }
    if (datanames.empty()) {
      ATOOLS::msg.Error()<<"Potential ERROR in GetObservable for Multis:"<<std::endl 
			 <<"   No datanames specified, use default 'Data' instead"
			 <<" and hope for the best."<<std::endl;
      datanames.push_back("Data");
    }
    if (inputfile!=("")) obs = new Class(inputfile,datanames,list,outputfile);
    obs->SetTitle(title);
    obs->SetLegend(legend,xmin,ymin,xmax,ymax);
    return obs;
  }

  ATOOLS::msg.Error()<<"ERROR in GetObservable for Multis:"<<std::endl
		     <<"   Did not understand the analyses input:"<<std::endl;
  for (size_t i=0;i<parameters.size();++i) {
    for (size_t j=0;j<parameters[i].size();++j) ATOOLS::msg.Error()<<parameters[i][j]<<" ";
    ATOOLS::msg.Error()<<std::endl;
  }						  
  return NULL;
}									

#define DEFINE_GETTER_METHOD(CLASS,NAME)				\
  Primitive_Observable_Base *						\
  NAME::operator()(const String_Matrix &parameters) const		\
  {   return GetObservable<CLASS>(parameters); }

#define DEFINE_PRINT_METHOD(NAME)					\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"min max bins Lin|LinErr|Log|LogErr [list]"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,String_Matrix);	\
  DEFINE_GETTER_METHOD(CLASS,NAME)					\
  DEFINE_PRINT_METHOD(NAME)

#include "MathTools.H"
#include "Particle_Qualifier.H"

using namespace ATOOLS;

DEFINE_OBSERVABLE_GETTER(Multiplicity,Multiplicity_Getter,"Multi")




Multiplicity::Multiplicity(const std::string & infilename,
			   std::vector<std::string> & datanames,
			   const std::string & listname,
			   const std::string & outfilename) :
  Primitive_Observable_Base(infilename,outfilename,datanames)
{
  //std::cout<<METHOD<<" :"<<infilename<<" / "<<datanames[0]<<" ("<<datanames.size()<<")"
  //	   <<listname<<"/"<<outfilename<<std::endl;
  m_listname = listname;
  if (outfilename=="") m_name = listname+"_multi";
                  else m_name = outfilename;
}


Multiplicity::Multiplicity(int type,double xmin,double xmax,int nbins,
			   const std::string & listname) :
  Primitive_Observable_Base(type,xmin,xmax,nbins)
{
  if (listname=="") {
    msg.Error()<<"ERROR in Multiplcity::Multilicity("<<type<<", ... ,"<<listname<<") : "<<std::endl
	       <<"   will abort."<<std::endl;
    abort();
  }
  m_listname = listname;
  m_name     = listname+"_multi.dat";
}

Multiplicity::Multiplicity(const Multiplicity * old) :
  Primitive_Observable_Base(*old)
{
  m_listname = old->m_listname;
  m_name     = old->m_name;
}

Primitive_Observable_Base * Multiplicity::Copy() const {
  //std::cout<<"In "<<METHOD<<std::endl;
  return new Multiplicity(this);
}

void Multiplicity::Evaluate(const ATOOLS::Particle_List & pl,
			    double weight, int ncount)
{
  //std::cout<<"In "<<METHOD<<" "<<pl.size()<<std::endl;
  p_histo->Insert(pl.size(),weight,ncount); 
  //std::cout<<"out of "<<METHOD<<std::endl;
}





