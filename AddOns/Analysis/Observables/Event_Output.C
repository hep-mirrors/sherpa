#include "AddOns/Analysis/Observables/Primitive_Observable_Base.H"
#include "AddOns/Analysis/Main/Analysis_Handler.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/Run_Parameter.H"

#include "SHERPA/Tools/Output_Base.H"
#include "SHERPA/Single_Events/Event_Handler.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"

#include <limits>

using namespace ANALYSIS;
using namespace ATOOLS;
using namespace SHERPA;

namespace ANALYSIS {

  class Event_Output : public Primitive_Observable_Base {
    Output_Vector m_outputs;
    std::string m_inlist;
    double m_n, m_sum, m_sumsqr;
    size_t m_wit;

    inline double TotalXS() const { return m_sum/m_n; }
    inline double TotalVar() const
    { return (m_sumsqr-m_sum*m_sum/m_n)/(m_n-1); }
    inline double TotalErr() const {
      if (m_n==1) return TotalXS();
      else if (ATOOLS::IsEqual
	       (m_sumsqr*m_n,m_sum*m_sum,1.0e-6)) return 0.0;
      else return sqrt(TotalVar()/m_n);
    }

  public:

    Event_Output(Data_Reader *reader, const std::string& inlist) :
      Primitive_Observable_Base(), m_inlist(inlist),
      m_n(0.0), m_sum(0.0), m_sumsqr(0.0),
      m_wit(std::numeric_limits<size_t>::max())
    {
      m_splitt_flag=false;
      std::string outpath=reader->GetValue<std::string>("EVT_FILE_PATH",".");
      std::string format=reader->GetValue<std::string>("EVENT_FORMAT","None");
      std::vector<std::string> outputs;
      Data_Reader readline(",",";","#","");
      readline.SetString(format);
      readline.VectorFromString(outputs);
      for (size_t i=0; i<outputs.size(); ++i) {
	if (outputs[i]=="None") continue;
	std::string outfile;
	size_t bpos(outputs[i].find('[')), epos(outputs[i].rfind(']'));
	if (bpos!=std::string::npos && epos!=std::string::npos) {
	  outfile=outputs[i].substr(bpos+1,epos-bpos-1);
	  outputs[i]=outputs[i].substr(0,bpos);
	}
	std::string libname(outputs[i]);
	if (libname.find('_')) libname=libname.substr(0,libname.find('_'));
	Output_Base* out=Output_Base::Getter_Function::GetObject
	  (outputs[i],Output_Arguments(outpath,outfile,reader));
	if (out==NULL) {
	  if (!s_loader->LoadLibrary("Sherpa"+libname+"Output")) 
	    THROW(missing_module,"Cannot load output library Sherpa"+libname+"Output.");
	  out=Output_Base::Getter_Function::GetObject
	    (outputs[i],Output_Arguments(outpath,outfile,reader));
	}
	if (out==NULL) THROW(fatal_error,"Cannot initialize "+outputs[i]+" output");
	m_outputs.push_back(out);
	out->Header();
      }
      double wit;
      if (reader->ReadFromFile(wit,"FILE_SIZE")) {
	if (wit<1.0) {
	  if (wit*rpa->gen.NumberOfEvents()>1.0)
	    m_wit=(size_t)(wit*rpa->gen.NumberOfEvents());
	}
	else m_wit=(size_t)(wit);
	msg_Info()<<METHOD<<"(): Set output interval "<<m_wit<<" events.\n";
      }
    }

    ~Event_Output()
    {
      while (m_outputs.size()>0) {
	delete m_outputs.back();
	m_outputs.pop_back();
      }
    }


    void Evaluate(const ATOOLS::Blob_List & blobs, double weight, double ncount)
    {
      if (m_outputs.empty()) return;
      Particle_List * pl=p_ana->GetParticleList(m_inlist);
      m_n+=ncount;
      if (pl->empty()) return;
      m_sum+=weight;
      m_sumsqr+=sqr(weight);
      Blob_List* blobsptr=const_cast<Blob_List*>(&blobs);
      if (blobs.size())
	for (Output_Vector::iterator it=m_outputs.begin(); it!=m_outputs.end(); ++it) {
	  (*it)->SetXS(p_ana->AnalysisHandler()->EventHandler()->TotalXS(),
		       p_ana->AnalysisHandler()->EventHandler()->TotalErr());
	  (*it)->Output((Blob_List*)&blobs,weight);
	}
      if (rpa->gen.NumberOfGeneratedEvents()>0 &&
	  rpa->gen.NumberOfGeneratedEvents()%m_wit==0 &&
	  rpa->gen.NumberOfGeneratedEvents()<rpa->gen.NumberOfEvents()) 
	for (Output_Vector::iterator it=m_outputs.begin(); 
	     it!=m_outputs.end(); ++it)
	  (*it)->ChangeFile();
    }


    void EndEvaluation(double scale=1.)
    {
      if (m_sum==0.0) return;
      for (Output_Vector::iterator it=m_outputs.begin();
	   it!=m_outputs.end(); ++it) (*it)->Footer();
      PRINT_FUNC("");
      double xs(TotalXS()), err(TotalErr());
      msg_Info()<<om::bold<<"Triggered XS"<<om::reset<<" is "
                <<om::blue<<om::bold<<xs<<" pb"<<om::reset<<" +- ( "
                <<om::red<<err<<" pb = "<<((int(err/xs*10000))/100.0)
                <<" %"<<om::reset<<" )";
    }


    Primitive_Observable_Base * Copy() const
    {
      // don't duplicate event output
      return new Event_Output(NULL, "");
    }

  };// end of class Event_Output

}



DECLARE_GETTER(Event_Output,"Event_Output",
	       Primitive_Observable_Base,Argument_Matrix);

Primitive_Observable_Base *ATOOLS::Getter<Primitive_Observable_Base,Argument_Matrix,Event_Output>::operator()(const Argument_Matrix &parameters) const
{
  Data_Reader reader(" ",";","!","=");
  reader.AddComment("#");
  reader.AddWordSeparator("\t");
  std::string inlist="";
  for (size_t i=0; i<parameters.size(); ++i) {
    std::string line = "";
    for (size_t j=0; j<parameters[i].size(); ++j) {
      line += parameters[i][j]+" ";
    }
    reader.AddFileContent(line);
    if (parameters[i][0]=="InList" && parameters[i].size()>1) inlist=parameters[i][1];
  }

  if (inlist=="") {
    THROW(fatal_error, "You didn't specify an InList for Event_Output");
  }
  return new Event_Output(&reader, inlist);
}

void ATOOLS::Getter<Primitive_Observable_Base,Argument_Matrix,Event_Output>::PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"{\n"
     <<std::setw(width+7)<<" "<<"InList  <triggeroutlist>\n"
     <<std::setw(width+7)<<" "<<"# event output settings cf. manual, e.g.:\n"
     <<std::setw(width+7)<<" "<<"HEPMC2_GENEVENT_OUTPUT <filename>\n"
     <<std::setw(width+7)<<" "<<"FILE_SIZE <n>\n"
     <<std::setw(width+7)<<" "<<"EVT_FILE_PATH <path>\n"
     <<std::setw(width+4)<<" "<<"}";
}
