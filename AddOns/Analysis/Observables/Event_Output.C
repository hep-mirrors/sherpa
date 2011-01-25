#include "AddOns/Analysis/Observables/Primitive_Observable_Base.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Exception.H"

#include "SHERPA/Tools/Input_Output_Handler.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"

using namespace ANALYSIS;
using namespace ATOOLS;
using namespace SHERPA;

namespace ANALYSIS {

  class Event_Output : public Primitive_Observable_Base {
    Input_Output_Handler* p_io;
    std::string m_inlist;
    double m_n, m_sum, m_sumsqr;

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

    Event_Output(Input_Output_Handler* io, const std::string& inlist) :
      Primitive_Observable_Base(), p_io(io), m_inlist(inlist),
      m_n(0.0), m_sum(0.0), m_sumsqr(0.0)
    {
      m_splitt_flag=false;
    }

    ~Event_Output()
    {
      if (p_io) delete p_io; p_io=NULL;
    }


    void Evaluate(const ATOOLS::Blob_List & blobs, double weight, double ncount)
    {
      if (!p_io) return;
      Particle_List * pl=p_ana->GetParticleList(m_inlist);
      m_n+=ncount;
      if (pl->empty()) return;
      m_sum+=weight;
      m_sumsqr+=sqr(weight);
      Blob_List* blobsptr=const_cast<Blob_List*>(&blobs);
      p_io->OutputToFormat(blobsptr);
    }


    void EndEvaluation(double scale=1.)
    {
      if (m_sum==0.0) return;
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



DECLARE_GETTER(Event_Output_Getter,"Event_Output",
	       Primitive_Observable_Base,Argument_Matrix);

Primitive_Observable_Base *
Event_Output_Getter::operator()(const Argument_Matrix &parameters) const
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
  return new Event_Output(new Input_Output_Handler(&reader), inlist);
}

void Event_Output_Getter::PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"{\n"
     <<std::setw(width+7)<<" "<<"InList  <triggeroutlist>\n"
     <<std::setw(width+7)<<" "<<"# event output settings cf. manual, e.g.:\n"
     <<std::setw(width+7)<<" "<<"HEPMC2_GENEVENT_OUTPUT <filename>\n"
     <<std::setw(width+7)<<" "<<"FILE_SIZE <n>\n"
     <<std::setw(width+7)<<" "<<"EVT_FILE_PATH <path>\n"
     <<std::setw(width+4)<<" "<<"}";
}
