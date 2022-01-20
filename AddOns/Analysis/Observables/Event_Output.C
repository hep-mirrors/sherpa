#include "AddOns/Analysis/Observables/Primitive_Observable_Base.H"
#include "AddOns/Analysis/Main/Analysis_Handler.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Scoped_Settings.H"

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

    Event_Output() :
      Primitive_Observable_Base(),
      m_inlist(""),
      m_n(0.0),
      m_sum(0.0), m_sumsqr(0.0),
      m_wit(std::numeric_limits<size_t>::max())
    {
      m_splitt_flag=false;
    }

    Event_Output(Scoped_Settings s) :
      Primitive_Observable_Base(),
      m_n(0.0), m_sum(0.0), m_sumsqr(0.0),
      m_wit(std::numeric_limits<size_t>::max())
    {
      m_inlist = s["InList"].SetDefault("").Get<std::string>();
      if (m_inlist.empty()) {
        THROW(fatal_error, "You didn't specify an InList for Event_Output");
      }
      m_splitt_flag=false;
      std::string outpath(
          s["EVENT_FILE_PATH"].SetDefault(".").Get<std::string>());
      auto outputs = s["EVENT_FORMAT"]
        .SetDefault<std::string>({})
        .UseNoneReplacements()
        .GetVector<std::string>();
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
	  (outputs[i],Output_Arguments(outpath,outfile));
	if (out==NULL) {
	  if (!s_loader->LoadLibrary("Sherpa"+libname+"Output")) 
	    THROW(missing_module,"Cannot load output library Sherpa"+libname+"Output.");
	  out=Output_Base::Getter_Function::GetObject
	    (outputs[i],Output_Arguments(outpath,outfile));
	}
	if (out==NULL) THROW(fatal_error,"Cannot initialize "+outputs[i]+" output");
	m_outputs.push_back(out);
	out->Header();
      }
      if (s["FILE_SIZE"].IsCustomised()) {
        const double filesize{
          s["FILE_SIZE"].SetDefault(0.0).Get<double>() };
	if (filesize<1.0) {
	  if (filesize*rpa->gen.NumberOfEvents()>1.0)
	    m_wit=(size_t)(filesize*rpa->gen.NumberOfEvents());
	}
	else m_wit=(size_t)(filesize);
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
      if (blobs.size())
	for (Output_Vector::iterator it=m_outputs.begin(); it!=m_outputs.end(); ++it) {
	  (*it)->SetXS(p_ana->AnalysisHandler()->EventHandler()->TotalXS().Nominal(),
		       p_ana->AnalysisHandler()->EventHandler()->TotalErr().Nominal());
	  (*it)->Output((Blob_List*)&blobs);
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
      return new Event_Output();
    }

  };// end of class Event_Output

}



DECLARE_GETTER(Event_Output, "Event_Output",
	       Primitive_Observable_Base, Analysis_Key);

Primitive_Observable_Base*
ATOOLS::Getter<Primitive_Observable_Base,
               Analysis_Key,
               Event_Output>::operator()(const Analysis_Key& key) const
{
  return new Event_Output(key.m_settings);
}

void ATOOLS::Getter<Primitive_Observable_Base,
                    Analysis_Key,
                    Event_Output>::PrintInfo(
                        std::ostream &str, const size_t width) const
{
  str<<"{\n"
     <<std::setw(width+7)<<" "<<"InList: <triggeroutlist>,\n"
     <<std::setw(width+7)<<" "<<"# event output settings; cf. manual, e.g.:\n"
     <<std::setw(width+7)<<" "<<"EVENT_FORMAT: HepMC_GenEvent[<filename>],\n"
     <<std::setw(width+7)<<" "<<"FILE_SIZE: <n>\n"
     <<std::setw(width+4)<<" "<<"}";
}
