#ifndef Grid_Creator_C
#define Grid_Creator_C

#include "Phase_Space_Integrator.H"

namespace AMISIC {

  template <class Argument_Type,class Result_Type>
  Grid_Creator<Argument_Type,Result_Type>::Grid_Creator(GridHandlerVector _p_gridhandler,
							EXTRAXS::Simple_XS *_p_processes):
    GridCreatorBaseType(),
    p_gridhandler(_p_gridhandler),
    p_processes(_p_processes),
    m_storemax(true),
    m_xsextension("_xs.dat"),
    m_maxextension("_max.dat")
  {
    m_maxpoints[0]=2500000;
    if (p_gridhandler.size()<2) {
      ATOOLS::msg.Out()<<"Grid_Creator::Grid_Creator("<<&_p_gridhandler<<","<<_p_processes<<"): "
		       <<"Constructor called with one grid handler only. "<<std::endl
		       <<"   Maxima will not be stored."<<std::endl;
      m_storemax=false;
      if (p_gridhandler.size()==0) {
	ATOOLS::msg.Error()<<"Grid_Creator_Base::Grid_Creator_Base("<<&_p_gridhandler<<"): "
			   <<"Grid handler is not initialized! "<<std::endl
			   <<"   Run cannot continue."<<std::endl;
	abort();
      }
    }
    if (p_processes==NULL) {
      ATOOLS::msg.Error()<<"Grid_Creator::Grid_Creator("<<&_p_gridhandler<<","<<_p_processes<<"): "
			 <<"Process handler is not initialized! "<<std::endl
			 <<"   Run cannot continue."<<std::endl;
      abort();
    }
    m_maxpoints[1]=PHASIC::Phase_Space_Integrator::MaxPoints();
    PHASIC::Phase_Space_Integrator::SetMaxPoints(m_maxpoints[0]);
    p_gridhandler[0]->Grid()->SetMonotony(GridFunctionType::None);
    p_gridhandler[1]->Grid()->SetMonotony(GridFunctionType::None);
  }
  
  template <class Argument_Type,class Result_Type>
  Grid_Creator<Argument_Type,Result_Type>::~Grid_Creator()
  {
    PHASIC::Phase_Space_Integrator::SetMaxPoints(m_maxpoints[1]);
  }

  template <class Argument_Type,class Result_Type>
  bool Grid_Creator<Argument_Type,Result_Type>::
  ReadInArguments(std::string tempifile,std::string tempipath)
  {
    bool success=true;
    if (tempipath!=ATOOLS::nullstring) this->SetInputPath(tempipath);
    if (tempifile!=ATOOLS::nullstring) this->SetInputFile(tempifile);
    if (this->InputFile()==ATOOLS::nullstring) {
      ATOOLS::msg.Error()<<"Grid_Creator_Base::Arguments("<<tempifile<<"): "
			 <<"No input file specified! Abort."<<std::endl;
      abort();
    }
    success=success&&ReadSingleArguments(p_gridhandler[0]);
    success=success&&ReadSingleArguments(p_gridhandler[1]);
    return success;
  }
  
  template <class Argument_Type,class Result_Type>
  bool Grid_Creator<Argument_Type,Result_Type>::InitializeCalculation()
  {
    p_xaxis=p_gridhandler[0]->Grid()->XAxis();
    p_yaxis=p_gridhandler[0]->Grid()->YAxis();
    m_criterion=ATOOLS::Variable::TypeToSelectorID(p_xaxis->Variable().Type());
    m_initialdata=p_processes->SelectorData()->RemoveData(m_criterion);
    return true;
  }
  
  template <class Argument_Type,class Result_Type>
  unsigned int Grid_Creator<Argument_Type,Result_Type>::
  CreateSinglePoint(GridArgumentType *boundary,bool newpoint,bool force)
  {
    unsigned int success=1;
    GridArgumentType lower, upper;
    GridResultType newxs, newmax;
    lower=ATOOLS::Max((*p_xaxis)[(*p_xaxis)(boundary[1])*(GridArgumentType)(3.0/4.0)
				 +(*p_xaxis)(boundary[2])*(GridArgumentType)(1.0/4.0)],(*p_xaxis)[GridXMin()]);
    upper=ATOOLS::Min((*p_xaxis)[(*p_xaxis)(boundary[1])*(GridArgumentType)(1.0/4.0)
				 +(*p_xaxis)(boundary[2])*(GridArgumentType)(3.0/4.0)],(*p_xaxis)[GridXMax()]);
    if (lower==upper) return 0;
    p_processes->Reset();
    p_processes->SelectorData()->SetData(m_criterion,m_initialdata.flavs,m_initialdata.help,lower,upper);
    p_processes->ResetSelector(p_processes->SelectorData());
    p_processes->SetMax(0.0,1);
    p_processes->CalculateTotalXSec("");
    p_processes->SetMax(0.0,0);
    newxs=(GridResultType)p_processes->TotalXS()/(upper-lower);
    newmax=(GridResultType)p_processes->Max();
    if (((GridResultType)dabs((*p_yaxis)(newxs)-(*p_yaxis)(m_lastxs))>GridDeltaYMin())||(force)) {
      if (!newpoint) {
	if (!p_gridhandler[0]->Grid()->DeleteXPoint(boundary[0])) success=0;
	if (m_storemax) if (!p_gridhandler[1]->Grid()->DeleteXPoint(boundary[0])) success=0;
      }
      if (p_gridhandler[0]->Grid()->AddPoint(boundary[0],newxs)) {
	m_lastxs=newxs;
      }
      else {
	ATOOLS::msg.Error()<<"Grid_Creator_Base::CreateInitialGrid(): Could not add last cross section! "
			   <<"Ignored it instead."<<std::endl
			   <<"   Please do either reduce the grid point distance "<<std::endl
			   <<"   or select higher precision for the integration step."<<std::endl;
	success=0;
      }
      if (m_storemax) {
	if (p_gridhandler[1]->Grid()->AddPoint(boundary[0],newmax)) {
	  m_lastmax=newmax;
	}
	else {
	  ATOOLS::msg.Error()<<"Grid_Creator_Base::CreateInitialGrid(): Could not add last maximum! "
			     <<"Ignored it instead."<<std::endl
			     <<"   Please do either reduce the grid point distance "<<std::endl
			     <<"   or select higher precision for the integration step."<<std::endl;
	  success=0;
	}
      }
    }
    ATOOLS::msg.Out()<<"Grid_Creator::CalculateSingleValue(): Got value for "<<boundary[0]<<" GeV"<<std::endl
		     <<"   Calculation for "<<lower<<" GeV < "<<p_xaxis->Variable().Name()
		     <<" < "<<upper<<" GeV yielded "<<newxs*rpa.Picobarn()<<" pb/GeV ( max = "
		     <<newmax*rpa.Picobarn()<<" pb/GeV )"<<std::endl;
    return success;
  }
  
  template <class Argument_Type,class Result_Type>
  bool Grid_Creator<Argument_Type,Result_Type>::WriteOutGrid(std::vector<std::string> addcomments,
							     std::string tempofile,std::string tempopath)
  {
    if (tempopath!=ATOOLS::nullstring) this->SetOutputPath(tempopath);
    if (tempofile!=ATOOLS::nullstring) this->SetOutputFile(tempofile);
    if (this->OutputFile()==ATOOLS::nullstring) {
      ATOOLS::msg.Error()<<"Grid_Creator_Base::WriteOutGrid("<<&addcomments<<","<<tempofile<<","<<tempopath<<"): "
			 <<"No output file specified!"<<std::endl
			 <<"   Writing grid to 'output_"<<m_maxextension<<"'"<<std::endl;
      this->SetOutputFile("output_");
    }
    bool success=true;
    tempofile=this->OutputFile();
    this->SetOutputFile(tempofile+m_xsextension);
    addcomments.push_back("--------------------");
    addcomments.push_back(std::string("max file : ")+tempofile+m_maxextension);
    success=success&&WriteSingleGrid(p_gridhandler[0],addcomments);
    if (m_storemax) {
      p_gridhandler[1]->Grid()->XAxis()->
	SetVariable(p_gridhandler[0]->Grid()->XAxis()->Variable().Name());
      p_gridhandler[1]->Grid()->YAxis()->
	SetVariable(std::string("\\frac{\\partial \\sigma}{\\partial \\Omega}_{max}"));
      this->SetOutputFile(tempofile+m_maxextension);
      addcomments[addcomments.size()-1]=std::string("xs file : ")+tempofile+m_xsextension;
      success=success&&WriteSingleGrid(p_gridhandler[1],addcomments);
    }
    this->SetOutputFile(tempofile);
    return success;
  }

} // end of namespace AMISIC

#endif
