#ifndef Grid_Creator_Base_C
#define Grid_Creator_Base_C

#include "Grid_Creator_Base.H"
#include "Data_Reader.H"
#include "Data_Writer.H"
#include "Type.H"

namespace AMISIC {

  template <class Argument_Type,class Result_Type>
  Grid_Creator_Arguments<Argument_Type,Result_Type>::Grid_Creator_Arguments():
    m_useymax(false),
    m_useymin(false),
    m_relativexmin(false),
    m_relativedeltaxmin(false),
    m_relativeymin(false),
    m_relativedeltaymin(false) {}

  template <class Argument_Type,class Result_Type>
  Grid_Creator_Base<Argument_Type,Result_Type>::Grid_Creator_Base(GridHandlerType *_p_gridhandler):
    p_gridhandler(_p_gridhandler)
  {
    if (p_gridhandler==NULL) {
      ATOOLS::msg.Error()<<"Grid_Creator_Base::Grid_Creator_Base("<<_p_gridhandler<<"): "
			 <<"Grid handler is not initialized! Abort."<<std::endl;
      abort();
    }
    SetInputFile("");
    SetOutputFile("");
  }

  template <class Argument_Type,class Result_Type>
  Grid_Creator_Base<Argument_Type,Result_Type>::~Grid_Creator_Base() {}

  template <class Argument_Type,class Result_Type>
  bool Grid_Creator_Base<Argument_Type,Result_Type>::InitializeCalculation()
  {
    ATOOLS::msg.Error()<<"Grid_Creator_Base::InitializeCalculation(): "
		       <<"Virtual method called!"<<std::endl;
    return false;
  }

  template <class Argument_Type,class Result_Type>
  Result_Type Grid_Creator_Base<Argument_Type,Result_Type>::CalculateSingleValue(GridArgumentType nextleft,
								     GridArgumentType nextright)
  {
    ATOOLS::msg.Error()<<"Grid_Creator_Base::CalculateSingleValue("<<nextleft<<","<<nextright<<"): "
		       <<"Virtual method called!"<<std::endl;
    return (GridResultType)0.0;
  }

  template <class Argument_Type,class Result_Type>
  bool Grid_Creator_Base<Argument_Type,Result_Type>::SetArguments(ArgumentsType _arguments)
  {
    SetGridXMin(_arguments.GridXMin());
    SetGridXMax(_arguments.GridXMax());
    SetGridDeltaXMin(_arguments.GridDeltaXMin());
    SetGridDeltaXMax(_arguments.GridDeltaXMax());
    SetGridYMin(_arguments.GridYMin());
    SetGridYMax(_arguments.GridYMax());
    SetGridDeltaYMin(_arguments.GridDeltaYMin());
    SetGridDeltaYMax(_arguments.GridDeltaYMax());
    m_usedeltaymax=_arguments.m_usedeltaymax;
    m_usedeltaymin=_arguments.m_usedeltaymin;
    m_outputlevel=_arguments.m_outputlevel;
    m_gridxvariable=_arguments.m_gridxvariable;
    m_gridyvariable=_arguments.m_gridyvariable;
    p_gridhandler->Grid()->XAxis()->SetVariable(m_gridxvariable);
    p_gridhandler->Grid()->YAxis()->SetVariable(m_gridyvariable);
    p_gridhandler->SetScaling(m_gridxscaling,m_gridyscaling);
    return true;
  }

  template <class Argument_Type,class Result_Type>
  bool Grid_Creator_Base<Argument_Type,Result_Type>::ReadInArguments(std::string tempifile)
  {
    if (tempifile!=ATOOLS::nullstring) SetInputFile(tempifile);
    if (m_inputfile==ATOOLS::nullstring) {
      ATOOLS::msg.Error()<<"Grid_Creator_Base::Arguments("<<tempifile<<"): "
			 <<"No input file specified! Abort."<<std::endl;
      abort();
    }
    std::vector<std::string> temp;
    ATOOLS::Data_Reader *reader = new ATOOLS::Data_Reader("=",";","!");
    reader->SetFileName(m_inputfile);
    if (!reader->VectorFromFile(temp,"X_VAR",ATOOLS::noinputtag,reader->VHorizontal)) {
      m_gridxvariable=std::string("");
    }
    else {
      m_gridxvariable=temp[0];
      for (unsigned int i=1;i<temp.size();++i) m_gridxvariable+=std::string(" ")+temp[i];
    }
    if (!reader->ReadFromFile(m_gridxscaling,"X_SCALE")) {
      m_gridxscaling=std::string("id");
    }
    if (!reader->ReadFromFile(m_gridxmax,"V_X_MAX")) SetGridXMax((GridArgumentType)0.0);
    if (reader->VectorFromFile(temp,"V_X_MIN",ATOOLS::noinputtag,reader->VHorizontal)) { 
      reader->ReadFromString(m_gridxmin,ATOOLS::nullstring,temp[0]);
      if (temp.size()>1) if (temp[1]==std::string("relative")) m_relativexmin=true;
    }
    else {
      m_gridxmin=(GridArgumentType)1.0e-6;
      m_relativexmin=true;
    }
    if (!reader->ReadFromFile(m_griddeltaxmax,"D_X_MAX")) {
      m_griddeltaxmax=(GridXMax()-GridXMin())/(GridArgumentType)2.0;
    }
    if (reader->VectorFromFile(temp,"D_X_MIN",ATOOLS::noinputtag,reader->VHorizontal)) {
      reader->ReadFromString(m_griddeltaxmin,ATOOLS::nullstring,temp[0]); 
      if (temp.size()>1) if (temp[1]==std::string("relative")) m_relativedeltaxmin=true;
    }
    else {
      m_griddeltaxmin=(GridResultType)1.0e-6;
      m_relativedeltaxmin=true;
    }
    if (!reader->VectorFromFile(temp,"Y_VAR",ATOOLS::noinputtag,reader->VHorizontal)) {
      m_gridyvariable=std::string("");
    }
    else {
      m_gridyvariable=temp[0];
      for (unsigned int i=1;i<temp.size();++i) m_gridyvariable+=std::string(" ")+temp[i];
    }
    if (!reader->ReadFromFile(m_gridyscaling,"Y_SCALE")) 
      m_gridyscaling=std::string("id");
    if (reader->ReadFromFile(m_gridymax,"V_Y_MAX")) m_useymax=true;
    if (reader->VectorFromFile(temp,"V_Y_MIN",ATOOLS::noinputtag,reader->VHorizontal)) { 
      reader->ReadFromString(m_gridymin,ATOOLS::nullstring,temp[0]);
      if (temp.size()>1) if (temp[1]==std::string("relative")) m_relativeymin=true;
      m_useymin=true;
    }
    if (!reader->ReadFromFile(m_griddeltaymax,"D_Y_MAX")) {
      m_griddeltaymax=(GridYMax()-GridYMin())/(GridResultType)2.0;
    }
    if (reader->VectorFromFile(temp,"D_Y_MIN",ATOOLS::noinputtag,reader->VHorizontal)) { 
      reader->ReadFromString(m_griddeltaymin,ATOOLS::nullstring,temp[0]);
      if (temp.size()>1) if (temp[1]==std::string("relative")) m_relativedeltaymin=true;
    }
    else {
      m_griddeltaymin=(GridResultType)1.0e-6;
      m_relativedeltaymin=true;
    }
    if (!reader->ReadFromFile(m_outputlevel,"OUTPUT")) m_outputlevel=ATOOLS::msg.Level();
    delete reader;
    p_gridhandler->Grid()->XAxis()->SetVariable(m_gridxvariable);
    p_gridhandler->Grid()->YAxis()->SetVariable(m_gridyvariable);
    p_gridhandler->Grid()->XAxis()->SetScaling(m_gridxscaling);
    p_gridhandler->Grid()->YAxis()->SetScaling(m_gridyscaling);
    return true;
  }
  
  template <class Argument_Type,class Result_Type>
  bool Grid_Creator_Base<Argument_Type,Result_Type>::CreateGrid(std::string tempofile)
  {
    bool success=true;
    int formerlevel=ATOOLS::msg.Level();
    ATOOLS::msg.SetLevel(m_outputlevel);
    if (!CreateInitialGrid()) {
      ATOOLS::msg.Error()<<"Grid_Creator_Base::CreateGrid(..): Initial grid creation failed! Abort."<<std::endl;
      return false;
    }
    if (!CreateOptimizedGrid()) {
      ATOOLS::msg.Out()<<"Grid_Creator_Base::CreateGrid(..): Sorry, grid cannot be optimized."<<std::endl;
      success=false;
    }
    if (tempofile!=ATOOLS::nullstring) { 
      SetOutputFile(tempofile);
      if (!WriteOutGrid()) {
	ATOOLS::msg.Out()<<"Grid_Creator_Base::CreateGrid(..): "
			 <<"Grid cannot be written to '"<<m_outputfile<<"'"<<std::endl;
	success=false;
      }
    }
    ATOOLS::msg.Out()<<"Grid_Creator_Base: Wrote grid to '"<<m_outputfile<<"'"<<std::endl;
    ATOOLS::msg.SetLevel(formerlevel);
    return success;
  }

  template <class Argument_Type,class Result_Type>
  bool Grid_Creator_Base<Argument_Type,Result_Type>::CreateInitialGrid()
  {
    if ((GridXMin()>=GridXMax())||
	((int)((GridXMax()-GridXMin())/GridDeltaXMax())<2)) {
      ATOOLS::msg.Error()<<"Grid_Creator_Base::CreateInitialGrid("<<m_inputfile<<"): "
			 <<"Argument boundaries improperly specified!"<<std::endl
			 <<"   Abort grid creation."<<std::endl;
      return false;
    }
    GridArgumentType left, right, middle;
    GridResultType oldvalue, newvalue;
    ATOOLS::Axis<GridArgumentType> *xaxis=p_gridhandler->Grid()->XAxis();
    ATOOLS::Axis<GridResultType> *yaxis=p_gridhandler->Grid()->YAxis();
    int npoints=(int)((GridXMax()-GridXMin())/GridDeltaXMax());
    SetGridDeltaXMax((GridXMax()-GridXMin())/(GridArgumentType)npoints);
    InitializeCalculation();
    for (int i=npoints;i>=0;--i) {
      left=(*xaxis)[GridXMin()+((double)i-1.0)*GridDeltaXMax()];
      middle=(*xaxis)[GridXMin()+(double)i*GridDeltaXMax()];
      right=(*xaxis)[GridXMin()+((double)i+1.0)*GridDeltaXMax()];
      ATOOLS::msg.Out()<<"Grid_Creator_Base::CreateInitialGrid(): "
		       <<"Calculation for "<<npoints-i<<"th of "<<npoints<<" points in progress."<<std::endl;
      newvalue=CalculateSingleValue(left,right);
      if (((GridResultType)dabs((*yaxis)(newvalue)-(*yaxis)(oldvalue))>GridDeltaYMin())||
	  (i==npoints)||(i==0)) {
	if (p_gridhandler->Grid()->AddPoint(middle,newvalue)) {
	  oldvalue=newvalue;
	}
	else {
	  ATOOLS::msg.Error()<<"Grid_Creator_Base::CreateInitialGrid(): Tried to add a value twice! "
			     <<"Ignored last value."<<std::endl
			     <<"   Please do either reduce the grid point distance "<<std::endl
			     <<"   or select higher precision for the integration step."<<std::endl;
	}
      }
    }
    return true;
  }

  template <class Argument_Type,class Result_Type>
  bool Grid_Creator_Base<Argument_Type,Result_Type>::CreateOptimizedGrid()
  {
    if ((GridYMin()>GridYMax())||(GridDeltaYMax()<=(GridArgumentType)0.0)) {
      ATOOLS::msg.Error()<<"Grid_Creator_Base::CreateOptimizedGrid(): "
			 <<"Grid result boundaries improperly specified! "<<std::endl
			 <<"   Abort optimization."<<std::endl;
      return false;
    }
    unsigned int errorcounter=0;
    GridArgumentType left, right, middle, nextleft, nextright;
    GridResultType min, max, deltaymax;
    ATOOLS::Data_To_Function<GridArgumentType,GridResultType> *grid=p_gridhandler->Grid();
    ATOOLS::Axis<GridArgumentType> *xaxis=grid->XAxis();
    ATOOLS::Axis<GridResultType> *yaxis=grid->YAxis();
    InitializeCalculation();
    yaxis->SetScalingMode(yaxis->Identical);
    for (;(deltaymax=grid->DeltaYMax(left,right))>GridDeltaYMax();yaxis->SetScalingMode(yaxis->Identical)) {
      if (((*xaxis)(right)-(*xaxis)(left))<GridDeltaXMin()) break;
      min=ATOOLS::Min(grid->Y(left,grid->Data),grid->Y(right,grid->Data));
      max=ATOOLS::Max(grid->Y(left,grid->Data),grid->Y(right,grid->Data));
      if (!m_useymin) SetGridYMin(grid->YMin());
      if (!m_useymax) SetGridYMax(grid->YMax());
      yaxis->SetScalingMode(yaxis->Reference);
      if ((min>=GridYMin())&&(max<=GridYMax())) {
	middle=(*xaxis)[((*xaxis)(left)+(*xaxis)(right))/(GridArgumentType)2.0];
	ATOOLS::msg.Out()<<"Grid_Creator_Base::CreateOptimizedGrid(): "
			 <<"Calculation for new grid point in progress."<<std::endl
			 <<"   Currently \\Delta y_{max} = "<<deltaymax
			 <<" vs. \\Delta y_{limit} = "<<GridDeltaYMax()<<std::endl;
	if (!grid->AddPoint(middle,CalculateSingleValue(left,right))) {
	  ATOOLS::msg.Error()<<"Grid_Creator_Base::CreateOptimizedGrid(): Tried to add a value twice! "
			     <<"Ignored last value."<<std::endl
			     <<"   Please do either reduce the grid point distance "<<std::endl
			     <<"   or choose higher precision for the integration step."<<std::endl;
	  if (++errorcounter>grid->XDataSize()/10) {
	    ATOOLS::msg.Error()<<"Grid_Creator_Base::CreateOptimizedGrid(): Too many errors occured! "
			       <<"Abort optimization."<<std::endl;
	    return false;
	  }
	}
	else {
	  nextleft=(*xaxis)[(*xaxis)(left)-((*xaxis)(right)-(*xaxis)(left))/(GridArgumentType)2.0 ];
	  ATOOLS::msg.Out()<<"Grid_Creator_Base::CreateOptimizedGrid(): "
			   <<"Calculation for corrected left point in progress."<<std::endl
			   <<"   Currently \\Delta y_{max} = "<<deltaymax
			   <<" vs. \\Delta y_{limit} = "<<GridDeltaYMax()<<std::endl;
	  if (!grid->ReplaceXPoint(left,CalculateSingleValue(nextleft,middle))) {
	    ATOOLS::msg.Error()<<"Grid_Creator_Base::CreateOptimizedGrid(): Cannot delete obsolete grid point at "
			       <<left<<"! Abort optimization."<<std::endl
			       <<"   Warning! Grid values will be unreliable."<<std::endl;
	    return false;
	  }
	  nextright=(*xaxis)[(*xaxis)(right)+((*xaxis)(right)-(*xaxis)(left))/(GridArgumentType)2.0 ];
	  ATOOLS::msg.Out()<<"Grid_Creator_Base::CreateOptimizedGrid(): "
			   <<"Calculation for corrected right point in progress."<<std::endl
			   <<"   Currently \\Delta y_{max} = "<<deltaymax
			   <<" vs. \\Delta y_{limit} = "<<GridDeltaYMax()<<std::endl;
	  if (!grid->ReplaceXPoint(right,CalculateSingleValue(middle,nextright))) {
	    ATOOLS::msg.Error()<<"Grid_Creator_Base::CreateOptimizedGrid(): Cannot delete obsolete grid point at "
			       <<right<<"! Abort optimization."<<std::endl
			       <<"   Warning! Grid values will be unreliable."<<std::endl;
	    return false;
	  }
	}
      }
      else {
	if (min<GridYMin()) grid->DeleteYPoint((*yaxis)[min]);
	else if (max>GridYMax()) grid->DeleteYPoint((*yaxis)[max]);
      }
    } 
    yaxis->SetScalingMode(yaxis->Reference);
    return true;
  }

  template <class Argument_Type,class Result_Type>
  bool Grid_Creator_Base<Argument_Type,Result_Type>::WriteOutGrid(std::string tempofile,
							     std::vector<std::string> addcomments)
  {
    if (tempofile!=ATOOLS::nullstring) SetOutputFile(tempofile);
    if (m_outputfile==ATOOLS::nullstring) {
      ATOOLS::msg.Error()<<"Grid_Creator_Base::WriteOutGrid("<<tempofile<<"): "
			 <<"No output file specified! Abort()."<<std::endl;
      abort();
    }
    std::vector<std::string> comments;
    std::string xvar, yvar;
    ATOOLS::Data_To_Function<GridArgumentType,GridResultType> *grid=p_gridhandler->Grid(); 
    xvar=grid->XAxis()->Variable().Name();
    yvar=grid->YAxis()->Variable().Name();
    comments.push_back(std::string("x : ")+xvar);
    comments.push_back(std::string("y : ")+yvar);
    comments.push_back("--------------------");
    comments.push_back(std::string("x scale : ")+
		       p_gridhandler->Grid()->XAxis()->Scaling()->Name());
    comments.push_back(std::string("y scale : ")+
		       p_gridhandler->Grid()->YAxis()->Scaling()->Name());
    comments.push_back("--------------------");
    comments.push_back("boundary conditions ");
    comments.push_back("--------------------");
    comments.push_back(std::string("x_{min} = ")+ATOOLS::ToString(GridXMin()));
    comments.push_back(std::string("x_{max} = ")+ATOOLS::ToString(GridXMax()));
    comments.push_back(std::string("y_{min} = ")+ATOOLS::ToString(GridYMin()));
    comments.push_back(std::string("y_{max} = ")+ATOOLS::ToString(GridYMax()));
    comments.push_back("--------------------");
    comments.push_back(std::string("\\Delta x_{max} = ")+ATOOLS::ToString(GridDeltaXMax()));
    comments.push_back(std::string("\\Delta x_{min} = ")+ATOOLS::ToString(GridDeltaXMin()));
    comments.push_back(std::string("\\Delta y_{max} = ")+ATOOLS::ToString(GridDeltaYMax()));
    comments.push_back(std::string("\\Delta y_{min} = ")+ATOOLS::ToString(GridDeltaYMin()));
    if (addcomments.size()!=0) {
      comments.push_back("--------------------");
      for (unsigned int i=0;i<addcomments.size();++i) comments.push_back(addcomments[i]);
    }
    comments.push_back("--------------------");
    comments.push_back("  Data Set follows  ");
    comments.push_back("--------------------");
    return p_gridhandler->WriteOut(ATOOLS::Type::TFStream,m_outputfile,comments);
  }

} // end of namespace AMISIC

#endif

