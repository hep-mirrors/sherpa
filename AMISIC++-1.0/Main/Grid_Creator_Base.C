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
    m_optimize(true),
    m_relativexmin(false),
    m_relativedeltaxmin(false),
    m_relativeymin(false),
    m_relativedeltaymin(false) {}

  template <class Argument_Type,class Result_Type>
  Grid_Creator_Base<Argument_Type,Result_Type>::Grid_Creator_Base() {}

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
  bool Grid_Creator_Base<Argument_Type,Result_Type>::
  SetArguments(GridHandlerType *grid,ArgumentsType _arguments)
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
    grid->Grid()->XAxis()->SetVariable(m_gridxvariable);
    grid->Grid()->YAxis()->SetVariable(m_gridyvariable);
    grid->Grid()->XAxis()->SetScaling(m_gridxscaling);
    grid->Grid()->YAxis()->SetScaling(m_gridyscaling);
    return true;
  }

  template <class Argument_Type,class Result_Type>
  bool Grid_Creator_Base<Argument_Type,Result_Type>::
  ReadInArguments(std::string tempifile,std::string tempipath)
  {
    ATOOLS::msg.Error()<<"Grid_Creator_Base::ReadInArguments(): "
		       <<"Virtual method called!"<<std::endl;
    return false;
  }

  template <class Argument_Type,class Result_Type>
  bool Grid_Creator_Base<Argument_Type,Result_Type>::ReadSingleArguments(GridHandlerType *grid)
  {
    if (!CheckInputFile()) return false;
    std::vector<std::string> temp;
    ATOOLS::Data_Reader *reader = new ATOOLS::Data_Reader("=",";","!");
    reader->SetInputPath(InputPath());
    reader->SetInputFile(InputFile());
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
      if (m_useymax&&m_useymin) m_griddeltaymax=(GridYMax()-GridYMin())/(GridResultType)2.0;
      else m_optimize=false;
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
    grid->Grid()->XAxis()->SetVariable(m_gridxvariable);
    grid->Grid()->YAxis()->SetVariable(m_gridyvariable);
    grid->Grid()->XAxis()->SetScaling(m_gridxscaling);
    grid->Grid()->YAxis()->SetScaling(m_gridyscaling);
    return true;
  }
  
  template <class Argument_Type,class Result_Type>
  bool Grid_Creator_Base<Argument_Type,Result_Type>::CreateGrid()
  {
    bool success=true;
    int formerlevel=ATOOLS::msg.Level();
    ATOOLS::msg.SetLevel(m_outputlevel);
    if (!CreateInitialGrid()) {
      ATOOLS::msg.Error()<<"Grid_Creator_Base::CreateGrid(..): "
			 <<"Initial grid creation failed! Abort."<<std::endl;
      return false;
    }
    if (!CreateOptimizedGrid()) {
      ATOOLS::msg.Out()<<"Grid_Creator_Base::CreateGrid(..): "
		       <<"Sorry, grid cannot be optimized."<<std::endl;
      success=false;
    }
    if (!WriteOutGrid()) {
      ATOOLS::msg.Out()<<"Grid_Creator_Base::CreateGrid(..): "
		       <<"Sorry, grid cannot be written to '"<<OutputFile()<<"'"<<std::endl;
      success=false;
    }
    ATOOLS::msg.Out()<<"Grid_Creator_Base: Wrote grid to '"<<OutputFile()<<"'"<<std::endl;
    ATOOLS::msg.SetLevel(formerlevel);
    return success;
  }

  template <class Argument_Type,class Result_Type>
  unsigned int Grid_Creator_Base<Argument_Type,Result_Type>::
  CreateSinglePoint(GridArgumentType *boundary,bool newpoint,bool force)
  {
    ATOOLS::msg.Error()<<"Grid_Creator_Base::CreateSinglePoint(): "
		       <<"Virtual method called!"<<std::endl;
    return false;
  }

  template <class Argument_Type,class Result_Type>
  bool Grid_Creator_Base<Argument_Type,Result_Type>::CreateInitialGrid()
  {
    ATOOLS::msg.Error()<<"Grid_Creator_Base::CreateInitialGrid(): "
		       <<"Virtual method called!"<<std::endl;
    return false;
  }

  template <class Argument_Type,class Result_Type>
  bool Grid_Creator_Base<Argument_Type,Result_Type>::InitializeSingleGrid(GridHandlerType *grid)
  {
    if ((GridXMin()>=GridXMax())||
	((int)((GridXMax()-GridXMin())/GridDeltaXMax())<2)) {
      ATOOLS::msg.Error()<<"Grid_Creator_Base::CreateInitialGrid("<<InputFile()<<"): "
			 <<"Argument boundaries improperly specified!"<<std::endl
			 <<"   Abort grid creation."<<std::endl;
      return false;
    }
    GridArgumentType *boundary = new GridArgumentType[3];
    ATOOLS::Axis<GridArgumentType> *xaxis=grid->Grid()->XAxis();
    int npoints=(int)((GridXMax()-GridXMin())/GridDeltaXMax());
    SetGridDeltaXMax((GridXMax()-GridXMin())/(GridArgumentType)npoints);
    InitializeCalculation();
    for (int i=npoints;i>=0;--i) {
      boundary[0]=(*xaxis)[GridXMin()+(double)i*GridDeltaXMax()];
      boundary[1]=(*xaxis)[GridXMin()+((double)i-1.0)*GridDeltaXMax()];
      boundary[2]=(*xaxis)[GridXMin()+((double)i+1.0)*GridDeltaXMax()];
      ATOOLS::msg.Out()<<"Grid_Creator_Base::CreateInitialGrid(): "
		       <<"Calculation for "<<npoints-i<<"th of "<<npoints<<" points in progress."<<std::endl;
      if (!CreateSinglePoint(boundary,true,(i==0)||(i==npoints))) {
	ATOOLS::msg.Error()<<"Grid_Creator_Base::CreateInitialGrid(): "
			   <<"Cannot determine point at "<<m_gridxvariable<<" = "<<boundary[0]<<"!"<<std::endl
			   <<"   Abort attempt and continue."<<std::endl;
      }
    }
    delete [] boundary;
    return true;
  }

  template <class Argument_Type,class Result_Type>
  bool Grid_Creator_Base<Argument_Type,Result_Type>::OptimizeSingleGrid(GridHandlerType *grid)
  {
    if (!m_optimize) return true;
    if ((GridYMin()>GridYMax())||(GridDeltaYMax()<=(GridArgumentType)0.0)) {
      ATOOLS::msg.Error()<<"Grid_Creator_Base::CreateOptimizedGrid(): "
			 <<"Grid result boundaries improperly specified! "<<std::endl
			 <<"   Abort optimization."<<std::endl;
      return false;
    }
    unsigned int errorcounter=0, correct;
    GridArgumentType left, right, middle, delta, *boundary = new GridArgumentType[3];
    GridResultType deltaymax;
    if (grid==NULL) return false;
    ATOOLS::Axis<GridArgumentType> *xaxis=grid->Grid()->XAxis();
    ATOOLS::Axis<GridResultType> *yaxis=grid->Grid()->YAxis();
    InitializeCalculation();
    yaxis->SetScalingMode(yaxis->Identical);
    for (;(deltaymax=grid->Grid()->DeltaYMax(left,right))>GridDeltaYMax();yaxis->SetScalingMode(yaxis->Identical)) {
      if (((*xaxis)(right)-(*xaxis)(left))<GridDeltaXMin()) break;
      if (!m_useymin) SetGridYMin(grid->Grid()->YMin());
      if (!m_useymax) SetGridYMax(grid->Grid()->YMax());
      if ((grid->Grid()->YMin()>=GridYMin())&&(grid->Grid()->YMax()<=GridYMax())) {
	yaxis->SetScalingMode(yaxis->Reference);
	middle=(*xaxis)[((*xaxis)(left)+(*xaxis)(right))/(GridArgumentType)2.0];
	boundary[0]=middle;
	boundary[1]=left;
	boundary[2]=right;
	ATOOLS::msg.Out()<<"Grid_Creator_Base::OptimizeSingleGrid(): "
			 <<"Calculation for new grid point in progress."<<std::endl
			 <<"   Currently \\Delta y_{max} = "<<deltaymax
			 <<" vs. \\Delta y_{limit} = "<<GridDeltaYMax()<<std::endl;
	if ((correct=CreateSinglePoint(boundary))==0) {
	  ATOOLS::msg.Error()<<"Grid_Creator_Base::OptimizeSingleGrid(): "
			     <<"Cannot determine point at "<<m_gridxvariable<<" = "<<boundary[0]<<"!"<<std::endl
			     <<"   Abort attempt and continue."<<std::endl;
	  if (++errorcounter>grid->Grid()->XDataSize()/10) {
	    ATOOLS::msg.Error()<<"Grid_Creator_Base::OptimizeSingleGrid(): Too many errors occured! "
			       <<"Abort optimization."<<std::endl;
	    return false;
	  }
	}
	else {
	  delta=((*xaxis)(right)-(*xaxis)(left))/(GridArgumentType)2.0;
	  for (int i=0;i<(int)correct;++i) {
	    boundary[0]=(*xaxis)[(*xaxis)(left)-delta*(GridArgumentType)i];
	    boundary[1]=(*xaxis)[(*xaxis)(left)-delta*(GridArgumentType)(i+1)];
	    boundary[2]=(*xaxis)[(*xaxis)(left)-delta*(GridArgumentType)(i-1)];
	    if ((boundary[0]>=(*xaxis)[GridXMin()])&&(boundary[0]<=(*xaxis)[GridXMax()])) {
	      ATOOLS::msg.Out()<<"Grid_Creator_Base::OptimizeSingleGrid(): "
			       <<"Calculation for corrected "<<i+1<<"th left point in progress."<<std::endl
			       <<"   Currently \\Delta y_{max} = "<<deltaymax
			       <<" vs. \\Delta y_{limit} = "<<GridDeltaYMax()<<std::endl;
	      if (!CreateSinglePoint(boundary,false)) {
		ATOOLS::msg.Error()<<"Grid_Creator_Base::OptimizeSingleGrid(): Cannot delete obsolete grid point at "
				   <<boundary[0]<<"! Abort optimization."<<std::endl
				   <<"   Warning! Grid values will be unreliable."<<std::endl;
		return false;
	      }
	    }
	    boundary[0]=(*xaxis)[(*xaxis)(right)+delta*(GridArgumentType)i];
	    boundary[1]=(*xaxis)[(*xaxis)(right)+delta*(GridArgumentType)(i-1)];
	    boundary[2]=(*xaxis)[(*xaxis)(right)+delta*(GridArgumentType)(i+1)];
	    if ((boundary[0]>=(*xaxis)[GridXMin()])&&(boundary[0]<=(*xaxis)[GridXMax()])) {
	      ATOOLS::msg.Out()<<"Grid_Creator_Base::OptimizeSingleGrid(): "
			       <<"Calculation for corrected "<<i+1<<"th right point in progress."<<std::endl
			       <<"   Currently \\Delta y_{max} = "<<deltaymax
			       <<" vs. \\Delta y_{limit} = "<<GridDeltaYMax()<<std::endl;
	      if (!CreateSinglePoint(boundary,false)) {
		ATOOLS::msg.Error()<<"Grid_Creator_Base::OptimizeSingleGrid(): Cannot delete obsolete grid point at "
				   <<boundary[0]<<"! Abort optimization."<<std::endl
				   <<"   Warning! Grid values will be unreliable."<<std::endl;
		return false;
	      }
	    }
	  }
	}
      }
      else {
	if (grid->Grid()->YMin()<GridYMin()) grid->Grid()->DeleteYPoint(grid->Grid()->YMin());
	else if (grid->Grid()->YMax()>GridYMax()) grid->Grid()->DeleteYPoint(grid->Grid()->YMax());
	yaxis->SetScalingMode(yaxis->Reference);
      }
    } 
    delete boundary;
    yaxis->SetScalingMode(yaxis->Reference);
    return true;
  }

  template <class Argument_Type,class Result_Type>
  bool Grid_Creator_Base<Argument_Type,Result_Type>::CreateOptimizedGrid()
  {
    ATOOLS::msg.Error()<<"Grid_Creator_Base::CreateOptimizedGrid(): "
		       <<"Virtual method called!"<<std::endl;
    return false;
  }

  template <class Argument_Type,class Result_Type>
  bool Grid_Creator_Base<Argument_Type,Result_Type>::
  WriteOutGrid(std::vector<std::string> addcomments,std::string tempofile,std::string tempopath)
  {
    ATOOLS::msg.Error()<<"Grid_Creator_Base::WriteOutGrid(..): "
		       <<"Virtual method called!"<<std::endl;
    return false;
  }

  template <class Argument_Type,class Result_Type>
  bool Grid_Creator_Base<Argument_Type,Result_Type>::
  WriteSingleGrid(GridHandlerType *grid,std::vector<std::string> addcomments)
  {
    if (!CheckOutputFile()) return false;
    std::vector<std::string> comments;
    std::string xvar, yvar;
    xvar=grid->Grid()->XAxis()->Variable().Name();
    yvar=grid->Grid()->YAxis()->Variable().Name();
    comments.push_back(std::string("x : ")+xvar);
    comments.push_back(std::string("y : ")+yvar);
    comments.push_back("--------------------");
    comments.push_back(std::string("x scale : ")+grid->Grid()->XAxis()->Scaling()->Name());
    comments.push_back(std::string("y scale : ")+grid->Grid()->YAxis()->Scaling()->Name());
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
    return grid->WriteOut(ATOOLS::Type::TFStream,OutputPath()+OutputFile(),comments);
  }

} // end of namespace AMISIC

#endif

