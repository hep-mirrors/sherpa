#ifndef Grid_Handler_C
#define Grid_Handler_C

#include "Grid_Handler.H"
#include "Message.H"
#include "Data_Reader.H"
#include "Data_Writer.H"

namespace AMISIC {

  template <class Argument_Type,class Result_Type>
  Grid_Handler<Argument_Type,Result_Type>::Grid_Handler()
  { Init(); }
  
  template <class Argument_Type,class Result_Type>
  Grid_Handler<Argument_Type,Result_Type>::
  Grid_Handler(ATOOLS::Type::ID streamtype,std::string streamname)
  { Init(); SetStreamType(streamtype); SetStreamName(streamname); }
  
  template <class Argument_Type,class Result_Type>
  Grid_Handler<Argument_Type,Result_Type>::
  Grid_Handler(GridFunctionType *grid)
  { Init(); Grid()->Import(*grid); }
  
  template <class Argument_Type,class Result_Type>
  Grid_Handler<Argument_Type,Result_Type>::~Grid_Handler() 
  { delete p_grid; }

  template <class Argument_Type,class Result_Type>
  void Grid_Handler<Argument_Type,Result_Type>::Init()
  {
    m_streamtype = ATOOLS::Type::TUnknown; 
    m_streamname = std::string();
    m_datatag = std::string("[x,y] =");
    p_grid = new GridFunctionType();
  }
  
  template <class Argument_Type,class Result_Type>
  bool Grid_Handler<Argument_Type,Result_Type>::
  ReadIn(Type::ID temptype,std::string tempname)
  {
    if (temptype==Type::TUnknown) temptype=m_streamtype;
    switch (temptype) {
    case ATOOLS::Type::TFStream:
      return ReadFromFile(tempname);
      break;
    case ATOOLS::Type::TString:
      // return ReadFromString(tempname);
      break;
    case ATOOLS::Type::TUnknown:
    default:
      ATOOLS::msg.Error()<<"Grid_Handler::ReadIn(..): "
			 <<"No stream type specified."
			 <<"   Abort Reading."<<std::endl;
      return false;
      break;
    }
    return false;
  }
    
  template <class Argument_Type,class Result_Type>
  bool Grid_Handler<Argument_Type,Result_Type>::
  WriteOut(Type::ID temptype,std::string tempname,
	   std::vector<std::string> comments)
  {
    if (temptype==Type::TUnknown) temptype=m_streamtype;
    switch (temptype) {
    case ATOOLS::Type::TFStream:
      return WriteToFile(tempname,comments);
      break;
    case ATOOLS::Type::TString:
      // return WriteToString(tempname,comments);
      break;
    case ATOOLS::Type::TUnknown:
    default:
      ATOOLS::msg.Error()<<"Grid_Handler::WriteOut(..): "
			 <<"No stream type specified."
			 <<"   Abort writing."<<std::endl;
      return false;
      break;
    }
    return false;
  }
  
  template <class Argument_Type,class Result_Type>
  bool Grid_Handler<Argument_Type,Result_Type>::
  ReadFromFile(std::string tempname)
  {
    if (tempname!=std::string("")) m_streamname=tempname;
    if (m_streamname==std::string("")) {
      ATOOLS::msg.Error()<<"Grid_Handler::ReadFromFile(..): "
			 <<"No filename specified."<<std::endl
			 <<"   Abort reading."<<std::endl;
      return false;
    }
    std::string gridxscaling, gridyscaling, gridxvariable, gridyvariable;
    ATOOLS::Data_Reader *reader = new ATOOLS::Data_Reader("=",";","#");
    reader->SetInputFile(m_streamname);
    reader->AddIgnore("!");
    reader->SetVectorType(reader->VHorizontal);
    if (!reader->ReadFromFile(gridxscaling,"x scale :")) {
      msg_Tracking()<<"Grid_Handler::ReadFromFile("<<tempname
		    <<"): Aborted reading.\n"
		    <<"   No x scaling information in "<<tempname
		    <<"! "<<std::endl;
      return false;
    }
    if (!reader->ReadFromFile(gridyscaling,"y scale :")) {
      msg_Tracking()<<"Grid_Handler::ReadFromFile("<<tempname
		    <<"): Aborted reading.\n"
		    <<"   No y scaling information in "<<tempname
		    <<"! "<<std::endl;
      return false;
    }
    std::vector<std::string> temp;
    if (!reader->VectorFromFile(temp,"x :")) {
      gridxvariable=std::string("Unknown");
    }
    else {
      gridxvariable=temp[0];
      for (unsigned int i=1;i<temp.size();++i) 
	gridxvariable+=std::string(" ")+temp[i];
    }
    if (!reader->VectorFromFile(temp,"y :")) {
      gridyvariable=std::string("Unknown");
    }
    else {
      gridyvariable=temp[0];
      for (unsigned int i=1;i<temp.size();++i) 
	gridyvariable+=std::string(" ")+temp[i];
    }
    reader->SetComment("!");
    reader->SetIgnore(":");
    std::vector<std::vector<GridArgumentType> > xydata;
    reader->MatrixFromFile(xydata,m_datatag);
    if (xydata.size()<2) return false;
    std::vector<GridArgumentType> xdata = xydata[0];
    std::vector<GridResultType> ydata = 
      std::vector<GridResultType>(xydata[1].size());
    for (unsigned int i=0;i<xydata[1].size();++i) ydata[i]=xydata[1][i];
    delete reader;
    if (ydata.size()!=xdata.size()) {
      xydata[1].clear();
      xdata.clear();
      ydata.clear();
      return false;
    }
    p_grid->XAxis()->SetScaling(gridxscaling);
    p_grid->YAxis()->SetScaling(gridyscaling);
    p_grid->XAxis()->SetVariable(ATOOLS::Variable(gridxvariable));
    p_grid->YAxis()->SetVariable(ATOOLS::Variable(gridyvariable));    
    typename ATOOLS::Axis<GridArgumentType>::ScalingModeID 
      xscalingmode=p_grid->XAxis()->ScalingMode();
    typename ATOOLS::Axis<GridResultType>::ScalingModeID 
      yscalingmode=p_grid->YAxis()->ScalingMode();
    p_grid->XAxis()->SetScalingMode(p_grid->XAxis()->Identical);
    p_grid->YAxis()->SetScalingMode(p_grid->YAxis()->Identical);
    p_grid->Import(&xdata,&ydata);
    p_grid->XAxis()->SetScalingMode(xscalingmode);
    p_grid->YAxis()->SetScalingMode(yscalingmode);
    xydata[1].clear();
    xdata.clear();
    ydata.clear();
    return true;
  }

  template <class Argument_Type,class Result_Type>
  bool Grid_Handler<Argument_Type,Result_Type>::
  WriteToFile(std::string tempname,std::vector<std::string>& comments)
  {
    if (tempname!=std::string("")) m_streamname=tempname;
    if (m_streamname==std::string("")) {
      msg_Tracking()<<"Grid_Handler::WriteToFile(..): "
		    <<"No filename specified.\n"
		    <<"   Writing to 'last_grid.dat'."<<std::endl;
      m_streamname=std::string("last_grid.dat");
    }
    ATOOLS::Data_Writer *writer = new ATOOLS::Data_Writer(":",";","!");
    writer->SetOutputFile(m_streamname);
    writer->SetBlank(ATOOLS::defaultblank);
    writer->WriteComment("===================="); 
    writer->WriteComment(" AMISIC++ grid file "); 
    writer->WriteComment("===================="); 
    writer->WriteComment(comments);
    writer->SetBlank(ATOOLS::defaulttab);
    std::vector<GridArgumentType> oxdata;
    std::vector<GridResultType> ydata;
    typename ATOOLS::Axis<GridArgumentType>::ScalingModeID 
      xscalingmode=p_grid->XAxis()->ScalingMode();
    typename ATOOLS::Axis<GridResultType>::ScalingModeID 
      yscalingmode=p_grid->YAxis()->ScalingMode();
    p_grid->XAxis()->SetScalingMode(p_grid->XAxis()->Identical);
    p_grid->YAxis()->SetScalingMode(p_grid->YAxis()->Identical);
    p_grid->Export(&oxdata,&ydata);
    p_grid->XAxis()->SetScalingMode(xscalingmode);
    p_grid->YAxis()->SetScalingMode(yscalingmode);
    std::vector<GridResultType> xdata = 
      std::vector<GridResultType>(oxdata.size());
    for (unsigned int i=0;i<oxdata.size();++i) xdata[i]=oxdata[i];
    oxdata.clear();
    std::vector<std::vector<GridArgumentType> > xydata;
    xydata.push_back(xdata);
    xydata.push_back(ydata);
    writer->MatrixToFile(xydata,m_datatag,true,
			 ATOOLS::nullstring,writer->MNormal,12);
    delete writer;
    xdata.clear();
    ydata.clear();
    xydata.clear();
    return true;
  }

} // end of namespace AMISIC

#endif
