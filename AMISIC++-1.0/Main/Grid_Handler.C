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
