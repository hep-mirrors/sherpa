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
  Grid_Handler(ATOOLS::Type::ID _m_streamtype,std::string _m_streamname)
  { Init(); SetStreamType(_m_streamtype); SetStreamName(_m_streamname); }
  
  template <class Argument_Type,class Result_Type>
  Grid_Handler<Argument_Type,Result_Type>::Grid_Handler(GridFunctionType *_p_grid)
  { Init(); Grid()->Import(*_p_grid); }
  
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
//       return ReadFromString(tempname);
      break;
    case ATOOLS::Type::TUnknown:
    default:
      ATOOLS::msg.Error()<<"Grid_Handler::ReadIn(..): No stream type specified."
			 <<"   Abort Reading."<<std::endl;
      return false;
      break;
    }
    return false;
  }
    
  template <class Argument_Type,class Result_Type>
  bool Grid_Handler<Argument_Type,Result_Type>::
  WriteOut(Type::ID temptype,std::string tempname,std::vector<std::string> comments)
  {
    if (temptype==Type::TUnknown) temptype=m_streamtype;
    switch (temptype) {
    case ATOOLS::Type::TFStream:
      return WriteToFile(tempname,comments);
      break;
    case ATOOLS::Type::TString:
//       return WriteToString(tempname,comments);
      break;
    case ATOOLS::Type::TUnknown:
    default:
      ATOOLS::msg.Error()<<"Grid_Handler::WriteOut(..): No stream type specified."
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
      ATOOLS::msg.Error()<<"Grid_Handler::ReadFromFile(..): No filename specified."<<std::endl
			 <<"   Abort reading."<<std::endl;
      return false;
    }
    std::string gridxscaling, gridyscaling;
    ATOOLS::Data_Reader *reader = new ATOOLS::Data_Reader(":",";","#");
    reader->SetFileName(m_streamname);
    reader->AddIgnore("!");
    if (!reader->ReadFromFile(gridxscaling,"x scale")) {
      ATOOLS::msg.Error()<<"Grid_Handler::ReadFromFile("<<tempname<<"): Aborted reading."<<std::endl
			 <<"   No x scaling information in "<<tempname<<"! "<<std::endl;
      return false;
    }
    if (!reader->ReadFromFile(gridyscaling,"y scale")) {
      ATOOLS::msg.Error()<<"Grid_Handler::ReadFromFile("<<tempname<<"): Aborted reading."<<std::endl
			 <<"   No y scaling information in "<<tempname<<"! "<<std::endl;
      return false;
    }
    reader->SetComment("!");
    std::vector<std::vector<GridArgumentType> > _m_xydata;
    reader->ArrayFromFile(_m_xydata,m_datatag);
    if (_m_xydata.size()<2) return false;
    std::vector<GridArgumentType> _m_xdata = _m_xydata[0];
    std::vector<GridResultType> _m_ydata = std::vector<GridResultType>(_m_xydata[1].size());
    for (unsigned int i=0;i<_m_xydata[1].size();++i) _m_ydata[i]=_m_xydata[1][i];
    delete reader;
    if (_m_ydata.size()!=_m_xdata.size()) {
      _m_xydata[1].clear();
      _m_xdata.clear();
      _m_ydata.clear();
      return false;
    }
    SetScaling(gridxscaling,gridyscaling);
    typename ATOOLS::Axis<GridArgumentType>::ScalingMode xscalingmode=p_grid->XAxis()->GetScalingMode();
    typename ATOOLS::Axis<GridResultType>::ScalingMode yscalingmode=p_grid->YAxis()->GetScalingMode();
    p_grid->XAxis()->SetScalingMode(p_grid->XAxis()->Identical);
    p_grid->YAxis()->SetScalingMode(p_grid->YAxis()->Identical);
    p_grid->Import(&_m_xdata,&_m_ydata);
    p_grid->XAxis()->SetScalingMode(xscalingmode);
    p_grid->YAxis()->SetScalingMode(yscalingmode);
    _m_xydata[1].clear();
    _m_xdata.clear();
    _m_ydata.clear();
    return true;
  }

  template <class Argument_Type,class Result_Type>
  bool Grid_Handler<Argument_Type,Result_Type>::
  WriteToFile(std::string tempname,std::vector<std::string>& comments)
  {
    if (tempname!=std::string("")) m_streamname=tempname;
    if (m_streamname==std::string("")) {
      ATOOLS::msg.Error()<<"Grid_Handler::WriteToFile(..): No filename specified."<<std::endl
			 <<"   Writing to last_grid.dat ."<<std::endl;
      m_streamname=std::string("last_grid.dat");
    }
    ATOOLS::Data_Writer *writer = new ATOOLS::Data_Writer(":",";","!");
    writer->SetFileName(m_streamname);
    writer->SetBlank(ATOOLS::defaultblank);
    writer->WriteComment("===================="); 
    writer->WriteComment(" AMISIC++ grid file "); 
    writer->WriteComment("===================="); 
    writer->WriteComment(comments);
    writer->SetBlank(ATOOLS::defaulttab);
    std::vector<GridArgumentType> __m_xdata;
    std::vector<GridResultType> _m_ydata;
    typename ATOOLS::Axis<GridArgumentType>::ScalingMode xscalingmode=p_grid->XAxis()->GetScalingMode();
    typename ATOOLS::Axis<GridResultType>::ScalingMode yscalingmode=p_grid->YAxis()->GetScalingMode();
    p_grid->XAxis()->SetScalingMode(p_grid->XAxis()->Identical);
    p_grid->YAxis()->SetScalingMode(p_grid->YAxis()->Identical);
    p_grid->Export(&__m_xdata,&_m_ydata);
    p_grid->XAxis()->SetScalingMode(xscalingmode);
    p_grid->YAxis()->SetScalingMode(yscalingmode);
    std::vector<GridResultType> _m_xdata = std::vector<GridResultType>(__m_xdata.size());
    for (unsigned int i=0;i<_m_xdata.size();++i) _m_xdata[i]=(GridResultType)__m_xdata[i];
    __m_xdata.clear();
    std::vector<std::vector<GridArgumentType> > _m_xydata;
    _m_xydata.push_back(_m_xdata);
    _m_xydata.push_back(_m_ydata);
    writer->ArrayToFile(_m_xydata,m_datatag,true,ATOOLS::nullstring,writer->MNormal,12);
    delete writer;
    _m_xdata.clear();
    _m_ydata.clear();
    _m_xydata.clear();
    return true;
  }

  template <class Argument_Type,class Result_Type>
  void Grid_Handler<Argument_Type,Result_Type>::SetScaling(std::string gridxscaling,
							   std::string gridyscaling)
  {
    GridArgumentType argx;
    GridResultType argy;
    ATOOLS::Data_Reader *reader = new ATOOLS::Data_Reader();
    reader->SetString(gridxscaling);
    if (gridxscaling==std::string("Log")) {
      p_grid->XAxis()->SetScaling(new ATOOLS::Log_Scaling<GridArgumentType>()); 
    }
    else if (gridxscaling==std::string("Sqr")) {
      p_grid->XAxis()->SetScaling(new ATOOLS::Sqr_Scaling<GridArgumentType>()); 
    }
    else if (reader->ReadFromString(argx,"Log_B_")) { 
      p_grid->XAxis()->SetScaling(new ATOOLS::Log_B_Scaling<GridArgumentType>(argx)); 
    }
    else if (reader->ReadFromString(argx,"B_To_X_")) { 
      p_grid->XAxis()->SetScaling(new ATOOLS::B_To_X_Scaling<GridArgumentType>(argx)); 
    }
    else if (reader->ReadFromString(argx,"X_To_P_")) { 
      p_grid->XAxis()->SetScaling(new ATOOLS::X_To_P_Scaling<GridArgumentType>(argx)); 
    }
    reader->SetString(gridyscaling);
    if (gridyscaling==std::string("Log")) {
      p_grid->YAxis()->SetScaling(new ATOOLS::Log_Scaling<GridResultType>()); 
    }
    else if (gridyscaling==std::string("Sqr")) {
      p_grid->YAxis()->SetScaling(new ATOOLS::Sqr_Scaling<GridResultType>()); 
    }
    else if (reader->ReadFromString(argy,"Log_B_")) { 
      p_grid->YAxis()->SetScaling(new ATOOLS::Log_B_Scaling<GridResultType>(argy)); 
    }
    else if (reader->ReadFromString(argy,"B_To_X_")) { 
      p_grid->YAxis()->SetScaling(new ATOOLS::B_To_X_Scaling<GridResultType>(argy)); 
    }
    else if (reader->ReadFromString(argy,"X_To_P_")) { 
      p_grid->YAxis()->SetScaling(new ATOOLS::X_To_P_Scaling<GridResultType>(argy)); 
    }
    delete reader;
  }

} // end of namespace AMISIC

#endif
