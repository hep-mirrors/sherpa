#ifndef Data_To_Function_C
#define Data_To_Function_C

#include "Data_To_Function.H"
#include "MathTools.H"

#ifdef DEBUG__Data_To_Function
#include "My_IO_Stream.H"
#endif

namespace ATOOLS {

  template <class Argument_Type,class Result_Type>
  Data_To_Function<Argument_Type,Result_Type>::Data_To_Function()
  { Init(); }

  template <class Argument_Type,class Result_Type>
  Data_To_Function<Argument_Type,Result_Type>::
  Data_To_Function(DataToFunctionType &reference)
  { Init(); Import(reference); }

  template <class Argument_Type,class Result_Type>
  Data_To_Function<Argument_Type,Result_Type>::
  Data_To_Function(ArgumentVector *_p_xdata,ResultVector *_p_ydata)
  { Init(); Import(_p_xdata,_p_ydata); }
  
  template <class Argument_Type,class Result_Type>
  Data_To_Function<Argument_Type,Result_Type>::
  Data_To_Function(XYMap *_p_xydata)
  { Init(); Import(_p_xydata); }

  template <class Argument_Type,class Result_Type>
  Data_To_Function<Argument_Type,Result_Type>::~Data_To_Function()
  {
    p_xydata->clear(); delete p_xydata;
    p_yxdata->clear(); delete p_yxdata; 
    delete p_xaxis;
    delete p_yaxis;
  }

  template <class Argument_Type,class Result_Type>
  void Data_To_Function<Argument_Type,Result_Type>::Init()
  {
    p_xaxis = new ATOOLS::Axis<ArgumentType>();
    p_yaxis = new ATOOLS::Axis<ResultType>();
    p_xydata = new XYVector();
    p_yxdata = new YXVector();
    m_acquisitionmode = Interpolation;
    m_interpolationmode = Linear;
  }
  
  template <class Argument_Type,class Result_Type>
  void Data_To_Function<Argument_Type,Result_Type>::
  Import(ArgumentVector *_p_xdata,ResultVector *_p_ydata,bool normal)
  { 
    if (_p_xdata->size()==_p_ydata->size()) {
      Resize(_p_xdata->size());
      if (normal) {
	for (unsigned int i=0; i<_p_xdata->size(); ++i) {
	  (*p_yxdata)[i].second=(*p_xydata)[i].first=(*p_xaxis)((*_p_xdata)[i]); 
	  (*p_xydata)[i].second=(*p_yxdata)[i].first=(*p_yaxis)((*_p_ydata)[i]);
	}
      }
      else {
	for (unsigned int i=0; i<_p_xdata->size(); ++i) {
	  (*p_yxdata)[i].second=(*p_xydata)[i].first=(*p_xaxis)((ArgumentType)(*_p_ydata)[i]); 
	  (*p_xydata)[i].second=(*p_yxdata)[i].first=(*p_yaxis)((ResultType)(*_p_xdata)[i]);
	}
      }
    }
#ifdef DEBUG__Data_To_Function
    std::cout<<"Data_To_Function::Import("<<_p_xdata<<","<<_p_ydata<<") :"<<std::endl;
#endif
    Sort();
  }
  
  template <class Argument_Type,class Result_Type>
  void Data_To_Function<Argument_Type,Result_Type>::
  Import(XYMap *_p_xydata,bool normal)
  { 
    unsigned int i=-1;
    Resize(_p_xydata->size());
    if (normal) {
      for (XYVectorIterator xyit=_p_xydata->begin();xyit!=_p_xydata->end();++xyit) {
	++i;
	(*p_yxdata)[i].second=(*p_xydata)[i].first=(*p_xaxis)((*xyit).first);
	(*p_xydata)[i].second=(*p_yxdata)[i].first=(*p_yaxis)((*xyit).second);
      }
    }
    else {
      for (XYVectorIterator xyit=_p_xydata->begin();xyit!=_p_xydata->end();++xyit) {
	++i;
	(*p_yxdata)[i].second=(*p_xydata)[i].first=(*p_xaxis)((ArgumentType)(*xyit).second);
	(*p_xydata)[i].second=(*p_yxdata)[i].first=(*p_yaxis)((ResultType)(*xyit).first);
      }
    }
#ifdef DEBUG__Data_To_Function
    std::cout<<"Data_To_Function::Import("<<_p_xydata<<") :"<<std::endl;
#endif
    Sort();
  }
  
  template <class Argument_Type,class Result_Type>
  inline void Data_To_Function<Argument_Type,Result_Type>::
  Import(DataToFunctionType &reference)
  { 
    Resize(reference.p_xydata->size());
    for (unsigned int i=0;i<reference.p_xydata->size();++i) {
      ArgumentType x=(*p_xaxis)((*reference.p_xaxis)[(*reference.p_xydata)[i].first]);
      ResultType y=(*p_yaxis)((*reference.p_yaxis)[(*reference.p_xydata)[i].second]);
      (*p_xydata)[i]=XYPair(x,y);
      (*p_yxdata)[i]=YXPair(y,x);
    }
    p_xaxis->SetVariable(reference.p_xaxis->Variable());
    p_yaxis->SetVariable(reference.p_yaxis->Variable());
#ifdef DEBUG__Data_To_Function
    std::cout<<"Data_To_Function::Import("<<&reference<<") :"<<std::endl;
#endif
    Sort();
  }

  template <class Argument_Type,class Result_Type>
  void Data_To_Function<Argument_Type,Result_Type>::
  Export(ArgumentVector *_p_xdata,ResultVector *_p_ydata,bool normal)
  { 
    _p_xdata->clear();
    _p_ydata->clear();
    _p_xdata->resize(p_xydata->size());
    _p_ydata->resize(p_yxdata->size());
    if (normal) {
      for (unsigned int i=0;i<p_xydata->size();++i) {
	(*_p_xdata)[i]=(*p_xaxis)[(*p_xydata)[i].first];   
	(*_p_ydata)[i]=(*p_yaxis)[(*p_xydata)[i].second];   
      }
    }
    else {
      for (unsigned int i=0;i<p_yxdata->size();++i) {
	(*_p_xdata)[i]=(ArgumentType)(*p_yaxis)[(*p_yxdata)[i].first];   
	(*_p_ydata)[i]=(ResultType)(*p_xaxis)[(*p_yxdata)[i].second];   
      }
    }
  }
  
  template <class Argument_Type,class Result_Type>
  void Data_To_Function<Argument_Type,Result_Type>::
  Export(XYMap *_p_xydata,bool normal)
  { 
    _p_xydata->clear();
    if (normal) {
      for (unsigned int i=0;i<p_xydata->size();++i) {
	(*_p_xydata)[(*p_xaxis)[(*p_xydata)[i].first]]=
                     (*p_yaxis)[(*p_xydata)[i].second];   
      }
    }
    else {
      for (unsigned int i=0;i<p_yxdata->size();++i) {
	(*_p_xydata)[(ArgumentType)(*p_yaxis)[(*p_yxdata)[i].first]]=
	               (ResultType)(*p_xaxis)[(*p_yxdata)[i].second];   
      }
    }
  }
  
  template <class Argument_Type,class Result_Type>
  void Data_To_Function<Argument_Type,Result_Type>::SortX()
  {
    XYPair xy;
    bool cont;
#ifdef DEBUG__Data_To_Function
    std::cout<<"Data_To_Function::SortX() :"<<std::endl;
    std::cout<<"   before sorting: p_xdata = [ "; 
    for(XYVectorIterator xyit=p_xydata->begin();xyit!=p_xydata->end();std::cout<<(xyit++)->first<<" ");
    std::cout<<"]"<<std::endl;
#endif
    do {
      cont=false;
      for (unsigned int i=1; i<p_xydata->size(); ++i) {
	if ((*p_xydata)[i].first<(*p_xydata)[i-1].first) {
	  xy=(*p_xydata)[i]; (*p_xydata)[i]=(*p_xydata)[i-1]; (*p_xydata)[i-1]=xy;
	  cont=true;
	}
      }
    } while (cont); 
#ifdef DEBUG__Data_To_Function
    std::cout<<"   after sorting : p_xdata = [ "; 
    for(XYVectorIterator xyit=p_xydata->begin();xyit!=p_xydata->end();std::cout<<(xyit++)->first<<" ");
    std::cout<<"]"<<std::endl;
#endif
  }
  
  template <class Argument_Type,class Result_Type>
  void Data_To_Function<Argument_Type,Result_Type>::SortY()
  {
    YXPair yx;
    bool cont;
#ifdef DEBUG__Data_To_Function
    std::cout<<"Data_To_Function::SortY() :"<<std::endl;
    std::cout<<"   before sorting: p_ydata = [ "; 
    for(YXVectorIterator yxit=p_yxdata->begin();yxit!=p_yxdata->end();std::cout<<(yxit++)->first<<" ");
    std::cout<<"]"<<std::endl;
#endif
    do {
      cont=false;
      for (unsigned int i=1; i<p_yxdata->size(); ++i) {
	if ((*p_yxdata)[i].first<(*p_yxdata)[i-1].first) {
	  yx=(*p_yxdata)[i]; (*p_yxdata)[i]=(*p_yxdata)[i-1]; (*p_yxdata)[i-1]=yx;
	  cont=true;
	}
      }
    } while (cont); 
#ifdef DEBUG__Data_To_Function
    std::cout<<"   after sorting : p_ydata = [ "; 
    for(YXVectorIterator yxit=p_yxdata->begin();yxit!=p_yxdata->end();std::cout<<(yxit++)->first<<" ");
    std::cout<<"]"<<std::endl;
#endif
  }
  
  template <class Argument_Type,class Result_Type>
  unsigned int Data_To_Function<Argument_Type,Result_Type>::
  ClosestX(ArgumentType x,unsigned int &left,unsigned int &right)
  {
#ifdef DEBUG__Data_To_Function
    std::cout<<"Data_To_Function::ClosestX("<<x<<") :"<<std::endl;
    std::cout<<"   p_xdata = [ "; 
    for(XYVectorIterator xyit=p_xydata->begin();xyit!=p_xydata->end();std::cout<<(xyit++)->first<<" ");
    std::cout<<"]"<<std::endl;
#endif
    left=0; right=p_xydata->size()-1;
    unsigned int middle=(right-left)/2;
    if (x<(*p_xydata)[left].first) { 
#ifdef DEBUG__Data_To_Function
      std::cout<<"   value is out of range"<<std::endl;
      std::cout<<"   returning "<<left<<" => x value "<<(*p_xydata)[left].first<<std::endl;
#endif
      right=left+1;
      return left; 
    } 
    else if (x>(*p_xydata)[right].first) { 
#ifdef DEBUG__Data_To_Function
      std::cout<<"   value is out of range"<<std::endl;
      std::cout<<"   returning "<<right<<" => x value "<<(*p_xydata)[right].first<<std::endl;
#endif
      left=right-1;
      return right; 
    }
    do {
      if (x<(*p_xydata)[middle].first) { right=middle; middle=(middle+left)/2; }
      else { left=middle; middle=(right+middle)/2; } 
    } while (((right-middle)!=0)&&((middle-left)!=0));
    if ((*p_xydata)[right].first-x <= x-(*p_xydata)[left].first) { 
#ifdef DEBUG__Data_To_Function
      std::cout<<"   returning "<<right<<" / ["<<left<<","<<right<<"] => x values "
	       <<(*p_xydata)[right].first<<" ["<<(*p_xydata)[left].first
	       <<","<<(*p_xydata)[right].first<<"]"<<std::endl;
#endif
      return right;
    }
    else { 
#ifdef DEBUG__Data_To_Function
      std::cout<<"   returning "<<left<<" / ["<<left<<","<<right<<"] => x values "
	       <<(*p_xydata)[left].first<<" ["<<(*p_xydata)[left].first
	       <<","<<(*p_xydata)[right].first<<"]"<<std::endl;
#endif
      return left; 
    }
  }
  
  template <class Argument_Type,class Result_Type>
  unsigned int Data_To_Function<Argument_Type,Result_Type>::
  ClosestY(ResultType y,unsigned int &left,unsigned int &right)
  {
#ifdef DEBUG__Data_To_Function
    std::cout<<"Data_To_Function::ClosestY("<<y<<") :"<<std::endl;
    std::cout<<"   p_ydata = [ "; 
    for(YXVectorIterator yxit=p_yxdata->begin();yxit!=p_yxdata->end();std::cout<<(yxit++)->first<<" ");
    std::cout<<"]"<<std::endl;
#endif
    left=0; right=p_yxdata->size()-1;
    unsigned int middle=(right-left)/2;
    if (y<(*p_yxdata)[left].first) { 
#ifdef DEBUG__Data_To_Function
      std::cout<<"   value is out of range"<<std::endl;
      std::cout<<"   returning "<<left<<" => y value "<<(*p_yxdata)[left].first<<std::endl;
#endif
      right=left+1; 
      return left; 
    } 
    else if (y>(*p_yxdata)[right].first) { 
#ifdef DEBUG__Data_To_Function
      std::cout<<"   value is out of range"<<std::endl;
      std::cout<<"   returning "<<right<<" => y value "<<(*p_yxdata)[right].first<<std::endl;
#endif
      left=right-1; 
      return right; 
    }
    do {
      if (y<(*p_yxdata)[middle].first) { right=middle; middle=(middle+left)/2; }
      else { left=middle; middle=(right+middle)/2; } 
    } while (((right-middle)!=0)&&((middle-left)!=0));
    if ((*p_yxdata)[right].first-y <= y-(*p_yxdata)[left].first) { 
#ifdef DEBUG__Data_To_Function
      std::cout<<"   returning "<<right<<" / ["<<left<<","<<right<<"] => y values "
	       <<(*p_yxdata)[right].first<<" ["<<(*p_yxdata)[left].first
	       <<","<<(*p_yxdata)[right].first<<"]"<<std::endl;
#endif
      return right; 
    }
    else { 
#ifdef DEBUG__Data_To_Function
      std::cout<<"   returning "<<left<<" / ["<<left<<","<<right<<"] => y values "
	       <<(*p_yxdata)[left].first<<" ["<<(*p_yxdata)[left].first
	       <<","<<(*p_yxdata)[right].first<<"]"<<std::endl;
#endif
      return left; 
    }
  }
  
  template <class Argument_Type,class Result_Type>
  bool Data_To_Function<Argument_Type,Result_Type>::
  DeleteXPoint(ArgumentType _x)
  { 
#ifdef DEBUG__Data_To_Function
    std::cout<<"Data_To_Function::DeleteXPoint("<<_x<<") :"<<std::endl;
    std::cout<<"   before deletion: p_xdata = [ "; 
    for(XYVectorIterator xyit=p_xydata->begin();xyit!=p_xydata->end();std::cout<<(xyit++)->first<<" ");
    std::cout<<"]"<<std::endl;
#endif
    _x=(*p_xaxis)(_x);
    for (XYVectorIterator xyit=p_xydata->begin();xyit!=p_xydata->end();++xyit) 
      if (ATOOLS::IsZero(xyit->first-_x)) { 
	ArgumentType x=xyit->first;
	p_xydata->erase(xyit);
	for (YXVectorIterator yxit=p_yxdata->begin();yxit!=p_yxdata->end();++yxit) 
	  if (yxit->second==x) { 
	    p_yxdata->erase(yxit); 
	    break; 
	  }
#ifdef DEBUG__Data_To_Function
	std::cout<<"   deleted point "<<xyit->first<<std::endl;
	std::cout<<"   after deletion : p_xdata = [ "; 
	for(XYVectorIterator xyit=p_xydata->begin();xyit!=p_xydata->end();std::cout<<(xyit++)->first<<" ");
	std::cout<<"]"<<std::endl;
#endif
	return true;
      }
#ifdef DEBUG__Data_To_Function
    std::cout<<"   did not delete any point"<<std::endl;
#endif 
    return false;
  }
  
  template <class Argument_Type,class Result_Type>
  bool Data_To_Function<Argument_Type,Result_Type>::
  DeleteYPoint(ResultType _y)
  { 
#ifdef DEBUG__Data_To_Function
    std::cout<<"Data_To_Function::DeleteYPoint("<<_y<<") :"<<std::endl;
    std::cout<<"   before deletion: p_ydata = [ "; 
    for(YXVectorIterator yxit=p_yxdata->begin();yxit!=p_yxdata->end();std::cout<<(yxit++)->first<<" ");
    std::cout<<"]"<<std::endl;
#endif
    _y=(*p_yaxis)(_y);
    for (YXVectorIterator yxit=p_yxdata->begin();yxit!=p_yxdata->end();++yxit) 
      if (ATOOLS::IsZero(yxit->first-_y)) { 
	ResultType y=yxit->first;
	p_yxdata->erase(yxit);
	for (XYVectorIterator xyit=p_xydata->begin();xyit!=p_xydata->end();++xyit) 
	  if (xyit->second==y) { 
	    p_xydata->erase(xyit); 
	    break; 
	  }
#ifdef DEBUG__Data_To_Function
	std::cout<<"   deleted point "<<_y<<std::endl;
	std::cout<<"   after deletion : p_ydata = [ "; 
	for(YXVectorIterator yxit=p_yxdata->begin();yxit!=p_yxdata->end();std::cout<<(yxit++)->first<<" ");
	std::cout<<"]"<<std::endl;
#endif
	return true;
      }
#ifdef DEBUG__Data_To_Function
    std::cout<<"   did not delete any point"<<std::endl; 
#endif
    return false;
  }

  template <class Argument_Type,class Result_Type>
  Result_Type Data_To_Function<Argument_Type,Result_Type>::
  InterpolateY(ArgumentType x)
  { 
    switch (m_interpolationmode) {
    case Linear:
      unsigned int left, right;
      ClosestX(x,left,right);
      return LinearY(x,left,right); 
      break;
    case IUnknown:
      break;
    }
    return (ResultType)0;
  }
  
  template <class Argument_Type,class Result_Type>
  Argument_Type Data_To_Function<Argument_Type,Result_Type>::
  InterpolateX(ResultType y)
  { 
    switch (m_interpolationmode) {
    case Linear:
      unsigned int left, right;
      ClosestY(y,left,right);
      return LinearX(y,left,right); 
      break;
    case IUnknown:
      break;
    }
    return (ArgumentType)0;
  }

  template <class Argument_Type,class Result_Type>
  Result_Type Data_To_Function<Argument_Type,Result_Type>::
  Y(ArgumentType x,AcquisitionModeID tempmode)
  { 
    if (p_xydata->size()<2) {
      ATOOLS::msg.Debugging()<<"Data_To_Function::Y("<<x<<","<<tempmode<<"): "
			     <<"Less than 2 data points available."<<std::endl
			     <<"   Returning (ResultType)0."<<std::endl;
      return (ResultType)0.0;
    }
    if (tempmode==AUnknown) tempmode=m_acquisitionmode;
    switch(tempmode) {
    case Interpolation:
      return (*p_yaxis)[InterpolateY((*p_xaxis)(x))]; 
      break;
    case Data:
      return (*p_yaxis)[DataY((*p_xaxis)(x))];
      break;
    case AUnknown:
      break;
    }
    return (ResultType)0;
  }
  
  template <class Argument_Type,class Result_Type>
  Argument_Type Data_To_Function<Argument_Type,Result_Type>::
  X(ResultType y,AcquisitionModeID tempmode)
  { 
    if (p_xydata->size()<2) {
      ATOOLS::msg.Debugging()<<"Data_To_Function::X("<<y<<","<<tempmode<<"): "
			     <<"Less than 2 data points available."<<std::endl
			     <<"   Returning (ArgumentType)0."<<std::endl;
      return (ArgumentType)0.0;
    }
    if (tempmode==AUnknown) tempmode=m_acquisitionmode;
    switch(tempmode) {
    case Interpolation:
      return (*p_xaxis)[InterpolateX((*p_yaxis)(y))]; 
      break;
    case Data:
      return (*p_xaxis)[DataX((*p_yaxis)(y))];
      break;
    case AUnknown:
      break;
    }
    return (ArgumentType)0;
  }
  
  template <class Argument_Type,class Result_Type>
  Result_Type Data_To_Function<Argument_Type,Result_Type>::
  IntegrateY(ArgumentType xmin,ArgumentType xmax)
  { 
    ResultType integral=(ResultType)0.0;
    if (p_xydata->size()<2) return integral;
#ifdef DEBUG__Data_To_Function
    std::cout<<"Data_To_Function::IntegrateY("<<xmin<<","<<xmax<<"): starting integration"<<std::endl;
#endif
    xmin=(*p_xaxis)(xmin);
    xmax=(*p_xaxis)(xmax);
    unsigned int dummy, start, stop;
    ResultType yleft, yright;
    ClosestX(xmin,dummy,start);
    yleft=LinearY(xmin,dummy,start);
    if (xmin<=dummy) start=dummy;
    ClosestX(xmax,stop,dummy);
    yright=LinearY(xmax,stop,dummy);
    if (xmax>=dummy) stop=dummy;
    integral+=((*p_yaxis)[(*p_xydata)[start].second]+(*p_yaxis)[yleft])*
      ((*p_xaxis)[(*p_xydata)[start].first]-(*p_xaxis)[xmin])/(ResultType)2.0;
#ifdef DEBUG__Data_To_Function
    std::cout<<"   integral values are [ "<<integral<<" ";
#endif
    for (unsigned int i=start;i<stop;++i) {
      integral+=((*p_yaxis)[(*p_xydata)[i+1].second]+(*p_yaxis)[(*p_xydata)[i].second])*
	((*p_xaxis)[(*p_xydata)[i+1].first]-(*p_xaxis)[(*p_xydata)[i].first])/(ResultType)2.0;
#ifdef DEBUG__Data_To_Function
      std::cout<<integral<<" ";
#endif
    }
    integral+=((*p_yaxis)[yright]+(*p_yaxis)[(*p_xydata)[stop].second])*
      ((*p_xaxis)[xmax]-(*p_xaxis)[(*p_xydata)[stop].first])/(ResultType)2.0;
#ifdef DEBUG__Data_To_Function
    std::cout<<integral<<" ]"<<std::endl;
#endif
    return integral;
  }

  template <class Argument_Type,class Result_Type>
  Argument_Type Data_To_Function<Argument_Type,Result_Type>::
  IntegrateX(ResultType ymin,ResultType ymax)
  { 
    ArgumentType integral=(ArgumentType)0.0;
    if (p_yxdata->size()<2) return integral;
#ifdef DEBUG__Data_To_Function
    std::cout<<"Data_To_Function::IntegrateX("<<ymin<<","<<ymax<<"): starting integration"<<std::endl;
#endif
    ymin=(*p_yaxis)(ymin);
    ymax=(*p_yaxis)(ymax);
    unsigned int dummy, start, stop;
    ArgumentType xleft, xright;
    ClosestY(ymin,dummy,start);
    xleft=LinearX(ymin,dummy,start);
    if (ymin<=dummy) start=dummy;
    ClosestY(ymax,stop,dummy);
    xright=LinearX(ymax,stop,dummy);
    if (ymax>=dummy) stop=dummy;
    integral+=((*p_xaxis)[(*p_yxdata)[start].second]+(*p_xaxis)[xleft])*
      ((*p_yaxis)[(*p_yxdata)[start].first]-(*p_yaxis)[ymin])/(ArgumentType)2.0;
#ifdef DEBUG__Data_To_Function
    std::cout<<"   integral values are [ "<<integral<<" ";
#endif
    for (unsigned int i=start;i<stop;++i) {
      integral+=((*p_xaxis)[(*p_yxdata)[i+1].second]+(*p_xaxis)[(*p_yxdata)[i].second])*
	((*p_yaxis)[(*p_yxdata)[i+1].first]-(*p_yaxis)[(*p_yxdata)[i].first])/(ArgumentType)2.0;
#ifdef DEBUG__Data_To_Function
      std::cout<<integral<<" ";
#endif
    }
    integral+=((*p_xaxis)[xright]+(*p_xaxis)[(*p_yxdata)[stop].second])*
      ((*p_yaxis)[ymax]-(*p_yaxis)[(*p_yxdata)[stop].first])/(ArgumentType)2.0;
#ifdef DEBUG__Data_To_Function
    std::cout<<integral<<" ]"<<std::endl;
#endif
    return integral;
  }

  template <class Argument_Type,class Result_Type>
  Data_To_Function<Argument_Type,Result_Type> *Data_To_Function<Argument_Type,Result_Type>::
  IntegralY(ArgumentType xmin,ArgumentType xmax)
  { 
    DataToFunctionType *integrated = NULL;
    ResultType integral=(ResultType)0.0;
    if (p_xydata->size()<2) return integrated;
    if (xmin==xmax) {
      xmin=(*p_xaxis)[(*p_xydata)[0].first];
      xmax=(*p_xaxis)[(*p_xydata)[p_xydata->size()-1].first];
    }
    integrated = new DataToFunctionType();
    integrated->XAxis()->SetVariable(ATOOLS::Variable(p_xaxis->Variable().Name()));
    integrated->YAxis()->SetVariable(ATOOLS::Variable(std::string("\\int d")+p_xaxis->Variable().Name()
						      +std::string(" ")+p_yaxis->Variable().Name()));
    integrated->XAxis()->SetScaling(p_xaxis->Scaling()->Name());
    integrated->YAxis()->SetScaling(p_yaxis->Scaling()->Name());
#ifdef DEBUG__Data_To_Function
    std::cout<<"Data_To_Function::IntegralY("<<xmin<<","<<xmax<<"): starting integration"<<std::endl;
#endif
    xmin=(*p_xaxis)(xmin);
    xmax=(*p_xaxis)(xmax);
    unsigned int dummy, start, stop;
    ResultType yleft, yright;
    ClosestX(xmin,dummy,start);
    yleft=LinearY(xmin,dummy,start);
    if (xmin<=dummy) start=dummy;
    ClosestX(xmax,stop,dummy);
    yright=LinearY(xmax,stop,dummy);
    if (xmax>=dummy) stop=dummy;
    integrated->AddPoint((*p_xaxis)[xmin],integral);
    integral+=((*p_yaxis)[(*p_xydata)[start].second]+(*p_yaxis)[yleft])*
      ((*p_xaxis)[(*p_xydata)[start].first]-(*p_xaxis)[xmin])/(ResultType)2.0;
#ifdef DEBUG__Data_To_Function
    std::cout<<"   integral value for first step is ["<<integral<<"]"<<std::endl;
#endif
    for (unsigned int i=start;i<stop;++i) {
      integrated->AddPoint((*p_xaxis)[(*p_xydata)[i].first],integral);
      integral+=((*p_yaxis)[(*p_xydata)[i+1].second]+(*p_yaxis)[(*p_xydata)[i].second])*
	((*p_xaxis)[(*p_xydata)[i+1].first]-(*p_xaxis)[(*p_xydata)[i].first])/(ResultType)2.0;
#ifdef DEBUG__Data_To_Function
      std::cout<<"   integral value for step "<<i<<" is ["<<integral<<"]"<<std::endl;
#endif
    }
    integrated->AddPoint((*p_xaxis)[(*p_xydata)[stop].first],integral);
    integral+=((*p_yaxis)[yright]+(*p_yaxis)[(*p_xydata)[stop].second])*
      ((*p_xaxis)[xmax]-(*p_xaxis)[(*p_xydata)[stop].first])/(ResultType)2.0;
    integrated->AddPoint((*p_xaxis)[xmax],integral);
#ifdef DEBUG__Data_To_Function
    std::cout<<"   integral value for last step is ["<<integral<<"]"<<std::endl;
#endif
    return integrated;
  }

  template <class Argument_Type,class Result_Type>
  Data_To_Function<Argument_Type,Result_Type> *Data_To_Function<Argument_Type,Result_Type>::
  IntegralX(ResultType ymin,ResultType ymax)
  { 
    DataToFunctionType *integrated = NULL;
    ArgumentType integral=(ArgumentType)0.0;
    if (p_yxdata->size()<2) return integrated;
    if (ymin==ymax) {
      ymin=(*p_yaxis)[(*p_yxdata)[0].first];
      ymax=(*p_yaxis)[(*p_yxdata)[p_yxdata->size()-1].first];
    }
    integrated = new DataToFunctionType();
    integrated->XAxis()->SetVariable(ATOOLS::Variable(std::string("\\int d")+p_yaxis->Variable().Name()
						      +std::string(" ")+p_xaxis->Variable().Name()));
    integrated->YAxis()->SetVariable(ATOOLS::Variable(p_yaxis->Variable().Name()));
    integrated->XAxis()->SetScaling(p_xaxis->Scaling()->Name());
    integrated->YAxis()->SetScaling(p_yaxis->Scaling()->Name());
#ifdef DEBUG__Data_To_Function
    std::cout<<"Data_To_Function::IntegralX("<<ymin<<","<<ymax<<"): starting integration"<<std::endl;
#endif
    ymin=(*p_yaxis)(ymin);
    ymax=(*p_yaxis)(ymax);
    unsigned int dummy, start, stop;
    ArgumentType xleft, xright;
    ClosestY(ymin,dummy,start);
    xleft=LinearX(ymin,dummy,start);
    if (ymin<=dummy) start=dummy;
    ClosestY(ymax,stop,dummy);
    xright=LinearX(ymax,stop,dummy);
    if (ymax>=dummy) stop=dummy;
    integrated->AddPoint((*p_yaxis)[ymin],integral);
    integral+=((*p_xaxis)[(*p_yxdata)[start].second]+(*p_xaxis)[xleft])*
      ((*p_yaxis)[(*p_yxdata)[start].first]-(*p_yaxis)[ymin])/(ArgumentType)2.0;
#ifdef DEBUG__Data_To_Function
    std::cout<<"   integral value for first step is ["<<integral<<"]"<<std::endl;
#endif
    for (unsigned int i=start;i<stop;++i) {
      integrated->AddPoint((*p_yaxis)[(*p_yxdata)[i].first],integral);
      integral+=((*p_xaxis)[(*p_yxdata)[i+1].second]+(*p_xaxis)[(*p_yxdata)[i].second])*
	((*p_yaxis)[(*p_yxdata)[i+1].first]-(*p_yaxis)[(*p_yxdata)[i].first])/(ArgumentType)2.0;
      integrated->AddXPoint(integral);
#ifdef DEBUG__Data_To_Function
      std::cout<<"   integral value for step "<<i<<" is ["<<integral<<"]"<<std::endl;
#endif
    }
    integrated->AddPoint(integral,(*p_yaxis)[(*p_yxdata)[stop].first]);
    integral+=((*p_xaxis)[xright]+(*p_xaxis)[(*p_yxdata)[stop].second])*
      ((*p_yaxis)[ymax]-(*p_yaxis)[(*p_yxdata)[stop].first])/(ArgumentType)2.0;
    integrated->AddPoint(integral,(*p_yaxis)[ymax]);
#ifdef DEBUG__Data_To_Function
    std::cout<<"   integral value for last step is ["<<integral<<"]"<<std::endl;
#endif
    return integrated;
  }

  template <class Argument_Type,class Result_Type>
  void Data_To_Function<Argument_Type,Result_Type>::
  ScaleY(ResultType scalefactor)
  { 
    for (YXVectorIterator yxit=p_yxdata->begin();yxit!=p_yxdata->end();++yxit) 
      yxit->first=(*p_yaxis)((*p_yaxis)[yxit->first]*scalefactor);
    for (XYVectorIterator xyit=p_xydata->begin();xyit!=p_xydata->end();++xyit) 
      xyit->second=(*p_yaxis)((*p_yaxis)[xyit->second]*scalefactor);
  }

  template <class Argument_Type,class Result_Type>
  void Data_To_Function<Argument_Type,Result_Type>::
  ScaleX(ArgumentType scalefactor)
  { 
    for (XYVectorIterator xyit=p_xydata->begin();xyit!=p_xydata->end();++xyit) 
      xyit->first=(*p_xaxis)((*p_xaxis)[xyit->first]*scalefactor);
    for (YXVectorIterator yxit=p_yxdata->begin();yxit!=p_yxdata->end();++yxit) 
      yxit->second=(*p_xaxis)((*p_xaxis)[yxit->second]*scalefactor);
  }

  template <class Argument_Type,class Result_Type>
  void Data_To_Function<Argument_Type,Result_Type>::
  MoveY(ResultType distance)
  { 
    for (YXVectorIterator yxit=p_yxdata->begin();yxit!=p_yxdata->end();++yxit) 
      yxit->first=(*p_yaxis)((*p_yaxis)[yxit->first]+distance);
    for (XYVectorIterator xyit=p_xydata->begin();xyit!=p_xydata->end();++xyit) 
      xyit->second=(*p_yaxis)((*p_yaxis)[xyit->second]+distance);
  }

  template <class Argument_Type,class Result_Type>
  void Data_To_Function<Argument_Type,Result_Type>::
  MoveX(ArgumentType distance)
  { 
    for (XYVectorIterator xyit=p_xydata->begin();xyit!=p_xydata->end();++xyit) 
      xyit->first=(*p_xaxis)((*p_xaxis)[xyit->first]+distance);
    for (YXVectorIterator yxit=p_yxdata->begin();yxit!=p_yxdata->end();++yxit) 
      yxit->second=(*p_xaxis)((*p_xaxis)[yxit->second]+distance);
  }

  template <class Argument_Type,class Result_Type>
  Argument_Type  Data_To_Function<Argument_Type,Result_Type>::
  DeltaXMin(ResultType& left,ResultType& right)
  { 
    ArgumentType minimum=(ArgumentType)0.0, cur;
    if (p_xydata->size()>2) {
      XYVectorIterator xyit=p_xydata->begin()+1;
      minimum=(ArgumentType)dabs((*p_xaxis)[xyit->first]-(*p_xaxis)[(xyit-1)->first]);
      left=(*p_yaxis)[(xyit-1)->second];
      right=(*p_yaxis)[xyit->second];
      for (++xyit;xyit!=p_xydata->end();++xyit) {
	cur=(ArgumentType)dabs((*p_xaxis)[xyit->first]-(*p_xaxis)[(xyit-1)->first]);
	if (cur<minimum) {
	  minimum=cur;
	  left=(*p_yaxis)[(xyit-1)->second];
	  right=(*p_yaxis)[xyit->second];
	}
      }
    }
    if (left>right) {
      ResultType store=left;
      left=right;
      right=store;
    }
    return minimum;
  }

  template <class Argument_Type,class Result_Type>
  Argument_Type  Data_To_Function<Argument_Type,Result_Type>::
  DeltaXMax(ResultType& left,ResultType& right)
  { 
    ArgumentType maximum=(ArgumentType)0.0, cur;
    if (p_xydata->size()>2) {
      for (XYVectorIterator xyit=p_xydata->begin()+1;xyit!=p_xydata->end();++xyit) {
	cur=(ArgumentType)dabs((*p_xaxis)[xyit->first]-(*p_xaxis)[(xyit-1)->first]);
	if (cur>maximum) {
	  maximum=cur;
	  left=(*p_yaxis)[(xyit-1)->second];
	  right=(*p_yaxis)[xyit->second];
	}
      }
    }
    if (left>right) {
      ResultType store=left;
      left=right;
      right=store;
    }
    return maximum;
  }

  template <class Argument_Type,class Result_Type>
  Result_Type  Data_To_Function<Argument_Type,Result_Type>::
  DeltaYMin(ArgumentType& left,ArgumentType& right)
  { 
    ResultType minimum=(ResultType)0.0, cur;
    if (p_yxdata->size()>2) {
      YXVectorIterator yxit=p_yxdata->begin()+1;
      minimum=(ResultType)dabs((*p_yaxis)[yxit->first]-(*p_yaxis)[(yxit-1)->first]);
      left=(*p_xaxis)[(yxit-1)->second];
      right=(*p_xaxis)[yxit->second];
      for (++yxit;yxit!=p_yxdata->end();++yxit) {
	cur=(ResultType)dabs((*p_yaxis)[yxit->first]-(*p_yaxis)[(yxit-1)->first]);
	if (cur<minimum) {
	  minimum=cur;
	  left=(*p_xaxis)[(yxit-1)->second];
	  right=(*p_xaxis)[yxit->second];
	}
      }
    }
    if (left>right) {
      ArgumentType store=left;
      left=right;
      right=store;
    }
    return minimum;
  }

  template <class Argument_Type,class Result_Type>
  Result_Type  Data_To_Function<Argument_Type,Result_Type>::
  DeltaYMax(ArgumentType& left,ArgumentType& right)
  { 
    ResultType maximum=(ResultType)0.0, cur;
    if (p_yxdata->size()>2) {
      for (YXVectorIterator yxit=p_yxdata->begin()+1;yxit!=p_yxdata->end();++yxit) {
	cur=(ResultType)dabs((*p_yaxis)[yxit->first]-(*p_yaxis)[(yxit-1)->first]);
	if (cur>maximum) {
	  maximum=cur;
	  left=(*p_xaxis)[(yxit-1)->second];
	  right=(*p_xaxis)[yxit->second];
	}
      }
    }
    if (left>right) {
      ArgumentType store=left;
      left=right;
      right=store;
    }
    return maximum;
  }

} // end of namespace ATOOLS

#endif
