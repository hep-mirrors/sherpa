#ifndef Data_To_Function_C
#define Data_To_Function_C

#include "Data_To_Function.H"
#include "MathTools.H"

#ifdef DEBUG__Data_To_Function
#include <iostream>
#endif

#ifdef PROFILE__Data_To_Function
#include "prof.hh"
#endif

namespace ATOOLS {

  template <class Argument_Type,class Result_Type>
  Data_To_Function<Argument_Type,Result_Type>::Data_To_Function()
  { Init(); }

  template <class Argument_Type,class Result_Type>
  Data_To_Function<Argument_Type,Result_Type>::
  Data_To_Function(const DataToFunctionType &reference,bool adoptscaling)
  { Init(); Import(reference,adoptscaling); }

  template <class Argument_Type,class Result_Type>
  Data_To_Function<Argument_Type,Result_Type>::
  Data_To_Function(const ArgumentVector *_p_xdata,const ResultVector *_p_ydata)
  { Init(); Import(_p_xdata,_p_ydata); }
  
  template <class Argument_Type,class Result_Type>
  Data_To_Function<Argument_Type,Result_Type>::
  Data_To_Function(const XYMap *_p_xydata)
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
    m_monotony = MUnknown;
    m_acquisitionmode = Interpolation;
    m_interpolationmode = Linear;
  }
  
  template <class Argument_Type,class Result_Type>
  bool Data_To_Function<Argument_Type,Result_Type>::
  Import(const ArgumentVector *_p_xdata,const ResultVector *_p_ydata,bool normal)
  { 
#ifdef PROFILE__Data_To_Function
    PROFILE_HERE;
#endif
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
    return Sort();
  }
  
  template <class Argument_Type,class Result_Type>
  bool Data_To_Function<Argument_Type,Result_Type>::
  Import(const XYMap *_p_xydata,bool normal)
  { 
#ifdef PROFILE__Data_To_Function
    PROFILE_HERE;
#endif
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
    return Sort();
  }
  
  template <class Argument_Type,class Result_Type>
  bool Data_To_Function<Argument_Type,Result_Type>::
  Import(const DataToFunctionType &reference,bool adoptscaling)
  { 
#ifdef PROFILE__Data_To_Function
    PROFILE_HERE;
#endif
    if (adoptscaling) {
      p_xaxis->SetScaling(reference.p_xaxis->Scaling()->Name());
      p_yaxis->SetScaling(reference.p_yaxis->Scaling()->Name());
    }
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
    return Sort();
  }

  template <class Argument_Type,class Result_Type>
  void Data_To_Function<Argument_Type,Result_Type>::
  Export(ArgumentVector *_p_xdata,ResultVector *_p_ydata,bool normal)
  { 
#ifdef PROFILE__Data_To_Function
    PROFILE_HERE;
#endif
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
#ifdef PROFILE__Data_To_Function
    PROFILE_HERE;
#endif
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
  bool Data_To_Function<Argument_Type,Result_Type>::SortX()
  {
#ifdef PROFILE__Data_To_Function
    PROFILE_HERE;
#endif
    bool success=true;
#ifdef DEBUG__Data_To_Function
    std::cout<<"Data_To_Function::SortX() :"<<std::endl;
    std::cout<<"   before sorting: p_xdata = [ "; 
    for(XYVectorIterator xyit=p_xydata->begin();xyit!=p_xydata->end();std::cout<<(xyit++)->first<<" ");
    std::cout<<"]"<<std::endl;
#endif
    for (unsigned int i=1; i<p_xydata->size(); ++i) {
      XYPair xy=(*p_xydata)[i];
      unsigned int j;
      for (j=i; j>0;--j) {
	if (xy.first>=(*p_xydata)[j-1].first) break;
	(*p_xydata)[j]=(*p_xydata)[j-1]; 
      }
      if (j!=i) {
	(*p_xydata)[j]=xy;
	switch (m_monotony) {
	case Increasing:
	  if (((*p_xydata)[j+1].second-xy.second)<0) {
	    ATOOLS::msg.Error()<<"Data_To_Function::SortX(): Monotony violation! "<<std::endl
			       <<"   Removing point at x = "<<(*p_xaxis)[xy.first]<<std::endl;
	    DeleteXPoint((*p_xaxis)[xy.first]);
	    success=false;
	  }
	  break;
	case Decreasing:
	  if (((*p_xydata)[j+1].second-xy.second)>0) {
	    ATOOLS::msg.Error()<<"Data_To_Function::SortX(): Monotony violation! "<<std::endl
			       <<"   Removing point at x = "<<(*p_xaxis)[xy.first]<<std::endl;
	    DeleteXPoint((*p_xaxis)[xy.first]);
	    success=false;
	  }
	  break;
	case MUnknown:
	  if (((*p_xydata)[j+1].second-xy.second)>0) SetMonotony(Increasing);
	  else SetMonotony(Decreasing);
	case None:
	  break;
	}
      }
    }
#ifdef DEBUG__Data_To_Function
    std::cout<<"   after sorting : p_xdata = [ "; 
    for(XYVectorIterator xyit=p_xydata->begin();xyit!=p_xydata->end();std::cout<<(xyit++)->first<<" ");
    std::cout<<"]"<<std::endl;
#endif
    return success;
  }
  
  template <class Argument_Type,class Result_Type>
  bool Data_To_Function<Argument_Type,Result_Type>::SortY()
  {
#ifdef PROFILE__Data_To_Function
    PROFILE_HERE;
#endif
    bool success=true;
#ifdef DEBUG__Data_To_Function
    std::cout<<"Data_To_Function::SortY() :"<<std::endl;
    std::cout<<"   before sorting: p_ydata = [ "; 
    for(YXVectorIterator yxit=p_yxdata->begin();yxit!=p_yxdata->end();std::cout<<(yxit++)->first<<" ");
    std::cout<<"]"<<std::endl;
#endif
    for (unsigned int i=1; i<p_yxdata->size(); ++i) {
      YXPair yx=(*p_yxdata)[i];
      unsigned int j;
      for (j=i; j>0;--j) {
	if (yx.first>=(*p_yxdata)[j-1].first) break;
	(*p_yxdata)[j]=(*p_yxdata)[j-1]; 
      }
      if (j!=i) {
	(*p_yxdata)[j]=yx;
	switch (m_monotony) {
	case Increasing:
	  if (((*p_yxdata)[j+1].second-yx.second)<0) {
	    ATOOLS::msg.Error()<<"Data_To_Function::SortX(): Monotony violation! "<<std::endl
			       <<"   Removing point at y = "<<(*p_yaxis)[yx.first]<<std::endl;
	    DeleteXPoint((*p_yaxis)[yx.first]);
	    success=false;
	  }
	  break;
	case Decreasing:
	  if (((*p_yxdata)[j+1].second-yx.second)>0) {
	    ATOOLS::msg.Error()<<"Data_To_Function::SortX(): Monotony violation! "<<std::endl
			       <<"   Removing point at y = "<<(*p_yaxis)[yx.first]<<std::endl;
	    DeleteXPoint((*p_yaxis)[yx.first]);
	    success=false;
	  }
	  break;
	case MUnknown:
	  if (((*p_yxdata)[j+1].second-yx.second)>0) SetMonotony(Increasing);
	  else SetMonotony(Decreasing);
	case None:
	  break;
	}
      }
    }
#ifdef DEBUG__Data_To_Function
    std::cout<<"   after sorting : p_ydata = [ "; 
    for(YXVectorIterator yxit=p_yxdata->begin();yxit!=p_yxdata->end();std::cout<<(yxit++)->first<<" ");
    std::cout<<"]"<<std::endl;
#endif
    return success;
  }
  
  template <class Argument_Type,class Result_Type>
  unsigned int Data_To_Function<Argument_Type,Result_Type>::
  ClosestX(ArgumentType x,unsigned int &left,unsigned int &right)
  {
#ifdef PROFILE__Data_To_Function
    PROFILE_HERE;
#endif
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
#ifdef PROFILE__Data_To_Function
    PROFILE_HERE;
#endif
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
#ifdef PROFILE__Data_To_Function
    PROFILE_HERE;
#endif
#ifdef DEBUG__Data_To_Function
    std::cout<<"Data_To_Function::DeleteXPoint("<<_x<<") :"<<std::endl;
    std::cout<<"   before deletion: p_xdata = [ "; 
    for (XYVectorIterator xyit=p_xydata->begin();xyit!=p_xydata->end();std::cout<<(xyit++)->first<<" ");
    std::cout<<"]"<<std::endl;
#endif
    _x=(*p_xaxis)(_x);
    for (XYVectorIterator xyit=p_xydata->begin();xyit!=p_xydata->end();++xyit) 
      if (ATOOLS::IsEqual(xyit->first,_x)) { 
#ifdef DEBUG__Data_To_Function
	std::cout<<"   deleted point "<<xyit->first<<std::endl;
#endif
	ArgumentType x=xyit->first;
	p_xydata->erase(xyit);
	for (YXVectorIterator yxit=p_yxdata->begin();yxit!=p_yxdata->end();++yxit) 
	  if (yxit->second==x) { 
	    p_yxdata->erase(yxit); 
	    break; 
	  }
#ifdef DEBUG__Data_To_Function
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
#ifdef PROFILE__Data_To_Function
    PROFILE_HERE;
#endif
#ifdef DEBUG__Data_To_Function
    std::cout<<"Data_To_Function::DeleteYPoint("<<_y<<") :"<<std::endl;
    std::cout<<"   before deletion: p_ydata = [ "; 
    for (YXVectorIterator yxit=p_yxdata->begin();yxit!=p_yxdata->end();std::cout<<(yxit++)->first<<" ");
    std::cout<<"]"<<std::endl;
#endif
    _y=(*p_yaxis)(_y);
    for (YXVectorIterator yxit=p_yxdata->begin();yxit!=p_yxdata->end();++yxit) 
      if (ATOOLS::IsEqual(yxit->first,_y)) { 
#ifdef DEBUG__Data_To_Function
	std::cout<<"   deleted point "<<_y<<std::endl;
#endif
	ResultType y=yxit->first;
	p_yxdata->erase(yxit);
	for (XYVectorIterator xyit=p_xydata->begin();xyit!=p_xydata->end();++xyit) 
	  if (xyit->second==y) { 
	    p_xydata->erase(xyit); 
	    break; 
	  }
#ifdef DEBUG__Data_To_Function
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
  void Data_To_Function<Argument_Type,Result_Type>::
  RescaleX(const std::string scalename)
  { 
#ifdef PROFILE__Data_To_Function
    PROFILE_HERE;
#endif
#ifdef DEBUG__Data_To_Function
    std::cout<<"Data_To_Function::RescaleX("<<scalename<<") :"<<std::endl;
    std::cout<<"   before rescaling: p_xdata = [ "; 
    for (XYVectorIterator xyit=p_xydata->begin();xyit!=p_xydata->end();std::cout<<(xyit++)->first<<" ");
    std::cout<<"]"<<std::endl;
#endif
    Axis<ArgumentType> *xaxis = new Axis<ArgumentType>();
    xaxis->SetScaling(p_xaxis->Scaling()->Name());
    xaxis->SetScalingMode(p_xaxis->ScalingMode());
    p_xaxis->SetScaling(scalename);
    for (XYVectorIterator xyit=p_xydata->begin();xyit!=p_xydata->end();++xyit) {
      xyit->first=(*p_xaxis)((*xaxis)[xyit->first]);
    }
    for (YXVectorIterator yxit=p_yxdata->begin();yxit!=p_yxdata->end();++yxit) {
      yxit->second=(*p_xaxis)((*xaxis)[yxit->second]);
    }
#ifdef DEBUG__Data_To_Function
    std::cout<<"   after rescaling : p_xdata = [ "; 
    for(XYVectorIterator xyit=p_xydata->begin();xyit!=p_xydata->end();std::cout<<(xyit++)->first<<" ");
    std::cout<<"]"<<std::endl;
#endif
  }

  template <class Argument_Type,class Result_Type>
  void Data_To_Function<Argument_Type,Result_Type>::
  RescaleY(const std::string scalename)
  { 
#ifdef PROFILE__Data_To_Function
    PROFILE_HERE;
#endif
#ifdef DEBUG__Data_To_Function
    std::cout<<"Data_To_Function::RescaleY("<<scalename<<") :"<<std::endl;
    std::cout<<"   before rescaling: p_ydata = [ "; 
    for (YXVectorIterator yxit=p_yxdata->begin();yxit!=p_yxdata->end();std::cout<<(yxit++)->first<<" ");
    std::cout<<"]"<<std::endl;
#endif
    Axis<ResultType> *yaxis = new Axis<ResultType>();
    yaxis->SetScaling(p_yaxis->Scaling()->Name());
    yaxis->SetScalingMode(p_yaxis->ScalingMode());
    p_yaxis->SetScaling(scalename);
    for (YXVectorIterator yxit=p_yxdata->begin();yxit!=p_yxdata->end();++yxit) {
      yxit->first=(*p_yaxis)((*yaxis)[yxit->first]);
    }
    for (XYVectorIterator xyit=p_xydata->begin();xyit!=p_xydata->end();++xyit) {
      xyit->second=(*p_yaxis)((*yaxis)[xyit->second]);
    }
#ifdef DEBUG__Data_To_Function
    std::cout<<"   after rescaling : p_ydata = [ "; 
    for(YXVectorIterator yxit=p_yxdata->begin();yxit!=p_yxdata->end();std::cout<<(yxit++)->first<<" ");
    std::cout<<"]"<<std::endl;
#endif
  }

  template <class Argument_Type,class Result_Type>
  Result_Type Data_To_Function<Argument_Type,Result_Type>::
  InterpolateY(ArgumentType x)
  { 
#ifdef PROFILE__Data_To_Function
    PROFILE_HERE;
#endif
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
#ifdef PROFILE__Data_To_Function
    PROFILE_HERE;
#endif
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
#ifdef PROFILE__Data_To_Function
    PROFILE_HERE;
#endif
    if (p_xydata->size()<2) {
      return (ResultType)0.0;
    }
    if (tempmode==AUnknown) tempmode=m_acquisitionmode;
    switch(tempmode) {
    case Interpolation:
      return (*p_yaxis)[InterpolateY((*p_xaxis)(x))]; 
      break;
    case LowerData:
      return (*p_yaxis)[DataY((*p_xaxis)(x),LowerData)];
      break;
    case Data:
      return (*p_yaxis)[DataY((*p_xaxis)(x),Data)];
      break;
    case UpperData:
      return (*p_yaxis)[DataY((*p_xaxis)(x),UpperData)];
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
#ifdef PROFILE__Data_To_Function
    PROFILE_HERE;
#endif
    if (p_xydata->size()<2) {
      return (ArgumentType)0.0;
    }
    if (tempmode==AUnknown) tempmode=m_acquisitionmode;
    switch(tempmode) {
    case Interpolation:
      return (*p_xaxis)[InterpolateX((*p_yaxis)(y))]; 
      break;
    case LowerData:
      return (*p_xaxis)[DataX((*p_yaxis)(y),LowerData)];
      break;
    case Data:
      return (*p_xaxis)[DataX((*p_yaxis)(y),Data)];
      break;
    case UpperData:
      return (*p_xaxis)[DataX((*p_yaxis)(y),UpperData)];
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
#ifdef PROFILE__Data_To_Function
    PROFILE_HERE;
#endif
    ResultType integral=(ResultType)0.0;
    if (p_xydata->size()<2) return integral;
    if (xmin==xmax) {
      xmin=XMin();
      xmax=XMax();
    }
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
    std::cout<<"   integral value for first step is ["<<integral<<"]"<<std::endl;
#endif
    for (unsigned int i=start;i<stop;++i) {
      integral+=((*p_yaxis)[(*p_xydata)[i+1].second]+(*p_yaxis)[(*p_xydata)[i].second])*
	((*p_xaxis)[(*p_xydata)[i+1].first]-(*p_xaxis)[(*p_xydata)[i].first])/(ResultType)2.0;
#ifdef DEBUG__Data_To_Function
      std::cout<<"   integral value for step "<<i<<" is ["<<integral<<"]"<<std::endl;
#endif
    }
    integral+=((*p_yaxis)[yright]+(*p_yaxis)[(*p_xydata)[stop].second])*
      ((*p_xaxis)[xmax]-(*p_xaxis)[(*p_xydata)[stop].first])/(ResultType)2.0;
#ifdef DEBUG__Data_To_Function
    std::cout<<"   integral value for last step is ["<<integral<<"]"<<std::endl;
#endif
    return integral;
  }

  template <class Argument_Type,class Result_Type>
  Argument_Type Data_To_Function<Argument_Type,Result_Type>::
  IntegrateX(ResultType ymin,ResultType ymax)
  { 
#ifdef PROFILE__Data_To_Function
    PROFILE_HERE;
#endif
    ArgumentType integral=(ArgumentType)0.0;
    if (p_yxdata->size()<2) return integral;
    if (ymin==ymax) {
      ymin=YMin();
      ymax=YMax();
    }
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
    std::cout<<"   integral value for first step is ["<<integral<<"]"<<std::endl;
#endif
    for (unsigned int i=start;i<stop;++i) {
      integral+=((*p_xaxis)[(*p_yxdata)[i+1].second]+(*p_xaxis)[(*p_yxdata)[i].second])*
	((*p_yaxis)[(*p_yxdata)[i+1].first]-(*p_yaxis)[(*p_yxdata)[i].first])/(ArgumentType)2.0;
#ifdef DEBUG__Data_To_Function
      std::cout<<"   integral value for step "<<i<<" is ["<<integral<<"]"<<std::endl;
#endif
    }
    integral+=((*p_xaxis)[xright]+(*p_xaxis)[(*p_yxdata)[stop].second])*
      ((*p_yaxis)[ymax]-(*p_yaxis)[(*p_yxdata)[stop].first])/(ArgumentType)2.0;
#ifdef DEBUG__Data_To_Function
    std::cout<<"   integral value for last step is ["<<integral<<"]"<<std::endl;
#endif
    return integral;
  }

  template <class Argument_Type,class Result_Type>
  Data_To_Function<Argument_Type,Result_Type> *Data_To_Function<Argument_Type,Result_Type>::
  IntegralY(ArgumentType xmin,ArgumentType xmax,std::string xscaling,std::string yscaling,bool forward)
  { 
#ifdef PROFILE__Data_To_Function
    PROFILE_HERE;
#endif
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
    if (xscaling==nullstring) xscaling=p_xaxis->Scaling()->Name();
    if (yscaling==nullstring) yscaling=p_yaxis->Scaling()->Name();
    integrated->XAxis()->SetScaling(xscaling);
    integrated->YAxis()->SetScaling(yscaling);
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
    if (forward) {
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
    }
    else {
      integrated->AddPoint((*p_xaxis)[xmax],integral);
      integral+=((*p_yaxis)[(*p_xydata)[stop].second]+(*p_yaxis)[yright])*
	((*p_xaxis)[xmax]-(*p_xaxis)[(*p_xydata)[stop].first])/(ResultType)2.0;
#ifdef DEBUG__Data_To_Function
      std::cout<<"   integral value for first step is ["<<integral<<"]"<<std::endl;
#endif
      for (unsigned int i=stop;i>start;--i) {
	integrated->AddPoint((*p_xaxis)[(*p_xydata)[i].first],integral);
	integral+=((*p_yaxis)[(*p_xydata)[i].second]+(*p_yaxis)[(*p_xydata)[i-1].second])*
	  ((*p_xaxis)[(*p_xydata)[i].first]-(*p_xaxis)[(*p_xydata)[i-1].first])/(ResultType)2.0;
#ifdef DEBUG__Data_To_Function
	std::cout<<"   integral value for step "<<i<<" is ["<<integral<<"]"<<std::endl;
#endif
      }
      integrated->AddPoint((*p_xaxis)[(*p_xydata)[start].first],integral);
      integral+=((*p_yaxis)[yleft]+(*p_yaxis)[(*p_xydata)[start].second])*
	((*p_xaxis)[(*p_xydata)[start].first]-(*p_xaxis)[xmin])/(ResultType)2.0;
      integrated->AddPoint((*p_xaxis)[xmin],integral);
#ifdef DEBUG__Data_To_Function
      std::cout<<"   integral value for last step is ["<<integral<<"]"<<std::endl;
#endif
    }
    return integrated;
  }

  template <class Argument_Type,class Result_Type>
  Data_To_Function<Argument_Type,Result_Type> *Data_To_Function<Argument_Type,Result_Type>::
  IntegralX(ResultType ymin,ResultType ymax,std::string yscaling,std::string xscaling,bool forward)
  { 
#ifdef PROFILE__Data_To_Function
    PROFILE_HERE;
#endif
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
    if (yscaling==nullstring) yscaling=p_yaxis->Scaling()->Name();
    if (xscaling==nullstring) xscaling=p_xaxis->Scaling()->Name();
    integrated->YAxis()->SetScaling(yscaling);
    integrated->XAxis()->SetScaling(xscaling);
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
    if (forward) {
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
    }
    else {
      integrated->AddPoint((*p_yaxis)[ymax],integral);
      integral+=((*p_xaxis)[(*p_yxdata)[stop].second]+(*p_xaxis)[xright])*
	((*p_yaxis)[ymax]-(*p_yaxis)[(*p_yxdata)[stop].first])/(ArgumentType)2.0;
#ifdef DEBUG__Data_To_Function
      std::cout<<"   integral value for first step is ["<<integral<<"]"<<std::endl;
#endif
      for (unsigned int i=stop;i>start;--i) {
	integrated->AddPoint((*p_yaxis)[(*p_yxdata)[i].first],integral);
	integral+=((*p_xaxis)[(*p_yxdata)[i].second]+(*p_xaxis)[(*p_yxdata)[i-1].second])*
	  ((*p_yaxis)[(*p_yxdata)[i].first]-(*p_yaxis)[(*p_yxdata)[i-1].first])/(ArgumentType)2.0;
#ifdef DEBUG__Data_To_Function
	std::cout<<"   integral value for step "<<i<<" is ["<<integral<<"]"<<std::endl;
#endif
      }
      integrated->AddPoint(integral,(*p_yaxis)[(*p_yxdata)[start].first]);
      integral+=((*p_xaxis)[xleft]+(*p_xaxis)[(*p_yxdata)[start].second])*
	((*p_yaxis)[(*p_yxdata)[start].first]-(*p_yaxis)[ymin])/(ArgumentType)2.0;
      integrated->AddPoint(integral,(*p_yaxis)[ymin]);
#ifdef DEBUG__Data_To_Function
      std::cout<<"   integral value for last step is ["<<integral<<"]"<<std::endl;
#endif
    }
    return integrated;
  }

  template <class Argument_Type,class Result_Type>
  Result_Type Data_To_Function<Argument_Type,Result_Type>::DifferentialY(ArgumentType x)
  { 
#ifdef PROFILE__Data_To_Function
    PROFILE_HERE;
#endif
    unsigned int left, right;
    ClosestX(x,left,right);
    return ((*p_yaxis)[(*p_xydata)[right].second]-(*p_yaxis)[(*p_xydata)[left].second])/(ResultType)
      ((*p_xaxis)[(*p_xydata)[right].first]-(*p_xaxis)[(*p_xydata)[left].first]);	
  }
  
  template <class Argument_Type,class Result_Type>
  Argument_Type Data_To_Function<Argument_Type,Result_Type>::DifferentialX(ResultType y)
  { 
#ifdef PROFILE__Data_To_Function
    PROFILE_HERE;
#endif
    unsigned int left, right;
    ClosestY(y,left,right);
    return ((*p_xaxis)[(*p_yxdata)[right].second]-(*p_xaxis)[(*p_yxdata)[left].second])/(ResultType)
           ((*p_yaxis)[(*p_yxdata)[right].first]-(*p_yaxis)[(*p_yxdata)[left].first]);	
  }
  
  template <class Argument_Type,class Result_Type>
  Data_To_Function<Argument_Type,Result_Type> *Data_To_Function<Argument_Type,Result_Type>::
  DerivativeY(ArgumentType xmin,ArgumentType xmax,std::string xscaling,std::string yscaling)
  { 
#ifdef PROFILE__Data_To_Function
    PROFILE_HERE;
#endif
    DataToFunctionType *derivated = NULL;
    ResultType differential;
    if (p_xydata->size()<2) return derivated;
    if (xmin==xmax) {
      xmin=(*p_xaxis)[(*p_xydata)[0].first];
      xmax=(*p_xaxis)[(*p_xydata)[p_xydata->size()-1].first];
    }
    derivated = new DataToFunctionType();
    derivated->XAxis()->SetVariable(ATOOLS::Variable(p_xaxis->Variable().Name()));
    derivated->YAxis()->SetVariable(
      ATOOLS::Variable(std::string("\\frac {\\partial ")+p_yaxis->Variable().Name()
		       +std::string("}{\\partial ")+p_xaxis->Variable().Name()+std::string("}")));
    if (xscaling==nullstring) xscaling=p_xaxis->Scaling()->Name();
    if (yscaling==nullstring) yscaling=p_yaxis->Scaling()->Name();
    derivated->XAxis()->SetScaling(xscaling);
    derivated->YAxis()->SetScaling(yscaling);
#ifdef DEBUG__Data_To_Function
    std::cout<<"Data_To_Function::DerivativeY("<<xmin<<","<<xmax<<"): "<<std::endl;
#endif
    xmin=(*p_xaxis)(xmin);
    xmax=(*p_xaxis)(xmax);
    unsigned int dummy, start, stop, next;
    ResultType yleft, yright;
    ClosestX(xmin,dummy,start);
    yleft=LinearY(xmin,dummy,start);
    if (xmin<=dummy) start=dummy;
    ClosestX(xmax,stop,dummy);
    yright=LinearY(xmax,stop,next);
    if (xmax>=next) stop=next;
    if (xmin!=(*p_xydata)[start].first) {
      differential=((*p_yaxis)[(*p_xydata)[start].second]-(*p_yaxis)[yleft])/(ResultType)
	((*p_xaxis)[(*p_xydata)[start].first]-(*p_xaxis)[xmin]);
      derivated->AddPoint((*p_xaxis)[xmin],differential);
#ifdef DEBUG__Data_To_Function
      std::cout<<"   derivative for first step is ["<<differential<<"]"<<std::endl;
#endif
    }
    for (unsigned int i=start;i<stop;++i) {
      differential=((*p_yaxis)[(*p_xydata)[i+1].second]-(*p_yaxis)[(*p_xydata)[i].second])/(ResultType)
	((*p_xaxis)[(*p_xydata)[i+1].first]-(*p_xaxis)[(*p_xydata)[i].first]);
      derivated->AddPoint((*p_xaxis)[(*p_xydata)[i].first],differential);
#ifdef DEBUG__Data_To_Function
      std::cout<<"   derivative for step "<<i<<" is ["<<differential<<"]"<<std::endl;
#endif
    }
    if ((*p_xydata)[stop].first!=xmax) {
      differential=((*p_yaxis)[yright]-(*p_yaxis)[(*p_xydata)[stop].second])/(ResultType)
	((*p_xaxis)[xmax]-(*p_xaxis)[(*p_xydata)[stop].first]);
      derivated->AddPoint((*p_xaxis)[(*p_xydata)[stop].first],differential);
#ifdef DEBUG__Data_To_Function
      std::cout<<"   derivative for step "<<stop<<" is ["<<differential<<"]"<<std::endl;
#endif
    }
    if (next<p_xydata->size()) {
      if (xmax!=(*p_xydata)[next].first) {
	differential=((*p_yaxis)[(*p_xydata)[next].second]-(*p_yaxis)[yright])/(ResultType)
	  ((*p_xaxis)[(*p_xydata)[next].first]-(*p_xaxis)[xmax]);
	derivated->AddPoint((*p_xaxis)[xmax],differential);
#ifdef DEBUG__Data_To_Function
	std::cout<<"   derivative for last step is ["<<differential<<"]"<<std::endl;
#endif
      }
    }
    return derivated;
  }

  template <class Argument_Type,class Result_Type>
  Data_To_Function<Argument_Type,Result_Type> *Data_To_Function<Argument_Type,Result_Type>::
  DerivativeX(ResultType ymin,ResultType ymax,std::string yscaling,std::string xscaling)
  { 
#ifdef PROFILE__Data_To_Function
    PROFILE_HERE;
#endif
    DataToFunctionType *derivated = NULL;
    ArgumentType differential;
    if (p_yxdata->size()<2) return derivated;
    if (ymin==ymax) {
      ymin=(*p_yaxis)[(*p_yxdata)[0].first];
      ymax=(*p_yaxis)[(*p_yxdata)[p_yxdata->size()-1].first];
    }
    derivated = new DataToFunctionType();
    derivated->YAxis()->SetVariable(ATOOLS::Variable(p_yaxis->Variable().Name()));
    derivated->XAxis()->SetVariable(
      ATOOLS::Variable(std::string("\\frac {\\partial ")+p_xaxis->Variable().Name()
		       +std::string("}{\\partial ")+p_yaxis->Variable().Name()+std::string("}")));
    if (yscaling==nullstring) yscaling=p_yaxis->Scaling()->Name();
    if (xscaling==nullstring) xscaling=p_xaxis->Scaling()->Name();
    derivated->YAxis()->SetScaling(yscaling);
    derivated->XAxis()->SetScaling(xscaling);
#ifdef DEBUG__Data_To_Function
    std::cout<<"Data_To_Function::DerivativeY("<<ymin<<","<<ymax<<"): "<<std::endl;
#endif
    ymin=(*p_yaxis)(ymin);
    ymax=(*p_yaxis)(ymax);
    unsigned int dummy, start, stop, next;
    ArgumentType xleft, xright;
    ClosestY(ymin,dummy,start);
    xleft=LinearY(ymin,dummy,start);
    if (ymin<=dummy) start=dummy;
    ClosestY(ymax,stop,dummy);
    xright=LinearY(ymax,stop,next);
    if (ymax>=next) stop=next;
    if (ymin!=(*p_yxdata)[start].first) {
      differential=((*p_xaxis)[(*p_yxdata)[start].second]-(*p_xaxis)[xleft])/(ArgumentType)
	((*p_yaxis)[(*p_yxdata)[start].first]-(*p_yaxis)[ymin]);
      derivated->AddPoint(differential,(*p_yaxis)[ymin]);
#ifdef DEBUG__Data_To_Function
      std::cout<<"   derivative for first step is ["<<differential<<"]"<<std::endl;
#endif
    }
    for (unsigned int i=start;i<stop;++i) {
      differential=((*p_xaxis)[(*p_yxdata)[i+1].second]-(*p_xaxis)[(*p_yxdata)[i].second])/(ArgumentType)
	((*p_yaxis)[(*p_yxdata)[i+1].first]-(*p_yaxis)[(*p_yxdata)[i].first]);
      derivated->AddPoint(differential,(*p_yaxis)[(*p_yxdata)[i].first]);
#ifdef DEBUG__Data_To_Function
      std::cout<<"   derivative for step "<<i<<" is ["<<differential<<"]"<<std::endl;
#endif
    }
    if ((*p_yxdata)[stop].first!=ymax) {
      differential=((*p_xaxis)[xright]-(*p_xaxis)[(*p_yxdata)[stop].second])/(ArgumentType)
	((*p_yaxis)[ymax]-(*p_yaxis)[(*p_yxdata)[stop].first]);
      derivated->AddPoint(differential,(*p_yaxis)[(*p_yxdata)[stop].first]);
#ifdef DEBUG__Data_To_Function
      std::cout<<"   derivative for step "<<stop<<" is ["<<differential<<"]"<<std::endl;
#endif
    }
    if (next<p_yxdata->size()) {
      if (ymax!=(*p_yxdata)[next].first) {
	differential=((*p_xaxis)[(*p_yxdata)[next].second]-(*p_xaxis)[xright])/(ArgumentType)
	  ((*p_yaxis)[(*p_yxdata)[next].first]-(*p_yaxis)[ymax]);
	derivated->AddPoint(differential,(*p_yaxis)[ymax]);
#ifdef DEBUG__Data_To_Function
	std::cout<<"   derivative for last step is ["<<differential<<"]"<<std::endl;
#endif
      }
    }
    return derivated;
  }

  template <class Argument_Type,class Result_Type>
  void Data_To_Function<Argument_Type,Result_Type>::
  ScaleY(ResultType scalefactor)
  { 
#ifdef PROFILE__Data_To_Function
    PROFILE_HERE;
#endif
    for (YXVectorIterator yxit=p_yxdata->begin();yxit!=p_yxdata->end();++yxit) 
      yxit->first=(*p_yaxis)((*p_yaxis)[yxit->first]*scalefactor);
    for (XYVectorIterator xyit=p_xydata->begin();xyit!=p_xydata->end();++xyit) 
      xyit->second=(*p_yaxis)((*p_yaxis)[xyit->second]*scalefactor);
  }

  template <class Argument_Type,class Result_Type>
  void Data_To_Function<Argument_Type,Result_Type>::
  ScaleX(ArgumentType scalefactor)
  { 
#ifdef PROFILE__Data_To_Function
    PROFILE_HERE;
#endif
    for (XYVectorIterator xyit=p_xydata->begin();xyit!=p_xydata->end();++xyit) 
      xyit->first=(*p_xaxis)((*p_xaxis)[xyit->first]*scalefactor);
    for (YXVectorIterator yxit=p_yxdata->begin();yxit!=p_yxdata->end();++yxit) 
      yxit->second=(*p_xaxis)((*p_xaxis)[yxit->second]*scalefactor);
  }

  template <class Argument_Type,class Result_Type>
  void Data_To_Function<Argument_Type,Result_Type>::
  MoveY(ResultType distance)
  { 
#ifdef PROFILE__Data_To_Function
    PROFILE_HERE;
#endif
    for (YXVectorIterator yxit=p_yxdata->begin();yxit!=p_yxdata->end();++yxit) 
      yxit->first=(*p_yaxis)((*p_yaxis)[yxit->first]+distance);
    for (XYVectorIterator xyit=p_xydata->begin();xyit!=p_xydata->end();++xyit) 
      xyit->second=(*p_yaxis)((*p_yaxis)[xyit->second]+distance);
  }

  template <class Argument_Type,class Result_Type>
  void Data_To_Function<Argument_Type,Result_Type>::
  MoveX(ArgumentType distance)
  { 
#ifdef PROFILE__Data_To_Function
    PROFILE_HERE;
#endif
    for (XYVectorIterator xyit=p_xydata->begin();xyit!=p_xydata->end();++xyit) 
      xyit->first=(*p_xaxis)((*p_xaxis)[xyit->first]+distance);
    for (YXVectorIterator yxit=p_yxdata->begin();yxit!=p_yxdata->end();++yxit) 
      yxit->second=(*p_xaxis)((*p_xaxis)[yxit->second]+distance);
  }

  template <class Argument_Type,class Result_Type>
  Argument_Type  Data_To_Function<Argument_Type,Result_Type>::
  DeltaXMin(ResultType& left,ResultType& right)
  { 
#ifdef PROFILE__Data_To_Function
    PROFILE_HERE;
#endif
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
#ifdef PROFILE__Data_To_Function
    PROFILE_HERE;
#endif
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
#ifdef PROFILE__Data_To_Function
    PROFILE_HERE;
#endif
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
#ifdef PROFILE__Data_To_Function
    PROFILE_HERE;
#endif
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
