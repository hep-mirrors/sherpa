#include "Data_To_Function.H"

#include "MathTools.H"
#include "Message.H"

#ifdef PROFILE__Data_To_Function
#include "prof.hh"
#else 
#define PROFILE_HERE
#endif

using namespace ATOOLS;

template <class Argument_Type,class Result_Type>
Data_To_Function<Argument_Type,Result_Type>::Data_To_Function()
{ 
  Init(); 
}

template <class Argument_Type,class Result_Type>
Data_To_Function<Argument_Type,Result_Type>::
Data_To_Function(const DataToFunctionType &reference,bool adoptscaling)
{ 
  Init(); 
  Import(reference,adoptscaling); 
}

template <class Argument_Type,class Result_Type>
Data_To_Function<Argument_Type,Result_Type>::
Data_To_Function(const ArgumentVector *_p_xdata,const ResultVector *_p_ydata)
{ 
  Init(); 
  Import(_p_xdata,_p_ydata); 
}
  
template <class Argument_Type,class Result_Type>
Data_To_Function<Argument_Type,Result_Type>::~Data_To_Function()
{
  delete p_xydata;
  delete p_yxdata; 
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
	(*p_yxdata)[i].second=(*p_xydata)[i].first=(*p_xaxis)((*_p_ydata)[i]); 
	(*p_xydata)[i].second=(*p_yxdata)[i].first=(*p_yaxis)((*_p_xdata)[i]);
      }
    }
  }
  return Sort();
}
  
template <class Argument_Type,class Result_Type>
bool Data_To_Function<Argument_Type,Result_Type>::
Import(const DataToFunctionType &reference,bool adoptscaling)
{ 
  if (adoptscaling) {
    p_xaxis->SetScaling(reference.p_xaxis->Scaling()->Name());
    p_yaxis->SetScaling(reference.p_yaxis->Scaling()->Name());
  }
  Resize(reference.p_xydata->size());
  for (unsigned int i=0;i<reference.p_xydata->size();++i) {
    ArgumentType x=(*p_xaxis)((*reference.p_xaxis)
			      [(*reference.p_xydata)[i].first]);
    ResultType y=(*p_yaxis)((*reference.p_yaxis)
			    [(*reference.p_xydata)[i].second]);
    (*p_xydata)[i]=XYPair(x,y);
    (*p_yxdata)[i]=YXPair(y,x);
  }
  return Sort();
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
      (*_p_xdata)[i]=(*p_yaxis)[(*p_yxdata)[i].first];   
      (*_p_ydata)[i]=(*p_xaxis)[(*p_yxdata)[i].second];   
    }
  }
}

template <class Argument_Type,class Result_Type>
bool Data_To_Function<Argument_Type,Result_Type>::SortX()
{
  bool success=true;
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
	  ATOOLS::msg.Error()<<"Data_To_Function::SortX(): "
			     <<"Monotony violation! "<<std::endl
			     <<"   Removing point at x = "
			     <<(*p_xaxis)[xy.first]<<std::endl;
	  DeleteXPoint((*p_xaxis)[xy.first]);
	  success=false;
	}
	break;
      case Decreasing:
	if (((*p_xydata)[j+1].second-xy.second)>0) {
	  ATOOLS::msg.Error()<<"Data_To_Function::SortX(): "
			     <<"Monotony violation! "<<std::endl
			     <<"   Removing point at x = "
			     <<(*p_xaxis)[xy.first]<<std::endl;
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
  return success;
}
  
template <class Argument_Type,class Result_Type>
bool Data_To_Function<Argument_Type,Result_Type>::SortY()
{
  bool success=true;
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
	  ATOOLS::msg.Error()<<"Data_To_Function::SortY(): "
			     <<"Monotony violation! "<<std::endl
			     <<"   Removing point at y = "
			     <<(*p_yaxis)[yx.first]<<std::endl;
	  DeleteYPoint((*p_yaxis)[yx.first]);
	  success=false;
	}
	break;
      case Decreasing:
	if (((*p_yxdata)[j+1].second-yx.second)>0) {
	  ATOOLS::msg.Error()<<"Data_To_Function::SortY(): "
			     <<"Monotony violation! "<<std::endl
			     <<"   Removing point at y = "
			     <<(*p_yaxis)[yx.first]<<std::endl;
	    DeleteYPoint((*p_yaxis)[yx.first]);
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
  return success;
}
  
template <class Argument_Type,class Result_Type>
unsigned int Data_To_Function<Argument_Type,Result_Type>::
ClosestX(ArgumentType x,unsigned int &left,unsigned int &right)
{
  left=0; right=p_xydata->size()-1;
  unsigned int middle=(right-left)/2;
  if (x<(*p_xydata)[left].first) { 
    right=left+1;
    return left; 
  } 
  else if (x>(*p_xydata)[right].first) { 
    left=right-1;
    return right; 
  }
  do {
    if (x<(*p_xydata)[middle].first) { right=middle; middle=(middle+left)/2; }
    else { left=middle; middle=(right+middle)/2; } 
  } while (((right-middle)!=0)&&((middle-left)!=0));
  if ((*p_xydata)[right].first-x <= x-(*p_xydata)[left].first) return right;
  else return left; 
}
  
template <class Argument_Type,class Result_Type>
unsigned int Data_To_Function<Argument_Type,Result_Type>::
ClosestY(ResultType y,unsigned int &left,unsigned int &right)
{
  left=0; right=p_yxdata->size()-1;
  unsigned int middle=(right-left)/2;
  if (y<(*p_yxdata)[left].first) { 
    right=left+1; 
    return left; 
  } 
  else if (y>(*p_yxdata)[right].first) { 
    left=right-1; 
    return right; 
  }
  do {
    if (y<(*p_yxdata)[middle].first) { right=middle; middle=(middle+left)/2; }
    else { left=middle; middle=(right+middle)/2; } 
  } while (((right-middle)!=0)&&((middle-left)!=0));
  if ((*p_yxdata)[right].first-y <= y-(*p_yxdata)[left].first) return right; 
  else return left; 
}

template <class Argument_Type,class Result_Type>
bool Data_To_Function<Argument_Type,Result_Type>::
DeleteXPoint(ArgumentType _x)
{ 
  _x=(*p_xaxis)(_x);
  for (typename XYVector::iterator xyit=p_xydata->begin();
       xyit!=p_xydata->end();++xyit) 
    if (ATOOLS::IsEqual(xyit->first,_x)) { 
      ArgumentType x=xyit->first;
      p_xydata->erase(xyit);
      for (typename YXVector::iterator yxit=p_yxdata->begin();
	   yxit!=p_yxdata->end();++yxit) 
	if (yxit->second==x) { 
	  p_yxdata->erase(yxit); 
	  break; 
	}
      return true;
    }
  return false;
}

template <class Argument_Type,class Result_Type>
bool Data_To_Function<Argument_Type,Result_Type>::
DeleteYPoint(ResultType _y)
{ 
  _y=(*p_yaxis)(_y);
  for (typename YXVector::iterator yxit=p_yxdata->begin();
       yxit!=p_yxdata->end();++yxit) 
    if (ATOOLS::IsEqual(yxit->first,_y)) { 
      ResultType y=yxit->first;
      p_yxdata->erase(yxit);
      for (typename XYVector::iterator xyit=p_xydata->begin();
	   xyit!=p_xydata->end();++xyit) 
	if (xyit->second==y) { 
	  p_xydata->erase(xyit); 
	  break; 
	}
      return true;
    }
  return false;
}

template <class Argument_Type,class Result_Type>
void Data_To_Function<Argument_Type,Result_Type>::
RescaleX(const std::string scalename)
{ 
  Axis<ArgumentType> *xaxis = new Axis<ArgumentType>();
  xaxis->SetScaling(p_xaxis->Scaling()->Name());
  xaxis->SetScalingMode(p_xaxis->ScalingMode());
  p_xaxis->SetScaling(scalename);
  for (typename XYVector::iterator xyit=p_xydata->begin();
       xyit!=p_xydata->end();++xyit) {
    xyit->first=(*p_xaxis)((*xaxis)[xyit->first]);
  }
  for (typename YXVector::iterator yxit=p_yxdata->begin();
       yxit!=p_yxdata->end();++yxit) {
    yxit->second=(*p_xaxis)((*xaxis)[yxit->second]);
  }
}

template <class Argument_Type,class Result_Type>
void Data_To_Function<Argument_Type,Result_Type>::
RescaleY(const std::string scalename)
{ 
  Axis<ResultType> *yaxis = new Axis<ResultType>();
  yaxis->SetScaling(p_yaxis->Scaling()->Name());
  yaxis->SetScalingMode(p_yaxis->ScalingMode());
  p_yaxis->SetScaling(scalename);
  for (typename YXVector::iterator yxit=p_yxdata->begin();
       yxit!=p_yxdata->end();++yxit) {
    yxit->first=(*p_yaxis)((*yaxis)[yxit->first]);
  }
  for (typename XYVector::iterator xyit=p_xydata->begin();
       xyit!=p_xydata->end();++xyit) {
    xyit->second=(*p_yaxis)((*yaxis)[xyit->second]);
  }
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
  return 0;
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
  return 0;
}

template <class Argument_Type,class Result_Type>
Result_Type Data_To_Function<Argument_Type,Result_Type>::
Y(ArgumentType x,AcquisitionModeID tempmode)
{ 
  if (p_xydata->size()<2) {
    return 0.0;
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
  return 0;
}

template <class Argument_Type,class Result_Type>
Argument_Type Data_To_Function<Argument_Type,Result_Type>::
X(ResultType y,AcquisitionModeID tempmode)
{ 
  if (p_xydata->size()<2) {
    return 0.0;
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
  return 0;
}

template <class Argument_Type,class Result_Type>
Result_Type Data_To_Function<Argument_Type,Result_Type>::
IntegrateY(ArgumentType xmin,ArgumentType xmax)
{ 
  ResultType integral=0.0;
  if (p_xydata->size()<2) return integral;
  if (xmin==xmax) {
    xmin=XMin();
    xmax=XMax();
  }
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
  integral+=(*p_yaxis)[yleft]*
    ((*p_xaxis)[(*p_xydata)[start].first]-(*p_xaxis)[xmin]);
  for (unsigned int i=start;i<stop;++i) {
    integral+=(*p_yaxis)[(*p_xydata)[i].second]*
      ((*p_xaxis)[(*p_xydata)[i+1].first]-(*p_xaxis)[(*p_xydata)[i].first]);
  }
  integral+=(*p_yaxis)[(*p_xydata)[stop].second]*
    ((*p_xaxis)[xmax]-(*p_xaxis)[(*p_xydata)[stop].first]);
  return integral;
}

template <class Argument_Type,class Result_Type>
Argument_Type Data_To_Function<Argument_Type,Result_Type>::
IntegrateX(ResultType ymin,ResultType ymax)
{ 
  ArgumentType integral=0.0;
  if (p_yxdata->size()<2) return integral;
  if (ymin==ymax) {
    ymin=YMin();
    ymax=YMax();
  }
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
    ((*p_yaxis)[(*p_yxdata)[start].first]-(*p_yaxis)[ymin])/2.0;
  for (unsigned int i=start;i<stop;++i) {
    integral+=((*p_xaxis)[(*p_yxdata)[i+1].second]+
	       (*p_xaxis)[(*p_yxdata)[i].second])*
      ((*p_yaxis)[(*p_yxdata)[i+1].first]-
       (*p_yaxis)[(*p_yxdata)[i].first])/2.0;
  }
  integral+=((*p_xaxis)[xright]+(*p_xaxis)[(*p_yxdata)[stop].second])*
    ((*p_yaxis)[ymax]-(*p_yaxis)[(*p_yxdata)[stop].first])/2.0;
  return integral;
}

template <class Argument_Type,class Result_Type>
Data_To_Function<Argument_Type,Result_Type> *
Data_To_Function<Argument_Type,Result_Type>::
IntegralY(ArgumentType xmin,ArgumentType xmax,std::string xscaling,
	  std::string yscaling,bool forward,const MonotonyID monotony)
{ 
  DataToFunctionType *integrated = NULL;
  ResultType integral=0.0;
  if (p_xydata->size()<2) return integrated;
  if (xmin==xmax) {
    xmin=(*p_xaxis)[(*p_xydata)[0].first];
    xmax=(*p_xaxis)[(*p_xydata)[p_xydata->size()-1].first];
  }
  integrated = new DataToFunctionType();
  integrated->SetMonotony(monotony);
  if (xscaling=="") xscaling=p_xaxis->Scaling()->Name();
  if (yscaling=="") yscaling=p_yaxis->Scaling()->Name();
  integrated->XAxis()->SetScaling(xscaling);
  integrated->YAxis()->SetScaling(yscaling);
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
      ((*p_xaxis)[(*p_xydata)[start].first]-(*p_xaxis)[xmin])/2.0;
    for (unsigned int i=start;i<stop;++i) {
      integrated->AddPoint((*p_xaxis)[(*p_xydata)[i].first],integral);
      integral+=((*p_yaxis)[(*p_xydata)[i+1].second]+
		 (*p_yaxis)[(*p_xydata)[i].second])*
	((*p_xaxis)[(*p_xydata)[i+1].first]-
	 (*p_xaxis)[(*p_xydata)[i].first])/2.0;
    }
    integrated->AddPoint((*p_xaxis)[(*p_xydata)[stop].first],integral);
    integral+=((*p_yaxis)[yright]+(*p_yaxis)[(*p_xydata)[stop].second])*
      ((*p_xaxis)[xmax]-(*p_xaxis)[(*p_xydata)[stop].first])/2.0;
    integrated->AddPoint((*p_xaxis)[xmax],integral);
  }
  else {
    integrated->AddPoint((*p_xaxis)[xmax],integral);
    integral+=((*p_yaxis)[(*p_xydata)[stop].second]+(*p_yaxis)[yright])*
      ((*p_xaxis)[xmax]-(*p_xaxis)[(*p_xydata)[stop].first])/2.0;
    for (unsigned int i=stop;i>start;--i) {
      integrated->AddPoint((*p_xaxis)[(*p_xydata)[i].first],integral);
      integral+=((*p_yaxis)[(*p_xydata)[i].second]+
		 (*p_yaxis)[(*p_xydata)[i-1].second])*
	((*p_xaxis)[(*p_xydata)[i].first]-
	 (*p_xaxis)[(*p_xydata)[i-1].first])/2.0;
    }
    integrated->AddPoint((*p_xaxis)[(*p_xydata)[start].first],integral);
    integral+=((*p_yaxis)[yleft]+(*p_yaxis)[(*p_xydata)[start].second])*
      ((*p_xaxis)[(*p_xydata)[start].first]-(*p_xaxis)[xmin])/2.0;
    integrated->AddPoint((*p_xaxis)[xmin],integral);
  }
  return integrated;
}

template <class Argument_Type,class Result_Type>
Data_To_Function<Argument_Type,Result_Type> *
Data_To_Function<Argument_Type,Result_Type>::
IntegralX(ResultType ymin,ResultType ymax,std::string yscaling,
	  std::string xscaling,bool forward,const MonotonyID monotony)
{ 
  DataToFunctionType *integrated = NULL;
  ArgumentType integral=0.0;
  if (p_yxdata->size()<2) return integrated;
  if (ymin==ymax) {
    ymin=(*p_yaxis)[(*p_yxdata)[0].first];
    ymax=(*p_yaxis)[(*p_yxdata)[p_yxdata->size()-1].first];
  }
  integrated = new DataToFunctionType();
  integrated->SetMonotony(monotony);
  if (yscaling=="") yscaling=p_yaxis->Scaling()->Name();
  if (xscaling=="") xscaling=p_xaxis->Scaling()->Name();
  integrated->YAxis()->SetScaling(yscaling);
  integrated->XAxis()->SetScaling(xscaling);
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
      ((*p_yaxis)[(*p_yxdata)[start].first]-(*p_yaxis)[ymin])/2.0;
    for (unsigned int i=start;i<stop;++i) {
      integrated->AddPoint((*p_yaxis)[(*p_yxdata)[i].first],integral);
      integral+=((*p_xaxis)[(*p_yxdata)[i+1].second]+
		 (*p_xaxis)[(*p_yxdata)[i].second])*
	((*p_yaxis)[(*p_yxdata)[i+1].first]-
	 (*p_yaxis)[(*p_yxdata)[i].first])/2.0;
    }
    integrated->AddPoint(integral,(*p_yaxis)[(*p_yxdata)[stop].first]);
    integral+=((*p_xaxis)[xright]+(*p_xaxis)[(*p_yxdata)[stop].second])*
      ((*p_yaxis)[ymax]-(*p_yaxis)[(*p_yxdata)[stop].first])/2.0;
    integrated->AddPoint(integral,(*p_yaxis)[ymax]);
  }
  else {
    integrated->AddPoint((*p_yaxis)[ymax],integral);
    integral+=((*p_xaxis)[(*p_yxdata)[stop].second]+(*p_xaxis)[xright])*
	((*p_yaxis)[ymax]-(*p_yaxis)[(*p_yxdata)[stop].first])/2.0;
    for (unsigned int i=stop;i>start;--i) {
      integrated->AddPoint((*p_yaxis)[(*p_yxdata)[i].first],integral);
      integral+=((*p_xaxis)[(*p_yxdata)[i].second]+
		 (*p_xaxis)[(*p_yxdata)[i-1].second])*
	((*p_yaxis)[(*p_yxdata)[i].first]-
	 (*p_yaxis)[(*p_yxdata)[i-1].first])/2.0;
    }
    integrated->AddPoint(integral,(*p_yaxis)[(*p_yxdata)[start].first]);
    integral+=((*p_xaxis)[xleft]+(*p_xaxis)[(*p_yxdata)[start].second])*
	((*p_yaxis)[(*p_yxdata)[start].first]-(*p_yaxis)[ymin])/2.0;
    integrated->AddPoint(integral,(*p_yaxis)[ymin]);
  }
  return integrated;
}

template <class Argument_Type,class Result_Type>
Result_Type Data_To_Function<Argument_Type,Result_Type>::
DifferentialY(ArgumentType x)
{ 
  unsigned int left, right;
  ClosestX(x,left,right);
  return ((*p_yaxis)[(*p_xydata)[right].second]-
	  (*p_yaxis)[(*p_xydata)[left].second])/
    ((*p_xaxis)[(*p_xydata)[right].first]-
     (*p_xaxis)[(*p_xydata)[left].first]);	
}

template <class Argument_Type,class Result_Type>
Argument_Type Data_To_Function<Argument_Type,Result_Type>::
DifferentialX(ResultType y)
{ 
  unsigned int left, right;
  ClosestY(y,left,right);
  return ((*p_xaxis)[(*p_yxdata)[right].second]-
	  (*p_xaxis)[(*p_yxdata)[left].second])/
    ((*p_yaxis)[(*p_yxdata)[right].first]-
     (*p_yaxis)[(*p_yxdata)[left].first]);	
}
  
template <class Argument_Type,class Result_Type>
Data_To_Function<Argument_Type,Result_Type> *
Data_To_Function<Argument_Type,Result_Type>::
DerivativeY(ArgumentType xmin,ArgumentType xmax,
	    std::string xscaling,std::string yscaling)
{ 
  DataToFunctionType *derivated = NULL;
  ResultType differential;
  if (p_xydata->size()<2) return derivated;
  if (xmin==xmax) {
    xmin=(*p_xaxis)[(*p_xydata)[0].first];
    xmax=(*p_xaxis)[(*p_xydata)[p_xydata->size()-1].first];
  }
  derivated = new DataToFunctionType();
  if (xscaling=="") xscaling=p_xaxis->Scaling()->Name();
  if (yscaling=="") yscaling=p_yaxis->Scaling()->Name();
  derivated->XAxis()->SetScaling(xscaling);
  derivated->YAxis()->SetScaling(yscaling);
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
    differential=((*p_yaxis)[(*p_xydata)[start].second]-(*p_yaxis)[yleft])/
      ((*p_xaxis)[(*p_xydata)[start].first]-(*p_xaxis)[xmin]);
    derivated->AddPoint((*p_xaxis)[xmin],differential);
  }
  for (unsigned int i=start;i<stop;++i) {
    differential=((*p_yaxis)[(*p_xydata)[i+1].second]-
		  (*p_yaxis)[(*p_xydata)[i].second])/
      ((*p_xaxis)[(*p_xydata)[i+1].first]-(*p_xaxis)[(*p_xydata)[i].first]);
    derivated->AddPoint((*p_xaxis)[(*p_xydata)[i].first],differential);
  }
  if ((*p_xydata)[stop].first!=xmax) {
    differential=((*p_yaxis)[yright]-(*p_yaxis)[(*p_xydata)[stop].second])/
      ((*p_xaxis)[xmax]-(*p_xaxis)[(*p_xydata)[stop].first]);
    derivated->AddPoint((*p_xaxis)[(*p_xydata)[stop].first],differential);
  }
  if (next<p_xydata->size()) {
    if (xmax!=(*p_xydata)[next].first) {
      differential=((*p_yaxis)[(*p_xydata)[next].second]-(*p_yaxis)[yright])/
	((*p_xaxis)[(*p_xydata)[next].first]-(*p_xaxis)[xmax]);
      derivated->AddPoint((*p_xaxis)[xmax],differential);
    }
  }
  return derivated;
}

template <class Argument_Type,class Result_Type>
Data_To_Function<Argument_Type,Result_Type> *
Data_To_Function<Argument_Type,Result_Type>::
DerivativeX(ResultType ymin,ResultType ymax,
	    std::string yscaling,std::string xscaling)
{ 
  DataToFunctionType *derivated = NULL;
  ArgumentType differential;
  if (p_yxdata->size()<2) return derivated;
  if (ymin==ymax) {
    ymin=(*p_yaxis)[(*p_yxdata)[0].first];
    ymax=(*p_yaxis)[(*p_yxdata)[p_yxdata->size()-1].first];
  }
  derivated = new DataToFunctionType();
  if (yscaling=="") yscaling=p_yaxis->Scaling()->Name();
  if (xscaling=="") xscaling=p_xaxis->Scaling()->Name();
  derivated->YAxis()->SetScaling(yscaling);
  derivated->XAxis()->SetScaling(xscaling);
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
    differential=((*p_xaxis)[(*p_yxdata)[start].second]-(*p_xaxis)[xleft])/
      ((*p_yaxis)[(*p_yxdata)[start].first]-(*p_yaxis)[ymin]);
    derivated->AddPoint(differential,(*p_yaxis)[ymin]);
  }
  for (unsigned int i=start;i<stop;++i) {
    differential=((*p_xaxis)[(*p_yxdata)[i+1].second]-
		  (*p_xaxis)[(*p_yxdata)[i].second])/
      ((*p_yaxis)[(*p_yxdata)[i+1].first]-(*p_yaxis)[(*p_yxdata)[i].first]);
    derivated->AddPoint(differential,(*p_yaxis)[(*p_yxdata)[i].first]);
  }
  if ((*p_yxdata)[stop].first!=ymax) {
    differential=((*p_xaxis)[xright]-(*p_xaxis)[(*p_yxdata)[stop].second])/
      ((*p_yaxis)[ymax]-(*p_yaxis)[(*p_yxdata)[stop].first]);
    derivated->AddPoint(differential,(*p_yaxis)[(*p_yxdata)[stop].first]);
  }
  if (next<p_yxdata->size()) {
    if (ymax!=(*p_yxdata)[next].first) {
      differential=((*p_xaxis)[(*p_yxdata)[next].second]-(*p_xaxis)[xright])/
	((*p_yaxis)[(*p_yxdata)[next].first]-(*p_yaxis)[ymax]);
      derivated->AddPoint(differential,(*p_yaxis)[ymax]);
    }
  }
  return derivated;
}

template <class Argument_Type,class Result_Type>
void Data_To_Function<Argument_Type,Result_Type>::
ScaleY(ResultType scalefactor)
{ 
  for (typename YXVector::iterator yxit=p_yxdata->begin();
       yxit!=p_yxdata->end();++yxit) 
    yxit->first=(*p_yaxis)((*p_yaxis)[yxit->first]*scalefactor);
  for (typename XYVector::iterator xyit=p_xydata->begin();
       xyit!=p_xydata->end();++xyit) 
    xyit->second=(*p_yaxis)((*p_yaxis)[xyit->second]*scalefactor);
}

template <class Argument_Type,class Result_Type>
void Data_To_Function<Argument_Type,Result_Type>::
ScaleX(ArgumentType scalefactor)
{ 
  for (typename XYVector::iterator xyit=p_xydata->begin();
       xyit!=p_xydata->end();++xyit) 
    xyit->first=(*p_xaxis)((*p_xaxis)[xyit->first]*scalefactor);
  for (typename YXVector::iterator yxit=p_yxdata->begin();
       yxit!=p_yxdata->end();++yxit) 
    yxit->second=(*p_xaxis)((*p_xaxis)[yxit->second]*scalefactor);
}

template <class Argument_Type,class Result_Type>
void Data_To_Function<Argument_Type,Result_Type>::
MoveY(ResultType distance)
{ 
  for (typename YXVector::iterator yxit=p_yxdata->begin();
       yxit!=p_yxdata->end();++yxit) 
    yxit->first=(*p_yaxis)((*p_yaxis)[yxit->first]+distance);
  for (typename XYVector::iterator xyit=p_xydata->begin();
       xyit!=p_xydata->end();++xyit) 
    xyit->second=(*p_yaxis)((*p_yaxis)[xyit->second]+distance);
}

template <class Argument_Type,class Result_Type>
void Data_To_Function<Argument_Type,Result_Type>::
MoveX(ArgumentType distance)
{ 
  for (typename XYVector::iterator xyit=p_xydata->begin();
       xyit!=p_xydata->end();++xyit) 
    xyit->first=(*p_xaxis)((*p_xaxis)[xyit->first]+distance);
  for (typename YXVector::iterator yxit=p_yxdata->begin();
       yxit!=p_yxdata->end();++yxit) 
    yxit->second=(*p_xaxis)((*p_xaxis)[yxit->second]+distance);
}

template <class Argument_Type,class Result_Type>
Argument_Type  Data_To_Function<Argument_Type,Result_Type>::
DeltaXMin(ResultType& left,ResultType& right)
{ 
  ArgumentType minimum=0.0, cur;
  if (p_xydata->size()>2) {
    typename XYVector::iterator xyit=p_xydata->begin()+1;
    minimum=dabs((*p_xaxis)[xyit->first]-(*p_xaxis)[(xyit-1)->first]);
    left=(*p_yaxis)[(xyit-1)->second];
    right=(*p_yaxis)[xyit->second];
    for (++xyit;xyit!=p_xydata->end();++xyit) {
      cur=dabs((*p_xaxis)[xyit->first]-(*p_xaxis)[(xyit-1)->first]);
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
  ArgumentType maximum=0.0, cur;
  if (p_xydata->size()>2) {
    for (typename XYVector::iterator xyit=p_xydata->begin()+1;
	 xyit!=p_xydata->end();++xyit) {
      cur=dabs((*p_xaxis)[xyit->first]-
	       (*p_xaxis)[(xyit-1)->first]);
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
  ResultType minimum=0.0, cur;
  if (p_yxdata->size()>2) {
    typename YXVector::iterator yxit=p_yxdata->begin()+1;
    minimum=dabs((*p_yaxis)[yxit->first]-
		 (*p_yaxis)[(yxit-1)->first]);
    left=(*p_xaxis)[(yxit-1)->second];
    right=(*p_xaxis)[yxit->second];
    for (++yxit;yxit!=p_yxdata->end();++yxit) {
      cur=dabs((*p_yaxis)[yxit->first]-
	       (*p_yaxis)[(yxit-1)->first]);
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
  ResultType maximum=0.0, cur;
  if (p_yxdata->size()>2) {
    for (typename YXVector::iterator yxit=p_yxdata->begin()+1;
	 yxit!=p_yxdata->end();++yxit) {
      cur=dabs((*p_yaxis)[yxit->first]-
	       (*p_yaxis)[(yxit-1)->first]);
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

template <class Argument_Type,class Result_Type>
bool Data_To_Function<Argument_Type,Result_Type>::
AddPoint(Argument_Type x,Result_Type y)
{
  ArgumentType _x=(*p_xaxis)(x);
  ResultType _y=(*p_yaxis)(y);
  if (p_xydata->size()!=0) {
    if (((*p_xydata)[ClosestX(_x)].first==_x)||
	((*p_yxdata)[ClosestY(_y)].second==_y)) return false;
  }
  p_xydata->push_back(XYPair(_x,_y)); 
  p_yxdata->push_back(YXPair(_y,_x));
  return Sort();
}

template <class Argument_Type,class Result_Type>
Result_Type Data_To_Function<Argument_Type,Result_Type>::
LinearY(ArgumentType x,unsigned int &left,unsigned int &right)
{ 
  ArgumentType xleft=(*p_xydata)[left].first, xright=(*p_xydata)[right].first;
  ResultType yleft=(*p_xydata)[left].second, yright=(*p_xydata)[right].second;
  return yleft+(x-xleft)*((yright-yleft)/(xright-xleft)); 
}

template <class Argument_Type,class Result_Type>
Argument_Type Data_To_Function<Argument_Type,Result_Type>::
LinearX(ResultType y,unsigned int &left,unsigned int &right)
{ 
  ResultType yleft=(*p_yxdata)[left].first, yright=(*p_yxdata)[right].first;
  ArgumentType xleft=(*p_yxdata)[left].second, xright=(*p_yxdata)[right].second;
  return xleft+(y-yleft)*((xright-xleft)/(yright-yleft)); 
}

template <class Argument_Type,class Result_Type>
unsigned int Data_To_Function<Argument_Type,Result_Type>::
ClosestX(ArgumentType x)
{
  unsigned int dummyleft, dummyright;
  return ClosestX(x,dummyleft,dummyright);
}

template <class Argument_Type,class Result_Type>
unsigned int Data_To_Function<Argument_Type,Result_Type>::
ClosestY(ResultType y)
{
  unsigned int dummyleft, dummyright;
  return ClosestY(y,dummyleft,dummyright);
}

template <class Argument_Type,class Result_Type>
Result_Type Data_To_Function<Argument_Type,Result_Type>::
DataY(ArgumentType x,AcquisitionModeID tempmode)
{ 
  if (tempmode!=Data) {
    unsigned int left, right;
    ClosestX(x,left,right);
    if (tempmode==UpperData) return (*p_xydata)[right].second;
    if (tempmode==LowerData) return (*p_xydata)[left].second;
  }
  return (*p_xydata)[ClosestX(x)].second; 
}

template <class Argument_Type,class Result_Type>
Argument_Type Data_To_Function<Argument_Type,Result_Type>::
DataX(ResultType y,AcquisitionModeID tempmode)
{ 
  if (tempmode!=Data) {
    unsigned int left, right;
    ClosestY(y,left,right);
    if (tempmode==UpperData) return (*p_yxdata)[right].second;
    if (tempmode==LowerData) return (*p_yxdata)[left].second;
  }
  return (*p_yxdata)[ClosestY(y)].second; 
}

template <class Argument_Type,class Result_Type>
void Data_To_Function<Argument_Type,Result_Type>::
NormalizeY(ArgumentType xmin,ArgumentType xmax,ResultType norm)
{ 
  if (p_xydata->size()<2) return;
  if (xmin==xmax) ScaleY(norm/IntegrateY()); 
  ScaleY(norm/IntegrateY(xmin,xmax)); 
}

template <class Argument_Type,class Result_Type>
void Data_To_Function<Argument_Type,Result_Type>::
NormalizeX(ResultType ymin,ResultType ymax,ArgumentType norm)
{ 
  if (p_yxdata->size()<2) return;
  if (ymin==ymax) ScaleX(norm/IntegrateX()); 
  ScaleX(norm/IntegrateX(ymin,ymax)); 
}

template <class Argument_Type,class Result_Type>
int Data_To_Function<Argument_Type,Result_Type>::XPosition(ArgumentType x)
{ 
  int position=(int)ClosestX((*p_xaxis)(x)); 
  if ((*p_xydata)[position].first==(*p_xaxis)(x)) return position; 
  else return -1; 
}

template <class Argument_Type,class Result_Type>
int Data_To_Function<Argument_Type,Result_Type>::YPosition(ResultType y)
{ 
  int position=(int)ClosestY((*p_yaxis)(y)); 
  if ((*p_yxdata)[position].first==(*p_yaxis)(y)) return position; 
  else return -1; 
}

template <class Argument_Type,class Result_Type>
bool Data_To_Function<Argument_Type,Result_Type>::
ReplaceXPoint(ArgumentType x,ResultType y)
{ 
  bool deletepoint=DeleteXPoint(x);
  if (deletepoint) return AddPoint(x,y); 
  return false; 
}

template <class Argument_Type,class Result_Type>
bool Data_To_Function<Argument_Type,Result_Type>::
ReplaceYPoint(ArgumentType x,ResultType y)
{ 
  bool deletepoint=DeleteYPoint(y);
  if (deletepoint) return AddPoint(x,y);
  return false;
}

template <class Argument_Type,class Result_Type>
void Data_To_Function<Argument_Type,Result_Type>::
Resize(unsigned int newsize)
{ 
  p_xydata->resize(newsize); 
  p_yxdata->resize(newsize); 
}

template class Data_To_Function<double,double>;
