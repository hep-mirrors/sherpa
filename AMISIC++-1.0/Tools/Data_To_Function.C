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
  Data_To_Function(const Data_To_Function<ArgumentType,ResultType> &reference)
  { Init(); Import(reference.p_xydata); }

  template <class Argument_Type,class Result_Type>
  Data_To_Function<Argument_Type,Result_Type>::
  Data_To_Function(ArgumentVector *_p_xdata,ResultVector *_p_ydata)
  { Init(); Import(_p_xdata,_p_ydata); }
  
  template <class Argument_Type,class Result_Type>
  Data_To_Function<Argument_Type,Result_Type>::
  Data_To_Function(XYDataMap *_p_xydata)
  { Init(); Import(_p_xydata); }

  template <class Argument_Type,class Result_Type>
  Data_To_Function<Argument_Type,Result_Type>::~Data_To_Function()
  {
    p_xydata->clear(); delete p_xydata;
    p_yxdata->clear(); delete p_yxdata;
    p_xdata->clear(); delete p_xdata;
    p_ydata->clear(); delete p_ydata; 
    delete p_xaxis;
    delete p_yaxis;
  }

  template <class Argument_Type,class Result_Type>
  void Data_To_Function<Argument_Type,Result_Type>::Init()
  {
    p_xaxis = new ATOOLS::Axis<ArgumentType>();
    p_yaxis = new ATOOLS::Axis<ResultType>();
    p_xdata = new ArgumentVector();
    p_ydata = new ResultVector();
    p_xydata = new XYDataMap();
    p_yxdata = new YXDataMap();
    m_acquisitionmode = Interpolation;
    m_interpolationmode = Linear;
  }
  
  template <class Argument_Type,class Result_Type>
  void Data_To_Function<Argument_Type,Result_Type>::
  Import(ArgumentVector *_p_xdata,ResultVector *_p_ydata)
  { 
    if (_p_xdata->size()==_p_ydata->size()) {
      Resize(_p_xdata->size());
      ArgumentType x; ResultType y;
      for (unsigned int i=0; i<_p_xdata->size(); ++i) {
	x=(*p_xaxis)((*_p_xdata)[i]); 
	y=(*p_yaxis)((*_p_ydata)[i]);
	(*p_xdata)[i]=x; 
	(*p_ydata)[i]=y;
	(*p_xydata)[x]=y; 
	(*p_yxdata)[y]=x;
      }
    }
#ifdef DEBUG__Data_To_Function
    std::cout<<"Data_To_Function::Import("<<_p_xdata<<","<<_p_ydata<<") :"<<std::endl;
#endif
    Sort();
  }
  
  template <class Argument_Type,class Result_Type>
  void Data_To_Function<Argument_Type,Result_Type>::
  Import(XYDataMap *_p_xydata)
  { 
    unsigned int i=-1;
    Resize(_p_xydata->size());
    for (XYDataIterator xyit=_p_xydata->begin();xyit!=_p_xydata->end();++xyit) {
      ++i;
      (*p_xdata)[i]=(*p_xaxis)((*xyit).first);
      (*p_ydata)[i]=(*p_yaxis)((*xyit).second);
      (*p_xydata)[(*p_xaxis)((*xyit).first)]=(*p_yaxis)((*xyit).second);
      (*p_yxdata)[(*p_yaxis)((*xyit).second)]=(*p_xaxis)((*xyit).first);
    }
#ifdef DEBUG__Data_To_Function
    std::cout<<"Data_To_Function::Import("<<_p_xydata<<") :"<<std::endl;
#endif
    Sort();
  }
  
  template <class Argument_Type,class Result_Type>
  void Data_To_Function<Argument_Type,Result_Type>::
  Export(ArgumentVector *_p_xdata,ResultVector *_p_ydata)
  { 
    _p_xdata->clear();
    _p_ydata->clear();
    _p_xdata->resize(p_xdata->size());
    _p_ydata->resize(p_ydata->size());
    for (unsigned int i=0;i<p_xdata->size();++i) {
      (*_p_xdata)[i]=(*p_xaxis)[(*p_xdata)[i]];   
      (*_p_ydata)[i]=(*p_yaxis)[(*p_xydata)[(*p_xdata)[i]]];   
    }
  }
  
  template <class Argument_Type,class Result_Type>
  void Data_To_Function<Argument_Type,Result_Type>::
  Export(XYDataMap *_p_xydata)
  { 
    _p_xydata->clear();
    for (unsigned int i=0;i<p_xydata->size();++i) {
      (*_p_xydata)[(*p_xaxis)[(*p_xdata)[i]]]=(*p_yaxis)[(*p_xydata)[(*p_xdata)[i]]];   
    }
  }
  
  template <class Argument_Type,class Result_Type>
  void Data_To_Function<Argument_Type,Result_Type>::SortX()
  {
    ArgumentType x;
    bool cont;
#ifdef DEBUG__Data_To_Function
    std::cout<<"Data_To_Function::SortX() :"<<std::endl;
    std::cout<<"   before sorting: p_xdata = [ "; 
    for(ArgumentIterator xit=p_xdata->begin();xit!=p_xdata->end();std::cout<<*(xit++)<<" ");
    std::cout<<"]"<<std::endl;
#endif
    do {
      cont=false;
      for (unsigned int i=1; i<p_xdata->size(); ++i) {
	if ((*p_xdata)[i]<(*p_xdata)[i-1]) {
	  x=(*p_xdata)[i]; (*p_xdata)[i]=(*p_xdata)[i-1]; (*p_xdata)[i-1]=x;
	  cont=true;
	}
      }
    } while (cont); 
#ifdef DEBUG__Data_To_Function
    std::cout<<"   after sorting : p_xdata = [ "; 
    for(ArgumentIterator xit=p_xdata->begin();xit!=p_xdata->end();std::cout<<*(xit++)<<" ");
    std::cout<<"]"<<std::endl;
#endif
  }
  
  template <class Argument_Type,class Result_Type>
  void Data_To_Function<Argument_Type,Result_Type>::SortY()
  {
    ResultType y;
    bool cont;
#ifdef DEBUG__Data_To_Function
    std::cout<<"Data_To_Function::SortY() :"<<std::endl;
    std::cout<<"   before sorting: p_ydata = [ "; 
    for(ResultIterator yit=p_ydata->begin();yit!=p_ydata->end();std::cout<<*(yit++)<<" ");
    std::cout<<"]"<<std::endl;
#endif
    do {
      cont=false;
      for (unsigned int i=1; i<p_ydata->size(); ++i) {
	if ((*p_ydata)[i]<(*p_ydata)[i-1]) {
	  y=(*p_ydata)[i]; (*p_ydata)[i]=(*p_ydata)[i-1]; (*p_ydata)[i-1]=y;
	  cont=true;
	}
      }
    } while (cont); 
#ifdef DEBUG__Data_To_Function
    std::cout<<"   after sorting : p_ydata = [ "; 
    for(ResultIterator yit=p_ydata->begin();yit!=p_ydata->end();std::cout<<*(yit++)<<" ");
    std::cout<<"]"<<std::endl;
#endif
  }
  
  template <class Argument_Type,class Result_Type>
  unsigned int Data_To_Function<Argument_Type,Result_Type>::
  ClosestX(ArgumentType x,unsigned int &left,unsigned int &right)
  {
#ifdef DEBUG__Data_To_Function
    std::cout<<"Data_To_Function::GetClosestX("<<x<<") :"<<std::endl;
    std::cout<<"   p_xdata = [ "; 
    for(ArgumentIterator xit=p_xdata->begin();xit!=p_xdata->end();std::cout<<*(xit++)<<" ");
    std::cout<<"]"<<std::endl;
#endif
    left=0; right=p_xdata->size()-1;
    unsigned int middle=(right-left)/2;
    if (x<(*p_xdata)[left]) { 
#ifdef DEBUG__Data_To_Function
      std::cout<<"   value is out of range"<<std::endl;
      std::cout<<"   returning "<<left<<" => x value "<<(*p_xdata)[left]<<std::endl;
#endif
      right=left+1;
      return left; 
    } 
    else if (x>(*p_xdata)[right]) { 
#ifdef DEBUG__Data_To_Function
      std::cout<<"   value is out of range"<<std::endl;
      std::cout<<"   returning "<<right<<" => x value "<<(*p_xdata)[right]<<std::endl;
#endif
      left=right-1;
      return right; 
    }
    do {
      if (x<(*p_xdata)[middle]) { right=middle; middle=(middle+left)/2; }
      else { left=middle; middle=(right+middle)/2; } 
    } while (((right-middle)!=0)&&((middle-left)!=0));
    if ((*p_xdata)[right]-x <= x-(*p_xdata)[left]) { 
#ifdef DEBUG__Data_To_Function
      std::cout<<"   returning "<<right<<" / ["<<left<<","<<right<<"] => x values "
	       <<(*p_xdata)[right]<<" ["<<(*p_xdata)[left]<<","<<(*p_xdata)[right]<<"]"<<std::endl;
#endif
      return right;
    }
    else { 
#ifdef DEBUG__Data_To_Function
      std::cout<<"   returning "<<left<<" / ["<<left<<","<<right<<"] => x values "
	       <<(*p_xdata)[left]<<" ["<<(*p_xdata)[left]<<","<<(*p_xdata)[right]<<"]"<<std::endl;
#endif
      return left; 
    }
  }
  
  template <class Argument_Type,class Result_Type>
  unsigned int Data_To_Function<Argument_Type,Result_Type>::
  ClosestY(ResultType y,unsigned int &left,unsigned int &right)
  {
#ifdef DEBUG__Data_To_Function
    std::cout<<"Data_To_Function::GetClosestY("<<y<<") :"<<std::endl;
    std::cout<<"   p_ydata = [ "; 
    for(ResultIterator yit=p_ydata->begin();yit!=p_ydata->end();std::cout<<*(yit++)<<" ");
    std::cout<<"]"<<std::endl;
#endif
    left=0; right=p_ydata->size()-1;
    unsigned int middle=(right-left)/2;
    if (y<(*p_ydata)[left]) { 
#ifdef DEBUG__Data_To_Function
      std::cout<<"   value is out of range"<<std::endl;
      std::cout<<"   returning "<<left<<" => y value "<<(*p_ydata)[left]<<std::endl;
#endif
      right=left+1; 
      return left; 
    } 
    else if (y>(*p_ydata)[right]) { 
#ifdef DEBUG__Data_To_Function
      std::cout<<"   value is out of range"<<std::endl;
      std::cout<<"   returning "<<right<<" => y value "<<(*p_ydata)[right]<<std::endl;
#endif
      left=right-1; 
      return right; 
    }
    do {
      if (y<(*p_ydata)[middle]) { right=middle; middle=(middle+left)/2; }
      else { left=middle; middle=(right+middle)/2; } 
    } while (((right-middle)!=0)&&((middle-left)!=0));
    if ((*p_ydata)[right]-y <= y-(*p_ydata)[left]) { 
#ifdef DEBUG__Data_To_Function
      std::cout<<"   returning "<<right<<" / ["<<left<<","<<right<<"] => y values "
	       <<(*p_ydata)[right]<<" ["<<(*p_ydata)[left]<<","<<(*p_ydata)[right]<<"]"<<std::endl;
#endif
      return right; 
    }
    else { 
#ifdef DEBUG__Data_To_Function
      std::cout<<"   returning "<<left<<" / ["<<left<<","<<right<<"] => y values "
	       <<(*p_ydata)[left]<<" ["<<(*p_ydata)[left]<<","<<(*p_ydata)[right]<<"]"<<std::endl;
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
    for(ArgumentIterator xit=p_xdata->begin();xit!=p_xdata->end();std::cout<<*(xit++)<<" ");
    std::cout<<"]"<<std::endl;
#endif
    _x=(*p_xaxis)(_x);
    for (ArgumentIterator xit=p_xdata->begin();xit!=p_xdata->end();++xit) 
      if (ATOOLS::IsZero(*xit-_x)) { 
	ArgumentType x=*xit;
	p_xdata->erase(xit);
	ResultType y=(*p_xydata)[x];
	for (ResultIterator yit=p_ydata->begin();yit!=p_ydata->end();++yit) 
	  if (*yit==y) { 
	    p_ydata->erase(yit); 
	    break; 
	  }
	for (XYDataIterator xyit=p_xydata->begin();xyit!=p_xydata->end();++xyit) 
	  if ((*xyit).first==x) { 
	    p_xydata->erase(xyit);
	    break; 
	  }
	for (YXDataIterator yxit=p_yxdata->begin();yxit!=p_yxdata->end();++yxit) 
	  if ((*yxit).first==x) { 
	    p_yxdata->erase(yxit);
	    break; 
	  }
#ifdef DEBUG__Data_To_Function
	std::cout<<"   deleted point "<<x<<std::endl;
	std::cout<<"   after deletion : p_xdata = [ "; 
	for(ArgumentIterator xit=p_xdata->begin();xit!=p_xdata->end();std::cout<<*(xit++)<<" ");
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
    for(ResultIterator yit=p_ydata->begin();yit!=p_ydata->end();std::cout<<*(yit++)<<" ");
    std::cout<<"]"<<std::endl;
#endif
    _y=(*p_yaxis)(_y);
    for (ResultIterator yit=p_ydata->begin();yit!=p_ydata->end();++yit) 
      if (ATOOLS::IsZero(*yit-_y)) { 
	ResultType y=*yit;
	p_ydata->erase(yit);
	ArgumentType x=(*p_yxdata)[y];
	for (ArgumentIterator xit=p_xdata->begin();xit!=p_xdata->end();++xit) 
	  if (*xit==x) { 
	    p_xdata->erase(xit); 
	    break; 
	  }
	for (YXDataIterator yxit=p_yxdata->begin();yxit!=p_yxdata->end();++yxit) 
	  if ((*yxit).first==x) { 
	    p_yxdata->erase(yxit);
	    break; 
	  }
	for (XYDataIterator xyit=p_xydata->begin();xyit!=p_xydata->end();++xyit) 
	  if ((*xyit).first==x) { 
	    p_xydata->erase(xyit);
	    break; 
	  }
#ifdef DEBUG__Data_To_Function
	std::cout<<"   deleted point "<<y<<std::endl;
	std::cout<<"   after deletion : p_ydata = [ "; 
	for(ResultIterator yit=p_ydata->begin();yit!=p_ydata->end();std::cout<<*(yit++)<<" ");
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
  Y(ArgumentType x,AcquisitionMode tempmode)
  { 
    if (p_xdata->size()<2) {
      ATOOLS::msg.Error()<<"Data_To_Function::Y("<<x<<","<<tempmode<<"): "
			 <<"Error! Less than 2 data points available."<<std::endl
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
  X(ResultType y,AcquisitionMode tempmode)
  { 
    if (p_xdata->size()<2) {
      ATOOLS::msg.Error()<<"Data_To_Function::X("<<y<<","<<tempmode<<"): "
			 <<"Error! Less than 2 data points available."<<std::endl
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
    if (p_xdata->size()<2) return integral;
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
    integral+=((*p_yaxis)[(*p_xydata)[(*p_xdata)[start]]]+(*p_yaxis)[yleft])*
      ((*p_xaxis)[(*p_xdata)[start]]-(*p_xaxis)[xmin])/(ResultType)2.0;
#ifdef DEBUG__Data_To_Function
    std::cout<<"   integral values are [ "<<integral<<" ";
#endif
    for (unsigned int i=start;i<stop;++i) {
      integral+=((*p_yaxis)[(*p_xydata)[(*p_xdata)[i+1]]]+(*p_yaxis)[(*p_xydata)[(*p_xdata)[i]]])*
	((*p_xaxis)[(*p_xdata)[i+1]]-(*p_xaxis)[(*p_xdata)[i]])/(ResultType)2.0;
#ifdef DEBUG__Data_To_Function
      std::cout<<integral<<" ";
#endif
    }
    integral+=((*p_yaxis)[yright]+(*p_yaxis)[(*p_xydata)[(*p_xdata)[stop]]])*
      ((*p_xaxis)[xmax]-(*p_xaxis)[(*p_xdata)[stop]])/(ResultType)2.0;
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
    if (p_ydata->size()<2) return integral;
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
    integral+=((*p_xaxis)[(*p_yxdata)[(*p_ydata)[start]]]+(*p_xaxis)[xleft])*
      ((*p_yaxis)[(*p_ydata)[start]]-(*p_yaxis)[ymin])/(ArgumentType)2.0;
#ifdef DEBUG__Data_To_Function
    std::cout<<"   integral values are [ "<<integral<<" ";
#endif
    for (unsigned int i=start;i<stop;++i) {
      integral+=((*p_xaxis)[(*p_yxdata)[(*p_ydata)[i+1]]]+(*p_xaxis)[(*p_yxdata)[(*p_ydata)[i]]])*
	((*p_yaxis)[(*p_ydata)[i+1]]-(*p_yaxis)[(*p_ydata)[i]])/(ArgumentType)2.0;
#ifdef DEBUG__Data_To_Function
      std::cout<<integral<<" ";
#endif
    }
    integral+=((*p_xaxis)[xright]+(*p_xaxis)[(*p_yxdata)[(*p_ydata)[stop]]])*
      ((*p_yaxis)[ymax]-(*p_yaxis)[(*p_ydata)[stop]])/(ArgumentType)2.0;
#ifdef DEBUG__Data_To_Function
    std::cout<<integral<<" ]"<<std::endl;
#endif
    return integral;
  }

  template <class Argument_Type,class Result_Type>
  void Data_To_Function<Argument_Type,Result_Type>::
  ScaleY(ResultType scalefactor)
  { 
    ArgumentVector *_p_xdata = new ArgumentVector();
    ResultVector *_p_ydata = new ResultVector();
    Export(_p_xdata,_p_ydata);
    for (ResultIterator yit=_p_ydata->begin();yit!=_p_ydata->end();++yit) (*yit)*=scalefactor;
    Import(_p_xdata,_p_ydata);
    delete _p_xdata;
    delete _p_ydata;
  }

  template <class Argument_Type,class Result_Type>
  void Data_To_Function<Argument_Type,Result_Type>::
  ScaleX(ArgumentType scalefactor)
  { 
    ArgumentVector *_p_xdata = new ArgumentVector();
    ResultVector *_p_ydata = new ResultVector();
    Export(_p_xdata,_p_ydata);
    for (ArgumentIterator xit=_p_xdata->begin();xit!=_p_xdata->end();++xit) (*xit)*=scalefactor;
    Import(_p_xdata,_p_ydata);
    delete _p_xdata;
    delete _p_ydata;
  }

  template <class Argument_Type,class Result_Type>
  void Data_To_Function<Argument_Type,Result_Type>::
  MoveY(ResultType distance)
  { 
    ArgumentVector *_p_xdata = new ArgumentVector();
    ResultVector *_p_ydata = new ResultVector();
    Export(_p_xdata,_p_ydata);
    for (ResultIterator yit=_p_ydata->begin();yit!=_p_ydata->end();++yit) (*yit)+=distance;
    Import(_p_xdata,_p_ydata);
    delete _p_xdata;
    delete _p_ydata;
  }

  template <class Argument_Type,class Result_Type>
  void Data_To_Function<Argument_Type,Result_Type>::
  MoveX(ArgumentType distance)
  { 
    ArgumentVector *_p_xdata = new ArgumentVector();
    ResultVector *_p_ydata = new ResultVector();
    Export(_p_xdata,_p_ydata);
    for (ArgumentIterator xit=_p_xdata->begin();xit!=_p_xdata->end();++xit) (*xit)+=distance;
    Import(_p_xdata,_p_ydata);
    delete _p_xdata;
    delete _p_ydata;
  }

  template <class Argument_Type,class Result_Type>
  Argument_Type  Data_To_Function<Argument_Type,Result_Type>::
  GetDeltaXMin(ResultType& left,ResultType& right)
  { 
    ArgumentType minimum=(ArgumentType)0.0, cur;
    if (p_xdata->size()>2) {
      ArgumentIterator xit=p_xdata->begin()+1;
      minimum=(ArgumentType)dabs((*p_xaxis)[*xit]-(*p_xaxis)[*(xit-1)]);
      left=(*p_yaxis)[(*p_xydata)[*(xit-1)]];
      right=(*p_yaxis)[(*p_xydata)[*xit]];
      for (++xit;xit!=p_xdata->end();++xit) {
	cur=(ArgumentType)dabs((*p_xaxis)[*xit]-(*p_xaxis)[*(xit-1)]);
	if (cur<minimum) {
	  minimum=cur;
	  left=(*p_yaxis)[(*p_xydata)[*(xit-1)]];
	  right=(*p_yaxis)[(*p_xydata)[*xit]];
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
  GetDeltaXMax(ResultType& left,ResultType& right)
  { 
    ArgumentType maximum=(ArgumentType)0.0, cur;
    if (p_xdata->size()>2) {
      for (ArgumentIterator xit=p_xdata->begin()+1;xit!=p_xdata->end();++xit) {
	cur=(ArgumentType)dabs((*p_xaxis)[*xit]-(*p_xaxis)[*(xit-1)]);
	if (cur>maximum) {
	  maximum=cur;
	  left=(*p_yaxis)[(*p_xydata)[*(xit-1)]];
	  right=(*p_yaxis)[(*p_xydata)[*xit]];
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
  GetDeltaYMin(ArgumentType& left,ArgumentType& right)
  { 
    ResultType minimum=(ResultType)0.0, cur;
    if (p_ydata->size()>2) {
      ResultIterator yit=p_ydata->begin()+1;
      minimum=(ResultType)dabs((*p_yaxis)[*yit]-(*p_yaxis)[*(yit-1)]);
      left=(*p_xaxis)[(*p_yxdata)[*(yit-1)]];
      right=(*p_xaxis)[(*p_yxdata)[*yit]];
      for (++yit;yit!=p_ydata->end();++yit) {
	cur=(ResultType)dabs((*p_yaxis)[*yit]-(*p_yaxis)[*(yit-1)]);
	if (cur<minimum) {
	  minimum=cur;
	  left=(*p_xaxis)[(*p_yxdata)[*(yit-1)]];
	  right=(*p_xaxis)[(*p_yxdata)[*yit]];
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
  GetDeltaYMax(ArgumentType& left,ArgumentType& right)
  { 
    ResultType maximum=(ResultType)0.0, cur;
    if (p_ydata->size()>2) {
      for (ResultIterator yit=p_ydata->begin()+1;yit!=p_ydata->end();++yit) {
	cur=(ResultType)dabs((*p_yaxis)[*yit]-(*p_yaxis)[*(yit-1)]);
	if (cur>maximum) {
	  maximum=cur;
	  left=(*p_xaxis)[(*p_yxdata)[*(yit-1)]];
	  right=(*p_xaxis)[(*p_yxdata)[*yit]];
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
