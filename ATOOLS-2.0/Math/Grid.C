#include "Grid.H"

#include "MathTools.H"
#include "Message.H"
#include "MyStrStream.H"
#ifdef DEBUG__Grid
#include "Blob.H"
#endif

using namespace ATOOLS;

Grid::Grid(const size_t dim,const size_t points):
  m_dim(dim),
  m_points(points),
  m_order(3)
{
  p_root = new Node<double>(0.0,true);
  Allocate(p_root,1);
}

Grid::~Grid()
{
  delete p_root;
}

void Grid::Allocate(Node<double> *current,const size_t node)
{
#ifdef DEBUG__Grid
  msg.Debugging()<<"Grid<"<<m_dim<<","<<m_points
		 <<">::Allocate("<<current<<","<<node<<") {\n";
  {
    msg_Indent();
#endif
    (*current).resize(m_points);
    if (node<=m_dim) (*current)().resize(m_points);
    for (size_t i=0;i<m_points;++i) {
      (*current)[i]=0.0;
      if (node<=m_dim) {
	(*current)()[i] = new Node<double>(0.0,node+1<=m_dim);
	if (node<m_dim) Allocate((*current)()[i],node+1);
      }
    }
#ifdef DEBUG__Grid
  }
  msg.Debugging()<<"}"<<std::endl;
#endif
}

double Grid::Y(const std::vector<size_t> &position,
	       Node<double> *current,const size_t node) const
{
  if (node>position.size()+1 || 
      (node<=position.size() && position[node-1]>m_points)) {
    msg.Error()<<"Grid<"<<m_dim<<","<<m_points<<">::Y(..): "
	       <<"Position specification failed."<<std::endl;
    return 0.0;
  }
#ifdef DEBUG__Grid
  msg.Debugging()<<"Grid<"<<m_dim<<","<<m_points<<">::Y("
		 <<&position<<","<<current<<","<<node
		 <<") ["<<(node>position.size()?0:position[node-1])<<"]\n";
  msg_Indent();
#endif
  if (current->operator->()==NULL) return (*current)[0];
  return Y(position,(*current)()[position[node-1]],node+1);
}

void Grid::SetY(const std::vector<size_t> &position,
		const double value,
		Node<double> *current,const size_t node) 
{
  if (node>position.size()+1 || 
      (node<=position.size() && position[node-1]>m_points)) {
    msg.Error()<<"Grid<"<<m_dim<<","<<m_points<<">::SetY(..): "
	       <<"Position specification failed."<<std::endl;
    return;
  }
#ifdef DEBUG__Grid
  msg.Debugging()<<"Grid<"<<m_dim<<","<<m_points<<">::SetY("
		 <<&position<<","<<value<<","<<current<<","<<node
		 <<") ["<<(node>position.size()?0:position[node-1])<<"]\n";
  msg_Indent();
#endif
  if (current->operator->()==NULL) (*current)[0]=value;
  else SetY(position,value,(*current)()[position[node-1]],node+1);
}

void Grid::SetX(const size_t position,const std::vector<double> &x,
		Node<double> *current,const size_t node) 
{
  if (position>m_dim) {
    msg.Error()<<"Grid<"<<m_dim<<","<<m_points<<">::SetX(..): "
	       <<"Position specification failed."<<std::endl;
    return;
  }
#ifdef DEBUG__Grid
  msg.Debugging()<<"Grid<"<<m_dim<<","<<m_points<<">::SetX("
		 <<position<<","<<x<<","<<current<<","<<node
		 <<")\n";
  msg_Indent();
#endif
  if (node-1==position) {
    if (x.size()!=current->size()) {
      msg.Error()<<"Grid<"<<m_dim<<","<<m_points<<">::SetX(..): "
		 <<"Incorrect point number."<<std::endl;
      return;
    }
    (*(std::vector<double>*)current)=x;
  }
  else {
    for (size_t j=0;j<(*current)().size();++j)
      SetX(position,x,(*current)()[j],node+1);
  }
}

double Grid::Interpolate(const std::vector<double> &xgrid,
			 const std::vector<double> &ygrid,
			 const double &x) const
{
  // This is NR::polin from "Numerical Recipes in C++" ISBN 0 521 75033 4
#ifdef DEBUG__Grid
  msg_Debugging()<<"Grid::Interpolate("<<xgrid<<","
		 <<ygrid<<","<<x<<"):"<<std::endl;
#endif
  int i, m, ns=0;
  double den, dif, dift, ho, hp, w, y, dy;
  int n=xgrid.size();
  std::vector<double> c(n), d(n);
  dif=dabs(x-xgrid[0]);
  for (i=0;i<n;++i) {
    if ((dift=dabs(x-xgrid[i]))<dif) {
      ns=i;
      dif=dift;
    }
    c[i]=d[i]=ygrid[i];
  }
  y=ygrid[ns--];
  for (m=1;m<n;++m) {
    for (i=0;i<n-m;++i) {
      ho=xgrid[i]-x;
      hp=xgrid[i+m]-x;
      w=c[i+1]-d[i];
      if ((den=ho-hp)==0.0) 
	msg.Error()<<"Grid::Interpolate(..): Identical x values.\n"
		   <<"   Result will be unreliable."<<std::endl;
      den=w/den;
      d[i]=hp*den;
      c[i]=ho*den;
    }
    y+=dy=(2*(ns+1)<(n-m)?c[ns+1]:d[ns--]);
  }
  return y;
}

void Grid::Interpolate(const std::vector<double> &x,const size_t i,
		       Node<double> *const node,std::vector<double> &y,
		       const size_t j) const
{
  // This is NR::polin from "Numerical Recipes in C++" ISBN 0 521 75033 4
#ifdef DEBUG__Grid
  msg_Debugging()<<"Grid::Interpolate("<<x<<","<<i<<","
		 <<node<<","<<y<<","<<j<<"):"<<std::endl;
#endif
  if (node->operator->()==NULL) y[j]=(*node)[0];
  else {
    msg_Indent();
    size_t l,r;
    Find(x[i],node,l,r);
    l=Max(l-m_order/2,(size_t)0);
    l=Min(l,(*node)().size()-2);
    std::vector<double> yt(Min((*node)().size()-l,m_order)), xt(yt.size());
    for (size_t k=0;k<yt.size();++k) {
      Interpolate(x,i+1,(*node)()[l+k],yt,k);
      xt[k]=(*node)[l+k];
    }
    y[j]=Interpolate(xt,yt,x[i]);
  }
}

double Grid::Interpolate(const std::vector<double> &x) const
{
  // This is NR::polin from "Numerical Recipes in C++" ISBN 0 521 75033 4
  double cur;
#ifdef DEBUG__Grid
  msg_Debugging()<<"Grid::Interpolate("<<x<<"): {"<<std::endl;
  {
    msg_Indent();
#endif
    Grid *grid=(Grid*)this;
    size_t l, r;
    Find(x[0],p_root,l,r);
    l=Max(l-m_order/2,(size_t)0);
    l=Min(l,(*p_root)().size()-2);
    std::vector<double> yt(Min((*p_root)().size()-l,m_order)), xt(yt.size());
    for (size_t k=0;k<yt.size();++k) {
      Interpolate(x,1,(*grid->p_root)()[l+k],yt,k);
      xt[k]=(*grid->p_root)[l+k];
    }
    cur=Interpolate(xt,yt,x[0]);
#ifdef DEBUG__Grid
    msg_Debugging()<<"return "<<cur<<"\n";
  }
  msg_Debugging()<<"}"<<std::endl;
#endif
  return cur;
}

bool Grid::WriteOut(const std::string &filename)
{
  p_file = new std::fstream((m_file=filename).c_str(),std::ios::out);
  if (p_file->bad()) {
    delete p_file;
    return false;
  }
  p_file->precision(14);
  (*p_file)<<m_dim<<" "<<m_points<<" {"<<std::endl;
  bool result=WriteOut(p_root,0);
  (*p_file)<<"}"<<std::endl;
  delete p_file;
  return result;
}

bool Grid::WriteOut(Node<double> *current,const size_t node)
{
  bool result=true;
  for (size_t i=0;i<current->size();++i) {
    (*p_file)<<(node!=m_dim?std::string(node,' '):" \t")<<(*current)[i]<<" ";
    if (node<m_dim-1) (*p_file)<<std::endl;
    if (current->operator->()!=NULL) 
      if (!WriteOut((*current)()[i],node+1)) result=false;
    if (node==m_dim-1) (*p_file)<<std::endl;
  }
  return result;
}

bool Grid::ReadIn(const std::string &filename)
{
  msg_Debugging()<<"Grid::ReadIn(\""<<filename<<"\"):"<<std::endl;
  p_file = new std::fstream((m_file=filename).c_str(),std::ios::in);
  if (!p_file->good()) {
    msg_Info()<<"Grid::ReadIn(\""<<filename<<"\"): "
	      <<"Cannot find file."<<std::endl;
    delete p_file;
    return false;
  }
  std::string dummy;
  size_t dim, points;
  p_file->precision(14);
  (*p_file)>>dim>>points>>dummy;
  if (dim!=m_dim || points!=m_points || dummy!="{") {
    msg.Error()<<"Grid::ReadIn(\""<<filename<<"\"): "
	       <<"Inconsistent grid dimension."<<std::endl;
    delete p_file;
    return false;
  }
  bool result=ReadIn(p_root,0);
  if (!p_file->eof()) (*p_file)>>dummy;
  if (dummy!="}") result=false;
  delete p_file;
  if (!result) msg.Error()<<"Grid::ReadIn(\""<<filename<<"\"): "
			  <<"Data error."<<std::endl;
  return result;
}

bool Grid::ReadIn(Node<double> *current,const size_t node)
{
  for (size_t i=0;i<current->size();++i) {
    if (p_file->eof()) {
      msg.Error()<<"Grid::ReadIn(..): "
		 <<"EOF at "<<i<<" ";
      return false;
    }
    std::string dummy;
    (*p_file)>>dummy;
    (*current)[i]=ATOOLS::ToType<double>(dummy);
    if (dummy=="nan") (*current)[i]=0.0;
    if (ATOOLS::msg.LevelIsDebugging()) {
      msg_Debugging()<<std::string(node,' ')<<(*current)[i]<<" ";
      if (node!=m_dim-1) msg_Debugging()<<std::endl;
    }
    if (current->operator->()!=NULL) {
      if (!ReadIn((*current)()[i],node+1)) {
        if (node==0) msg.Error()<<std::endl;
	else msg.Error()<<i<<" ";
	return false;
      }
    }
  }
  return true;
}

double Grid::Data(const std::vector<double> &x) const
{
  return Data(x,p_root,0);
}

double Grid::Data(const std::vector<double> &x,
		  ATOOLS::Node<double> *const node,const size_t i) const
{
  if (node->operator->()==NULL) return (*node)[0];
  size_t l, r;
  l=Find(x[i],node,l,r);
  return Data(x,(*node)()[l],i+1);
}

size_t Grid::Find(const double &x,ATOOLS::Node<double> *const node,
		  size_t &l,size_t &r) const
{
  l=0, r=node->size()-1;
  size_t j=(l+r)/2;
  double xj=(*node)[j];
  while (r-l>1) {
    if (x<xj) r=j;
    else l=j;
    j=(l+r)/2;
    xj=(*node)[j];
  }
  return dabs(x-(*node)[l])<
    dabs(x-(*node)[Min(r,node->size()-1)])?l:r;
}

