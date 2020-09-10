class NLOHisto1D : public YODA::Histo1D {

YODA::Histo1D* _tmphist;
int _current_event_number;

void _syncHists() {
  for (size_t i=0; i<_tmphist->bins().size(); ++i) {
    if (_tmphist->bin(i).area()) YODA::Histo1D::fillBin(i, _tmphist->bin(i).area());
  }
    if (_tmphist->overflow().sumW())  YODA::Histo1D::overflow()+=_tmphist->overflow();
    if (_tmphist->underflow().sumW()) YODA::Histo1D::underflow()+=_tmphist->underflow();
    _tmphist->reset();
  }

public:

  NLOHisto1D(size_t nbins, double lower, double upper, const string& path) :
    YODA::Histo1D(nbins, lower, upper, path),
    _current_event_number(-1)
  {
    _tmphist = new Histo1D(nbins, lower, upper, path+"_tmp");
  }
 
  NLOHisto1D(const vector<double>& binedges, const string& path) :
    YODA::Histo1D(binedges, path),
    _current_event_number(-1)
  {
    _tmphist = new Histo1D(binedges, path+"_tmp");
  }

  ~NLOHisto1D()
  {
    delete _tmphist;
  }
   
  void fill(double x, const Event& event)
  {
    if (_current_event_number==-1)
      _current_event_number = event.genEvent()->event_number();

    if (event.genEvent()->event_number()!=_current_event_number) {
      _syncHists();
      _current_event_number = event.genEvent()->event_number();
    }

    _tmphist->fill(x, event.weight());
  }

  void fillBin(size_t i, const Event& event, const double& fac)
  {
    if (_current_event_number==-1)
      _current_event_number = event.genEvent()->event_number();

    if (event.genEvent()->event_number()!=_current_event_number) {
      _syncHists();
      _current_event_number = event.genEvent()->event_number();
    }

    _tmphist->fillBin(i, event.weight()*fac);
  }

  void finalize()
  {
    _syncHists();
  }
};

typedef shared_ptr<NLOHisto1D> NLOHisto1DPtr;


NLOHisto1DPtr bookNLOHisto1D(const string& hname,
                                   size_t nbins, double lower, double upper)
{
  NLOHisto1DPtr hist( new NLOHisto1D(nbins, lower, upper, histoPath(hname)) );
  addAnalysisObject(hist);
  return hist;
}

NLOHisto1DPtr bookNLOHisto1D(const string& hname, const vector<double>& binedges)
{
  NLOHisto1DPtr hist( new NLOHisto1D(binedges, histoPath(hname)) );
  addAnalysisObject(hist);
  return hist;
}
