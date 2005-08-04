#ifndef ResonanceFlavour_h
#define ResonanceFlavour_h

namespace HADRONS {

  class ResonanceFlavour {

	std::string name;
	double mass, width;

	public:
	ResonanceFlavour() {};
	ResonanceFlavour( std::string _name, double _mass, double width );
	ResonanceFlavour( ATOOLS::kf::code _kfc, double _mass, double width );
	~ResonanceFlavour() {}

	void Set( ATOOLS::kf::code _kfc, double _mass, double width );
	 
	double Mass() { return mass; }
	double Width() { return width; }
	std::string Name() { return name; }

  };

}

#endif
  
