#include "TLorentzVector.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "TObject.h" 
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"

//utility function prototypes
double deltaphi(double phi1, double phi2);
double correct_phi(double phi);
double delta_R(double phi,double eta);
double Theta(double eta);

class CrystalInfo{
 public:
  // methods
  CrystalInfo();
  ~CrystalInfo();
  
  // variables
  //DetId id;
  uint32_t rawId;
  int ieta;
  int iphi;
  int ix;
  int iy;
  double energy;
  double time;
};

CrystalInfo::CrystalInfo()
{}
CrystalInfo::~CrystalInfo()
{}  
