#include "TFile.h"
#include "TBrowser.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TTree.h"
#include "/afs/cern.ch/user/s/sandhya/scratch0/CMSSW_3_1_2/src/Analysis/Analyzer/interface/CrystalInfo.h"
#include "FWCore/Framework/interface/TriggerNames.h"
#include <string>

class Analyzer : public edm::EDAnalyzer {
 public:
  explicit Analyzer(const edm::ParameterSet&);
  ~Analyzer();


 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  int nevents;
  int is_nunu_event;
  int is_elec_event;
  int is_gamjet_event;
  int is_diphoton_event;
  int is_mu_event;
  int is_tau_event;
  int is_Z_event;
  int n_diphoton_events;
  int n_gamJet_events;
  int ngenphotons;
  int nrecPhotons;
  int nmuons;
  int ntaus;
  int nelectrons;
  int nTracks;
  int number_of_jets;
  //HLT

  edm::TriggerNames triggerNames_;  // TriggerNames class
  std::vector<std::string>  hlNames_;  // name of each HLT algorithm
 
  int is_HLT_MET50_event;
  int is_HLT_MET75_event;
  int is_HLT_Photon25_event;
  std::map<std::string,int> HLT_chosen;
  std::map<std::string,int> L1_chosen;

  std::vector<std::string> JET_CORR;

  edm::InputTag eleLabel_;
  edm::InputTag muoLabel_;
  edm::InputTag jetLabel_;
  edm::InputTag tauLabel_;
  edm::InputTag metLabel_;
  edm::InputTag phoLabel_;
  edm::InputTag rechitBLabel_;
  edm::InputTag rechitELabel_;
  edm::InputTag hlTriggerResults_;  // Input tag for TriggerResults
  edm::InputTag Tracks_;
  edm::InputTag Vertices_;
  std::string outFile_;
  bool rungenParticleCandidates_;
  bool runphotons_;
  bool runmet_;
  bool runelectrons_;
  bool runmuons_;
  bool runjets_;
  bool runtaus_;
  bool runHLT_;
  bool runtracks_;
  bool runrechit_;
  bool runvertex_;
  bool init_;

  //output root file
  TFile *f;
  
  // output root tree
  TTree* myEvent;
  //variables to be filled in the tree


  //vertex variables
  float vx[100];
  float vy[100];
  float vz[100];
  float chi2[100];

  //track variables
  double trk_pt[250];
  double trk_px[250];
  double trk_py[250];
  double trk_pz[250];
  double trk_eta[250];
  double trk_phi[250];

  //jet variables
  double jet_pt[250];
  double jet_px[250];
  double jet_py[250];
  double jet_pz[250];
  double jet_eta[250];
  double jet_phi[250];
  double jet_emEnergyFraction[250];
  double jet_energyFractionHadronic[250];

  //electron variables
  double electron_pt[250];
  double electron_px[250];
  double electron_py[250];
  double electron_pz[250];
  double electron_energy[250];
  double electron_charge[250];
  double electron_eta[250];
  double electron_phi[250];
  double electron_trkIso[250];

  //muon variables
  double muon_pt[250];
  double muon_px[250];
  double muon_py[250];
  double muon_pz[250];
  double muon_energy[250];
  double muon_charge[250];
  double muon_eta[250];
  double muon_phi[250];

  //tau variables
  double tau_pt[250];
  double tau_px[250];
  double tau_py[250];
  double tau_pz[250];
  double tau_energy[250];
  double tau_charge[250];
  double tau_eta[250];
  double tau_phi[250];

  //gen level variables
  double gen_pho_pt[250];
  double gen_pho_px[250];
  double gen_pho_py[250];
  double gen_pho_pz[250];
  double gen_pho_phi[250];
  double gen_pho_eta[250];
  double gen_pho_E[250];
  
  double gen_Hpho_pt[2];
  double gen_Hpho_px[2];
  double gen_Hpho_py[2];
  double gen_Hpho_pz[2];
  double gen_Hpho_phi[2];
  double gen_Hpho_eta[2];
  double gen_Hpho_E[2];
  int gen_Hpho_SerialNumber[2];

  
  //gen level variables for Zdaughter
  double gen_Zdaughter_pt[2];
  double gen_Zdaughter_px[2];
  double gen_Zdaughter_py[2];
  double gen_Zdaughter_pz[2];
  double gen_Zdaughter_phi[2];
  double gen_Zdaughter_eta[2];
  double gen_Zdaughter_theta[2];
  double gen_Zdaughter_E[2];

  //gen level variables for Z
  double gen_Zboson_pt;
  double gen_Zboson_px;
  double gen_Zboson_py;
  double gen_Zboson_pz;
  double gen_Zboson_phi;
  double gen_Zboson_eta;
  double gen_Zboson_theta;
  double gen_Zboson_E;
  
  double pho_E[250];
  double pho_pt[250];
  double pho_eta[250];
  double pho_phi[250];
  double pho_px[250];
  double pho_py[250];
  double pho_pz[250];
  double pho_r9[250];
  double pho_theta[250];
  double pho_et[250];
  int    pho_isConverted[250];
  
  //isolation variables
  double pho_ecalRecHitIso[250];
  double pho_hcalRecHitIso[250];
  double pho_HollowTrackConeIso[250];
  double pho_HoE[250];

  //SC variables
  double pho_sc_energy[250];
  int    pho_size[250];
  double pho_sc_eta[250];
  double pho_sc_phi[250];
  double pho_sc_etaWidth[250];
  double pho_sc_phiWidth[250];
  double pho_sc_et[250];

  //gen matched photon 
  double matchpho_E[250];
  double matchpho_pt[250];
  double matchpho_eta[250];
  double matchpho_phi[250];
  double matchpho_px[250];
  double matchpho_py[250];
  double matchpho_pz[250];
  int    ismatchedpho[250];
  
  //converted photon variabes
  unsigned int pho_nTracks[250];
  double pho_pairInvariantMass[250];
  double pho_pairCotThetaSeparation[250];
  double pho_pairMomentum_x[250];
  double pho_pairMomentum_y[250];
  double pho_pairMomentum_z[250];
  double pho_EoverP[250];
  double pho_vertex_x[250];
  double pho_vertex_y[250];
  double pho_vertex_z[250];
  double pho_zOfPrimaryVertex[250];

  //rechit information
  double pho_timing_xtalEB[250][100];
  double pho_timingavg_xtalEB[250];
  double pho_energy_xtalEB[250][100];
  double pho_ieta_xtalEB[250][100];
  double pho_iphi_xtalEB[250][100];

  //calomet variables
  double CaloMetSig;
  double CaloMetCorr;
  double CaloMetEt;
  double CaloMetEx;
  double CaloMetEy;
  double CaloMetPhi;
  double CaloEtFractionHadronic;
  double CaloEmEtFraction;
  double CaloHadEtInHB;
  double CaloHadEtInHO;
  double CaloHadEtInHE;
  double CaloHadEtInHF;
  double CaloEmEtInEB;
  double CaloEmEtInEE;
  double CaloEmEtInHF;
  double CaloMetEz;
  double CaloMaxEtInEmTowers;
  double CaloSumEt;
  double CaloMaxEtInHadTowers;
  
  double Delta_phi;
};

