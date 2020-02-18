// System include files
#include <memory>
#include <iostream>
#include <vector>
#include <utility>


// CMSSW include files
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Scouting/interface/ScoutingCaloJet.h"
#include "DataFormats/Scouting/interface/ScoutingParticle.h"
#include "DataFormats/Scouting/interface/ScoutingVertex.h"
#include "DataFormats/Scouting/interface/ScoutingMuon.h"


// Root include files
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"

///////////////////////////////////////

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"

#include <vector>
#include <string>
#include <iostream>


// Standard C++ includes                                                                                                            
                                                                                                                                  
#include <memory>
#include <vector>
#include <string>
#include <iostream>
#include <map>

// ROOT includes                                                                                                                    
                                                                                                                                   
#include <TPRegexp.h>

// CMSSW framework includes                                                                                                         
                                                                                                                                    
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// CMSSW data formats                                                                                                               
                                                                                                                                    
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "DataFormats/Scouting/interface/ScoutingMuon.h"
#include "DataFormats/Scouting/interface/ScoutingVertex.h"

// Other relevant CMSSW includes                                                                                                    
                                                                                                                                    
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"


#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionData.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionEvaluator.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionParser.h"


#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"

// Beam Spot includes

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "TMath.h"

///////////////////////////////////////

class ScoutingNtuplizer : public edm::EDAnalyzer {
   public:
      explicit ScoutingNtuplizer(const edm::ParameterSet&);
      ~ScoutingNtuplizer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      virtual int GetCollections(const edm::Event&);
      
      void ResetVariables();


      //   private:

      //HLTPrescaleProvider hltPSProv_;
      //std::string hltProcess_; //name of HLT process, usually "HLT"                                                                      

      //   public:
      //explicit L1MenuExample(const edm::ParameterSet& iConfig);
      //~L1MenuExample(){}

      //private:
      //      virtual void beginJob(){}
      virtual void beginRun(const edm::Run& run,const edm::EventSetup& iSetup);
      //virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
      //virtual void endJob(){}

      //edm::InputTag l1AlgoTag;
      //edm::EDGetTokenT<BXVector<GlobalAlgBlk> > l1AlgoToken;


    // ----------member data ---------------------------


      //L1 trigger

      edm::InputTag l1AlgoTag;
      edm::EDGetTokenT<BXVector<GlobalAlgBlk> > l1AlgoToken;

      std::vector<int> triggerfilter;
      std::vector<int> prescaleL1;

      // Jets
      edm::EDGetTokenT<ScoutingCaloJetCollection> token_jets;
      edm::Handle<ScoutingCaloJetCollection> jets;

      int jet_num;
      std::vector<float> jet_pt;
      std::vector<float> jet_eta;
      std::vector<float> jet_phi;
      std::vector<float> jet_m;
      
      std::vector<float> jet_Area;
      std::vector<float> jet_maxEInEmTowers;
      std::vector<float> jet_maxEInHadTowers;
      std::vector<float> jet_hadEnergyInHB;
      std::vector<float> jet_hadEnergyInHE;
      std::vector<float> jet_hadEnergyInHF;
      std::vector<float> jet_emEnergyInEB;
      std::vector<float> jet_emEnergyInEE;
      std::vector<float> jet_emEnergyInHF;
      std::vector<float> jet_towersArea;
      std::vector<float> jet_mvaDiscriminator;
      std::vector<float> jet_btagDiscriminator;      
     
      double HT;

  
 
    // Muon Data    
      edm::EDGetTokenT<ScoutingMuonCollection> token_muons;
      edm::Handle<ScoutingMuonCollection> muons;
  
      int muon_num;

      std::vector<int> muon_q;
      std::vector<float> muon_pt;
      std::vector<float> muon_eta;
      std::vector<float> muon_phi;
      
      std::vector<int> muon_isGlobalMuon;
      std::vector<int> muon_isTrackerMuon;
      
      std::vector<float> muon_dxy;
      std::vector<float> muon_edxy;
      std::vector<float> muon_dz;
      std::vector<float> muon_edz;
      std::vector<float> muon_chisquared;
      std::vector<float> muon_ndof;
      std::vector<float> muon_trackIso;
      std::vector<int> muon_nmatchedstations;
      std::vector<int> muon_ntrklayersmeasurement;
      std::vector<int> muon_nmuonhits;
      std::vector<int> muon_npixelhits;
      std::vector<int> muon_nstriphits;

      std::vector<std::vector<int>> muon_vtxindex;


    // Event Data
      edm::EDGetTokenT<double> token_rho;
      edm::EDGetTokenT<double> token_MET_pt;
      edm::EDGetTokenT<double> token_MET_phi;
      
      edm::Handle<double> handle_rho;
      edm::Handle<double> handle_MET_pt;
      edm::Handle<double> handle_MET_phi;
      
      edm::EDGetTokenT<ScoutingVertexCollection> token_dispvertices;
      edm::Handle<ScoutingVertexCollection> dispvertices;

      int dispvertex_num;   
      
      std::vector<float> dispvtx_x;
      std::vector<float> dispvtx_y;
      std::vector<float> dispvtx_z;
      
      std::vector<float> dispvtx_ex;
      std::vector<float> dispvtx_ey;
      std::vector<float> dispvtx_ez;

      std::vector<float> dispvtx_chisquared;
      std::vector<float> dispvtx_ndof;
      std::vector<float> dispvtx_tracksize;

      edm::EDGetTokenT<ScoutingVertexCollection> token_privertices;
      edm::Handle<ScoutingVertexCollection> privertices;

      int privertex_num;

      std::vector<float> privtx_x;
      std::vector<float> privtx_y;
      std::vector<float> privtx_z;
      
      std::vector<float> privtx_ex;
      std::vector<float> privtx_ey;
      std::vector<float> privtx_ez;

      
      float rho;
      float MET_pt;
      float MET_phi;
      
      int run;
      int lumi;
      double event;
      
      std::string file_name;
      TFile *file;
      TTree *tree;
      
      int min_muons;

      //  trigger information

      HLTPrescaleProvider hltPSProv_;
      std::string hltProcess_; //name of HLT process, usually "HLT"

      edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
      edm::Handle<edm::TriggerResults> triggerBits;

      std::vector<std::string> hltseedsvector;
      std::vector<pair<string,int>> hltbitmap;
      //      std::vector<pair<string,int>> hltprescalemap;

      std::vector<std::string> l1seedsvector;
      std::vector<pair<string,int>> l1bitmap;
      std::vector<pair<string,int>> l1prescalemap;

      // Beamspot

      edm::EDGetTokenT<reco::BeamSpot> beamSpotToken;
      edm::Handle<reco::BeamSpot> beamSpotH;

      float BS_x;
      float BS_y;
      float BS_z;

      edm::EDGetTokenT<std::vector<int>> token_muonexpectedhits;
      edm::Handle<std::vector<int> > handle_muonexpectedhits;

      std::vector<int> muon_nexpectedhitsmultiple;

};
