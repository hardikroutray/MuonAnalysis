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

//
// class declaration
//

class TheoryMiniTuplizer : public edm::EDAnalyzer {
public:
    explicit TheoryMiniTuplizer(const edm::ParameterSet&);
    ~TheoryMiniTuplizer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;
    virtual int GetCollections(const edm::Event&);
    
    virtual void addObject(TLorentzVector, double, double, double, double,
                           double, double, double, double, double, double);
    virtual void addObject(TLorentzVector, double, double, double, double,
                           double, double, double, double, double, double,
                           double);
    virtual void addObject(TLorentzVector, double, double);
    virtual void addMET(double, double, double);
    virtual void addScalar(double, double);
    virtual void addEvtNum(double, double, double);

    // ----------member data ---------------------------
    edm::EDGetTokenT<ScoutingCaloJetCollection> token_jets;
    edm::EDGetTokenT<double> token_rho;
    edm::EDGetTokenT<ScoutingMuonCollection> token_muons;
    edm::EDGetTokenT<double> token_MET_pt;
    edm::EDGetTokenT<double> token_MET_phi;

    edm::Handle<ScoutingCaloJetCollection> jets;
    edm::Handle<double> handle_rho;
    edm::Handle<ScoutingMuonCollection> muons;
    edm::Handle<double> MET_pt;
    edm::Handle<double> MET_phi;

    int nobject;
    double weight;
    
    std::vector<double>* m_object_type;
    std::vector<double>* m_object_btag;
    std::vector<double>* m_object_dum1;
    std::vector<double>* m_object_dum2;
    std::vector<double>* m_object_dum3;
    std::vector<double>* m_object_dum4;
    std::vector<double>* m_object_dum5;
    std::vector<double>* m_object_dum6;
    std::vector<double>* m_object_dum7;
    std::vector<double>* m_object_dum8;
    std::vector<double>* m_object_scalar;
    std::vector<double>* m_evtnum;

    TClonesArray *m_object;

    Double_t m_theoryID_muon;
    Double_t m_theoryID_tau_isotrack;
    Double_t m_theoryID_tau_1prongpi0s;
    Double_t m_theoryID_jet;
    Double_t m_theoryID_pfmet;
    Double_t m_theoryID_tcmet;
    Double_t m_theoryID_calomet;
    Double_t m_theoryID_st;
    Double_t m_theoryID_lt;
    Double_t m_theoryID_ht;
    Double_t m_theoryID_dymass;
    Double_t m_theoryID_wtmass;
    
    std::string file_name;
    TFile *file;
    TTree *tree;

};
