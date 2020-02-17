#include "MuonAnalysis/Scouting/plugins/TheoryMiniTuplizer.h"
using namespace std;
using namespace edm;


//
// constructors and destructor
//
TheoryMiniTuplizer::TheoryMiniTuplizer(const edm::ParameterSet& iConfig):
    token_jets(consumes<ScoutingCaloJetCollection>(
                   iConfig.getParameter<InputTag>("jet_collection"))),
    token_rho(consumes<double>(
                  iConfig.getParameter<InputTag>("rho"))),
    token_muons(consumes<ScoutingMuonCollection>(
                    iConfig.getParameter<InputTag>("muon_collection"))),
    token_MET_pt(consumes<double>(iConfig.getParameter<InputTag>("MET_pt"))),
    token_MET_phi(consumes<double>(iConfig.getParameter<InputTag>("MET_phi"))),
    file_name(iConfig.getParameter<string>("output_file_name"))
{
   //now do what ever initialization is needed
   file = new TFile(file_name.c_str(), "RECREATE");
   tree = new TTree("theoryminituplizer", "Tree for scouting data");
   
   // Initialize the ID Codes for the theory Ntuple
   m_theoryID_muon = 2.0;
   m_theoryID_tau_isotrack = 3.0;
   m_theoryID_tau_1prongpi0s = 3.1;
   m_theoryID_jet = 4.0;
   m_theoryID_pfmet = 6.0;
   m_theoryID_tcmet = 6.2;
   m_theoryID_calomet = 6.1;
   m_theoryID_st = 7.0;
   m_theoryID_lt = 8.0;
   m_theoryID_ht = 9.0;
   m_theoryID_dymass = 10.0;
   m_theoryID_wtmass = 11.0;

   nobject = 0;
   weight = 1.0;
    
   m_object = new TClonesArray("TLorentzVector");

   m_object_type = new vector<double>;
   m_object_btag = new vector<double>;
   m_object_dum1 = new vector<double>;
   m_object_dum2 = new vector<double>;
   m_object_dum3 = new vector<double>;
   m_object_dum4 = new vector<double>;
   m_object_dum5 = new vector<double>;
   m_object_dum6 = new vector<double>;
   m_object_dum7 = new vector<double>;
   m_object_dum8 = new vector<double>;
   m_object_scalar = new vector<double>;
   m_evtnum = new vector<double>;
   
   tree->Branch("weight", &weight, "weight/D");
   tree->Branch("object", "TClonesArray", &m_object, 32000, 0);
   tree->Branch("nobject", &nobject, "nobject/I");
   tree->Branch("object_type", "vector<double>", &m_object_type);
   tree->Branch("object_btag", "vector<double>", &m_object_btag);
   tree->Branch("object_dum1", "vector<double>", &m_object_dum1);
   tree->Branch("object_dum2", "vector<double>", &m_object_dum2);
   tree->Branch("object_dum3", "vector<double>", &m_object_dum3);
   tree->Branch("object_dum4", "vector<double>", &m_object_dum4);
   tree->Branch("object_dum5", "vector<double>", &m_object_dum5);
   tree->Branch("object_dum6", "vector<double>", &m_object_dum6);
   tree->Branch("object_dum7", "vector<double>", &m_object_dum7);
   tree->Branch("object_dum8", "vector<double>", &m_object_dum8);
   tree->Branch("object_scalar", "vector<double>", &m_object_scalar);
   tree->Branch("eventnum", "vector<double>", &m_evtnum);
  

}


TheoryMiniTuplizer::~TheoryMiniTuplizer() {
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   file->cd();
   tree->Write();
   file->Close();

   delete m_object;
   delete m_object_type;
   delete m_object_btag;
   delete m_object_dum1;
   delete m_object_dum2;
   delete m_object_dum3;
   delete m_object_dum4;
   delete m_object_dum5;
   delete m_object_dum6;
   delete m_object_dum7;
   delete m_object_dum8;
   delete m_object_scalar;
   delete m_evtnum;
}


//
// member functions
//

// ------------ method called for each event  ------------
void TheoryMiniTuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
   using namespace edm;
   int getCollectionsResult = GetCollections(iEvent);
    if (getCollectionsResult)
	return;

    double Ht = 0.0;
    double Lt = 0.0;
    double St = 0.0;

    addEvtNum(iEvent.id().run(), iEvent.id().event(),
              iEvent.id().luminosityBlock());

    addMET(m_theoryID_pfmet, *MET_pt, *MET_phi);
    St += *MET_pt;
    
//~    double rho = *handle_rho;

    for (auto &jet: *jets) {
        double type = m_theoryID_jet;
        double btag = 0.0;
		double dumb = 0.0;

        TLorentzVector tmp_vector;
        tmp_vector.SetPtEtaPhiM(jet.pt(), jet.eta(), jet.phi(), jet.m());

        addObject(tmp_vector, type, btag, dumb, dumb, dumb, dumb, dumb, dumb, dumb, dumb);
        Ht += jet.pt();
    }		
	
    
    for (auto &muon: *muons) {
        double type = -m_theoryID_muon*muon.charge();
        double btag = -999.0;

        TLorentzVector tmp_vector;
        tmp_vector.SetPtEtaPhiM(muon.pt(), muon.eta(), muon.phi(), muon.m());

        addObject(tmp_vector, type, btag);
        Lt += muon.pt();
    }

    St = Ht + Lt;
    addScalar(m_theoryID_ht, Ht);
    addScalar(m_theoryID_lt, Lt);
    addScalar(m_theoryID_st, St);

    tree->Fill();

    nobject = 0;

    m_object_type->clear();
    m_object_btag->clear();
    m_object_dum1->clear();
    m_object_dum2->clear();
    m_object_dum3->clear();
    m_object_dum4->clear();
    m_object_dum5->clear();
    m_object_dum6->clear();
    m_object_dum7->clear();
    m_object_dum8->clear();
    m_object_scalar->clear();
    m_evtnum->clear();

    m_object->Clear();

    return;
   
}


// ------------ method called once each job just before starting event loop  ------------
void TheoryMiniTuplizer::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void TheoryMiniTuplizer::endJob() {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TheoryMiniTuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

int TheoryMiniTuplizer::GetCollections(const edm::Event& iEvent) {
    // Get collections from ntuple
    // Returns nonzero if there is a problem getting a collection

    // Get jets
    iEvent.getByToken(token_jets, jets);
    if (!jets.isValid()) {
        throw edm::Exception(edm::errors::ProductNotFound)
            << "Could not find ScoutingCaloJetCollection." << endl;
        return 1;
    }
    
    // Get rho
    iEvent.getByToken(token_rho, handle_rho);
    if (!handle_rho.isValid()) {
        throw edm::Exception(edm::errors::ProductNotFound)
            << "Could not find rho." << endl;
        return 1;
    }

    // Get MET
    iEvent.getByToken(token_MET_pt, MET_pt);
    if (!MET_pt.isValid()) {
        throw edm::Exception(edm::errors::ProductNotFound)
            << "Could not find MET pt." << endl;
        return 1;
    }

    iEvent.getByToken(token_MET_phi, MET_phi);
    if (!MET_phi.isValid()) {
        throw edm::Exception(edm::errors::ProductNotFound)
            << "Could not find MET phi." << endl;
        return 1;
    }

    // Get muons
    iEvent.getByToken(token_muons, muons);
    if (!muons.isValid()) {
        throw edm::Exception(edm::errors::ProductNotFound)
            << "Could not find ScoutingMuonCollection." << endl;
        return 1;
    }


	return 0;
}



// Extra user defined functions below //

void TheoryMiniTuplizer::addObject(TLorentzVector fourvec, double type, double btag) {
    Double_t dum = -999.0;

    addObject(fourvec, type, btag, dum, dum, dum, dum, dum, dum, dum, dum);

    return;
}

void TheoryMiniTuplizer::addObject(TLorentzVector fourvec, double type,
                                   double btag, double dum1, double dum2,
                                   double dum3, double dum4, double dum5,
                                   double dum6, double dum7, double dum8) 
{
									   
    addObject(fourvec, type, btag, dum1, dum2, dum3, dum4, dum5, dum6, dum7,
              dum8, -999);

    return;
}

void TheoryMiniTuplizer::addObject(TLorentzVector fourvec, double type,
                                   double btag, double dum1, double dum2,
                                   double dum3, double dum4, double dum5,
                                   double dum6, double dum7, double dum8,
                                   double scalar)
{
    m_object_type->push_back(type);
    m_object_btag->push_back(btag);
    m_object_dum1->push_back(dum1);
    m_object_dum2->push_back(dum2);
    m_object_dum3->push_back(dum3);
    m_object_dum4->push_back(dum4);
    m_object_dum5->push_back(dum5);
    m_object_dum6->push_back(dum6);
    m_object_dum7->push_back(dum7);
    m_object_dum8->push_back(dum8);
    m_object_scalar->push_back(scalar);

    new ((*m_object)[nobject]) TLorentzVector(fourvec);

    ++nobject;

    return;
}


void TheoryMiniTuplizer::addEvtNum(double run, double evt, double lumi) {
    m_evtnum->push_back(run);
    m_evtnum->push_back(evt);
    m_evtnum->push_back(lumi);

    return;
}

void TheoryMiniTuplizer::addMET(double type, double pt, double phi) {
    double e = pt;
    double px = pt*cos(phi);
    double py = pt*sin(phi);
    double pz = 0;

    Double_t btag = 0;
    Double_t dum = 0;
    Double_t scalar = 0;

    TLorentzVector fourvec(px, py, pz, e);

    addObject(fourvec, type, btag, dum, dum, dum, dum, dum, dum, dum, dum,
              scalar);

    return;
}

void TheoryMiniTuplizer::addScalar(double type, double value) {
    Double_t dum = -999.0;

    addObject(TLorentzVector(0, 0, 0, 0), type, 0, dum, dum, dum, dum, dum, dum,
              dum, dum, value);

    return;
}

//define this as a plug-in
DEFINE_FWK_MODULE(TheoryMiniTuplizer);
