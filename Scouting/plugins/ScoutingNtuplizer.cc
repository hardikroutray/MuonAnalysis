#include "MuonAnalysis/Scouting/plugins/ScoutingNtuplizer.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
//#include "PrescaleProvider.h"
using namespace std;
using namespace edm;



//*****************************//
// constructors and destructor //
//*****************************//
ScoutingNtuplizer::ScoutingNtuplizer(const edm::ParameterSet& iConfig):
    token_jets(consumes<ScoutingCaloJetCollection>(iConfig.getParameter<InputTag>("jet_collection"))),
    token_muons(consumes<ScoutingMuonCollection>(iConfig.getParameter<InputTag>("muon_collection"))),
    //    token_trgResults(consumes<edm::TriggerResults>(iConfig.getParameter<InputTag>("triggerResults"))),
    token_rho(consumes<double>(iConfig.getParameter<InputTag>("rho"))),
    token_MET_pt(consumes<double>(iConfig.getParameter<InputTag>("MET_pt"))),
    token_MET_phi(consumes<double>(iConfig.getParameter<InputTag>("MET_phi"))),
    token_dispvertices(consumes<ScoutingVertexCollection>(iConfig.getParameter<InputTag>("displaced_vertex_collection"))),
    token_privertices(consumes<ScoutingVertexCollection>(iConfig.getParameter<InputTag>("primary_vertex_collection"))),
    file_name(iConfig.getParameter<string>("output_file_name")),
    min_muons(iConfig.getParameter<int>("mu_min")),
    hltPSProv_(iConfig,consumesCollector(),*this), //it needs a referernce to the calling module for some reason, hence the *this   
    hltProcess_(iConfig.getParameter<std::string>("hltProcess")),
    triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
    hltseedsvector(iConfig.getParameter<std::vector<std::string>>("hltseeds")),
    l1seedsvector(iConfig.getParameter<std::vector<std::string>>("l1seeds")),
    beamSpotToken(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("bsCollection")))
  
{
   //now do what ever initialization is needed
   file = new TFile(file_name.c_str(), "RECREATE");
   tree = new TTree("scoutingntuplizer", "Tree for scouting data");

   tree->Branch("rho", &rho, "rho/F");
   tree->Branch("Run", &run, "Run/I");
   tree->Branch("Lumi", &lumi, "Lumi/I");
   tree->Branch("Event", &event, "Event/D");

   tree->Branch("BS_x", &BS_x, "BS_x/F");
   tree->Branch("BS_y", &BS_y, "BS_y/F");
   tree->Branch("BS_z", &BS_z, "BS_z/F");
   
   tree->Branch("HT", &HT);
   
   tree->Branch("jet_num", &jet_num, "jet_num/I");
   tree->Branch("jet_pt",  &jet_pt);
   tree->Branch("jet_eta", &jet_eta);
   tree->Branch("jet_phi", &jet_phi);
   tree->Branch("jet_m", &jet_m);
      
   tree->Branch("jet_Area", &jet_Area);
   tree->Branch("jet_maxEInEmTowers", &jet_maxEInEmTowers);
   tree->Branch("jet_maxEInHadTowers", &jet_maxEInHadTowers);
   tree->Branch("jet_hadEnergyInHB", &jet_hadEnergyInHB);
   tree->Branch("jet_hadEnergyInHE", &jet_hadEnergyInHE);
   tree->Branch("jet_hadEnergyInHF", &jet_hadEnergyInHF);
   tree->Branch("jet_emEnergyInEB", &jet_emEnergyInEB);
   tree->Branch("jet_emEnergyInEE", &jet_emEnergyInEE);
   tree->Branch("jet_emEnergyInHF", &jet_emEnergyInHF);
   tree->Branch("jet_towersArea", &jet_towersArea);
   tree->Branch("jet_mvaDiscriminator", &jet_mvaDiscriminator);
   tree->Branch("jet_btagDiscriminator", &jet_btagDiscriminator);

   
   tree->Branch("muon_num", &muon_num, "muon_num/I");
   tree->Branch("muon_q", &muon_q);
   tree->Branch("muon_pt",  &muon_pt);
   tree->Branch("muon_eta", &muon_eta);
   tree->Branch("muon_phi", &muon_phi);
   tree->Branch("isGlobalMuon", &muon_isGlobalMuon);
   tree->Branch("isTrackerMuon", &muon_isTrackerMuon);
   tree->Branch("muon_dxy", &muon_dxy);
   tree->Branch("muon_edxy", &muon_edxy);
   tree->Branch("muon_dz", &muon_dz);
   tree->Branch("muon_edz", &muon_edz);
   tree->Branch("muon_chi2", &muon_chisquared);
   tree->Branch("muon_ndof", &muon_ndof);
   tree->Branch("muon_trackIso", &muon_trackIso);
   tree->Branch("nmatchedstations", &muon_nmatchedstations);
   tree->Branch("ntrackerlayerswithmeasurement", &muon_ntrklayersmeasurement);
   tree->Branch("nvalidmuonhits", &muon_nmuonhits);
   tree->Branch("nvalidpixelhits", &muon_npixelhits);
   tree->Branch("nvalidstriphits", &muon_nstriphits);
   tree->Branch("vertex_index", &muon_vtxindex);
   
   
   tree->Branch("dispvertex_num", &dispvertex_num, "dispvertex_num/I");
   tree->Branch("dispvertex_tracksize", &dispvtx_tracksize);
   tree->Branch("dispvertex_chi2", &dispvtx_chisquared);
   tree->Branch("dispvertex_ndof", &dispvtx_ndof);
   tree->Branch("dispvertex_x", &dispvtx_x);
   tree->Branch("dispvertex_y", &dispvtx_y);
   tree->Branch("dispvertex_z", &dispvtx_z);
   tree->Branch("dispvertex_ex", &dispvtx_ex);
   tree->Branch("dispvertex_ey", &dispvtx_ey);
   tree->Branch("dispvertex_ez", &dispvtx_ez);
   

   tree->Branch("privertex_num", &privertex_num, "privertex_num/I");
   tree->Branch("privertex_x", &privtx_x);
   tree->Branch("privertex_y", &privtx_y);
   tree->Branch("privertex_z", &privtx_z);
   tree->Branch("privertex_ex", &privtx_ex);
   tree->Branch("privertex_ey", &privtx_ey);
   tree->Branch("privertex_ez", &privtx_ez);


   tree->Branch("MET_pt",  &MET_pt);
   tree->Branch("MET_phi", &MET_phi);

   tree->Branch("hltbitmap", &hltbitmap);
   //   tree->Branch("hltprescalemap", &hltprescalemap);

   tree->Branch("l1bitmap", &l1bitmap);
   tree->Branch("l1prescalemap", &l1prescalemap);

   //   tree->Branch("triggerPass", &triggerfilter);
   //   tree->Branch("prescaleL1", &prescaleL1);

}


ScoutingNtuplizer::~ScoutingNtuplizer() {
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)
    file->cd();
    tree->Write();
    file->Close();
}


//
// member functions
//

//we need to initalise the menu each run (menu can and will change on run boundaries)                                                
void ScoutingNtuplizer::beginRun(const edm::Run& run,const edm::EventSetup& setup)
{
  bool changed=false;
  hltPSProv_.init(run,setup,hltProcess_,changed);
  const l1t::L1TGlobalUtil& l1GtUtils = hltPSProv_.l1tGlobalUtil();
  std::cout <<"l1 menu "<<l1GtUtils.gtTriggerMenuName()<<" version "<<l1GtUtils.gtTriggerMenuVersion()<<" comment "<<std::endl;
  std::cout <<"hlt name "<<hltPSProv_.hltConfigProvider().tableName()<<std::endl;
}


// ------------ method called for each event  ------------
void ScoutingNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {


   using namespace edm;
   int getCollectionsResult = GetCollections(iEvent);
    if (getCollectionsResult)
	return;
    
    if (muons->size() < (unsigned int)min_muons)
      return;
    
    	ResetVariables();

    // int HLTpass = 0;



    /*
    const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
    //std::cout << "\n == TRIGGER PATHS= " << std::endl;
    for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
      //std::cout << "Trigger " << names.triggerName(i) <<                                                                         
      // ", prescale " << triggerPrescales->getPrescaleForIndex(i) <<                                                               	  //": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)")                                                             
      if (triggerBits->accept(i)){
	//std::cout << "HLT Triggers passed" << names.triggerName(i) << std::endl;
      
      //string HLTstring = (names.triggerName(i)).substr(0, 34);
	if( (names.triggerName(i).find("DST_DoubleMu3_noVtx_CaloScouting_v") != std::string::npos)){
	  HLTpass = 1;
	}
      }
    }
 
    
    //std::cout<<"HLTpass"<<HLTpass<<std::endl;

    if (HLTpass != 1)
      return;      
    
    //std::cout<<"HLTpass"<<HLTpass<<std::endl;
    */

	
    rho = *handle_rho;
    run = iEvent.id().run();
    lumi = iEvent.id().luminosityBlock();
    event = iEvent.id().event();

    //    PrescaleProvider psProv("../hltJsons/triggerData2017");

    iEvent.getByToken(triggerBits_, triggerBits);

    const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits); 

    for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {                                                          
      const std::string& hltbitName = names.triggerName(i);
      bool hltpassFinal = triggerBits->accept(i);

      for(size_t i = 0; i < hltseedsvector.size(); i++){
	std::string hltName = hltseedsvector[i];
	std::string hltpathName = hltbitName;
        if(hltbitName.find(hltName) != std::string::npos){
          hltbitmap.push_back(std::make_pair(hltseedsvector[i],hltpassFinal));
	  //	  hltprescalemap.push_back(std::make_pair(hltseedsvector[i],psProv.hltPrescale(names.triggerName(i),run,lumi)));
        }
      }

      if (triggerBits->accept(i)){
	std::cout << "HLT Triggers passed" << names.triggerName(i) << std::endl;
      }

    }




    //Beamspot
    
    iEvent.getByToken(beamSpotToken, beamSpotH);
    BS_x = beamSpotH->position().x();
    BS_y = beamSpotH->position().y();
    BS_z = beamSpotH->position().z();
 
    
    //Jets

    for (auto &j: *jets) {
		jet_pt.push_back(j.pt());
		jet_eta.push_back(j.eta());
		jet_phi.push_back(j.phi());
		jet_m.push_back(j.m());
		
		jet_Area.push_back(j.jetArea());
		jet_maxEInEmTowers.push_back(j.maxEInEmTowers());
		jet_maxEInHadTowers.push_back(j.maxEInHadTowers());
		jet_hadEnergyInHB.push_back(j.hadEnergyInHB());
		jet_hadEnergyInHE.push_back(j.hadEnergyInHE());
		jet_hadEnergyInHF.push_back(j.hadEnergyInHF());
		jet_emEnergyInEB.push_back(j.emEnergyInEB());
		jet_emEnergyInEE.push_back(j.emEnergyInEE());
		jet_emEnergyInHF.push_back(j.emEnergyInHF());
		jet_towersArea.push_back(j.towersArea());
		jet_mvaDiscriminator.push_back(j.mvaDiscriminator());
		jet_btagDiscriminator.push_back(j.btagDiscriminator());

		
		HT += j.pt();
		jet_num += 1;
	}


    int length = 0;
    std::vector<int> Detach(int len,std::vector<int> vertex);


    //Muons
    for (auto &m: *muons) {
	   muon_q.push_back(m.charge());
	   muon_pt.push_back(m.pt());
	   muon_eta.push_back(m.eta()); 
	   muon_phi.push_back(m.phi()); 
	   muon_isGlobalMuon.push_back(m.isGlobalMuon());
	   muon_isTrackerMuon.push_back(m.isTrackerMuon());
	   muon_dxy.push_back(m.dxy());
	   muon_edxy.push_back(m.dxyError());
	   muon_dz.push_back(m.dz());
	   muon_edz.push_back(m.dzError());
	   muon_chisquared.push_back(m.chi2());
	   muon_ndof.push_back(m.ndof());
	   muon_trackIso.push_back(m.trackIso());
	   muon_nmatchedstations.push_back(m.nMatchedStations());
	   muon_ntrklayersmeasurement.push_back(m.nTrackerLayersWithMeasurement());
	   muon_nmuonhits.push_back(m.nValidMuonHits());
	   muon_npixelhits.push_back(m.nValidPixelHits());
	   muon_nstriphits.push_back(m.nValidStripHits());
	   //muon_vtxindex.push_back(m.vtxIndx());
	   if (muon_num == 0){
	     muon_vtxindex.push_back(m.vtxIndx());
	     length = 0;
           }else{
	     muon_vtxindex.push_back(Detach(length, m.vtxIndx()));

           }length = (m.vtxIndx().size());
	   muon_num += 1;
	}

    muon_num = muons->size();
    
    //MET
    MET_pt = *handle_MET_pt;
    MET_phi = *handle_MET_phi;
    
    //Vertices
    dispvertex_num = dispvertices->size();
    
    for (auto &v: *dispvertices) {
		dispvtx_x.push_back(v.x());
		dispvtx_y.push_back(v.y());
		dispvtx_z.push_back(v.z());
		dispvtx_ex.push_back(v.xError());
		dispvtx_ey.push_back(v.yError());
		dispvtx_ez.push_back(v.zError());
		dispvtx_chisquared.push_back(v.chi2());
		dispvtx_ndof.push_back(v.ndof());
		dispvtx_tracksize.push_back(v.tracksSize());

	}

    privertex_num = privertices->size();

    for (auto &pv: *privertices) {
                privtx_x.push_back(pv.x());
                privtx_y.push_back(pv.y());
                privtx_z.push_back(pv.z());
		privtx_ex.push_back(pv.xError());
		privtx_ey.push_back(pv.yError());
		privtx_ez.push_back(pv.zError());

    }


    /*
    std::vector<int> triggerEval = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    std::vector<string> triggers{"L1_DoubleMu_12_5",
        "L1_DoubleMu_12_8",
        "L1_DoubleMu_13_6",
        "L1_DoubleMu_15_5",
        "L1_DoubleMu_15_7",
        "L1_DoubleMu18er2p1",
        "L1_DoubleMu22er2p1",
        "L1_TripleMu_4_4_4",
        "L1_TripleMu_5_0_0",
        "L1_TripleMu_5_3_3",
        "L1_TripleMu_5_5_3",
        "L1_QuadMu0",
        "L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4",
        "L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18",
        "L1_DoubleMu4_SQ_OS_dR_Max1p2",
        "L1_DoubleMu5_SQ_OS_Mass7to18",
        "L1_DoubleMu_20_2_SQ_Mass_Max20",
        "L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4",
        "L1_DoubleMu4p5_SQ_OS_dR_Max1p2",
        "L1_DoubleMu6_SQ_OS",
	"L1_DoubleMu0er1p5_SQ_dR_Max1p4",
	"L1_DoubleMu0er2_SQ_dR_Max1p4",
	"L1_DoubleMu0_SQ"};

    std::vector<int> prescaler = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    */


    //    for (size_t i = 0; i < l1seedsvector.size(); i++) {
    //      l1map.push_back(std::make_pair(l1seedsvector[i],-1));
    //    }






      ///////////////////////////////////////////////////////

  //I seem to recall this function being slow so perhaps cache for a given lumi                                                     
  //(it only changes on lumi boundaries)                                                                                            
    int psColumn = hltPSProv_.prescaleSet(iEvent,iSetup);
    //  std::cout <<"PS column "<<psColumn<<std::endl;
  if(psColumn==0 && iEvent.isRealData()){
    std::cout <<"PS column zero detected for data, this is unlikely (almost all triggers are disabled in normal menus here) and its more likely that you've not loaded the correct global tag in "<<std::endl;
  }

  //note to the reader, what I'm doing is extremely dangerious (a const cast), never do this!                                       
  //however in this narrow case, it fixes a bug in l1t::L1TGlobalUtil (the method should be const)                                  
  //and it is safe for this specific instance                                                                                       
  //l1t::L1TGlobalUtil& l1GtUtils = const<l1t::L1TGlobalUtil&> (hltPSProv_.l1tGlobalUtil());                                        

  l1t::L1TGlobalUtil& l1GtUtils = const_cast<l1t::L1TGlobalUtil&> (hltPSProv_.l1tGlobalUtil());
  //l1t::L1TGlobalUtil& l1GtUtils->retrieveL1(iEvent, iSetup, l1AlgoToken);                                                        
  //  std::cout <<"l1 menu: name decisions prescale "<<std::endl;
  for(size_t bitNr=0;bitNr<l1GtUtils.decisionsFinal().size();bitNr++){

    const std::string& bitName = l1GtUtils.decisionsFinal()[bitNr].first; 
// l1GtUtils.decisionsFinal() is of type std::vector<std::pair<std::string,bool> >
                                                          
//    bool passInitial = l1GtUtils.decisionsInitial()[bitNr].second; //before masks and prescales, so if we have a 15 GeV electron passing L1_SingleEG10, it will show up as true but will likely not cause a L1 acccept due to the seeds high prescale                   
//    bool passInterm = l1GtUtils.decisionsInterm()[bitNr].second; //after mask (?, unsure what this is)                              
    bool passFinal = l1GtUtils.decisionsFinal()[bitNr].second; //after masks & prescales, true means it gives a L1 accept to the HLT
    int prescale = l1GtUtils.prescales()[bitNr].second;    

    //    if (passFinal != 0){
    //    std::cout <<"   "<<bitNr<<" "<<bitName<<" "<<passInitial<<" "<<passInterm<<" "<<passFinal<<" "<<prescale<<std::endl;       
    //    }
    for(size_t i = 0; i < l1seedsvector.size(); i++){
      //TPRegexp pattern(l1seedsvector[i]);
      std::string l1Name = l1seedsvector[i];
      std::string pathName = bitName;
      if(bitName.compare(l1Name) == 0){
	l1bitmap.push_back(std::make_pair(l1seedsvector[i],passFinal));
	l1prescalemap.push_back(std::make_pair(l1seedsvector[i],prescale));
      }
    }
    



 
    /*
    
    for(int i=0; (unsigned int)i<triggerEval.size();i++){
      string trig = triggers[i];
      if(trig.compare(bitName) == 0){
	//	std::cout << prescale;
	//	prescaler[i] = prescale;
	//	std::cout << prescaler[i];
	if(passFinal != 0){
	  triggerEval[i] = 1;
	}
	triggerfilter.push_back(triggerEval[i]);
	std::cout << triggerEval[i];
	std::cout <<"   "<<bitNr<<" "<<bitName<<" "<<passInitial<<" "<<passInterm<<" "<<passFinal<<" "<<prescale<<std::endl;
      }
      //      prescaleL1.push_back(prescaler[i]);
    }

    */




  }

  tree->Fill();


}

// ------------ method called to detach vertex vectors from previous ones  ------------                                              
std::vector<int> Detach(int len,std::vector<int> vertex){
  if (vertex.size() != 0 and len != 0 and vertex.size() >= (unsigned)len){
    std::vector<int> newvertex;
    for(int j=len; (unsigned int)j<vertex.size(); j++){
      newvertex.push_back(vertex[j]);
    }return newvertex;
  }else{
    return vertex;
  }
}



// ------------ method called once each job just before starting event loop  ------------
void ScoutingNtuplizer::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void ScoutingNtuplizer::endJob() {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ScoutingNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

int ScoutingNtuplizer::GetCollections(const edm::Event& iEvent) {
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
	
	// Get muons
    iEvent.getByToken(token_muons, muons);
    if (!muons.isValid()) {
        throw edm::Exception(edm::errors::ProductNotFound)
            << "Could not find ScoutingMuonCollection." << endl;
        return 1;
    }
    
    // Get MET
    iEvent.getByToken(token_MET_pt, handle_MET_pt);
    if (!handle_MET_pt.isValid()) {
        throw edm::Exception(edm::errors::ProductNotFound)
	    << "Could not find MET." << endl;
	return 1;
    }

    iEvent.getByToken(token_MET_phi, handle_MET_phi);
    if (!handle_MET_phi.isValid()) {
        throw edm::Exception(edm::errors::ProductNotFound)
	    << "Could not find MET_phi." << endl;
	return 1;
    }
   
    // Get Vertices
    iEvent.getByToken(token_dispvertices, dispvertices);
    if (!dispvertices.isValid()) {
        throw edm::Exception(edm::errors::ProductNotFound)
	    << "Could not find DisplacedScoutingVertexCollection." << endl;
	return 1;
    }
    
    iEvent.getByToken(token_privertices, privertices);
    if (!privertices.isValid()) {
      throw edm::Exception(edm::errors::ProductNotFound)
	<< "Could not find PrimaryScoutingVertexCollection." << endl;
      return 1;
    }
	
	return 0;
}

void ScoutingNtuplizer::ResetVariables() {
    rho = 0.0;
    run = 0;
    lumi = 0;
    event = 0;

    BS_x = 0.0;
    BS_y = 0.0;
    BS_z = 0.0;

    // Reset Jets
    jet_num = 0;
    jet_pt.clear();
    jet_eta.clear();
    jet_phi.clear();
    jet_m.clear();
    
    jet_Area.clear();
    jet_maxEInEmTowers.clear();
    jet_maxEInHadTowers.clear();
    jet_hadEnergyInHB.clear();
    jet_hadEnergyInHE.clear();
    jet_hadEnergyInHF.clear();
    jet_emEnergyInEB.clear();
    jet_emEnergyInEE.clear();
    jet_emEnergyInHF.clear();
    jet_towersArea.clear();
    jet_mvaDiscriminator.clear();
    jet_btagDiscriminator.clear(); 
    
    HT = 0;
    
    // Reset Muons
    muon_num = 0;
    muon_q.clear();
    muon_pt.clear();
    muon_eta.clear();
    muon_phi.clear();
    muon_isGlobalMuon.clear();
    muon_isTrackerMuon.clear();
    muon_dxy.clear();
    muon_edxy.clear();
    muon_dz.clear();
    muon_edz.clear();
    muon_chisquared.clear();
    muon_ndof.clear();
    muon_trackIso.clear();
    muon_nmatchedstations.clear();
    muon_ntrklayersmeasurement.clear();
    muon_nmuonhits.clear();
    muon_npixelhits.clear();
    muon_nstriphits.clear();
    muon_vtxindex.clear();


    
    // Reset MET
    MET_phi = 0;
    MET_pt = 0;
    
    // Reset Vertices
    dispvertex_num = 0;
    dispvtx_x.clear();
    dispvtx_y.clear();
    dispvtx_z.clear();
    
    dispvtx_ex.clear();
    dispvtx_ey.clear();
    dispvtx_ez.clear();

    dispvtx_chisquared.clear();
    dispvtx_ndof.clear();
    dispvtx_tracksize.clear();

    privertex_num = 0;
    privtx_x.clear();
    privtx_y.clear();
    privtx_z.clear();

    privtx_ex.clear();
    privtx_ey.clear();
    privtx_ez.clear();

    // Reset Triggers                                                                                                               
    hltbitmap.clear();
    //    hltprescalemap.clear();
    //    triggerfilter.clear();
    //    prescaleL1.clear();
    l1bitmap.clear();
    l1prescalemap.clear();
}

//define this as a plug-in
DEFINE_FWK_MODULE(ScoutingNtuplizer);
