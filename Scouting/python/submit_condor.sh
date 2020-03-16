#!/bin/bash                                                                                                             
            
### NEED to get CMSSW software working with condor                                                                      
            
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc630
export CMSSW_GIT_REFERENCE=/cvmfs/cms.cern.ch/cmssw.git
cd /cms/routray/muon_ntuples_sl7/CMSSW_10_2_5/src/
eval `scramv1 runtime -sh`
export HOME=/users/h2/routray/
###                                                                                                                     

cd MuonAnalysis/Scouting/python/

echo $1 

#cmsRun ScoutingNtuplizer_condor_cfg.py condInputFolder="/cms/scoutingmuon/hardik/Samples_Production_sl7/CMSSW_9_4_7/src/BtoLLP_m_2_ct_50_6000211_lumi/*" condInputNum=$1 is2017MC=True

cmsRun ScoutingNtuplizer_condor_cfg.py condInputFolder="/cms/scoutingmuon/hardik/Samples_Production_sl7/CMSSW_9_4_0_patch1/src/ggHmumu_AOD_Hmass2_Hctau0_5/*" condInputNum=$1 is2017MC=True

#cmsRun ScoutingNtuplizer_condor_cfg.py condInputtxt="BToPhi_params_mphi2_ctau20mm_RAWSIM_v0.txt" condInputNum=$1 is2017MC=True     
                                                                          
echo condInputNum

#mv scouting_ntuple_$1.root /cms/scoutingmuon/hardik/condor_output/BPhi_ntuples/BPhi_m2_ct50 
mv scouting_ntuple_$1.root /cms/scoutingmuon/hardik/condor_output/ggPhi_ntuples/ggPhi_m2_ct0p5
