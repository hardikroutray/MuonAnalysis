#!/bin/bash

### NEED to get CMSSW software working with condor
#source /osg/osg3.2/osg-wn-client/setup.sh
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc630
export CMSSW_GIT_REFERENCE=/cvmfs/cms.cern.ch/cmssw.git
cd /cms/routray/muon_ntuples_sl7/CMSSW_10_2_5/src/
eval `scramv1 runtime -sh`
export HOME=/users/h2/routray/
###

cd MuonAnalysis/genproductions/python/ThirteenTeV/B2Phi

echo $(($1+0))

echo "Running on GENSIM"

cmsRun step0_v0.py $(($1+0))

echo "Running on STEP1"

cmsRun step1.py $(($1+0))

echo "Running on STEP2"

cmsRun step2.py $(($1+0))

