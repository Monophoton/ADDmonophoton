
#!/bin/bash                                                                                                                                            
## variables from arguments string in jdl

#                            
#To Actually RUN this:
#condor_submit blah.jdl                                                                                                                                      
# 1: condor cluster number                                                                                                   
# 2: condor process number   
# 3: CMSSW_DIR
# 4: RUN_DIR
# 5: SOURCE_MACRO
# 6: NUM_EVENTS_PER_JOB

CONDOR_CLUSTER=$1
CONDOR_PROCESS=$2
CMSSW_DIR=$3
RUN_DIR=$4
SOURCE_MACRO=$5
NUM_EVENTS_PER_JOB=$6
let CONDOR_PROCESS=$CONDOR_PROCESS+0
#                          
# header                            
#                                                                                                                                  
echo ""
echo "ROOT Macro on Condor"
echo ""

START_TIME=`/bin/date`
echo "started at $START_TIME"

echo ""
echo "parameter set:"
echo "CONDOR_CLUSTER: $CONDOR_CLUSTER"
echo "CONDOR_PROCESS: $CONDOR_PROCESS"
echo "CMSSW_DIR: $CMSSW_DIR"
echo "RUN_DIR: $RUN_DIR"
echo "SOURCE_MACRO: $SOURCE_MACRO"
echo "NUM_EVENTS_PER_JOB: $NUM_EVENTS_PER_JOB"

# setup software environment at FNAL for the given CMSSW release                                                                                       
source /uscmst1/prod/sw/cms/shrc prod
#export SCRAM_ARCH=slc5_ia32_gcc434
export SCRAM_ARCH=slc5_amd64_gcc434
cd $CMSSW_DIR
echo "we setup cmssw here:" `pwd`
eval `scramv1 runtime -sh`

cd $RUN_DIR
#
# modify parameter-set                                                                         
#

#FINAL_SOURCE_MACRO_NAME=`echo batch_${CONDOR_CLUSTER}_${CONDOR_PROCESS}`
#FINAL_SOURCE_MACRO=`echo $FINAL_SOURCE_MACRO_NAME.py`

#FINAL_FILENAME=`echo $FINAL_SOURCE_MACRO_NAME.root`

#INPUTFILE="batch_118359_${CONDOR_PROCESS}_REC.root"

echo $INPUTFILE is the inputfile name
echo ""
#echo "Writing final parameter-set: $FINAL_SOURCE_MACRO to RUN_DIR: $RUN_DIR"
echo ""

#let "skip = $CONDOR_PROCESS * NUM_EVENTS_PER_JOB"

#cat $SOURCE_MACRO | sed -e s/CONDOR_INPUTFILENAME/$INPUTFILE/ | sed -e s/CONDOR_SKIPEVENTS/$skip/ | sed -e s/CONDOR_OUTPUTFILENAME/$FINAL_FILENAME/ > $FINAL_SOURCE_MACRO



#
# run cmssw 
#
cd ${_CONDOR_SCRATCH_DIR}
#cp ${RUN_DIR}/${SOURCE_MACRO} .
echo "ls gives"
echo `ls`

#for i in {1..2}
#do
FINAL_LOG=`echo  batch_${CONDOR_CLUSTER}_${CONDOR_PROCESS}.log`
#FINAL_SOURCE_MACRO=`echo \'${SOURCE_MACRO}\(${i}\)\'`
echo $FINAL_LOG
#echo $FINAL_SOURCE_MACRO
#echo "run: time root -b -l -q " $FINAL_SOURCE_MACRO
#root -b -l -q $FINAL_SOURCE_MACRO >> $FINAL_LOG 2>&1
#echo $i
i=`echo ${CONDOR_PROCESS}`
./$SOURCE_MACRO $i >> $FINAL_LOG 2>&1
exitcode=$?
echo "exitcode=" $exitcode
echo `pwd`
#done
echo "ls gives"
echo `ls`

#cp *.root ${RUN_DIR}/
#                                                                                                                           
# end run                                                                                                                  
#                                                                                                                          

echo ""
END_TIME=`/bin/date`
echo "finished at $END_TIME"
exit $exitcode
