/*
./compile runMonoPhodata.C runMonoPhodata
*/




#include <TList.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>
#include "TRint.h"
#include "TROOT.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TMath.h"
#include <iostream>
#include <TDCacheFile.h>
#include <vector>
#include <math.h>
#include <sstream>
#include "makeHistogramsMonoPho_8.C"
using namespace std;
void getHistos(int f){
  TChain *ap = new TChain("myEvent");
if(f==0){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_100_1_zzk.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_101_1_9Hc.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_102_1_ElJ.root");
}else if(f==1){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_103_1_pFa.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_104_1_Qok.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_105_1_53O.root");
}else if(f==2){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_106_1_MZ7.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_107_1_iQM.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_108_1_cyT.root");
}else if(f==3){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_109_1_O8f.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_10_1_99w.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_110_1_hTw.root");
}else if(f==4){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_111_1_mm4.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_112_1_BGB.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_113_1_CJF.root");
}else if(f==5){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_114_1_Yvw.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_115_1_uGO.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_116_1_knY.root");
}else if(f==6){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_117_1_4gy.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_118_1_53G.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_119_1_xSS.root");
}else if(f==7){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_11_1_zze.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_120_1_b1k.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_121_1_FOs.root");
}else if(f==8){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_122_1_Xob.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_123_1_GP5.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_124_1_Dtd.root");
}else if(f==9){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_125_1_2zK.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_126_1_1Un.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_127_1_xCP.root");
}else if(f==10){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_128_1_yFG.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_129_1_oR2.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_12_1_SJb.root");
}else if(f==11){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_130_1_6XC.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_131_1_IJt.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_132_1_bX4.root");
}else if(f==12){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_133_1_NW1.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_134_1_5fQ.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_135_1_Dbd.root");
}else if(f==13){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_136_1_wty.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_137_1_5FK.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_138_1_Aja.root");
}else if(f==14){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_139_1_ysM.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_13_1_lpb.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_140_1_yHC.root");
}else if(f==15){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_141_1_9An.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_142_1_BSh.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_143_1_i2p.root");
}else if(f==16){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_144_1_fXb.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_145_1_XBb.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_146_1_3hb.root");
}else if(f==17){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_147_1_tBC.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_148_1_Ave.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_149_1_CTS.root");
}else if(f==18){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_14_1_Y60.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_150_1_aAv.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_151_1_PVj.root");
}else if(f==19){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_152_1_jx9.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_153_1_NGQ.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_154_1_aDj.root");
}else if(f==20){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_155_1_oZL.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_156_1_Ot4.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_157_1_BRG.root");
}else if(f==21){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_158_1_77O.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_159_1_z4n.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_15_1_0Uz.root");
}else if(f==22){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_160_1_AWC.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_161_1_iKs.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_162_1_Gok.root");
}else if(f==23){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_163_1_eTw.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_164_1_y89.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_165_1_Kub.root");
}else if(f==24){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_166_1_Ibx.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_167_1_mSd.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_168_1_6JJ.root");
}else if(f==25){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_169_1_qfT.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_16_1_4T1.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_170_1_Ve2.root");
}else if(f==26){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_171_1_4QU.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_172_1_oYA.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_173_1_YPT.root");
}else if(f==27){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_174_1_hCl.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_175_1_Gj3.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_176_1_MHp.root");
}else if(f==28){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_177_1_eZ7.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_178_1_Fs5.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_179_1_SwQ.root");
}else if(f==29){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_17_1_VM0.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_180_1_znM.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_181_1_27X.root");
}else if(f==30){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_182_1_fBS.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_183_1_ZFv.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_184_1_wB1.root");
}else if(f==31){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_185_1_F42.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_186_1_Ez6.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_187_1_NpH.root");
}else if(f==32){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_188_1_wuT.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_189_1_d07.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_18_1_XyJ.root");
}else if(f==33){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_190_1_BWq.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_191_1_F8U.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_192_1_pD9.root");
}else if(f==34){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_193_1_mjH.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_194_1_TNr.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_195_1_ajz.root");
}else if(f==35){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_196_1_uct.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_197_1_FTP.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_198_1_oQR.root");
}else if(f==36){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_199_1_a10.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_19_1_NEi.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_1_1_tsV.root");
}else if(f==37){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_200_1_nJY.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_201_1_XGQ.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_202_1_IHQ.root");
}else if(f==38){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_203_1_kqr.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_204_1_TAp.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_205_1_MK6.root");
}else if(f==39){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_206_1_TiC.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_207_1_Wz4.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_208_1_vOB.root");
}else if(f==40){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_209_1_EOu.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_20_1_38W.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_210_1_u3o.root");
}else if(f==41){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_211_1_8za.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_212_1_LHu.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_213_1_z8W.root");
}else if(f==42){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_214_1_nDh.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_215_1_71q.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_216_1_U4i.root");
}else if(f==43){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_217_1_w8w.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_218_1_Mud.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_219_1_MFo.root");
}else if(f==44){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_21_1_lzZ.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_220_1_Hb8.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_221_1_xVZ.root");
}else if(f==45){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_222_1_Wnm.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_223_1_ReU.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_224_1_DZI.root");
}else if(f==46){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_225_1_Hmh.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_226_1_JDz.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_227_1_viQ.root");
}else if(f==47){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_228_1_pfs.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_229_1_b5d.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_22_1_1q1.root");
}else if(f==48){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_230_1_b3D.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_231_1_LjE.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_232_1_cEp.root");
}else if(f==49){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_233_1_XOu.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_234_1_Hj0.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_235_1_0FM.root");
}else if(f==50){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_236_1_Kdw.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_237_1_uYc.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_238_1_uIx.root");
}else if(f==51){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_239_1_rEe.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_23_1_FDJ.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_240_1_qjD.root");
}else if(f==52){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_241_1_KKm.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_242_1_HyH.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_243_1_a83.root");
}else if(f==53){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_244_1_BMc.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_245_1_78R.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_246_1_Eua.root");
}else if(f==54){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_247_1_9MJ.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_248_1_MIz.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_249_1_2zd.root");
}else if(f==55){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_24_1_vca.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_250_1_QDw.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_251_1_oWP.root");
}else if(f==56){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_252_1_3U4.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_253_1_vlp.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_254_1_EQJ.root");
}else if(f==57){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_255_1_eXP.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_256_1_Aae.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_257_1_sM5.root");
}else if(f==58){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_258_1_88l.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_259_1_4lA.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_25_1_U9Z.root");
}else if(f==59){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_260_1_G06.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_261_1_5U6.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_262_1_Amf.root");
}else if(f==60){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_263_1_WGe.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_264_1_CoQ.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_265_1_CWh.root");
}else if(f==61){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_266_1_xEk.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_267_1_SV3.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_268_1_CcR.root");
}else if(f==62){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_269_1_ADD.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_26_1_S3m.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_270_1_xpe.root");
}else if(f==63){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_271_1_M8f.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_272_1_5po.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_273_1_IKT.root");
}else if(f==64){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_274_1_Gz4.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_275_1_0N4.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_276_1_Emd.root");
}else if(f==65){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_277_1_0su.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_278_1_UXU.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_279_1_U1B.root");
}else if(f==66){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_27_1_Lz4.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_280_1_zA8.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_281_1_pht.root");
}else if(f==67){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_282_1_at0.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_283_1_Pc6.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_284_1_t88.root");
}else if(f==68){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_285_1_eGr.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_286_1_ncX.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_287_1_V8X.root");
}else if(f==69){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_288_1_U96.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_289_1_OI2.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_28_1_B9R.root");
}else if(f==70){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_290_1_F1c.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_291_1_40O.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_292_1_RCK.root");
}else if(f==71){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_293_1_EFb.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_294_1_LPE.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_295_1_oZH.root");
}else if(f==72){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_296_1_64A.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_297_1_8iR.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_298_1_kTP.root");
}else if(f==73){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_299_1_7O1.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_29_1_3My.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_2_1_oxN.root");
}else if(f==74){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_300_1_n1g.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_301_1_j8q.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_302_1_TNa.root");
}else if(f==75){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_303_1_ER7.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_304_1_scc.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_305_1_R6M.root");
}else if(f==76){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_306_1_1gc.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_307_1_i86.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_308_1_Yy1.root");
}else if(f==77){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_309_1_vS1.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_30_1_3iN.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_310_1_K4M.root");
}else if(f==78){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_311_1_3Ec.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_312_1_0IX.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_313_1_dci.root");
}else if(f==79){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_314_1_fYb.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_315_1_UP0.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_31_1_J4M.root");
}else if(f==80){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_32_1_7k4.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_33_1_eoU.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_34_1_mQb.root");
}else if(f==81){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_35_1_NGO.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_36_1_sEw.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_37_1_aX8.root");
}else if(f==82){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_38_1_36s.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_39_1_LNQ.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_3_1_RGq.root");
}else if(f==83){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_40_1_6ln.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_41_1_fzD.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_42_1_J7o.root");
}else if(f==84){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_43_1_M7u.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_44_1_JnE.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_45_1_HKp.root");
}else if(f==85){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_46_1_cDa.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_47_1_eoV.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_48_1_luR.root");
}else if(f==86){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_49_1_JTK.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_4_1_SCx.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_50_1_vrY.root");
}else if(f==87){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_51_1_z14.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_52_1_CXd.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_53_1_bmp.root");
}else if(f==88){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_54_1_DVr.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_55_1_FcX.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_56_1_bNq.root");
}else if(f==89){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_57_1_Fw9.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_58_1_yGe.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_59_1_BTA.root");
}else if(f==90){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_5_1_Ot0.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_60_1_YOJ.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_61_1_MkC.root");
}else if(f==91){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_62_1_j16.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_63_1_ydx.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_64_1_239.root");
}else if(f==92){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_65_1_2Lo.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_66_1_Zyu.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_67_1_FNe.root");
}else if(f==93){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_68_1_BKO.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_69_1_ROH.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_6_1_fO8.root");
}else if(f==94){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_70_1_zLH.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_71_1_kSD.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_72_1_liN.root");
}else if(f==95){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_73_1_Cmt.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_74_1_VvN.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_75_1_3sq.root");
}else if(f==96){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_76_1_Zvc.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_77_1_55d.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_78_1_ILq.root");
}else if(f==97){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_79_1_JOk.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_7_1_tsG.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_80_1_YfS.root");
}else if(f==98){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_81_1_tY7.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_82_1_wRi.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_83_1_Dnr.root");
}else if(f==99){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_84_1_tIg.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_85_1_lhY.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_86_1_BYn.root");
}else if(f==100){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_87_1_CKW.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_88_1_vJN.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_89_1_Czl.root");
}else if(f==101){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_8_1_fGD.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_90_1_38D.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_91_1_xHu.root");
}else if(f==102){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_92_1_6UJ.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_93_1_wIq.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_94_1_zxv.root");
}else if(f==103){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_95_1_Bui.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_96_1_FXJ.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_97_1_ZLd.root");
}else if(f==104){
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_98_1_Ldf.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_99_1_h4D.root");
ap->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/cms/WAX/11/store/user/sushil/Summer11/Data/May10thReReco_New/Histo_Data_AOD_9_1_02O.root");
}
cout<<"all loaded"<<endl;

makeHistogramsMonoPho_8 *h = new makeHistogramsMonoPho_8(ap);
stringstream outName_h;
outName_h<<"MonoPhodata2011_"<<f<<".root";
h->Loop('w',outName_h.str().data(),true,true,true,false);

makeHistogramsMonoPho_8 *h_ignoreMu = new makeHistogramsMonoPho_8(ap);
stringstream outName_h_ignoreMu;
outName_h_ignoreMu<<"MonoPhodata2011_ignoreMu_"<<f<<".root";
h_ignoreMu->Loop('w',outName_h_ignoreMu.str().data(),true,false,true,false);

makeHistogramsMonoPho_8 *h_1vertex = new makeHistogramsMonoPho_8(ap);
stringstream outName_h_1vertex;
outName_h_1vertex<<"MonoPhodata2011_1vertex_"<<f<<".root";
h_1vertex->Loop('w',outName_h_1vertex.str().data(),true,true,true,true);

makeHistogramsMonoPho_8 *h_ignoreMu_1vertex = new makeHistogramsMonoPho_8(ap);
stringstream outName_h_ignoreMu_1vertex;
outName_h_ignoreMu_1vertex<<"MonoPhodata2011_ignoreMu_1vertex"<<f<<".root";
h_ignoreMu_1vertex->Loop('w',outName_h_ignoreMu_1vertex.str().data(),true,false,true,true);

}

int main(int argc, char *argv[]){
    int inputs = atoi(argv[argc-1]);
    cout<<"you input "<<inputs<<endl;
    getHistos(inputs);
    return 0;
}

