//////////////////////////////////////////////////////////
// Original Author: Tia Miceli
// Wed Aug 18 15:35:09 2010
// Purpose:
// Provide common tools for removing non collision backgrounds.
//////////////////////////////////////////////////////////

#ifndef NonCollisionBG_cxx
#define NonCollisionBG_cxx
#include "NonCollisionBG.h"

ClassImp(NonCollisionBG)
//////////////////////////////////////////////////////////
// 
// Cosmic functions
//
//////////////////////////////////////////////////////////

// Example useage fisherCosmic3(Photon_Roundness[pho],Photon_Angle[pho],Photon_phHasConversionTracks[pho])
float NonCollisionBG::fisherCosmic3(float Photon_Roundness, float Photon_Angle, bool Photon_phHasConversionTracks){

  TH1F *sigRound          = (TH1F*)file_PHOTON->Get("h_Rho_MC");                                                                                                                     
  TH1F *sigAngle          = (TH1F*)file_PHOTON->Get("h_Angle_MC");
  TH1F *sighasconversion  = (TH1F*)file_PHOTON->Get("h_hasconversion_MC");                                                                                                  
  TH1F *bkgRound          = (TH1F*)file_COSMIC->Get("h_Rho_cos");
  TH1F *bkgAngle          = (TH1F*)file_COSMIC->Get("h_Angle_cos");
  TH1F *bkghasconversion  = (TH1F*)file_COSMIC->Get("h_hasconversion_cos");
  
  float sigProb1          = prob(sigRound, Photon_Roundness);
  float bkgProb1          = prob(bkgRound, Photon_Roundness);
  
  float sigProb2          = prob(sigAngle, Photon_Angle);
  float bkgProb2          = prob(bkgAngle, Photon_Angle);
    
  float sigProb3          = prob2(sighasconversion,Photon_phHasConversionTracks);
  float bkgProb3          = prob2(bkghasconversion,Photon_phHasConversionTracks);
  
  //cout<<"(signal like event if prob_Signal ~1 and bkg like event if prob_Signal ~0)"<<endl;
  return sigProb1*sigProb2*sigProb3/(sigProb1*sigProb2*sigProb3 + bkgProb1*bkgProb2*bkgProb3);
    
}

// Example useage fisherCosmic2(Photon_Roundness[pho],Photon_Angle[pho])
float NonCollisionBG::fisherCosmic2(float Photon_Roundness, float Photon_Angle){
//======= calculation fisher-discriminant  ===
  TH1F *sigRound          = (TH1F*)file_PHOTON->Get("h_Rho_MC");
  TH1F *sigAngle          = (TH1F*)file_PHOTON->Get("h_Angle_MC");
  TH1F *bkgRound          = (TH1F*)file_COSMIC->Get("h_Rho_cos");
  TH1F *bkgAngle          = (TH1F*)file_COSMIC->Get("h_Angle_cos");

  float sigProb1          = prob(sigRound, Photon_Roundness);
  float bkgProb1          = prob(bkgRound, Photon_Roundness);
  
  float sigProb2          = prob(sigAngle, Photon_Angle);
  float bkgProb2          = prob(bkgAngle, Photon_Angle);
    
  //probability of finding signal or background using shower shape variables(signal like event if prob_2 ~1 and bkg like event if prob_2 ~0)
  return sigProb1*sigProb2/(sigProb1*sigProb2 + bkgProb1*bkgProb2);
}

float NonCollisionBG::prob(TH1F* h, float value){
  int iBin = h->FindBin(value);
  float norm = h->Integral(1, h->GetNbinsX());
  return float(h->GetBinContent(iBin))/norm;
}

float NonCollisionBG::prob2(TH1F* h, bool Photon_phHasConversionTracks){
   float value=0.;
   if(Photon_phHasConversionTracks){value=1.;}
   else{value=0.;}
  int iBin = h->FindBin(value);
  float norm = h->Integral(1, h->GetNbinsX());
  return float(h->GetBinContent(iBin))/norm;
}

//example useage: float mySigmaRCosmic = bg->sigmaRCosmic(Photon_sc_x[i],  Photon_sc_y[i], Photon_sc_z[i],CosmicMuon_OuterTrack_InnerPoint_x,CosmicMuon_OuterTrack_InnerPoint_y,CosmicMuon_OuterTrack_InnerPoint_z, CosmicMuon_OuterTrack_InnerPoint_px,CosmicMuon_OuterTrack_InnerPoint_py,CosmicMuon_OuterTrack_InnerPoint_pz, CosmicMuon_n)
float NonCollisionBG::sigmaRCosmic(float Pho_sc_x, float Pho_sc_y, float Pho_sc_z,float CosmicMuon_OuterTrack_InnerPoint_x[],float CosmicMuon_OuterTrack_InnerPoint_y[],float CosmicMuon_OuterTrack_InnerPoint_z[], float  CosmicMuon_OuterTrack_InnerPoint_px[],float CosmicMuon_OuterTrack_InnerPoint_py[],float CosmicMuon_OuterTrack_InnerPoint_pz[], int CosmicMuon_n){
    float dr_min_cos = 99999;
    float dx_min_cos = 10000;
    float dy_min_cos = 10000;
    float dz_min_cos = 10000;
    
    vector< pair<bool,TVector3> > solutions;
    
    for(int iMuon=0 ; iMuon !=CosmicMuon_n; iMuon++){   
        bool singlesol = false;
        bool bothsol = false;
        
        float averageSolution_x1 = 0.;
        float averageSolution_y1 = 0.;
        float averageSolution_z1 = 0.;
        float averageSolution_x2 = 0.;
        float averageSolution_y2 = 0.;
        float averageSolution_z2 = 0.;
        TVector3 point(CosmicMuon_OuterTrack_InnerPoint_x[iMuon],CosmicMuon_OuterTrack_InnerPoint_y[iMuon],CosmicMuon_OuterTrack_InnerPoint_z[iMuon]);		
        TVector3 direction(CosmicMuon_OuterTrack_InnerPoint_px[iMuon],CosmicMuon_OuterTrack_InnerPoint_py[iMuon],CosmicMuon_OuterTrack_InnerPoint_pz[iMuon]);
        solutions= getEBIntersections(point, direction);
        
        if(solutions.size()==4){
            
            if(solutions[0].first==true && solutions[1].first==false){
                averageSolution_x1 = solutions[0].second.X();
                averageSolution_y1 = solutions[0].second.Y();
                averageSolution_z1 = solutions[0].second.Z();
            }else if(solutions[0].first==false && solutions[1].first==true){
                averageSolution_x1 = solutions[1].second.X();
                averageSolution_y1 = solutions[1].second.Y();
                averageSolution_z1 = solutions[1].second.Z();
            }else if(solutions[0].first==true && solutions[1].first==true){
                averageSolution_x1 = (solutions[0].second.X() + solutions[1].second.X())/2.;
                averageSolution_y1 = (solutions[0].second.Y() + solutions[1].second.Y())/2.;
                averageSolution_z1 = (solutions[0].second.Z() + solutions[1].second.Z())/2.;
            }
            
            if(solutions[2].first==true && solutions[3].first==false){
                averageSolution_x2 = solutions[2].second.X();
                averageSolution_y2 = solutions[2].second.Y();
                averageSolution_z2 = solutions[2].second.Z();
            }else if(solutions[2].first==false && solutions[3].first==true){
                averageSolution_x2 = solutions[3].second.X();
                averageSolution_y2 = solutions[3].second.Y();
                averageSolution_z2 = solutions[3].second.Z();
            }else if(solutions[2].first==true && solutions[3].first==true){
                averageSolution_x2 = (solutions[2].second.X() + solutions[3].second.X())/2.;
                averageSolution_y2 = (solutions[2].second.Y() + solutions[3].second.Y())/2.;
                averageSolution_z2 = (solutions[2].second.Z() + solutions[3].second.Z())/2.;
            }
            
        }else if(solutions.size()==2){
            
            if(solutions[0].first==true && solutions[1].first==false){
                averageSolution_x1 = solutions[0].second.X();
                averageSolution_y1 = solutions[0].second.Y();
                averageSolution_z1 = solutions[0].second.Z();
            }else if(solutions[0].first==false && solutions[1].first==true){
                averageSolution_x1 = solutions[1].second.X();
                averageSolution_y1 = solutions[1].second.Y();
                averageSolution_z1 = solutions[1].second.Z();
            }else if(solutions[0].first==true && solutions[1].first==true){
                averageSolution_x1 = (solutions[0].second.X() + solutions[1].second.X())/2.;
                averageSolution_y1 = (solutions[0].second.Y() + solutions[1].second.Y())/2.;
                averageSolution_z1 = (solutions[0].second.Z() + solutions[1].second.Z())/2.;
            }
            
        }else if(solutions.size()==1){
            if(solutions[0].first == true){
                averageSolution_x1 = solutions[0].second.X();
                averageSolution_y1 = solutions[0].second.Y();
                averageSolution_z1 = solutions[0].second.Z();
            }
        }
        //========distance of separation between first muon intersection in ECAL and SC position 
        
        //cout << "**dx1** "<< dx2_cos << " dy1: " << dy2_cos << "dz1: "<< dz2_cos << endl;
        
        //single solution means that the muon track only reported a DCA (solutions.size==1)
        // or that the muon grazed the ecal, never entering the tracking region
        // (solutions.size==2) where both of these solutions are intersections with the outer barrel only
        if(solutions.size()==2 || solutions.size()==1) singlesol = true;
        //both solutions means that the muon track passed through the tracking region so that
        // there is an intersection with the outer barrel and inner barrel on it's way in, and
        // there is an intersection with the inner barrel and outer barrel on it's way out hence
        // (solutions.size()==4)
        if(solutions.size()==4) bothsol = true;
        
        //initialize these
        float dr_cos = -1.;
        float dx_cos = 10000.;
        float dy_cos = 10000.;
        float dz_cos = 10000.;
        float dr2_cos = -1.;
        float dx2_cos = 10000.;
        float dy2_cos = 10000.;
        float dz2_cos = 10000.;
        
        if(bothsol == true){
            dx_cos = Pho_sc_x - averageSolution_x1;
            dy_cos = Pho_sc_y - averageSolution_y1;
            dz_cos = Pho_sc_z - averageSolution_z1;
            dx2_cos = Pho_sc_x - averageSolution_x2;
            dy2_cos = Pho_sc_y - averageSolution_y2;
            dz2_cos = Pho_sc_z - averageSolution_z2;
            dr_cos  = sqrt((dx_cos*dx_cos/4.4/4.4) + (dy_cos*dy_cos/4.1/4.1) + (dz_cos*dz_cos/4.6/4.6));
            dr2_cos  = sqrt((dx2_cos*dx2_cos/4.4/4.4) + (dy2_cos*dy2_cos/4.1/4.1) + (dz2_cos*dz2_cos/4.6/4.6));
            if( dr2_cos < dr_min_cos ){
                dr_min_cos = dr2_cos;
                dx_min_cos = dx2_cos;
                dy_min_cos = dy2_cos;
                dz_min_cos = dz2_cos;
            }
            if( dr_cos < dr_min_cos ){
                dr_min_cos = dr_cos;
                dx_min_cos = dx_cos;
                dy_min_cos = dy_cos;
                dz_min_cos = dz_cos;
            }
            
        }
        if (singlesol== true){
            dx_cos = Pho_sc_x - averageSolution_x1;
            dy_cos = Pho_sc_y - averageSolution_y1;
            dz_cos = Pho_sc_z - averageSolution_z1;
            dr_cos  = sqrt((dx_cos*dx_cos/5.2/5.2) + (dy_cos*dy_cos/4.8/4.8) + (dz_cos*dz_cos/9/9));
            
            if( dr_cos < dr_min_cos ){
                dr_min_cos = dr_cos;
                dx_min_cos = dx_cos;
                dy_min_cos = dy_cos;
                dz_min_cos = dz_cos;
            }
        }
        //cout << "mydrcos: "<< dr_min_cos << endl;
    }//loop muon
    
    return dr_min_cos; 
    
}

//example useage: float myDeltaRCosmic = bg->deltaRCosmic(Photon_sc_x[i],  Photon_sc_y[i], Photon_sc_z[i],CosmicMuon_OuterTrack_InnerPoint_x,CosmicMuon_OuterTrack_InnerPoint_y,CosmicMuon_OuterTrack_InnerPoint_z, CosmicMuon_OuterTrack_InnerPoint_px,CosmicMuon_OuterTrack_InnerPoint_py,CosmicMuon_OuterTrack_InnerPoint_pz, CosmicMuon_n)
float NonCollisionBG::deltaRCosmic(float Pho_sc_x, float Pho_sc_y, float Pho_sc_z,float CosmicMuon_OuterTrack_InnerPoint_x[],float CosmicMuon_OuterTrack_InnerPoint_y[],float CosmicMuon_OuterTrack_InnerPoint_z[], float  CosmicMuon_OuterTrack_InnerPoint_px[],float CosmicMuon_OuterTrack_InnerPoint_py[],float CosmicMuon_OuterTrack_InnerPoint_pz[], int CosmicMuon_n){
    float dr_min_cos = 99999;
    float dx_min_cos = 10000;
    float dy_min_cos = 10000;
    float dz_min_cos = 10000;
    
    vector< pair<bool,TVector3> > solutions;
    
    for(int iMuon=0 ; iMuon !=CosmicMuon_n; iMuon++){   
        bool singlesol = false;
        bool bothsol = false;
        
        float averageSolution_x1 = 0.;
        float averageSolution_y1 = 0.;
        float averageSolution_z1 = 0.;
        float averageSolution_x2 = 0.;
        float averageSolution_y2 = 0.;
        float averageSolution_z2 = 0.;
        TVector3 point(CosmicMuon_OuterTrack_InnerPoint_x[iMuon],CosmicMuon_OuterTrack_InnerPoint_y[iMuon],CosmicMuon_OuterTrack_InnerPoint_z[iMuon]);		
        TVector3 direction(CosmicMuon_OuterTrack_InnerPoint_px[iMuon],CosmicMuon_OuterTrack_InnerPoint_py[iMuon],CosmicMuon_OuterTrack_InnerPoint_pz[iMuon]);
        solutions= getEBIntersections(point, direction);
        
        if(solutions.size()==4){
            
            if(solutions[0].first==true && solutions[1].first==false){
                averageSolution_x1 = solutions[0].second.X();
                averageSolution_y1 = solutions[0].second.Y();
                averageSolution_z1 = solutions[0].second.Z();
            }else if(solutions[0].first==false && solutions[1].first==true){
                averageSolution_x1 = solutions[1].second.X();
                averageSolution_y1 = solutions[1].second.Y();
                averageSolution_z1 = solutions[1].second.Z();
            }else if(solutions[0].first==true && solutions[1].first==true){
                averageSolution_x1 = (solutions[0].second.X() + solutions[1].second.X())/2.;
                averageSolution_y1 = (solutions[0].second.Y() + solutions[1].second.Y())/2.;
                averageSolution_z1 = (solutions[0].second.Z() + solutions[1].second.Z())/2.;
            }
            
            if(solutions[2].first==true && solutions[3].first==false){
                averageSolution_x2 = solutions[2].second.X();
                averageSolution_y2 = solutions[2].second.Y();
                averageSolution_z2 = solutions[2].second.Z();
            }else if(solutions[2].first==false && solutions[3].first==true){
                averageSolution_x2 = solutions[3].second.X();
                averageSolution_y2 = solutions[3].second.Y();
                averageSolution_z2 = solutions[3].second.Z();
            }else if(solutions[2].first==true && solutions[3].first==true){
                averageSolution_x2 = (solutions[2].second.X() + solutions[3].second.X())/2.;
                averageSolution_y2 = (solutions[2].second.Y() + solutions[3].second.Y())/2.;
                averageSolution_z2 = (solutions[2].second.Z() + solutions[3].second.Z())/2.;
            }
            
        }else if(solutions.size()==2){
            
            if(solutions[0].first==true && solutions[1].first==false){
                averageSolution_x1 = solutions[0].second.X();
                averageSolution_y1 = solutions[0].second.Y();
                averageSolution_z1 = solutions[0].second.Z();
            }else if(solutions[0].first==false && solutions[1].first==true){
                averageSolution_x1 = solutions[1].second.X();
                averageSolution_y1 = solutions[1].second.Y();
                averageSolution_z1 = solutions[1].second.Z();
            }else if(solutions[0].first==true && solutions[1].first==true){
                averageSolution_x1 = (solutions[0].second.X() + solutions[1].second.X())/2.;
                averageSolution_y1 = (solutions[0].second.Y() + solutions[1].second.Y())/2.;
                averageSolution_z1 = (solutions[0].second.Z() + solutions[1].second.Z())/2.;
            }
            
        }else if(solutions.size()==1){
            if(solutions[0].first == true){
                averageSolution_x1 = solutions[0].second.X();
                averageSolution_y1 = solutions[0].second.Y();
                averageSolution_z1 = solutions[0].second.Z();
            }
        }
        //========distance of separation between first muon intersection in ECAL and SC position 
        
        //cout << "**dx1** "<< dx2_cos << " dy1: " << dy2_cos << "dz1: "<< dz2_cos << endl;
        
        //single solution means that the muon track only reported a DCA (solutions.size==1)
        // or that the muon grazed the ecal, never entering the tracking region
        // (solutions.size==2) where both of these solutions are intersections with the outer barrel only
        if(solutions.size()==2 || solutions.size()==1) singlesol = true;
        //both solutions means that the muon track passed through the tracking region so that
        // there is an intersection with the outer barrel and inner barrel on it's way in, and
        // there is an intersection with the inner barrel and outer barrel on it's way out hence
        // (solutions.size()==4)
        if(solutions.size()==4) bothsol = true;
        
        //initialize these
        float dr_cos = -1.;
        float dx_cos = 10000.;
        float dy_cos = 10000.;
        float dz_cos = 10000.;
        float dr2_cos = -1.;
        float dx2_cos = 10000.;
        float dy2_cos = 10000.;
        float dz2_cos = 10000.;
        
        if(bothsol == true){
            dx_cos = Pho_sc_x - averageSolution_x1;
            dy_cos = Pho_sc_y - averageSolution_y1;
            dz_cos = Pho_sc_z - averageSolution_z1;
            dx2_cos = Pho_sc_x - averageSolution_x2;
            dy2_cos = Pho_sc_y - averageSolution_y2;
            dz2_cos = Pho_sc_z - averageSolution_z2;
            dr_cos  = sqrt((dx_cos*dx_cos) + (dy_cos*dy_cos) + (dz_cos*dz_cos));
            dr2_cos  = sqrt((dx2_cos*dx2_cos) + (dy2_cos*dy2_cos) + (dz2_cos*dz2_cos));
            if( dr2_cos < dr_min_cos ){
                dr_min_cos = dr2_cos;
                dx_min_cos = dx2_cos;
                dy_min_cos = dy2_cos;
                dz_min_cos = dz2_cos;
            }
            if( dr_cos < dr_min_cos ){
                dr_min_cos = dr_cos;
                dx_min_cos = dx_cos;
                dy_min_cos = dy_cos;
                dz_min_cos = dz_cos;
            }
            
        }
        if (singlesol== true){
            dx_cos = Pho_sc_x - averageSolution_x1;
            dy_cos = Pho_sc_y - averageSolution_y1;
            dz_cos = Pho_sc_z - averageSolution_z1;
            dr_cos  = sqrt((dx_cos*dx_cos) + (dy_cos*dy_cos) + (dz_cos*dz_cos));
            
            if( dr_cos < dr_min_cos ){
                dr_min_cos = dr_cos;
                dx_min_cos = dx_cos;
                dy_min_cos = dy_cos;
                dz_min_cos = dz_cos;
            }
        }
        //cout << "mydrcos: "<< dr_min_cos << endl;
    }//loop muon
    
    return dr_min_cos; 
    
}

		
vector<pair <bool,TVector3> > NonCollisionBG::getEBIntersections(const TVector3& X, const TVector3& P){
  
  //Ecal Barrel parameters
  float OuterEBZ   = 291.12;
  float InnerEBZ   = 270.89;
  //note OuterEBRho is defined as a function because it varies with z
  float InnerEBRho = 129.;
  
  double t_closestZ = tDCA(X,P);  
  TVector3 X_closestZ = P*t_closestZ + X;
  double OuterEBRho_closestZ = OuterEBRho( X_closestZ.Z() );
  
  if( fabs(X_closestZ.Perp())>OuterEBRho_closestZ){
    TVector3 PointOnEBtemp(0.,0.,0.);
    pair<bool,TVector3> PointOnEB(false,PointOnEBtemp);
    if( fabs(X_closestZ.Z())< OuterEBZ){
      PointOnEB.first = true;
      double phi = X_closestZ.Phi();
      double x = OuterEBRho_closestZ*sin(phi);
      double y = OuterEBRho_closestZ*cos(phi);
      TVector3 temp(x,y,X_closestZ.Z());
      PointOnEB.second = temp;
    }
    vector<pair<bool,TVector3> > vectorOfEBPoints;
    vectorOfEBPoints.push_back(PointOnEB);
    return vectorOfEBPoints;
    
  }   
  
  
  else if( fabs(InnerEBRho)<fabs(X_closestZ.Perp()) && fabs(X_closestZ.Perp())<fabs(OuterEBRho_closestZ)){
    vector<double> outerTs = outerIntersectionTs(X,P,t_closestZ);
    
    //intersection = X0 + t P_X, where X = <x, y, z>
    TVector3 intersectionOuterEB1 = X + P * outerTs[0];
    TVector3 intersectionOuterEB2 = X + P * outerTs[1];
    //cout << "intersectionOuterEB1.Z()" << intersectionOuterEB1.Z()<< "OuterEBZ"<<OuterEBZ<< endl;
    //cout << "intersectionOuterEB2.Z()" << intersectionOuterEB2.Z()<< "OuterEBZ"<<OuterEBZ<< endl;
    vector<pair<bool,TVector3> > vectorOfEBPoints;
    
    // mytest z of 1st outer solution
    
    if(fabs(intersectionOuterEB1.Z())<OuterEBZ){ 
      vectorOfEBPoints.push_back(make_pair(true,intersectionOuterEB1));
    }else{
      TVector3 zero(0.,0.,0.);
      vectorOfEBPoints.push_back(make_pair(false,zero));
    }
	
    if(fabs(intersectionOuterEB2.Z())<OuterEBZ){
      vectorOfEBPoints.push_back(make_pair(true,intersectionOuterEB2));
    }else{
      TVector3 zero(0.,0.,0.);
      vectorOfEBPoints.push_back(make_pair(false,zero));
    }
    
    return vectorOfEBPoints;
    
  }else{ // rho_closestZ must be inside the Inner radius of the EB
    // we will report 4 intersections 2 with Ecal inner surface,
    // and 2 with the outer surface.
		
    vector<double> innerTs = innerIntersectionTs(129.,X,P);
    vector<double> outerTs = outerIntersectionTs(X,P,t_closestZ);
    
    TVector3 intersectionInnerEB1 = X +P * innerTs[0];
    TVector3 intersectionInnerEB2 = X + P * innerTs[1];
    TVector3 intersectionOuterEB1 = X + P * outerTs[0];
    TVector3 intersectionOuterEB2 = X + P * outerTs[1];
		
    vector< pair<bool,TVector3> > vectorOfEBPoints;
    if(fabs(intersectionOuterEB1.Z())<OuterEBZ){
      vectorOfEBPoints.push_back(make_pair(true,intersectionOuterEB1));
      
    }else{
      TVector3 zero(0.,0.,0.);
      vectorOfEBPoints.push_back(make_pair(false,zero));
    }
    
    if(fabs(intersectionInnerEB1.Z())<InnerEBZ){
      vectorOfEBPoints.push_back(make_pair(true,intersectionInnerEB1));
    }else{
      TVector3 zero(0.,0.,0.);
      vectorOfEBPoints.push_back(make_pair(false,zero));
    }
 
    if(fabs(intersectionInnerEB2.Z())<InnerEBZ){
      vectorOfEBPoints.push_back(make_pair(true,intersectionInnerEB2));
    }else{
      TVector3 zero(0.,0.,0.);
      vectorOfEBPoints.push_back(make_pair(false,zero));
    }
   
    // mytest z of 2nd outer solution
    if(fabs(intersectionOuterEB2.Z())<OuterEBZ){
      vectorOfEBPoints.push_back(make_pair(true,intersectionOuterEB2));
      //cout << "*2nd outer solution" <<intersectionOuterEB2.Z() <<  "OuterEBZ " << OuterEBZ <<endl;
    }else{
      TVector3 zero(0.,0.,0.);
      vectorOfEBPoints.push_back(make_pair(false,zero));
    }
    
    return vectorOfEBPoints;
  }
  
}

double NonCollisionBG::h(double z, double& length){
  // all barrel crystals point toward the tangent of a 6.67 cm circle around the IP
  // so they are semi-projective
  double r = 6.67;
  
  // z is the distance from the IP to the z of the inner surface of the crystal.
  // zprime is the distance from the IP to the point where the tangent line described
  // above crosses the z-axis in terms of z.
  // note: 129 cm is the radius of the inner barrel surface
  double zprime = -2.*r*r*z + sqrt(4*r*r*r*r*z*z + 4*r*r*(z*z+129.*129.)*(129.*129.-r*r));
  zprime = zprime/2./(129.*129. - r*r);
  
  // height is the projection of the crystal length onto the xy-plane.
  // note 23 cm is the length of all of the crystals
  double height = 23.*r/zprime;
  //cout << "height " << height <<endl;
  // length is the projection of the crystal length onto the z-axis
  length = sqrt(23.*23. - height*height);
  
  return height;
  
}

// iterations = 5 by default in the .h file
double NonCollisionBG::OuterEBRho(double zprimeprime, int iterations ){
  // zprimeprime is just a funny name for the z location where you would
  // like to know the radius of the outer envelope of the barrel Ecal.
  // This is computed numerically below.
  double lengthOfTheCrystal = 23.0;

  double minDZ = 1000;
  double bestH = 1000;
  for(int i = 0; i < iterations; ++i) {
    // newZ is a trial Z point we are testing
    double newZ = zprimeprime - lengthOfTheCrystal*(1. - (1./iterations)*i);

    // get length and height of this crystal whose front face is at newZ
    double length = 0;
    double height = h(newZ, length); // length is passed by reference!
    // face z postion, newZ, we must vary newZ until
    // zprimeprime = newZ + length
    double dz = fabs(zprimeprime - (newZ + length));
    if( dz < minDZ ){
      minDZ = dz;
      bestH = height;
      ////////cout << "bestH" << bestH << endl;
      //       //////cout << "z\'\' =  " << zprimeprime 
      //    << "\t new z = " << newZ << "\t length = " << length 
      //    << "\t approx value (z + length) =  " << newZ + length 
      //    << "\t height = " << height
      //    << endl;
    }
  }

  return bestH+129.0; // return the radius of the outer EB surface
}

double NonCollisionBG::tDCA(const TVector3& X, const TVector3& P){
 
  // get position, and momentum components
  double x0 = X.X();
  double y0 = X.Y();
  double px = P.X();
  double py = P.Y();

  // solve for t of closest approach to z-axis:
  double t_closestZ = -(y0*py + x0*px)/(px*px + py*py);
	
  return t_closestZ;
}

vector <double> NonCollisionBG::outerIntersectionTs(const TVector3& X, const TVector3& P, double timeDCA){
  // get limits of t inside which the outer intersection solutions must exist
  double x0 = X.X();
  double y0 = X.Y();
  double z0 = X.Z();
  double px = P.X();
  double py = P.Y();
  double pz = P.Z();
  
  //vector <double> tLimits = innerIntersectionTs(155.,X,P);
  vector <double> tLimits = innerIntersectionTs(155.,X,P);
  // make steps of 1cm
  TVector3 startPoint = tLimits[0]*P + X;
  TVector3 endPoint   = tLimits[1]*P + X;
  
  double pathLength = sqrt( (endPoint.X()-startPoint.X())*(endPoint.X()-startPoint.X())
			    + (endPoint.Y()-startPoint.Y())*(endPoint.Y()-startPoint.Y())
			    + (endPoint.Z()-startPoint.Z())*(endPoint.Z()-startPoint.Z()) );
  
  double timeStep = (tLimits[1]-tLimits[0])/pathLength;
  double testTime = tLimits[0];
  double t1 = testTime;
  double smallestDiff1 = 1000.;
  double t2 = testTime;
  double smallestDiff2 = 1000.;
  
  while(1){
    testTime+=timeStep;
    if(testTime>tLimits[1]) break;
    // look at quadratic form for t described in innerIntersectionTs
    // here R is no longer a constant, but varies with z. We now scan
    // t to make this equation as close to 0 as possible.
    double testZ = pz*testTime + z0;
    //cout << "testZ " << testZ << "z0 "<< z0<< "y0" << y0 << "z0 "<< z0 <<endl;
    double testR = OuterEBRho(testZ);
    double testDiff = fabs(  testTime*testTime * (px*px + py*py)
			     + testTime          * 2.*(px*x0 + py*y0)
			     + x0*x0 + y0*y0 - testR*testR
			     );
    
    if(testTime<timeDCA){
      if(testDiff<smallestDiff1){
	smallestDiff1 = testDiff;
	
	t1 = testTime;
	
      }
    }else{ 
      if(testDiff<smallestDiff2){
	smallestDiff2 = testDiff;
	t2 = testTime;
	
      }
    }//else timeDCA
  }//while
  
  vector<double> t;
  t.push_back(t1);
  t.push_back(t2);
  return t;
}

vector <double> NonCollisionBG::innerIntersectionTs(double R, const TVector3& X, const TVector3& P){  
  double x0 = X.X();
  double y0 = X.Y();
  double px = P.X();
  double py = P.Y();
  
  // this b_over_2 = b/2 from eq 3.
  double b_over_2 = x0 * px + y0 * py ;
  
  // this thingUnderSqrt = (b/2)^2 - a c from eq 3.
  double thingUnderSqrt = b_over_2 * b_over_2 - (px*px + py*py) * (x0*x0 + y0*y0 - R*R) ;
  vector <double> t;
  t.push_back(   (- x0*px - y0*py - sqrt(thingUnderSqrt)) / (px*px + py*py)  );
  t.push_back(   (- x0*px - y0*py + sqrt(thingUnderSqrt)) / (px*px + py*py)  );
  return t;
  
}



//////////////////////////////////////////////////////////
// 
// Halo functions
// 
//////////////////////////////////////////////////////////

// Example useage isHEHalo(Photon_sc_phi[x], HERecHit_subset_n, HERecHit_subset_x, HERecHit_subset_y, HERecHit_subset_energy, HERecHit_subset_time, false)
// useTime is defualt false (specified in .h file)
bool NonCollisionBG::isHEHalo(float photonSCphi, int nAllHERecHits, float HERecHitX[], float HERecHitY[], float HERecHitEnergy[], float HERecHitTime[], bool useTime){

  float HERecHitEnergy_MIN_CUT = 1.;
  float deltaPhi_he_MAX_CUT = 0.2;
  float rho_he_MIN_CUT = 115.;
  float rho_he_MAX_CUT = 130.;
  float HERecHitTime_MAX_CUT = 0.;
    
  for(int hehit = 0; hehit < nAllHERecHits; hehit++){
    if (HERecHitEnergy[hehit]>HERecHitEnergy_MIN_CUT){
      bool deltaPhiHE_isHalo = false;
      bool rhoHE_isHalo = false;
      bool timeHE_isHalo = false;
	    
      float HERecHitPhi = TMath::ATan2(HERecHitY[hehit],HERecHitX[hehit]);
      float deltaPhi_he = absDeltaPhi(fixPhi(photonSCphi),fixPhi(HERecHitPhi));

      float rho_he = sqrt( HERecHitX[hehit]*HERecHitX[hehit] + HERecHitY[hehit]*HERecHitY[hehit]);
	    
      if(deltaPhi_he<deltaPhi_he_MAX_CUT) deltaPhiHE_isHalo = true;
      if(rho_he_MIN_CUT < rho_he && rho_he < rho_he_MAX_CUT) rhoHE_isHalo = true;
      if(HERecHitTime[hehit]<HERecHitTime_MAX_CUT) timeHE_isHalo = true;
	    
      if(useTime){
        if(deltaPhiHE_isHalo && rhoHE_isHalo && timeHE_isHalo) return true;
      }else{ // don't useTime
        if(deltaPhiHE_isHalo && rhoHE_isHalo) return true;
      }
    }// if HE hit meets energy requirements
	
  }// loop over HE hits
  return false;
}

//Example usage isTrackHalo(Photon_sc_phi[x], CosmicMuon_n, CosmicMuon_OuterTrack_InnerPoint_x, CosmicMuon_OuterTrack_InnerPoint_y)
// where Photon_sc_phi [0,2pi]
bool NonCollisionBG::isTrackHalo(float photonSCphi, int nCosMu, float CosTrackX[], float CosTrackY[]){
  float deltaPhi_cosTrack_MAX_CUT = 0.2;
  float rho_cosTrack_MIN_CUT = 115.;
  float rho_cosTrack_MAX_CUT = 170.;
    
  for(int CosTrack = 0; CosTrack < nCosMu; CosTrack++){
    bool deltaPhiCosTrack0_isHalo = false;
    bool rhoCosTrack0_isHalo = false;

    float CosTrackPhi = TMath::ATan2(CosTrackY[CosTrack],CosTrackX[CosTrack]);
    float deltaPhi0 = absDeltaPhi(fixPhi(photonSCphi),fixPhi(CosTrackPhi));

    float rho0 = sqrt(CosTrackX[CosTrack]*CosTrackX[CosTrack] + CosTrackY[CosTrack]*CosTrackY[CosTrack]);

    if(deltaPhi0 <= deltaPhi_cosTrack_MAX_CUT) deltaPhiCosTrack0_isHalo = true;
    if(rho_cosTrack_MIN_CUT < rho0 && rho0 < rho_cosTrack_MAX_CUT) rhoCosTrack0_isHalo = true;
    if( deltaPhiCosTrack0_isHalo && rhoCosTrack0_isHalo) return true;
  }// loop over CosTracks
    
  return false;
}

// Example useage isCSCHalo(Photon_sc_phi[x], CSCseg_n, CSCseg_x, CSCseg_y, CSCseg_time, false)
// useTime is defualt false (specified in .h file)
bool NonCollisionBG::isCSCHalo(float photonSCphi, int nAllCSCSegments, float CSCSegmentX[], float CSCSegmentY[], float CSCSegmentTime[], bool useTime){
  
  float deltaPhi_csc_MAX_CUT = 0.2;
  float CSCSegmentTime_MAX_CUT = 0.;
  
  for(int cscseg = 0; cscseg < nAllCSCSegments; cscseg++){
    bool deltaPhiCSC_isHalo = false;
    //bool rhoCSC_isHalo = false;
    bool timeCSC_isHalo = false;
    
    float CSCSegmentPhi = TMath::ATan2(CSCSegmentY[cscseg],CSCSegmentX[cscseg]);
    float deltaPhi_csc = absDeltaPhi(fixPhi(photonSCphi),fixPhi(CSCSegmentPhi));
    
    float rho_csc = sqrt( CSCSegmentX[cscseg]*CSCSegmentX[cscseg] + CSCSegmentY[cscseg]*CSCSegmentY[cscseg]);
    
    if(deltaPhi_csc<deltaPhi_csc_MAX_CUT) deltaPhiCSC_isHalo = true;
    if(CSCSegmentTime[cscseg]<CSCSegmentTime_MAX_CUT) timeCSC_isHalo = true;
    
    if(useTime){
      if(deltaPhiCSC_isHalo && timeCSC_isHalo) return true;
    }else{ // don't useTime
      if(deltaPhiCSC_isHalo) return true;
    }
  }// loop over CSC segments
  return false;
}

// Example useage deltaPhiCSCHalo(Photon_sc_phi[x], CSCseg_x[x], CSCseg_y[x])
float NonCollisionBG::deltaPhiCSCHalo(float photonSCphi, float CSCSegmentX, float CSCSegmentY){

  float CSCSegmentPhi = TMath::ATan2(CSCSegmentY,CSCSegmentX);
  float deltaPhi_csc = absDeltaPhi(fixPhi(photonSCphi),fixPhi(CSCSegmentPhi));
  //cout<<"CSCSegmentY: "<<CSCSegmentY<<"CSCSegmentX: "<<CSCSegmentX<<"PhotonSCphi:"<< fixPhi(photonSCphi)<<"CSCSegmentPhi: "<<fixPhi(CSCSegmentPhi)<<"deltaPhi_csc: " <<deltaPhi_csc<<endl; 
  return deltaPhi_csc;
}

//before using this, make sure phi1 and phi2 are defined on the same range
float NonCollisionBG::absDeltaPhi(float phi1, float phi2){
  float dPhi = fabs(phi1 - phi2);
  if(dPhi > TMath::Pi()) dPhi = 2.*TMath::Pi() - dPhi;
  return dPhi;
}

float NonCollisionBG::fixPhi(float phi){
  if (phi < 0.) phi+=2.0*TMath::Pi();
  return phi;
}

//////////////////////////////////////////////////////////
// 
// Spike functions
// 
//////////////////////////////////////////////////////////

// Example usage simpleBarrelE2E9(Photon_ncrys[x], thisPho_ietaRH, thisPho_iphiRH, thisPho_eRH)
// arrayLimit is default to 100 because that's what is so for our ntuples (specified in .h file)
float NonCollisionBG::simpleBarrelE2E9(int nRH, vector<int> &ietaRH, vector<int> &iphiRH, vector<float> &eRH, int arrayLimit){
  float e9 = 0;
  float NeighborHiE = -10.;
    
  for(int cry=0;cry<arrayLimit && cry<nRH;cry++){

    //only consider barrel hits because only they have ieta & iphi !
    if(abs(ietaRH[cry])>85) continue;

    //get rid of no ieta=0 problem :P
    int ieta_pho_cry = ietaRH[cry];
    if(ieta_pho_cry<0) ieta_pho_cry++;

    int ieta_pho_cry_MaxE = ietaRH[0];
    if(ieta_pho_cry_MaxE<0) ieta_pho_cry_MaxE++;

    int delIEta = ieta_pho_cry_MaxE - ieta_pho_cry;

    //take care of phi wrapping
    int delIPhi = iphiRH[0] - iphiRH[cry];
    if(delIPhi==-359) delIPhi= 1;
    if(delIPhi== 359) delIPhi=-1;
	
    if(abs(delIPhi)<=1 && abs(delIEta)<=1){
      e9+=eRH[cry];
      if(eRH[cry] > NeighborHiE && cry!=0) //find highest neighbor to seed (that is not the seed itself)
        NeighborHiE = eRH[cry];
    }
  }
  return (NeighborHiE+eRH[0])/e9;
}

// Example usage isE2E9Spike(Photon_ncrys[x], thisPho_ietaRH, thisPho_iphiRH, thisPho_eRH)
// arrayLimit is default to 100 because that's what is so for our ntuples (specified in .h file)
bool NonCollisionBG::isE2E9Spike(int nRH, vector<int> &ietaRH, vector<int> &iphiRH, vector<float> &eRH, int arrayLimit){
  if(simpleBarrelE2E9(nRH,ietaRH,iphiRH,eRH,arrayLimit)>0.95)
    return true;
  else
    return false;
}

// Example usage isTimeSpike(Photon_ncrys[x],thisPho_tRH)
bool NonCollisionBG::isTimeSpike(vector<float> &tRH){
  if(tRH[0]<-3.5)
    return true;
  else
    return false;
}

#endif 
