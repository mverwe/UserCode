#include "TCanvas.h"
#include "TFile.h"
#include "TH3F.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <locale>

#include "/Users/mverweij/wrk/macros/plotUtils.C"
#include "/Users/mverweij/wrk/macros/style.C"

#include "/Users/mverweij/wrk/pPb/PlotMacros/readFile.C"

using namespace std;

void compareTwoJECs(TString strRaw = "/Users/mverweij/wrk/pPb/JEC/calo/SigOnly/ForestJEC/v9/JECOnTheFlyRawOnly/mergeRawOnly/avgResponseRawOnlyVsRefPt/unf/MV_Unf_L2Relative_ak4Calo_RawOnlyVsRefPt.txt", TString strRef = "/Users/mverweij/wrk/pPb/JEC/calo/SigOnly/ForestJEC/v3/merge/avgResponse/unfRefLessRebin/MV_Unf_L2Relative_ak4Calo_RefLessBins.txt", TString strLeg1 = "ConvMult", TString strLeg2 = "Conv") {

  const int nEtaBins = 82;

  std::vector<TF1*> vabscorRaw;
  std::vector<TF1*> vabscorRef;
  std::vector<TF1*> vabscorUnf;
  std::vector<double> vetamin;
  std::vector<double> vetamax;

  //Raw / ConvMult : [0]+[1]/(pow(log10(x),[2])+[3]))*([4]+[5]/(pow(log10(([0]+[1]/(pow(log10(x),[2])+[3]))*x),[6])+[7]))
  // ([0]+[1]/(pow(log10(x),[2])+[3]))*([4]+[5]/(pow(log10(([0]+[1]/(pow(log10(x),[2])+[3]))*x),[6])+[7]))

  //Ref / Conv : [0]+[1]/(pow(log10(([4]+[5]/(pow(log10(x),[6])+[7]))*x),[2])+[3])
  
  for(int etaBin = 0; etaBin<nEtaBins; ++etaBin) {
    std::string lineRaw = readFile(strRaw.Data(),etaBin+2);
    float parRaw[13];
    sscanf(lineRaw.c_str(),"%f %f %f %f %f %f %f %f %f %f %f %f %f",&parRaw[0],&parRaw[1],&parRaw[2],&parRaw[3],&parRaw[4],&parRaw[5],&parRaw[6],&parRaw[7],&parRaw[8],&parRaw[9],&parRaw[10],&parRaw[11],&parRaw[12]);

    std::string lineRef = readFile(strRef.Data(),etaBin+2);
    float parRef[13];
    sscanf(lineRef.c_str(),"%f %f %f %f %f %f %f %f %f %f %f %f %f",&parRef[0],&parRef[1],&parRef[2],&parRef[3],&parRef[4],&parRef[5],&parRef[6],&parRef[7],&parRef[8],&parRef[9],&parRef[10],&parRef[11],&parRef[12]);
    
    TF1 *fCorrV1Raw = new TF1(Form("fCorrV1Raw_EtaBin%d",etaBin),"([0]+[1]/(pow(log10(x),[2])+[3]))*([4]+[5]/(pow(log10(([0]+[1]/(pow(log10(x),[2])+[3]))*x),[6])+[7]))",4.,1000.); //calo 
    fCorrV1Raw->SetParameters(parRaw[5],parRaw[6],parRaw[7],parRaw[8],parRaw[9],parRaw[10],parRaw[11],parRaw[12]);

    TF1 *fCorrV1Ref = new TF1(Form("fCorrV1Ref_EtaBin%d",etaBin),"[0]+[1]/(pow(log10(([4]+[5]/(pow(log10(x),[6])+[7]))*x),[2])+[3])",4.,1000.); //calo 
    // fCorrV1Ref->SetParameters(parRef[5],parRef[6],parRef[7],parRef[8]);
    fCorrV1Ref->SetParameters(parRef[5],parRef[6],parRef[7],parRef[8],parRef[9],parRef[10],parRef[11],parRef[12]);
    
    vetamin.push_back(parRaw[0]);
    vetamax.push_back(parRaw[1]);
    
    TCanvas *c6 = new TCanvas(Form("c6_%d",etaBin),Form("c6_%d: JEC",etaBin),500,400);
    //vcanv.push_back(c6);
    double ymax = 5.;//2.5;//1.3;
    double ymin = 0.;
    TH1F *fr6 = DrawFrame(10.,1000.,ymin,ymax,"#it{p}_{T,raw/ref} (GeV/#it{c})","1/#LT(#it{p}_{T,raw}/#it{p}_{T,gen})#GT");
    fr6->GetXaxis()->SetMoreLogLabels(kTRUE);
    gPad->SetLogx();

    fCorrV1Raw->SetLineColor(4);
    fCorrV1Ref->SetLineColor(kGreen+3);
    
    fCorrV1Raw->Draw("same");
    fCorrV1Ref->Draw("same");

    DrawLatex(0.45,0.88,Form("%.2f<#eta<%.2f",parRaw[0],parRaw[1]),0.05);

    DrawLatex(0.6,0.8,strLeg1.Data(),0.05,fCorrV1Raw->GetLineColor());
    DrawLatex(0.6,0.74,strLeg2.Data(),0.05,fCorrV1Ref->GetLineColor());
        
    c6->SaveAs(Form("compareJECsEtaBin%dUnf.png",etaBin));
    //c6->SaveAs(Form("macros/ak4CaloCorrectionEtaBin%d.C",etaBin));
  }


    
}
