//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldExample.cxx 279 2011-02-11 18:23:44Z T.J.Adye $
//
// Description:
//      Simple example usage of the RooUnfold package using toy MC.
//
// Authors: Tim Adye <T.J.Adye@rl.ac.uk> and Fergus Wilson <fwilson@slac.stanford.edu>
//
//==============================================================================
// C/C++


void SetHistoUnf(TH1* h, Color_t color);

void InputDataToUnfold(TH1* h) {
  TFile *fInFitMassRes = TFile::Open(Form("source/data/fit/mass/%s.root",namFitMassRes.Data()), "read");
  fInFitMassRes->GetObject("hFitMassRes;1",h);
}

//==============================================================================
// Example Unfolding
//==============================================================================

void RooUnfolding(Int_t nIter = 10){

  TTree *tInput;

  TH1F *hPt2Rec, *hPt2Gen, *hPt2Meas, *hUnfolded, *hMeas, *hPt2GenUnl;

  Int_t nEvents(0), iEv, xcorr, ycorr;
  Double_t tempPt,tempPtGen, tempRelAv, tempAvErr;
  Double_t mtxRes [nIter][binrec_num];
  Double_t unfYield[bingen_num], unfErrSta[bingen_num];
  Double_t genYield[bingen_num], genErrSta[bingen_num];
  Double_t genUnlYield[20]     , genUnlErrSta[20];
  Double_t recYield[bingen_num], recErrSta[bingen_num];
  Double_t pt;

  TMatrixD mtxCorr(bingen_num,bingen_num);
  TMatrixD mtxCov(bingen_num,bingen_num);
  TMatrixD mtxInvCov(bingen_num,bingen_num);
  TMatrixD mtxTemp(bingen_num,bingen_num);
  TMatrixD mtxUnf(bingen_num,bingen_num);


  // open input file and set the tree branches
//  LoadMCAndTree();

  tInput = tRec;// easy switch between trees

  /*
   * Unfolding has 2 components - training and testing
   * training:
   *    - uses MC generated and MC reconstructed
   *    - binning is 2:1 in favor of reconstructed MC
   *    - difference between adjascent MC generated bins should be constant to make regularization of response matrix easier
   *    - used to fill response matrix:
   *                    - RooUnfoldResponse (Binning_of_reconstructed, Binning_of_generated)
   *                    - RooUnfoldResponse.Fill(x_reconstructed, x_generated)
   *                    - RooUnfoldResponse.Miss(x_generated) for case when x_reconstructed doesn't exist for x_generated (particle not catched in a detector)
   * testing:
   *    - uses distribution you want to fill
   *    - can be data OR MC reconstructed events, which were NOT used for training
   *    - MC reconstructed is used, ratio 9:1 (training:testing) is fine
   *
   * results:
   *    - 3 methods can be chosen - Bayesian, SVD, bin-by-bin
   *    Bayesian:
   *        - needs to regularize a priory singular response matrix
   *        - done iteratively
   *        - RooUnfoldBayes (pointer_to_RooUnfoldResponse_object, measured_distribution(what_we_want_to_unfold),number_of_iterations) does the job
   *        - RooUnfoldBayes.Hreco() returns unfolded_measured_distribution
   *        - RooUnfoldBayes.Hmeasured() returns measured_distribution
   */

  // Prepare Output files
  TFile* fUnfOut = new TFile(Form("%s/source/unfolding/%s_%s_UnfRes.root",namPath.Data(),namInMCcoh.Data(),namBinNo.Data()),"recreate");
  TTree* tUnfOut = new TTree("fUnfTree","tree of unfolding results");
  tUnfOut->Branch("fMtxCorr",        &mtxCorr);
  tUnfOut->Branch("fMtxInvCov",      &mtxInvCov);
  tUnfOut->Branch("fMtxCov",         &mtxCov);
  tUnfOut->Branch("fUnfYield",       &unfYield,Form("fUnfYield[%i]/D",bingen_num));
  tUnfOut->Branch("fUnfErrSta",      &unfErrSta,Form("fUnfErrSta[%i]/D",bingen_num));
  tUnfOut->Branch("fGenYield",       &genYield,Form("fGenYield[%i]/D",bingen_num));
  tUnfOut->Branch("fGenErrSta",      &genErrSta,Form("fGenErrSta[%i]/D",bingen_num));
  tUnfOut->Branch("fGenUnlYield",    &genUnlYield,Form("fGenUnlYield[%i]/D",20));
  tUnfOut->Branch("fGenUnlErrSta",   &genUnlErrSta,Form("fGenUnlErrSta[%i]/D",20));
  tUnfOut->Branch("fRecYield",       &recYield,Form("fRecYield[%i]/D",bingen_num));
  tUnfOut->Branch("fRecErrSta",      &recErrSta,Form("fRecErrSta[%i]/D",bingen_num));
  TFile* fGenOut = new TFile(Form("%s/source/MC/%s_GeneratedPt.root",namPath.Data(),namInMCcoh.Data()),"recreate");
  TTree* tGenOut = new TTree("fTree","tree of generated pt");
  tGenOut->Branch("fPt",&pt);

  hPt2Gen     = new TH1F("hPt2Gen ",Form("Generated p^{2}_{t} for #mu %s;p^{2}_{t}^{GEN} [GeV^{2}/c^{2}];entries [-]",namCoh.Data()),bingen_num,bingen_bnd);
  hPt2Rec     = new TH1F("hPt2Rec ",Form("Reconstructed p^{2}_{t} for #mu %;p^{2}_{t}^{REC} [GeV^{2}/c^{2}];entries [-]",namCoh.Data()),bingen_num,bingen_bnd);
  hPt2Meas    = new TH1F("hPt2Meas",Form("Measured p^{2}_{t} for #mu %;p^{2}_{t}^{REC} [GeV^{2}/c^{2}];entries [-]",namCoh.Data()),bingen_num,bingen_bnd);
  hPt2GenUnl  = new TH1F("hPt2GenUnl ",Form("Generated p^{2}_{t} for #mu %;p^{2}_{t}^{GEN} [GeV^{2}/c^{2}];entries [-]",namCoh.Data()),20,0.,0.02);


  std::cout << "==================================== TRAIN ====================================" << std::endl;
  RooUnfoldResponse response (hPt2Rec,hPt2Gen);

  Int_t nCut(0),nRec(0),nGen(0);
    // loop over reconstructed
    tInput = tRec;
    for (Int_t irec(0); irec<0.8*tInput->GetEntries();irec++){
      if (irec%100000==0) std::cout<<irec<<std::endl;
      tInput->GetEntry(irec);
      if (fPt*fPt>binrec_bnd[binrec_num] || 3.2<fM) continue;   // protect against wrongly reconstructed events
      nCut++;
      if (fPt!=0.) {
        response.Fill (fPt*fPt, fPtGen*fPtGen); nRec++;}
//      else {
//        response.Miss (fPtGen*fPtGen); nGen++; }
    }

    TMatrixD mtx = response.HresponseNoOverflow();
    TPaveText *labCov = new TPaveText(0.35, 0.9, 0.65, 1.0, "brNDC");
    labCov->AddText("Response matrix");
    TCanvas *cNeco = new TCanvas("cNeco","cNeco",1600,900);
    cNeco->cd();
      mtx.Draw("colzTEXT");
      labCov->Draw("same");


    // loop over generated
    tInput = tGen;
    for (Int_t igen(0); igen<tInput->GetEntries();igen++){
      if (igen%100000==0) std::cout<<igen<<std::endl;
      tInput->GetEntry(igen);
      hPt2GenUnl->Fill(fPtGen*fPtGen);
      if (igen<1000000) {pt = fPtGen; tGenOut->Fill();}
    }
  // ---------------------------------------------------------------------------------------------------------------------
  // ------------------------------ PREPARE MC FOR TESTING - FILL HISTOGRAMS ---------------------------------------------
  // ---------------------------------------------------------------------------------------------------------------------
  std::cout << "==================================== TEST =====================================" << std::endl;
  // loop over reconstructed
  tInput = tRec;
  iEv = 0.8*tInput->GetEntries()+1;
  while (nEvents<nTestedEvents){
    if (iEv%100==0) Printf("Unfolding - MC test sample event %i",iEv);
    tInput->GetEntry(iEv);
    iEv++;
    if (fPt*fPt>binrec_bnd[binrec_num] || fM<2.8 || 3.2<fM || fPt==0) continue;   // protect against wrongly reconstructed events
    nEvents++;
    hPt2Gen     ->Fill(fPtGen*fPtGen);  // Inserted generated spectrum (MC truth)
    hPt2Rec     ->Fill(fPt*fPt);        // Inserted reconstructed spectrum (MC smeared - after detector effects)
  }

  // ---------------------------------------------------------------------------------------------------------------------
  // ------------------------------ PREPARE DATA - LOAD YIELD AND CORRECT ------------------------------------------------
  // ---------------------------------------------------------------------------------------------------------------------

  std::cout << "==================================== DATA =====================================" << std::endl;

  InputDataToUnfold(hPt2Meas);

  for (Int_t ibin(1); ibin<=nBins; ibin++){
    hPt2Meas->SetBinContent(ibin,CorrectYield(hPt2Meas->GetBinContent(ibin),ibin-1));
    Printf("Unfolding - DATA yield after corrections: bin %i; content %.2f",ibin,hPt2Meas->GetBinContent(ibin));
  }

  // ---------------------------------------------------------------------------------------------------------------------
  // ------------------------------ PERFORM UNFOLDING ON DATA ------------------------------------------------------------
  // ---------------------------------------------------------------------------------------------------------------------

  std::cout << "==================================== UNFOLD ===================================" << std::endl;
  // *********************************
  // TEXT OUTPUT START
  // *********************************
  std::ofstream textFile, textError, textResults, textStatErr;
  textFile.   open(Form("%s/txt/unfolding/LatexTable_%s_%s.txt",namPath.Data(),namBinNo.Data(),namCoh.Data()));
  textResults.open(Form("%s/txt/unfolding/Results_%s_%s.txt",namPath.Data(),namBinNo.Data(),namCoh.Data()));
  textStatErr.open(Form("%s/txt/unfolding/Errors_%s_%s.txt",namPath.Data(),namBinNo.Data(),namCoh.Data()));
  textError.  open(Form("%s/txt/unfolding/LatexError_%s_%s.txt",namPath.Data(),namBinNo.Data(),namCoh.Data()));
  textError << "Aver. Err [$\\%$]";

  // *********************************
  // LOOP OVER ITERATIONS
  // *********************************
  for (Int_t iter(1);iter<nIter+1;iter++){
    Printf("================================== Iteration %i ===============================\n",iter);
    RooUnfoldBayes   unfold (&response, hPt2Meas, iter);    // OR
//    RooUnfoldSvd     unfold (&response, hPt2Meas, iter);   // OR
//    RooUnfoldTUnfold unfold (&response, hPt2Meas);

    Printf("Regularization parameter: %.5f\n",unfold.GetRegParm());

    TCanvas *cRes = new TCanvas("cRes","cRes",1600,900);
    TCanvas *cCov = new TCanvas("cCov","cCov",1600,900);
    TCanvas *cCor = new TCanvas("cCor","cCor",1600,900);
    TCanvas *cUnf = new TCanvas("cUnf","cUnf",1600,900);

    hUnfolded = (TH1F*) unfold.Hreco(); // Unfolded spectrum (corrected data)
    hMeas = (TH1F*) unfold.Hmeasured(); // Inserted measured spectrum (data)
    mtxUnf = unfold.UnfoldingMatrix();  // Unfolding matrix

    Printf("Covariance matrix:");
    unfold.Ereco(RooUnfold::kCovariance).Print();   // Covariance matrix of unfolded spectrum
    TPaveText *labCov = new TPaveText(0.35, 0.9, 0.65, 1.0, "brNDC");
    labCov->AddText("Covariance matrix");
    cCov->cd();
    unfold.Ereco(RooUnfold::kCovariance).Draw("colzTEXT");   // Covariance matrix of unfolded spectrum
    mtxCov = unfold.Ereco(RooUnfold::kCovariance);
    Printf("Inverted Covariance matrix:");
    unfold.Ereco(RooUnfold::kCovariance).Invert().Print();   // Inverted covariance matrix of unfolded spectrum
    mtxInvCov = unfold.Ereco(RooUnfold::kCovariance).Invert();
    labCov->Draw("SAME");

    xcorr=0; ycorr=0;
    for(Int_t iArray(0);iArray<unfold.Ereco(RooUnfold::kCovariance).GetNoElements();iArray++){
//      std::cout << "GetMatrixArray: " << unfold.Ereco(RooUnfold::kCovariance).GetMatrixArray()[iArray] << std::endl;
      if (iArray%bingen_num==0 && iArray!=0) {xcorr=0;ycorr++;}
      mtxTemp(xcorr,ycorr)=unfold.Ereco(RooUnfold::kCovariance)(xcorr,ycorr)/(TMath::Sqrt(unfold.Ereco(RooUnfold::kCovariance)(xcorr,xcorr))*TMath::Sqrt(unfold.Ereco(RooUnfold::kCovariance)(ycorr,ycorr)));
      xcorr++;
//      std::cout << "Corr mtx: " << mtxCorr.GetMatrixArray()[iArray] << std::endl;
    }
    mtxCorr = mtxTemp;

    Printf("Correlation matrix:");
    mtxCorr.Print();   // Covariance matrix of unfolded spectrum
    TPaveText *labCor = new TPaveText(0.35, 0.9, 0.65, 1.0, "brNDC");
    labCor->AddText("Correlation matrix");
    cCor->cd();
    mtxCorr.Draw("colzTEXT");   // Covariance matrix of unfolded spectrum
    labCor->Draw("SAME");

    gStyle->SetOptStat(0);

    SetHistoUnf(hUnfolded,kGreen+3);
    SetHistoUnf(hMeas,kBlue);
    SetHistoUnf(hPt2Gen,kRed);
    hMeas->Sumw2();

//    unfold.PrintTable (cout, hPt2Gen);
    cRes->cd();
    hUnfolded->SetTitle(Form("Unfolded p_{t}^{2} spectrum of %s: %i iterations",namCoh.Data(),iter));
    hUnfolded->SetMarkerStyle(kCircle);
    hUnfolded->Draw("SAME");
    hMeas->Draw("SAME");
//    hPt2Gen->Draw("SAME");

    TLegend *legRecGen = new TLegend(0.55,0.6,0.95,0.9);
    legRecGen->SetBorderSize(0);
    legRecGen->SetFillColor(0);
    legRecGen->SetFillStyle(0);
    legRecGen->AddEntry(hUnfolded,"Unfolded","l");
    legRecGen->AddEntry(hMeas,"Measured","l");
//    legRecGen->AddEntry(hPt2Gen,"Truth","l");
    legRecGen->Draw("SAME");

    Printf("Unfolding matrix:");
    mtxUnf.Print();
    TPaveText *labUnf = new TPaveText(0.35, 0.9, 0.65, 1.0, "brNDC");
    cUnf->cd();
    labUnf->AddText("Unfolding matrix");
    mtxUnf.Draw("colzTEXT");   // Covariance matrix of unfolded spectrum
    labUnf->Draw("SAME");

    if  (iter==1)      cRes->SaveAs(Form("%s/fig/unfolding/%s_%s.pdf(",namPath.Data(),namOutputFile.Data(),namCoh.Data()));
    else               cRes->SaveAs(Form("%s/fig/unfolding/%s_%s.pdf" ,namPath.Data(),namOutputFile.Data(),namCoh.Data()));
                       cUnf->SaveAs(Form("%s/fig/unfolding/%s_%s.pdf" ,namPath.Data(),namOutputFile.Data(),namCoh.Data()));
    if  (iter!=nIter)  cCor->SaveAs(Form("%s/fig/unfolding/%s_%s.pdf" ,namPath.Data(),namOutputFile.Data(),namCoh.Data()));
    else               cCor->SaveAs(Form("%s/fig/unfolding/%s_%s.pdf)",namPath.Data(),namOutputFile.Data(),namCoh.Data()));


    // PRINT RESULTS FOR FITTING
    for(Int_t ibin(1);ibin<bingen_num+1;ibin++){
      textResults << Form("%.0f",hUnfolded->GetBinContent(ibin)) << std::endl; unfYield [ibin-1]=hUnfolded->GetBinContent(ibin);
      textStatErr << Form("%.0f",hUnfolded->GetBinError(ibin)) << std::endl;   unfErrSta[ibin-1]=hUnfolded->GetBinError(ibin);
      genYield [ibin-1]=hPt2Gen->GetBinContent(ibin);  genErrSta[ibin-1]=hPt2Gen->GetBinError(ibin);
      recYield [ibin-1]=hPt2Meas->GetBinContent(ibin); recErrSta[ibin-1]=hPt2Meas->GetBinError(ibin);
    }
    for(Int_t ibin(1);ibin<21;ibin++){
      genUnlYield [ibin-1]=hPt2GenUnl->GetBinContent(ibin);  genUnlErrSta[ibin-1]=hPt2GenUnl->GetBinError(ibin);
//      Printf("?????????Check out unlimited: %f",genUnlYield[ibin-1]);
    }

    Printf("nCut: %i, nGen: %i, nRec %i, nGen/nRec %.2f",nCut,nGen,nRec,(Double_t)nGen/nRec);

    // PRINT RESULTS TO ROOT FILE FOR FITTING
    tUnfOut->Fill();

    // LATEX TABLE OUTPUT
    tempAvErr = 0.;
  // Print header
    textFile  << "\\begin{table}"<< std::endl;
    textFile  << "\\begin{tabular}{l c c}"<< std::endl;
    // print title
    textFile  << "\\toprule"<< std::endl;
    textFile  <<  iter <<" iterations/bin & \\textbf{$\\Delta$Error (Rec to Unf)} & \\textbf{$\\Delta$Result (Gen to Unf)} \\\\" << std::endl;
    for(Int_t ibin(1);ibin<binrec_num+1;ibin++) {
      //printf("Difference of errors for bin %i is %f\n",ibin,hPt2Meas->GetBinError(ibin)-hUnfolded->GetBinError(ibin));
      // print result to local variable
      mtxRes[iter-1][ibin-1] = hUnfolded->GetBinContent(ibin);
      if (hUnfolded->GetBinContent(ibin)!=0) tempAvErr = tempAvErr + 100*hUnfolded->GetBinError(ibin)/hUnfolded->GetBinContent(ibin);
    // print table containment
    textFile  << "Bin number "<< ibin <<" & "<<  hPt2Meas->GetBinError(ibin)-hUnfolded->GetBinError(ibin)
           <<" & "<< hPt2Gen->GetBinContent(ibin)-hUnfolded->GetBinContent(ibin)  <<" \\\\" << std::endl;
    }// end loop over bins
    // Print bottom
    textFile  << "\\bottomrule"<< std::endl;
    textFile << "Aver. Err [$\\%$] & " << tempAvErr/bingen_num << "&" << "\\\\" << std::endl;
    textFile  << "\\end{tabular}"<< std::endl;
    textFile  << "\\end{table}"<< std::endl;
    textFile  << "  "<< std::endl;

    textError << " & " << tempAvErr/bingen_num;


  }// end iter

  // *********************************
  // AFTERLOOP STUFFS
  // *********************************

  // LATEX OUTPUT END
  textError << "\\\\" << std::endl;
  textFile.close();
  textError.close();
  textResults.close();
  textStatErr.close();

  // print result to local variable and .txt file
  std::ofstream textIter;
  textIter. open(Form("%s/txt/unfolding/iterations_%s_%s.txt",namPath.Data(),namBinNo.Data(),namCoh.Data()));
//  Printf("Relative differences between iterations:");
  for(Int_t iIter(0);iIter<nIter-2;iIter++){
//    Printf("\n After %i iterations:\n",iIter+1);
    tempRelAv = 0.;
    for(Int_t ibin(0);ibin<bingen_num;ibin++) {
//      Printf("Bin %i: %.1f %; ",ibin+1,TMath::Abs(100*(mtxRes[iIter][ibin]-mtxRes[iIter+1][ibin])/mtxRes[iIter+1][ibin]));
//      Printf("Relative uncertainty: %.1f %; ",TMath::Abs(100*(mtxRes[iIter+1][ibin]-mtxRes[iIter+2][ibin])/mtxRes[iIter+1][ibin]));
      tempRelAv = tempRelAv+TMath::Abs(100*(mtxRes[iIter][ibin]-mtxRes[iIter+1][ibin])/mtxRes[iIter+1][ibin]);
    }// end loop over bins
//    Printf("Average: %.1f %",tempRelAv/bingen_num);
    textIter << "No. of iterations: " << iIter+1 << "  Average: " << tempRelAv/bingen_num << "%" << std::endl;
  }// end loop over iterations
  textIter.close();

//  unfold.Ereco(RooUnfold::kCovariance).Print();   // Covariance matrix of unfolded spectrum

  //CLOSE ROOT FILE
  fUnfOut->cd();
  tUnfOut->Write();
  fUnfOut->Close();
  Printf("fUnfOut closed!");
  fGenOut->cd();
  tGenOut->Write();
  fGenOut->Close();
  Printf("fGenOut closed!");

}// end unfolding



