#include "TRExFitter/LimitEvent.h"

LimitEvent::LimitEvent(TTree *tree) : fChain(0) 
{
   Init(tree);
}

LimitEvent::~LimitEvent()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}
Int_t LimitEvent::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

void LimitEvent::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;

   fChain->SetBranchAddress("parameter", &parameter, &b_parameter);
   fChain->SetBranchAddress("CLb_med", &CLb_med, &b_CLb_med);
   fChain->SetBranchAddress("pb_med", &pb_med, &b_pb_med);
   fChain->SetBranchAddress("CLs_med", &CLs_med, &b_CLs_med);
   fChain->SetBranchAddress("CLsplusb_med", &CLsplusb_med, &b_CLsplusb_med);
   fChain->SetBranchAddress("CLb_obs", &CLb_obs, &b_CLb_obs);
   fChain->SetBranchAddress("pb_obs", &pb_obs, &b_pb_obs);
   fChain->SetBranchAddress("CLs_obs", &CLs_obs, &b_CLs_obs);
   fChain->SetBranchAddress("CLsplusb_obs", &CLsplusb_obs, &b_CLsplusb_obs);
   fChain->SetBranchAddress("obs_upperlimit", &obs_upperlimit, &b_obs_upperlimit);
   fChain->SetBranchAddress("inj_upperlimit", &inj_upperlimit, &b_inj_upperlimit);
   fChain->SetBranchAddress("exp_upperlimit", &exp_upperlimit, &b_exp_upperlimit);
   fChain->SetBranchAddress("exp_upperlimit_plus1", &exp_upperlimit_plus1, &b_exp_upperlimit_plus1);
   fChain->SetBranchAddress("exp_upperlimit_plus2", &exp_upperlimit_plus2, &b_exp_upperlimit_plus2);
   fChain->SetBranchAddress("exp_upperlimit_minus1", &exp_upperlimit_minus1, &b_exp_upperlimit_minus1);
   fChain->SetBranchAddress("exp_upperlimit_minus2", &exp_upperlimit_minus2, &b_exp_upperlimit_minus2);
   fChain->SetBranchAddress("fit_status", &fit_status, &b_fit_status);
   fChain->SetBranchAddress("mu_hat_obs", &mu_hat_obs, &b_mu_hat_obs);
   fChain->SetBranchAddress("mu_hat_exp", &mu_hat_exp, &b_mu_hat_exp);
}
