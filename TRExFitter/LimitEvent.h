#ifndef LimitEvent_h
#define LimitEvent_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class LimitEvent {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         parameter;
   Float_t         CLb_med;
   Float_t         pb_med;
   Float_t         CLs_med;
   Float_t         CLsplusb_med;
   Float_t         CLb_obs;
   Float_t         pb_obs;
   Float_t         CLs_obs;
   Float_t         CLsplusb_obs;
   Float_t         obs_upperlimit;
   Float_t         inj_upperlimit;
   Float_t         exp_upperlimit;
   Float_t         exp_upperlimit_plus1;
   Float_t         exp_upperlimit_plus2;
   Float_t         exp_upperlimit_minus1;
   Float_t         exp_upperlimit_minus2;
   Float_t         fit_status;
   Float_t         mu_hat_obs;
   Float_t         mu_hat_exp;

   // List of branches
   TBranch        *b_parameter;   //!
   TBranch        *b_CLb_med;   //!
   TBranch        *b_pb_med;   //!
   TBranch        *b_CLs_med;   //!
   TBranch        *b_CLsplusb_med;   //!
   TBranch        *b_CLb_obs;   //!
   TBranch        *b_pb_obs;   //!
   TBranch        *b_CLs_obs;   //!
   TBranch        *b_CLsplusb_obs;   //!
   TBranch        *b_obs_upperlimit;   //!
   TBranch        *b_inj_upperlimit;   //!
   TBranch        *b_exp_upperlimit;   //!
   TBranch        *b_exp_upperlimit_plus1;   //!
   TBranch        *b_exp_upperlimit_plus2;   //!
   TBranch        *b_exp_upperlimit_minus1;   //!
   TBranch        *b_exp_upperlimit_minus2;   //!
   TBranch        *b_fit_status;   //!
   TBranch        *b_mu_hat_obs;   //!
   TBranch        *b_mu_hat_exp;   //!
   
   LimitEvent(TTree *tree=0);
   virtual ~LimitEvent();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual void     Init(TTree *tree);
};

#endif
