/// Run this code with root -l -b MuComparisons.cxx

struct Job {std::string name; std::string title; bool multiple=false; std::string mu_extension="";};

float HexToFloat(const std::string& s);
void DrawPlot (const std::vector<Job>& jobs,
               const std::vector<double>& POI,
               const std::vector<double>& Up,
               const std::vector<double>& Down,
               const float& range,
               const std::string& POIname,
               const std::string& xAxisName,
               const std::string& outDir,
               const std::string& ATLASlabel,
               const std::string& otherLabel
               );
bool startsWith (string bigString, string smallString);
bool endsWith(const string &str, const std::string &suffix);

  
void MuComparisons () {
  
  // ADJUST THESE LINES

  //The basic path of all TRExFitter jobs
  const std::string path = "/home/d/dado2/TRExFitter/";

  // Name of the POT
  const std::string POI = "topWidth";

  // Title on X axis
  const std::string xAxisName = "#Delta#Gamma_{t} [GeV]";

  // Path to the output directory
  const std::string outDir = "test";

  // range (difference) to be shown on the plots
  const float range = 1.5;

  // ATLAS label for plots
  const std::string ATLASlabel = "Internal";

  // other label for the plot
  const std::string otherLabel = "#sqrt{#it{s}} = 13 TeV, dilepton";

  //The list of individual fits to be compared
  // Arguments:
  //  1. (name) Folder name relative to "path"
  //  2. (title) Title on the plot
  //  3. (multiple) Boolean to teel the code that there are multiple names with the same subset of strings.
  //        e.g. "topWidth_ejets" and "topWidth_mujets". Default is false.
  //  4. (mu_extension) String that follows the POI name. In the example this would be "ejets". Default is "".
  //  E.g: 
  //    jobs.push_back({.name="topWidth_ejets",  .title="e+jets (2-#mu)", .multiple=true, .mu_extension="ejets"});
  //    jobs.push_back({.name="topWidth_muejts", .title="#mu+jets (2-#mu)", .multiple=true, .mu_extension="mujets"});
  
  std::vector<Job> jobs;
  jobs.push_back({.name="TopWidth_Data_mc16a_emu_Hidden", .title="2015+2016"});
  jobs.push_back({.name="TopWidth_Data_mc16d_emu_Hidden", .title="2017"});
 
  // END OF CONFIGURATION
 
  //To save the individual values (***DO NOT PRINT THEM!!***)
  std::vector<double> POIval;
  std::vector<double> POIval_Up;
  std::vector<double> POIval_Down;

  //Loop over all jobs, read and save the values
  for (const auto& ijob : jobs) {
    const std::string file = path + ijob.name + "/Fits/" + ijob.name + ".txt"; 
    std::ifstream input (file.c_str());

    if (!input) {
      std::cerr << "Cannot open file " << file << std::endl;
      exit(EXIT_FAILURE);
    }

    //Read the file and stop when POI is found
    std::string mu;
    while (input >> mu) {
      if (ijob.multiple) {
	    if (startsWith(mu,POI) && endsWith(mu,ijob.mu_extension)) {
	      break;
        }
      } else {
	    if (startsWith(mu,POI))
	    break;
      }
    }
    
    //Read the value, transform and save it
    input >> mu;
    POIval.push_back(HexToFloat(mu));
    input >> mu;
    POIval_Up.push_back(std::atof(mu.c_str()));
    input >> mu;
    POIval_Down.push_back(std::atof(mu.c_str()));
  }
  
  std::cout << std::endl;
  //Compare all the values
  for (int i_job = 0; i_job < jobs.size(); i_job++) {
    for (int j_job = i_job+1; j_job < jobs.size(); j_job++) {
      std::cout << "Difference " << jobs.at(i_job).title << " and " << jobs.at(j_job).title << std::endl;
      std::cout << POIval.at(i_job)-POIval.at(j_job);
      std::cout << std::endl << std::endl;
    }
  }
  DrawPlot (jobs, POIval, POIval_Up, POIval_Down, range, xAxisName, POI, outDir, ATLASlabel, otherLabel); 
}

// code to translate the hidden format to float
// taken from the main TRExFItter code
float HexToFloat(const std::string& s){
  const std::string first = s.substr(0,s.find('.'));
  const std::string rest = s.substr(s.find('.')+1, s.length());
  const std::string zeros = rest.substr(0,rest.find('.'));
  const std::string second = rest.substr(rest.find('.')+1, rest.length());

  unsigned int i1, i2, n0;

  std::stringstream ss;
  ss << std::hex << first;
  ss >> i1;

  std::stringstream ss1;
  ss1 << std::hex << second;
  ss1 >> i2;

  std::stringstream ss2;
  ss2 << std::hex << zeros;
  ss2 >> n0;


  int signed1 = static_cast<int>(i1);
  // need to subtract the 1234 we added
  signed1-= 1234;
  // need to substract the 5678
  n0-= 5678;

  std::string result = std::to_string(signed1)+".";
  for (int i = 0; i < n0; i++) {
    result += "0";
  }
  result += std::to_string(i2);

  return std::stof(result);
}

// plot the difference
void DrawPlot (const std::vector<Job>& jobs,
               const std::vector<double>& POIval,
               const std::vector<double>& Up,
               const std::vector<double>& Down,
               const float& range,
               const std::string& xAxisName,
               const std::string& POIname,
               const std::string& outDir,
               const std::string& ATLASlabel,
               const std::string& otherLabel
               ) {
  TCanvas c{};
  c.SetGridy();
  c.SetBottomMargin(0.15);

  TLatex lab1;
  lab1.SetTextAlign(13);
  lab1.SetTextFont(72);
  lab1.SetTextSize(1.0*0.045);
  lab1.SetNDC();
  TLatex lab2;
  lab2.SetTextAlign(13);
  lab2.SetTextFont(42);
  lab2.SetTextSize(1.0*0.045);
  lab2.SetNDC();
  TLatex lab3;
  lab3.SetTextAlign(9);
  lab3.SetTextFont(42);
  lab3.SetTextSize(0.8*0.045);
  lab3.SetNDC();
  TLatex lab4;
  lab4.SetTextAlign(9);
  lab4.SetTextFont(42);
  lab4.SetTextSize(0.5*0.045);
  lab4.SetNDC();

  TLegend leg(0.65,0.7,0.9,0.8);
  leg.SetBorderSize(0);
  
  const size_t n = POIval.size();
  Double_t x[n], y[n];
  Double_t hi[n], lo[n];
  Double_t eyhi[n], eylo[n];
  for (int i = 0; i < n; i++) {
    x[i] = POIval.at(i) - POIval.at(0);
    y[i] = n-0.5-i;
    hi[i] = Up.at(i);
    lo[i] = -Down.at(i);
    eyhi[i] = 0;
    eylo[i] = 0;
  }
  TH2F h("","",1,-range,range,n,0,n);
  for (int i = 1; i <=n; i++) {
    h.GetYaxis()->SetBinLabel(i,jobs.at(n-i).title.c_str());
  }
  gStyle->SetOptStat(0);
  double max = -1000;
  double min = 1000;
  for (int i = 0; i < n; i++) {
    if (x[i] > max)
      max = x[i];
    if (x[i] < min)
      min = x[i];
  }

  h.GetYaxis()->SetRangeUser(0,n+3);
  
  TGraphAsymmErrors gr(n,x,y,lo,hi,eylo,eyhi);
  Double_t x_dum[1] = {0};
  Double_t y_dum[1] = {0};
  Double_t xerrhi_dum[1] = {hi[0]};
  Double_t xerrlo_dum[1] = {lo[0]};
  Double_t yerrhi_dum[1] = {100};
  Double_t yerrlo_dum[1] = {100};

  TGraphAsymmErrors gr_dum(1,x_dum,y_dum,xerrlo_dum,xerrhi_dum,yerrlo_dum,yerrhi_dum);
  const double lsize = h.GetYaxis()->GetLabelSize();
  h.GetYaxis()->SetLabelSize(1.5*lsize);
  h.GetXaxis()->SetTitle(xAxisName.c_str());
  const double xsize = h.GetXaxis()->GetTitleSize();
  h.GetXaxis()->SetTitleSize(1.3*xsize);
  gPad->SetLeftMargin(0.2);
 
  h.Draw();

  gr.SetMarkerStyle(8);
  gr_dum.SetFillColorAlpha(kGreen+1,0.5);
  
  gr_dum.Draw("2");
  gr.SetLineWidth(2);
  leg.AddEntry(&gr,"Uncertainty","lp");
  gr.Draw("P");
  
  const double shift_x = 0.49;
  const double shift_y = -0.02;
  lab1.DrawLatex(-0.13+0.3+shift_x, 0.88+shift_y, "ATLAS");
  lab2.DrawLatex(0.3+shift_x, 0.88+shift_y, ATLASlabel.c_str());
  lab3.DrawLatex(-0.13+0.3+shift_x, 0.79+shift_y, otherLabel.c_str());
  
  gSystem->mkdir(outDir.c_str());

  c.Print((outDir+"/"+POIname+".png").c_str());
  c.Print((outDir+"/"+POIname+".pdf").c_str());
  c.Print((outDir+"/"+POIname+".eps").c_str());
}

// helper functions
bool startsWith (string bigString, string smallString) {
  if ( smallString.size() <= bigString.size() && (bigString.compare(0, smallString.size(), smallString) == 0) ) {
    return true;
  }
  return false;
}

bool endsWith(const string &str, const std::string &suffix) {
    return str.size() >= suffix.size() &&
           str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}
