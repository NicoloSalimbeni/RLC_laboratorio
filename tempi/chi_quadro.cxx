#include "TFitResult.h"
#include "TMatrixD.h"
#include "fstream"
#include "iostream"
#include "math.h"
#include "vector"
using namespace std;

double f_chi_LC(double* x, double* par)
{
  double R = par[0];
  vector<double> t;
  vector<double> Vt;
  vector<double> Vt_err;
  ifstream tutti_in("dati/dati_chi_quadro.txt");
  vector<double> tutti;
  double k;
  while (tutti_in >> k) {
    tutti.push_back(k);
  }

  for (int i = 0; i < tutti.size(); i = i + 4) {
    t.push_back(tutti[i]);
    Vt.push_back(tutti[i + 1]);
    Vt_err.push_back(tutti[i + 3]);
  }

  //X[0]=C e X[1]=L
  double chi = 0;
  for (int i = 0; i < t.size(); i++) {

    double Ft = ((10 * exp(-1 * (R * t[i]) / (2 * x[1]))) / sqrt((1 / (x[0] * pow(10, -9) * x[1])) - R * R / (2 * x[1] * x[1]))) * (R / (2 * x[1])) * (sin(t[i] * sqrt((1 / (x[0] * pow(10, -9) * x[1])) - R * R / (2 * x[1] * x[1]))) + sqrt((4 * x[1]) / (x[0] * pow(10, -9) * R * R)) * cos(t[i] * sqrt((1 / (x[0] * pow(10, -9) * x[1])) - R * R / (2 * x[1] * x[1]))));
    double x_i = pow(((Vt[i] - Ft) / (Vt_err[i])), 2);
    chi = chi + x_i;
  }
  return chi;
}

double f_chi_RC(double* x, double* par)
{
  double L = par[0];
  vector<double> t;
  vector<double> Vt;
  vector<double> Vt_err;
  ifstream tutti_in("dati/dati_chi_quadro.txt");
  vector<double> tutti;
  double k;
  while (tutti_in >> k) {
    tutti.push_back(k);
  }

  for (int i = 0; i < tutti.size(); i = i + 4) {
    t.push_back(tutti[i]);
    Vt.push_back(tutti[i + 1]);
    Vt_err.push_back(tutti[i + 3]);
  }

  //X[0]=C e X[1]=R
  double chi = 0;
  for (int i = 0; i < t.size(); i++) {

    double Ft = ((10 * exp(-1 * (x[1] * t[i]) / (2 * L))) / sqrt((1 / (x[0] * pow(10, -9) * L)) - x[1] * x[1] / (2 * L * L))) * (x[1] / (2 * L)) * (sin(t[i] * sqrt((1 / (x[0] * pow(10, -9) * L)) - x[1] * x[1] / (2 * L * L))) + sqrt((4 * L) / (x[0] * pow(10, -9) * x[1] * x[1])) * cos(t[i] * sqrt((1 / (x[0] * pow(10, -9) * L)) - x[1] * x[1] / (2 * L * L))));
    double x_i = pow(((Vt[i] - Ft) / (Vt_err[i])), 2);
    chi = chi + x_i;
  }
  return chi;
}

double f_chi_RL(double* x, double* par)
{
  double C = par[0];
  vector<double> t;
  vector<double> Vt;
  vector<double> Vt_err;
  ifstream tutti_in("dati/dati_chi_quadro.txt");
  vector<double> tutti;
  double k;
  while (tutti_in >> k) {
    tutti.push_back(k);
  }

  for (int i = 0; i < tutti.size(); i = i + 4) {
    t.push_back(tutti[i]);
    Vt.push_back(tutti[i + 1]);
    Vt_err.push_back(tutti[i + 3]);
  }

  //X[0]=R e X[1]=L
  double chi = 0;
  for (int i = 0; i < t.size(); i++) {

    double Ft = ((10 * exp(-1 * (x[0] * t[i]) / (2 * x[1]))) / sqrt((1 / (C * pow(10, -9) * x[1])) - x[0] * x[0] / (2 * x[1] * x[1]))) * (x[0] / (2 * x[1])) * (sin(t[i] * sqrt((1 / (C * pow(10, -9) * x[1])) - x[0] * x[0] / (2 * x[1] * x[1]))) + sqrt((4 * x[1]) / (C * pow(10, -9) * x[0] * x[0])) * cos(t[i] * sqrt((1 / (C * pow(10, -9) * x[1])) - x[0] * x[0] / (2 * x[1] * x[1]))));
    double x_i = pow(((Vt[i] - Ft) / (Vt_err[i])), 2);
    chi = chi + x_i;
  }
  return chi;
}

double fit_tot(double* x, double* par)
{
  double R = par[0];
  double C = par[1];
  double L = par[2];

  double Ft = ((10 * exp(-1 * (R * x[0]) / (2 * L))) / sqrt((1 / (C * L)) - R * R / (2 * L * L))) * (R / (2 * L)) * (sin(x[0] * sqrt((1 / (C * L)) - R * R / (2 * L * L))) + sqrt((4 * L) / (C * R * R)) * cos(x[0] * sqrt((1 / (C * L)) - R * R / (2 * L * L))));
  return Ft;
}

void chi_quadro()
{
  TCanvas* c_fit = new TCanvas("fit", "fit", 0, 450, 1000, 600);
  c_fit->SetGrid();
  TF1* f_fit = new TF1("fit", fit_tot, 0.068E-03, 0.11E-03, 3);
  f_fit->SetParNames("R", "C", "L");
  f_fit->SetLineWidth(4);
  f_fit->SetParameter(0, 50);
  f_fit->SetParameter(2, 0.0033);
  f_fit->SetParameter(1, 3.4E-09);
  f_fit->SetNpx(20000);
  TGraphErrors* g = new TGraphErrors("dati/dati_chi_quadro.txt", "%lg %lg %lg %lg");
  g->SetTitle("Fit totale; t (s); V(t) (V)");
  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.5);
  g->GetXaxis()->SetLabelSize(0.047);
  g->GetYaxis()->SetLabelSize(0.047);
  g->GetXaxis()->SetTitleSize(0.052);
  g->GetYaxis()->SetTitleSize(0.052);
  g->GetYaxis()->SetTitleOffset(1);
  g->GetXaxis()->SetMaxDigits(3);
  g->Draw("AP");

  

    TFitResultPtr result_par = g->Fit(f_fit, "R"); 
    TMatrixD cov_par = result_par->GetCovarianceMatrix(); 
    cov_par.Print();
    Double_t cov_LC = cov_par[2][1];
    Double_t cov_RL = cov_par[2][0];
    Double_t cov_RC = cov_par[1][0];

  TLegend* legend_int_lin_T = new TLegend(0.47, 0.7, 0.83, 0.87, NULL, "brNDC");
  legend_int_lin_T->AddEntry(f_fit, "interpolazione", "l");
  legend_int_lin_T->AddEntry(g, "V oscilloscopio", "pe");
  legend_int_lin_T->Draw("SAME");

  //c_fit->SaveAs("plot/fit_totale_chi_quadro.png");
  Double_t R = f_fit->GetParameter(0);
  Double_t C = (f_fit->GetParameter(1)) * pow(10, 9); //è in nF con il 10^9
  Double_t C_out = (f_fit->GetParameter(1));          //è in nF con il 10^9
  Double_t L = f_fit->GetParameter(2);

  Double_t R_err = f_fit->GetParError(0);
  Double_t C_err = (f_fit->GetParError(1));
  Double_t L_err = f_fit->GetParError(2);

  ofstream output("risultati_fit_totale_continua.txt");
  cout << "resistenza R= " << R << "+-" << R_err << endl;
  cout << "capacità C= " << C << "+-" << C_err << endl;
  cout << "induttanza L= " << L << "+-" << L_err << endl;
  cout << "chi quadro chi= " << f_fit->GetChisquare() << " con D.O.F.= " << f_fit->GetNDF() << endl;
  output << "resistenza R= " << R << "+-" << R_err << endl;
  output << "capacità C= " << C << "+-" << C_err << endl;
  output << "induttanza L= " << L << "+-" << L_err << endl;
  output << "chi quadro chi= " << f_fit->GetChisquare() << " con D.O.F.= " << f_fit->GetNDF() << endl;
  
  TCanvas* c_LC = new TCanvas("c1", "c1", 200, 10, 1000, 600);
  c_LC->SetGrid();
  c_LC->SetLogz();
  c_LC->cd();

  TF2* chi_LC = new TF2("f2", f_chi_LC, 2.8, 3.8, 0.0025, 0.004, 1);
  chi_LC->SetParNames("Res");
  chi_LC->SetParameter(0, R);
  chi_LC->SetNpx(200);
  chi_LC->SetNpy(200);
  chi_LC->SetTitle("#chi^{2} vs (C,L) smorzate;C (nF); L (H)");
  gStyle->SetPalette(kRainBow);
  chi_LC->GetXaxis()->SetTitleSize(0.04);
  chi_LC->GetXaxis()->SetLabelSize(0.047);
  chi_LC->GetYaxis()->SetLabelSize(0.047);
  chi_LC->GetXaxis()->SetTitleSize(0.052);
  chi_LC->GetYaxis()->SetTitleSize(0.052);
  chi_LC->GetYaxis()->SetTitleOffset(0.95);
  chi_LC->GetXaxis()->SetMaxDigits(3);
  chi_LC->GetYaxis()->SetMaxDigits(1);
  chi_LC->Draw("COLZ2");
  c_LC->SaveAs("plot/chi_quadro(C,L).png");

  TCanvas* c_RC = new TCanvas("c1", "c1", 200, 10, 1000, 600);
  c_RC->SetGrid();
  c_RC->SetLogz();
  c_RC->cd();
  TF2* chi_RC = new TF2("f2", f_chi_RC, 3.25, 3.575, 80, 145, 1);
  chi_RC->SetParNames("Ind");
  chi_RC->SetParameter(0, L);
  chi_RC->SetNpx(200);
  chi_RC->SetNpy(200);
  chi_RC->SetTitle("#chi^{2} vs (C,R) smorzate;C (nF); R (#Omega)");
  gStyle->SetPalette(kRainBow);
  chi_RC->GetXaxis()->SetTitleSize(0.04);
  chi_RC->GetXaxis()->SetLabelSize(0.047);
  chi_RC->GetYaxis()->SetLabelSize(0.047);
  chi_RC->GetXaxis()->SetTitleSize(0.052);
  chi_RC->GetYaxis()->SetTitleSize(0.052);
  chi_RC->GetYaxis()->SetTitleOffset(0.95);
  chi_RC->GetXaxis()->SetMaxDigits(3);
  chi_RC->GetYaxis()->SetMaxDigits(3);
  chi_RC->Draw("COLZ2");
  c_RC->SaveAs("plot/chi_quadro(C,R).png");

  TCanvas* c_RL = new TCanvas("c2", "c2", 200, 10, 1000, 600);
  c_RL->SetGrid();
  c_RL->SetLogz();
  c_RL->cd();
  TF2* chi_RL = new TF2("f2", f_chi_RL, 100, 130, 0.0031, 0.0034, 1);
  chi_RL->SetParNames("Cap");
  chi_RL->SetParameter(0, C);
  chi_RL->SetNpx(200);
  chi_RL->SetNpy(200);
  chi_RL->SetTitle("#chi^{2} vs (R,L) smorzate;R (#Omega); L (H)");
  gStyle->SetPalette(kRainBow);
  chi_RL->GetXaxis()->SetTitleSize(0.04);
  chi_RL->GetXaxis()->SetLabelSize(0.047);
  chi_RL->GetYaxis()->SetLabelSize(0.047);
  chi_RL->GetXaxis()->SetTitleSize(0.052);
  chi_RL->GetYaxis()->SetTitleSize(0.052);
  chi_RL->GetYaxis()->SetTitleOffset(0.95);
  chi_RL->GetXaxis()->SetMaxDigits(3);
  chi_RL->GetYaxis()->SetMaxDigits(1);
  chi_RL->Draw("COLZ2");
  c_RL->SaveAs("plot/chi_quadro(R,L).png");
}
