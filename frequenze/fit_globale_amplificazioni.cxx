#include "TAxis.h"
#include "TCanvas.h"
#include "TFitResult.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TH1.h"
#include "TLegend.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TPad.h"
#include "TStyle.h"
#include "fstream"
#include "iostream"
#include "string"
#include "vector"

vector<double> txt_to_vector(string);
void TVector_convert(vector<double>& v, TVectorD& v_root);
double mean(vector<double>& v);
double var(vector<double>& v);

double Trasfer_Function_C(double* x, double* par)
{

  Double_t delta = par[0];
  Double_t Omega = par[1];
  Double_t A = par[2];

  double TC_th = (A * Omega * Omega) / sqrt(pow(2 * 3.141592 * x[0], 4) - 2 * pow(2 * 3.141592 * x[0], 2) * (Omega * Omega - 2 * delta * delta) + pow(Omega, 4));

  return TC_th;
}

double profiling(double* x, double* par)
{

  //Omega x[0]
  //Delta x[1]

  vector<double> freq;
  vector<double> A;
  vector<double> err_A;
  vector<double> tutto = txt_to_vector("dati/dati_amp_C.txt");

  int size_tutto = tutto.size();

  for (int i = 0; i < size_tutto; i = i + 3) {
    freq.push_back(tutto[i]);
    A.push_back(tutto[i + 1]);
    err_A.push_back(tutto[i + 2]);
  }

  double precision = par[0];
  double range_sx = par[1];
  double range_dx = par[2];

  double step = (range_dx - range_sx) / precision;

  double min_chi = 9.99e+99;

  for (int i = 0; i < precision; i++) {
    double k = range_sx + i * step;
    double chi = 0;

    for (int j = 0; j < freq.size(); j++) {
      double TC_th = (k * x[0] * x[0] * 1000000 * 1000000) / sqrt(pow(2 * 3.141 * freq[j], 4) - 2 * pow(2 * 3.141 * freq[j], 2) * (x[0] * x[0] * 1000000 * 1000000 - 2 * x[1] * x[1]) + pow(x[0] * 1000000, 4));
      double chi_j = pow(((TC_th - A[j]) / err_A[j]), 2);
      chi = chi + chi_j;
    }

    if (chi < min_chi) {
      min_chi = chi;
    }
  }

  return min_chi;
}

void fit_globale_amplificazioni()
{
  //plot
  vector<double> freq;
  vector<double> err_freq_fake;
  vector<double> A;
  vector<double> err_A;
  vector<double> tutto = txt_to_vector("dati/dati_amp_C.txt");

  int size_tutto = tutto.size();

  for (int i = 0; i < size_tutto; i = i + 3) {
    freq.push_back(tutto[i]);
    A.push_back(tutto[i + 1]);
    err_A.push_back(tutto[i + 2]);
    err_freq_fake.push_back(tutto[i]*10e-06);
  }

  TVectorD freq_ROOT(size_tutto / 3);
  TVectorD err_freq_fake_ROOT(size_tutto / 3);
  TVectorD A_ROOT(size_tutto / 3);
  TVectorD err_A_ROOT(size_tutto / 3);
  TVector_convert(freq, freq_ROOT);
  TVector_convert(err_freq_fake, err_freq_fake_ROOT);
  TVector_convert(A, A_ROOT);
  TVector_convert(err_A, err_A_ROOT);

  ofstream output("risultati_fit_globale_amplificazioni.txt");
  output << "----------- RISULTATI FIT GLOBALE ---------------" << endl;

  TCanvas* c = new TCanvas("c", "c", 1000, 450, 1300, 650);
  c->SetGrid();
  c->SetBottomMargin(0.12591304);

  TGraphErrors* TC_plot = new TGraphErrors(freq_ROOT, A_ROOT, err_freq_fake_ROOT, err_A_ROOT);
  TC_plot->SetTitle("Funzione di Trasferimento Condensatore (T_{C}); f (Hz); A=V_{out}/V_{in}");
  TC_plot->SetMarkerStyle(20);
  TC_plot->SetMarkerSize(0.5);
  TC_plot->GetXaxis()->SetLabelSize(0.047);
  TC_plot->GetYaxis()->SetLabelSize(0.047);
  TC_plot->GetXaxis()->SetTitleSize(0.052);
  TC_plot->GetYaxis()->SetTitleSize(0.052);
  TC_plot->GetYaxis()->SetTitleOffset(1);
  TC_plot->GetXaxis()->SetMaxDigits(3);
  TC_plot->Draw("AP");

  TF1* Func_Tras = new TF1("Func_Tras", Trasfer_Function_C, 30e+03, 80E+03, 3);
  Func_Tras->SetParNames("Delta", "Omega", "A");
  Func_Tras->SetParameter(0, 1.1e+04);
  Func_Tras->SetParameter(1, 3.0445e+05);
  Func_Tras->SetParameter(2, 1.02676e+00);
  Func_Tras->SetNpx(400);
  Func_Tras->SetLineWidth(1);
  TFitResultPtr result_par_TC = TC_plot->Fit(Func_Tras, "S");

  TLegend* legend_TC = new TLegend(0.7547, 0.76234783, 0.8747, 0.893043, NULL, "brNDC");
  legend_TC->AddEntry(TC_plot, "A_{i}", "e");
  legend_TC->AddEntry(Func_Tras, "T_{C} (Fit)", "l");
  legend_TC->Draw("SAME");

  TMatrixD cov_par_TC = result_par_TC->GetCovarianceMatrix();
  cov_par_TC.Print();
  Double_t cov_OD = cov_par_TC[1][0];
  Double_t chisquare = Func_Tras->GetChisquare();
  int DOF = Func_Tras->GetNDF();

  double omega = Func_Tras->GetParameter(1);
  double err_omega = Func_Tras->GetParError(1);
  double delta = Func_Tras->GetParameter(0);
  double err_delta = Func_Tras->GetParError(0);

  output << "Dal Fit si ottiene: " << endl;
  output << "Delta = " << delta << " +- " << err_delta << endl;
  output << "Omega = " << omega << " +- " << err_omega << endl;
  output << "A = " << Func_Tras->GetParameter(2) << " +- " << Func_Tras->GetParError(2) << endl;
  output << "Chi = " << chisquare << " con " << DOF << " gradi di libertà " << endl;
  output << "------------------ Stime Derivate ----------------"<<endl;

  double risonanza = omega/(2*3.141);
  double errore_risonanza = err_omega/(2*3.141);
  double Q_valore = omega/(2*delta);
  double Q_valore_err = sqrt( pow( ((1/(2*delta))*err_omega) ,2 ) + pow((-omega*err_delta)/(2*delta*delta),2) - 2*(1/(2*delta))*((omega*err_delta)/(2*delta*delta))*cov_OD);

  output << "Risonanza= " << risonanza <<" +- " << errore_risonanza << endl;
  output << "Q-Value = " << Q_valore <<" +- " << Q_valore_err <<endl;
  
  
  c->SaveAs("fit_globale_amplificazioni_C.png");
/*
  TCanvas* c_prof = new TCanvas("c_prof", "c_prof", 1000, 450, 1300, 650);
  c_prof->cd();
  c_prof->SetGrid();
  c_prof->SetBottomMargin(0.12591304);

  TF2* chi_profiling = new TF2("f2", profiling, 3.01e-01, 3.075e-01, 0.95e+04, 1.35e+04, 3);
  chi_profiling->SetTitle("#chi^{2} vs (#Omega,#delta) A profiled; #Omega (Mrad/s);#delta (s^{-1})");
  chi_profiling->SetParameter(0, 300);
  chi_profiling->SetParameter(1, 0);
  chi_profiling->SetParameter(2, 2);
  chi_profiling->SetNpx(200);
  chi_profiling->SetNpy(200);
  gStyle->SetPalette(kRainBow);
  chi_profiling->GetXaxis()->SetLabelSize(0.047);
  chi_profiling->GetYaxis()->SetLabelSize(0.047);
  chi_profiling->GetXaxis()->SetTitleSize(0.052);
  chi_profiling->GetYaxis()->SetTitleSize(0.052);
  chi_profiling->GetYaxis()->SetTitleOffset(0.95);
  chi_profiling->GetXaxis()->SetMaxDigits(3);
  chi_profiling->GetYaxis()->SetMaxDigits(3);
  chi_profiling->Draw("COLZ2");
  c_prof->SaveAs("profiling_fit_globale_TC.png"); //Miniimo_Chi_Quadro_A
  */
}

vector<double> txt_to_vector(string file_name)
{

  fstream input_data(file_name);

  double k;
  vector<double> ausiliario;

  while (input_data >> k) {
    ausiliario.push_back(k);
  }

  return ausiliario;
}

void TVector_convert(vector<double>& v, TVectorD& v_root)
{
  for (int i = 0; i < v.size(); i++) {
    v_root(i) = v[i];
  }
}

double mean(vector<double>& v)
{
  double sum = 0;
  for (auto c : v) {
    sum = sum + c;
  }
  return sum / v.size();
}

double var(vector<double>& v)
{
  double sum = 0;
  for (auto c : v) {
    sum = sum + c;
  }
  double mean = sum / v.size();
  double res = 0;
  for (auto c : v) {
    res = res + pow(c - mean, 2);
  }
  return res / (v.size() - 1);
}
