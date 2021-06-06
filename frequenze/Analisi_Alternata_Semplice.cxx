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
#include "math.h"
#include "string"
#include "vector"

vector<double> txt_to_vector(string);
void TVector_convert(vector<double>& v, TVectorD& v_root);
double mean(vector<double>& v);
double var(vector<double>& v);

void Analisi_Alternata_Semplice()
{

  //-------Massimo di T_C--------

  vector<double> freq_tc;
  vector<double> err_freq_tc_fake;
  vector<double> A_tc;
  vector<double> err_A_tc;
  vector<double> tutto_tc = txt_to_vector("dati/dati_amp_C.txt");
  double TC_Max = 0;
  double f_TC_Max = 0;
  int size_tutto_tc = tutto_tc.size();

  for (int i = 0; i < size_tutto_tc; i = i + 3) {
    freq_tc.push_back(tutto_tc[i]);
    A_tc.push_back(tutto_tc[i + 1]);
    if (tutto_tc[i + 1] > TC_Max) {
      TC_Max = tutto_tc[i + 1];
      f_TC_Max = tutto_tc[i];
    }
    err_A_tc.push_back(tutto_tc[i + 2]);
    err_freq_tc_fake.push_back(tutto_tc[i] * 10e-06);
  }

  TVectorD freq_tc_ROOT(size_tutto_tc / 3);
  TVectorD err_freq_tc_fake_ROOT(size_tutto_tc / 3);
  TVectorD A_tc_ROOT(size_tutto_tc / 3);
  TVectorD err_A_tc_ROOT(size_tutto_tc / 3);
  TVector_convert(freq_tc, freq_tc_ROOT);
  TVector_convert(err_freq_tc_fake, err_freq_tc_fake_ROOT);
  TVector_convert(A_tc, A_tc_ROOT);
  TVector_convert(err_A_tc, err_A_tc_ROOT);

  ofstream output("risultati_alternata_semplice.txt");
  output << "----------- Q-VALUE FROM T_C ---------------" << endl;

  TCanvas* c_tc = new TCanvas("c_tc", "c_tc", 1000, 450, 1000, 600);
  c_tc->SetGrid();
  c_tc->SetBottomMargin(0.12591304);

  TGraphErrors* TC_plot = new TGraphErrors(freq_tc_ROOT, A_tc_ROOT, err_freq_tc_fake_ROOT, err_A_tc_ROOT);
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

  TF1* TC_Max_Par = new TF1("TC_Max", "pol2");
  TFitResultPtr result_par_TC = TC_plot->Fit(TC_Max_Par, "RS", "", f_TC_Max - 1000, f_TC_Max + 1000);
  TMatrixD cov_par_TC = result_par_TC->GetCovarianceMatrix();
  cov_par_TC.Print();
  Double_t cov_ab_TC = cov_par_TC[2][1];
  Double_t cov_ac_TC = cov_par_TC[2][0];
  Double_t cov_bc_TC = cov_par_TC[1][0];

  Double_t a_max_TC = TC_Max_Par->GetParameter(2);
  Double_t b_max_TC = TC_Max_Par->GetParameter(1);
  Double_t c_max_TC = TC_Max_Par->GetParameter(0);
  Double_t a_max_TC_err = TC_Max_Par->GetParError(2);
  Double_t b_max_TC_err = TC_Max_Par->GetParError(1);
  Double_t c_max_TC_err = TC_Max_Par->GetParError(0);
  Double_t delta_TC = b_max_TC * b_max_TC - 4 * a_max_TC * c_max_TC;
  Double_t delta_TC_err = sqrt(pow((b_max_TC_err) * (2 * b_max_TC), 2) + pow((a_max_TC_err * 4 * c_max_TC), 2) + pow((c_max_TC_err * 4 * a_max_TC), 2) - 2 * (2 * b_max_TC) * (4 * c_max_TC) * cov_ab_TC - 2 * (2 * b_max_TC) * (4 * a_max_TC) * cov_bc_TC + (2 * 4 * c_max_TC * 4 * a_max_TC * cov_ac_TC));

  Double_t Q_TC = (-delta_TC / (4 * a_max_TC));
  Double_t Q_TC_err = sqrt(pow((delta_TC_err / (4 * a_max_TC)), 2) + pow(((a_max_TC_err * delta_TC) / (4 * a_max_TC * a_max_TC)), 2));

  output << "Q_TC = " << Q_TC << " +- " << Q_TC_err << endl;

  TLegend* legend_TC = new TLegend(0.7947, 0.81234783, 0.8747, 0.893043, NULL, "brNDC");
  legend_TC->AddEntry(TC_plot, "A_{i}", "e");
  legend_TC->Draw("SAME");

  double f0_C_val = (-b_max_TC / (2 * a_max_TC));
  TLine* f0_C = new TLine(f0_C_val, 0.80, f0_C_val, 15.36);
  f0_C->SetLineWidth(2);
  f0_C->SetLineStyle(2);
  f0_C->SetLineColor(kGray + 3);
  f0_C->Draw("Same");

  c_tc->SaveAs("stima_semplice_TC.png");

  //-------Massimo di T_R--------
  output << "----------- T_R ANALYSIS ---------------" << endl;
  vector<double> freq_tr;
  vector<double> err_freq_tr_fake;
  vector<double> A_tr;
  vector<double> err_A_tr;
  vector<double> tutto_tr = txt_to_vector("dati/dati_amp_R.txt");
  double TR_Max = 0;
  double f_TR_Max = 0;
  int size_tutto_tr = tutto_tr.size();

  for (int i = 0; i < size_tutto_tr; i = i + 3) {
    freq_tr.push_back(tutto_tr[i]);
    A_tr.push_back(tutto_tr[i + 1]);
    if (tutto_tr[i + 1] > TR_Max) {
      TR_Max = tutto_tr[i + 1];
      f_TR_Max = tutto_tr[i];
    }
    err_A_tr.push_back(tutto_tr[i + 2]);
    err_freq_tr_fake.push_back(tutto_tr[i] * 10e-09);
  }

  TVectorD freq_tr_ROOT(size_tutto_tr / 3);
  TVectorD err_freq_tr_fake_ROOT(size_tutto_tr / 3);
  TVectorD A_tr_ROOT(size_tutto_tr / 3);
  TVectorD err_A_tr_ROOT(size_tutto_tr / 3);
  TVector_convert(freq_tr, freq_tr_ROOT);
  TVector_convert(err_freq_tr_fake, err_freq_tr_fake_ROOT);
  TVector_convert(A_tr, A_tr_ROOT);
  TVector_convert(err_A_tr, err_A_tr_ROOT);

  TCanvas* c_tr = new TCanvas("c_tr", "c_tr", 1000, 450, 1000, 600);
  c_tr->SetGrid();
  c_tr->SetBottomMargin(0.12591304);
  c_tr->cd();

  TGraphErrors* TR_plot = new TGraphErrors(freq_tr_ROOT, A_tr_ROOT, err_freq_tr_fake_ROOT, err_A_tr_ROOT);
  TR_plot->SetTitle("Funzione di Trasferimento Resistenza (T_{R}); f (Hz); A=V_{out}/V_{in}");
  TR_plot->SetMarkerStyle(20);
  TR_plot->SetMarkerSize(0.5);
  TR_plot->GetXaxis()->SetLabelSize(0.047);
  TR_plot->GetYaxis()->SetLabelSize(0.047);
  TR_plot->GetXaxis()->SetTitleSize(0.052);
  TR_plot->GetYaxis()->SetTitleSize(0.052);
  TR_plot->GetYaxis()->SetTitleOffset(1);
  TR_plot->GetXaxis()->SetMaxDigits(3);
  TR_plot->Draw("AP");

  TLegend* legend_TR = new TLegend(0.7947, 0.81234783, 0.8747, 0.893043, NULL, "brNDC");
  legend_TR->AddEntry(TR_plot, "A_{i}", "e");
  legend_TR->Draw("SAME");

  TF1* TR_Max_Par = new TF1("TR_Max", "pol2");
  TFitResultPtr result_par_TR = TR_plot->Fit(TR_Max_Par, "RS", "", f_TR_Max - 1000, f_TR_Max + 1000);
  TMatrixD cov_par_TR = result_par_TR->GetCovarianceMatrix();
  cov_par_TR.Print();
  Double_t cov_ab_TR = cov_par_TR[2][1];
  Double_t cov_ac_TR = cov_par_TR[2][0];
  Double_t cov_bc_TR = cov_par_TR[1][0];

  TR_Max_Par->SetNpx(200);

  Double_t a_max_TR = TR_Max_Par->GetParameter(2);
  Double_t b_max_TR = TR_Max_Par->GetParameter(1);
  Double_t c_max_TR = TR_Max_Par->GetParameter(0);
  Double_t a_max_TR_err = TR_Max_Par->GetParError(2);
  Double_t b_max_TR_err = TR_Max_Par->GetParError(1);
  Double_t c_max_TR_err = TR_Max_Par->GetParError(0);
  Double_t Delta_TR = b_max_TR * b_max_TR - 4 * a_max_TR * c_max_TR;
  Double_t Delta_TR_err = sqrt(pow((b_max_TR_err) * (2 * b_max_TR), 2) + pow((a_max_TR_err * 4 * c_max_TR), 2) + pow((c_max_TR_err * 4 * a_max_TR), 2) - 2 * (2 * b_max_TR) * (4 * c_max_TR) * cov_ab_TR - 2 * (2 * b_max_TR) * (4 * a_max_TR) * cov_bc_TR + (2 * 4 * c_max_TR * 4 * a_max_TR * cov_ac_TR));

  Double_t f0_TR = (-b_max_TR / (2 * a_max_TR));
  Double_t f0_TR_err = sqrt(pow((-b_max_TR_err) / (2 * a_max_TR), 2) + pow((((a_max_TR_err * b_max_TR)) / (2 * a_max_TR * a_max_TR)), 2) - 2 * (1 / (2 * a_max_TR)) * ((b_max_TR) / (2 * a_max_TR * a_max_TR)) * cov_ab_TR);
  Double_t A_f0_TR = (-Delta_TR / (4 * a_max_TR));
  Double_t A_f0_TR_err = sqrt(pow(((Delta_TR_err) / (4 * a_max_TR)), 2) + pow((((a_max_TR_err * Delta_TR)) / (4 * a_max_TR * a_max_TR)), 2));

  output << "f0_TR = " << f0_TR << " +- " << f0_TR_err << endl;
  output << "A_f0_TR = " << A_f0_TR << " +- " << A_f0_TR_err << endl;

  TF1* TR_f1 = new TF1("TR_f1", "pol1");
  Double_t f1_sx = 44750;
  Double_t f1_dx = 46250;
  TR_plot->Fit(TR_f1, "+R", "", f1_sx, f1_dx);

  TF1* TR_f2 = new TF1("TR_f2", "pol1");
  Double_t f2_sx = 48750;
  Double_t f2_dx = 50250;
  TR_plot->Fit(TR_f2, "+R", "", f2_sx, f2_dx);

  Double_t m_TR_f1 = TR_f1->GetParameter(1);
  Double_t q_TR_f1 = TR_f1->GetParameter(0);
  Double_t m_TR_f2 = TR_f2->GetParameter(1);
  Double_t q_TR_f2 = TR_f2->GetParameter(0);
  Double_t m_TR_f1_err = TR_f1->GetParError(1);
  Double_t q_TR_f1_err = TR_f1->GetParError(0);
  Double_t m_TR_f2_err = TR_f2->GetParError(1);
  Double_t q_TR_f2_err = TR_f2->GetParError(0);

  vector<double> x_rho_f1;
  vector<double> x_rho_f2;

  for (int i = 0; i < freq_tr.size(); i++) {
    double k = freq_tr_ROOT[i];
    if (k >= f1_sx and k <= f1_dx) {
      x_rho_f1.push_back(k);
    }
    if (k >= f2_sx and k <= f2_dx) {
      x_rho_f2.push_back(k);
    }
  }

  Double_t rho_f1 = -1 * mean(x_rho_f1) / sqrt(var(x_rho_f1) + pow(mean(x_rho_f1), 2));
  Double_t rho_f2 = -1 * mean(x_rho_f2) / sqrt(var(x_rho_f2) + pow(mean(x_rho_f2), 2));

  Double_t f1_TR = (A_f0_TR * (1 / sqrt(2)) - q_TR_f1) / m_TR_f1;
  Double_t f2_TR = (A_f0_TR * (1 / sqrt(2)) - q_TR_f2) / m_TR_f2;
  Double_t f1_TR_err = sqrt(pow((1 / (sqrt(2) * m_TR_f1)) * A_f0_TR_err, 2) + pow((-1 / m_TR_f1) * q_TR_f1_err, 2) + pow((((A_f0_TR / sqrt(2)) - q_TR_f1) / (m_TR_f1 * m_TR_f1)) * m_TR_f1_err, 2) + 2 * (((A_f0_TR / sqrt(2)) - q_TR_f1) / (m_TR_f1 * m_TR_f1)) * (1 / m_TR_f1) * rho_f1 * m_TR_f1_err * q_TR_f1_err);
  Double_t f2_TR_err = sqrt(pow((1 / (sqrt(2) * m_TR_f2)) * A_f0_TR_err, 2) + pow((-1 / m_TR_f2) * q_TR_f2_err, 2) + pow((((A_f0_TR / sqrt(2)) - q_TR_f2) / (m_TR_f2 * m_TR_f2)) * m_TR_f2_err, 2) + 2 * (((A_f0_TR / sqrt(2)) - q_TR_f2) / (m_TR_f2 * m_TR_f2)) * (1 / m_TR_f2) * rho_f2 * m_TR_f2_err * q_TR_f2_err);

  TF1* A_Taglio_TR = new TF1("TR_f2", "0.092937166", 0, 100000);
  A_Taglio_TR->SetLineWidth(2);
  A_Taglio_TR->SetLineStyle(2);
  A_Taglio_TR->SetLineColor(kGray + 3);
  A_Taglio_TR->Draw("Same");

  TLine* f0 = new TLine(f0_TR, 0.0095, f0_TR, 0.1465);
  f0->SetLineWidth(2);
  f0->SetLineStyle(2);
  f0->SetLineColor(kGray + 3);
  f0->Draw("Same");
  TLine* f1 = new TLine(f1_TR, 0.0095, f1_TR, 0.1465);
  f1->SetLineWidth(2);
  f1->SetLineStyle(2);
  f1->SetLineColor(kGray + 3);
  f1->Draw("Same");
  TLine* f2 = new TLine(f2_TR, 0.0095, f2_TR, 0.1465);
  f2->SetLineWidth(2);
  f2->SetLineStyle(2);
  f2->SetLineColor(kGray + 3);
  f2->Draw("Same");

  Double_t delta_f_TC = f2_TR - f1_TR;
  Double_t delta_f_TC_err = sqrt(f1_TR_err * f1_TR_err + f2_TR_err * f2_TR_err);
  Double_t Q_TR = f0_TR / delta_f_TC;
  Double_t Q_TR_err = sqrt(pow((f0_TR_err / delta_f_TC), 2) + pow(((f0_TR) / (delta_f_TC * delta_f_TC)) * delta_f_TC_err, 2));

  output << "f1_TR= " << f1_TR << " +- " << f1_TR_err << endl;
  output << "f2_TR= " << f2_TR << " +- " << f2_TR_err << endl;
  output << "Delta_f_TC= " << delta_f_TC << " +- " << delta_f_TC_err << endl;
  output << "Q_TR= " << Q_TR << " +- " << Q_TR_err << endl;

  c_tr->SaveAs("stima_semplice_TR.png");
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
