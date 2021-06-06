#include "/home/nicolo/Documenti/universita_2/programmazione/lib_analisi_dati.cxx"
#include "TGraphErrors.h"
#include "fstream"
#include "iostream"
#include "math.h"
#include "string"
#include "vector"
using namespace std;

double Puls(double w, double d);
double ErrP(double a, double b, double c, double d);
double Indu(double w, double C);
double ErrI(double a, double b, double c, double d);

void analisi()
{
  /////////////////////////RIEMPO TUTTI I VETTORI CON TUTTI I DATI/////////////////////////////
  vector<double> x0d;
  vector<double> y0d;
  vector<double> y_err0d;
  dati_in_x_y_yerr_txt("dati/dati_t-.txt", x0d, y0d, y_err0d, 3, 1, 2, 3);

  vector<double> x0c;
  vector<double> y0c;
  vector<double> y_err0c;
  dati_in_x_y_yerr_txt("dati/dati_t+.txt", x0c, y0c, y_err0c, 3, 1, 2, 3);

  vector<double> xmax;
  vector<double> ymax;
  vector<double> y_errmax;
  dati_in_x_y_yerr_txt("dati/dati_t_max.txt", xmax, ymax, y_errmax, 5, 1, 2, 3);

  vector<double> xmin;
  vector<double> ymin;
  vector<double> y_errmin;
  dati_in_x_y_yerr_txt("dati/dati_t_min.txt", xmin, ymin, y_errmin, 5, 1, 2, 3);

  TVectorD x0d_plot(x0d.size());
  TVectorD y0d_plot(x0d.size());
  TVectorD y_err0d_plot(x0d.size());
  TVectorD x0c_plot(x0c.size());
  TVectorD y0c_plot(x0c.size());
  TVectorD y_err0c_plot(x0c.size());
  TVectorD xmax_plot(xmax.size());
  TVectorD ymax_plot(xmax.size());
  TVectorD y_errmax_plot(xmax.size());
  TVectorD xmin_plot(xmin.size());
  TVectorD ymin_plot(xmin.size());
  TVectorD y_errmin_plot(xmin.size());
  TVectorD x_err_plot(x0d.size());

  TVector_convert(x0d, x0d_plot);
  TVector_convert(y0d, y0d_plot);
  TVector_convert(y_err0d, y_err0d_plot);
  TVector_convert(x0c, x0c_plot);
  TVector_convert(y0c, y0c_plot);
  TVector_convert(y_err0c, y_err0c_plot);
  TVector_convert(xmax, xmax_plot);
  TVector_convert(ymax, ymax_plot);
  TVector_convert(y_errmax, y_errmax_plot);
  TVector_convert(xmin, xmin_plot);
  TVector_convert(ymin, ymin_plot);
  TVector_convert(y_errmin, y_errmin_plot);
  TVector_fill_0(x_err_plot, x0d.size());

  for (int i = 0; i < x0d.size(); i++) {
    cout << x0d_plot(i) << " " << x_err_plot(i) << " " << y0d_plot(i) << " " << y_err0d_plot(i) << endl;
  }
  for (int i = 0; i < x0d.size(); i++) {
    cout << x0c_plot(i) << " " << x_err_plot(i) << " " << y0c_plot(i) << " " << y_err0c_plot(i) << endl;
  }
  for (int i = 0; i < x0d.size(); i++) {
    cout << xmin_plot(i) << " " << x_err_plot(i) << " " << ymin_plot(i) << " " << y_errmin_plot(i) << endl;
  }
  for (int i = 0; i < x0d.size(); i++) {
    cout << xmax_plot(i) << " " << x_err_plot(i) << " " << ymax_plot(i) << " " << y_errmax_plot(i) << endl;
  }
  //  ////////////////////////////////////////INIZIO FIT LINEARI/////////////////////////////////////////
  TCanvas* c_fit_T = new TCanvas("residui", "residui", 0, 0, 2000, 600);
  c_fit_T->SetRightMargin(0.07505);
  c_fit_T->SetBottomMargin(0.12591304);
  c_fit_T->Divide(2, 1);

  c_fit_T->cd(1);
  gPad->SetGrid();
  TGraphErrors* g_0d = new TGraphErrors(x0d_plot, y0d_plot, x_err_plot, y_err0d_plot);
  g_0d->SetMarkerStyle(24);
  g_0d->SetMarkerColor(kBlue);
  g_0d->SetTitle("Interpolazione lineare (n,t);n ;t (s)");
  g_0d->GetXaxis()->SetTitleSize(0.04);
  g_0d->GetXaxis()->SetLabelSize(0.047);
  g_0d->GetYaxis()->SetLabelSize(0.047);
  g_0d->SetMarkerSize(1.2);
  g_0d->GetXaxis()->SetLimits(0, 8.5);
  g_0d->GetYaxis()->SetRangeUser(0, 0.00015);
  g_0d->GetXaxis()->SetTitleSize(0.052);
  g_0d->GetYaxis()->SetTitleSize(0.052);
  g_0d->GetYaxis()->SetTitleOffset(0.95);
  g_0d->GetXaxis()->SetMaxDigits(3);
  g_0d->GetXaxis()->SetNdivisions(525, kTRUE);
  g_0d->Draw("AP");

  //faccio stampare tutti i punti
  TGraphErrors* g_0c = new TGraphErrors(x0c_plot, y0c_plot, x_err_plot, y_err0c_plot);
  g_0c->SetMarkerStyle(24);
  g_0c->SetMarkerColor(kRed);
  g_0c->Draw("P SAME");

  TGraphErrors* g_max = new TGraphErrors(xmax_plot, ymax_plot, x_err_plot, y_errmax_plot);
  g_max->SetMarkerStyle(20);
  g_max->SetMarkerColor(kRed);

  g_max->Draw("P SAME");

  TGraphErrors* g_min = new TGraphErrors(xmin_plot, ymin_plot, x_err_plot, y_errmin_plot);
  g_min->SetMarkerStyle(20);
  g_min->SetMarkerColor(kBlue);

  g_min->Draw("P SAME");

  //inizio i fit lineari
  TF1* f_0d = new TF1("f_0c", "[0]+x*[1]", 0, 8);
  f_0d->SetLineStyle(7);
  f_0d->SetLineColor(kBlue);
  TF1* f_0c = new TF1("f_0c", "[0]+x*[1]", 0, 8);
  f_0c->SetLineStyle(7);
  f_0c->SetLineColor(kRed);
  TF1* f_max = new TF1("f_0c", "[0]+x*[1]", 0, 8);
  f_max->SetLineColor(kRed);
  TF1* f_min = new TF1("f_0c", "[0]+x*[1]", 0, 8);
  f_min->SetLineColor(kBlue);

  g_0c->Fit(f_0c, "+R");
  g_0d->Fit(f_0d, "+R");
  g_min->Fit(f_min, "+R");
  g_max->Fit(f_max, "+R");

  //cout di tutti i risultati con i vari errori
  ofstream output("risultati.txt");

  cout << "-------------CHI QUADRO---------------------\n";
  cout << "chi quadro t0- chi= " << f_0d->GetChisquare() << "\t con D.O.F.= " << f_0d->GetNDF() << endl;
  cout << "chi quadro t0+ chi= " << f_0c->GetChisquare() << f_0c->GetNDF() << endl;
  cout << "chi quadro min chi= " << f_min->GetChisquare() << f_min->GetNDF() << endl;
  cout << "chi quadro max chi= " << f_max->GetChisquare() << f_max->GetNDF() << endl;
  cout << "-------------------------------------------\n";
  output << "-------------CHI QUADRO---------------------\n";
  output << "chi quadro t0- chi= " << f_0d->GetChisquare() << "\t con D.O.F.= " << f_0d->GetNDF() << endl;
  output << "chi quadro t0+ chi= " << f_0c->GetChisquare() << f_0c->GetNDF() << endl;
  output << "chi quadro min chi= " << f_min->GetChisquare() << f_min->GetNDF() << endl;
  output << "chi quadro max chi= " << f_max->GetChisquare() << f_max->GetNDF() << endl;
  output << "--------------------------------------------\n";

  Double_t T_0d = f_0d->GetParameter(1);
  Double_t T_0c = f_0c->GetParameter(1);
  Double_t T_min = f_min->GetParameter(1);
  Double_t T_max = f_max->GetParameter(1);

  Double_t T_0d_err = f_0d->GetParError(1);
  Double_t T_0c_err = f_0c->GetParError(1);
  Double_t T_min_err = f_min->GetParError(1);
  Double_t T_max_err = f_max->GetParError(1);
  //calcolo media ponderata tra le T compatibili
  vector<double> T_0 = { T_0d, T_0c, T_max };
  vector<double> T_0_err = { T_0d_err, T_0c_err, T_max_err };
  double T_fin = media_ponderata(T_0, T_0_err);
  double T_err_fin = media_ponderata_err(T_0, T_0_err);

  cout << "-------------PERIODO------------------------\n";
  cout << "periodo t0- T= " << T_0d << "+-" << T_0d_err << endl;
  cout << "periodo t0+ T= " << T_0c << "+-" << T_0c_err << endl;
  cout << "periodo min T= " << T_min << "+-" << T_min_err << endl;
  cout << "periodo max T= " << T_max << "+-" << T_max_err << endl;
  cout << "la media ponderata tra la t0- e t0+ vale T= " << T_fin << "+-" << T_err_fin << endl;
  cout << "--------------------------------------------\n";
  output << "-------------PERIODO------------------------\n";
  output << "periodo t0- T= " << T_0d << "+-" << T_0d_err << endl;
  output << "periodo t0+ T= " << T_0c << "+-" << T_0c_err << endl;
  output << "periodo min T= " << T_min << "+-" << T_min_err << endl;
  output << "periodo max T= " << T_max << "+-" << T_max_err << endl;
  output << "la media ponderata tra la t0- e t0+ vale T= " << T_fin << "+-" << T_err_fin << endl;
  output << "--------------------------------------------\n";

  Double_t omegad_0d = 2 * M_PI / T_0d;
  Double_t omegad_0c = 2 * M_PI / T_0c;
  Double_t omegad_min = 2 * M_PI / T_min;
  Double_t omegad_max = 2 * M_PI / T_max;

  Double_t omegad_0d_err = 2 * M_PI * T_0d_err / pow(T_0d, 2);
  Double_t omegad_0c_err = 2 * M_PI * T_0c_err / pow(T_0c, 2);
  Double_t omegad_min_err = 2 * M_PI * T_min_err / pow(T_min, 2);
  Double_t omegad_max_err = 2 * M_PI * T_max_err / pow(T_max, 2);
  //media ponderata per le omegad compatibili
  vector<double> omegad_vec = { omegad_0d, omegad_0c, omegad_max };
  vector<double> omegad_err_vec = { omegad_0d_err, omegad_0c_err, omegad_max_err };
  double omegad = media_ponderata(omegad_vec, omegad_err_vec);
  double omegad_err = media_ponderata_err(omegad_vec, omegad_err_vec);

  cout << "-------------OMEGA_d------------------------\n";
  cout << "t0- omega_d= " << omegad_0d << "+-" << omegad_0d_err << endl;
  cout << "t0+ omega_d= " << omegad_0c << "+-" << omegad_0c_err << endl;
  cout << "min omega_d= " << omegad_min << "+-" << omegad_min_err << endl;
  cout << "max omega_d= " << omegad_max << "+-" << omegad_max_err << endl;
  cout << "la media ponderata tra la t0- e t0+ vale omega_d= " << omegad << "+-" << omegad_err << endl;
  cout << "--------------------------------------------\n";
  output << "-------------OMEGA_d------------------------\n";
  output << "t0- omega_d= " << omegad_0d << "+-" << omegad_0d_err << endl;
  output << "t0+ omega_d= " << omegad_0c << "+-" << omegad_0c_err << endl;
  output << "min omega_d= " << omegad_min << "+-" << omegad_min_err << endl;
  output << "max omega_d= " << omegad_max << "+-" << omegad_max_err << endl;
  output << "la media ponderata tra la t0- e t0+ vale omega_d= " << omegad << "+-" << omegad_err << endl;
  output << "--------------------------------------------\n";

  //legenda
  TLegend* legend_int_lin_T = new TLegend(0.1147, 0.55234783, 0.2547, 0.873043, NULL, "brNDC");
  legend_int_lin_T->AddEntry(g_0d, "t^{-}", "pe");
  legend_int_lin_T->AddEntry(g_0c, "t^{+}", "pe");
  legend_int_lin_T->AddEntry(g_min, "t_{min}", "pe");
  legend_int_lin_T->AddEntry(g_max, "t_{max}", "pe");
  legend_int_lin_T->Draw("SAME");

  //stampo i residui
  TVectorD res_0d(x0d.size());
  TVectorD res_0c(x0c.size());
  TVectorD res_min(xmin.size());
  TVectorD res_max(xmax.size());
  fill_residui_lin(T_0d, f_0d->GetParameter(0), res_0d, x0d_plot, y0d_plot);
  fill_residui_lin(T_0c, f_0c->GetParameter(0), res_0c, x0c_plot, y0c_plot);
  fill_residui_lin(T_min, f_min->GetParameter(0), res_min, xmin_plot, ymin_plot);
  fill_residui_lin(T_max, f_max->GetParameter(0), res_max, xmax_plot, ymax_plot);

  c_fit_T->cd(2);
  gPad->SetGrid();
  TGraphErrors* g_0d_res = new TGraphErrors(x0d_plot, res_0d, x_err_plot, y_err0d_plot);
  g_0d_res->SetMarkerStyle(24);
  g_0d_res->SetMarkerColor(kBlue);
  g_0d_res->SetTitle("Residui;n ;misurato-fit");
  g_0d_res->GetXaxis()->SetTitleSize(0.04);
  g_0d_res->GetXaxis()->SetLabelSize(0.047);
  g_0d_res->GetYaxis()->SetLabelSize(0.047);
  g_0d_res->SetMarkerSize(1.2);
  g_0d_res->GetXaxis()->SetLimits(0, 8.5);
  g_0d_res->GetYaxis()->SetRangeUser(0, 0.00015);
  g_0d_res->GetXaxis()->SetTitleSize(0.052);
  g_0d_res->GetYaxis()->SetTitleSize(0.052);
  g_0d_res->GetYaxis()->SetTitleOffset(0.95);
  g_0d_res->GetXaxis()->SetMaxDigits(3);
  g_0d_res->GetXaxis()->SetNdivisions(525, kTRUE);
  g_0d_res->GetYaxis()->SetRangeUser(-0.0000006, 0.0000006);
  g_0d_res->Draw("AP");

  TGraphErrors* g_0c_res = new TGraphErrors(x0c_plot, res_0c, x_err_plot, y_err0c_plot);
  g_0c_res->SetMarkerStyle(24);
  g_0c_res->SetMarkerColor(kRed);
  g_0c_res->Draw("P SAME");

  TGraphErrors* g_min_res = new TGraphErrors(xmin_plot, res_min, x_err_plot, y_errmin_plot);
  g_min_res->SetMarkerStyle(20);
  g_min_res->SetMarkerColor(kBlue);
  g_min_res->Draw("P SAME");

  TGraphErrors* g_max_res = new TGraphErrors(xmax_plot, res_max, x_err_plot, y_errmax_plot);
  g_max_res->SetMarkerStyle(20);
  g_max_res->SetMarkerColor(kRed);
  g_max_res->Draw("P SAME");

  TF1* zero = new TF1("zero", "0", 0, 8.5);
  zero->SetLineColor(kRed);
  zero->SetLineStyle(7);
  zero->Draw("SAME");

  //Salvo il canvas
  c_fit_T->SaveAs("~/Documenti/universita_2/SF2/Laboratorio_di_Fisica_PD/relazione_6_RLC/tempi/plot/in_lin_T.png");

  ///////////////////////////////////////////////////////////FIT DECREMENTO LOGARITMICO////////////////////////////////////////////

  ifstream input_max("dati/dati_log_max.txt");
  ifstream input_min("dati/dati_log_min.txt");

  vector<double> nmax;
  vector<double> Lmax;
  vector<double> L_errmax;
  vector<double> residui_max;
  vector<double> nmin;
  vector<double> Lmin;
  vector<double> L_errmin;
  vector<double> residui_min;
  dati_in_x_y_yerr_txt("dati/dati_log_max.txt", nmax, Lmax, L_errmax, 3, 1, 2, 3);
  dati_in_x_y_yerr_txt("dati/dati_log_min.txt", nmin, Lmin, L_errmin, 3, 1, 2, 3);

  TVectorD nmin_plot(nmin.size());
  TVectorD Lmin_plot(Lmin.size());
  TVectorD L_errmin_plot(L_errmin.size());
  TVector_convert(nmin, nmin_plot);
  TVector_convert(Lmin, Lmin_plot);
  TVector_convert(L_errmin, L_errmin_plot);

  TVectorD nmax_plot(nmax.size());
  TVectorD Lmax_plot(Lmax.size());
  TVectorD L_errmax_plot(L_errmax.size());
  TVector_convert(nmax, nmax_plot);
  TVector_convert(Lmax, Lmax_plot);
  TVector_convert(L_errmax, L_errmax_plot);

  TVectorD n_err_plot(nmin.size());
  for (int i; i < nmin.size(); i++) {
    n_err_plot(i) = 0;
  }

  TCanvas* c2 = new TCanvas("c2", "c2", 200, 10, 2000, 600);
  c2->SetGrid();
  c2->Divide(2, 1);
  c2->cd(1);
  gPad->SetGrid();
  TGraphErrors* DATI_min = new TGraphErrors(nmin_plot, Lmin_plot, n_err_plot, L_errmin_plot);
  DATI_min->SetMarkerStyle(20);
  DATI_min->SetMarkerColor(kBlue);
  DATI_min->SetTitle("Interpolazione lineare coppie (n, log(V_{n}/1 V));n;log(V_{n}/1 V)");
  DATI_min->GetXaxis()->SetTitleSize(0.05);
  DATI_min->GetXaxis()->SetLabelSize(0.047);
  DATI_min->GetYaxis()->SetLabelSize(0.047);
  DATI_min->SetMarkerSize(1.2);
  DATI_min->GetXaxis()->SetLimits(0, 8.5);
  DATI_min->GetXaxis()->SetTitleSize(0.052);
  DATI_min->GetYaxis()->SetTitleSize(0.052);
  DATI_min->GetYaxis()->SetTitleOffset(0.95);
  DATI_min->GetXaxis()->SetMaxDigits(3);
  DATI_min->GetXaxis()->SetNdivisions(525, kTRUE);
  DATI_min->GetYaxis()->SetNdivisions(513, kTRUE);

  DATI_min->Draw("AP");

  TGraphErrors* DATI_max = new TGraphErrors(nmax_plot, Lmax_plot, n_err_plot, L_errmax_plot);
  DATI_max->SetMarkerStyle(20);
  DATI_max->SetMarkerColor(kRed);
  DATI_max->Draw("P SAME");

  TF1* fit_min = new TF1("f_2_min", "[0]*x+[1]", -100, 100);
  fit_min->SetParNames("a", "b");
  fit_min->SetLineColor(kBlue);

  fit_min->SetParameter(0, 1);
  fit_min->SetParameter(1, 1);

  DATI_min->Fit(fit_min, "R");

  TF1* fit_max = new TF1("f_2_max", "[0]*x+[1]", -100, 100);
  fit_max->SetParNames("a", "b");
  fit_min->SetLineColor(kBlue);

  fit_max->SetParameter(0, 1);
  fit_max->SetParameter(1, 1);

  DATI_max->Fit(fit_max, "R+");

  TLegend* legend = new TLegend(0.67, 0.6, 0.83, 0.87, NULL, "brNDC");
  legend->AddEntry(DATI_min, "Ln(V_{min})", "pe");
  legend->AddEntry(DATI_max, "Ln(V_{max})", "pe");
  legend->Draw("SAME");

  double m_min = fit_min->GetParameter(0);
  double q_min = fit_min->GetParameter(1);
  double m_min_err = fit_min->GetParError(0);
  double q_min_err = fit_min->GetParError(1);
  for (int i; i < Lmin.size(); i++) {
    double res = (Lmin[i] - (m_min * nmin[i] + q_min));
    residui_min.push_back(res);
  }

  TVectorD residui_min_plot(residui_min.size());
  TVector_convert(residui_min, residui_min_plot);

  double m_max = fit_max->GetParameter(0);
  double q_max = fit_max->GetParameter(1);
  double m_max_err = fit_max->GetParError(0);
  double q_max_err = fit_max->GetParError(1);
  for (int i = 0; i < Lmax.size(); i++) {
    double res = (Lmax[i] - (m_max * nmax[i] + q_max));
    residui_max.push_back(res);
  }
  cout << "--------------RISULTATI FIT LOGARITMICO--------------------\n";
  cout << "m_min= " << m_min << "+-" << m_min_err << endl;
  cout << "q_min= " << q_min << "+-" << q_min_err << endl;
  cout << "m_max= " << m_max << "+-" << m_max_err << endl;
  cout << "q_max= " << q_max << "+-" << q_max_err << endl;
  cout << "-----------------------------------------------------------\n";
  output << "--------------RISULTATI FIT LOGARITMICO--------------------\n";
  output << "m_min= " << m_min << "+-" << m_min_err << endl;
  output << "q_min= " << q_min << "+-" << q_min_err << endl;
  output << "m_max= " << m_max << "+-" << m_max_err << endl;
  output << "q_max= " << q_max << "+-" << q_max_err << endl;
  output << "-----------------------------------------------------------\n";

  TVectorD residui_max_plot(residui_max.size());
  TVector_convert(residui_max, residui_max_plot);

  c2->cd(2);
  gPad->SetGrid();
  TGraphErrors* g_res_min = new TGraphErrors(nmin_plot, residui_min_plot, n_err_plot, L_errmin_plot);
  g_res_min->SetTitle("Residui;n;misurato-fit");
  g_res_min->SetMarkerStyle(20);
  g_res_min->SetMarkerColor(kBlue);
  g_res_min->GetXaxis()->SetTitleSize(0.04);
  g_res_min->GetXaxis()->SetLabelSize(0.047);
  g_res_min->GetYaxis()->SetLabelSize(0.047);
  g_res_min->SetMarkerSize(1.2);
  g_res_min->GetXaxis()->SetLimits(0, 8.5);
  g_res_min->GetYaxis()->SetRangeUser(-0.06, 0.06);
  g_res_min->GetXaxis()->SetTitleSize(0.052);
  g_res_min->GetYaxis()->SetTitleSize(0.052);
  g_res_min->GetYaxis()->SetTitleOffset(0.95);
  g_res_min->GetXaxis()->SetMaxDigits(3);
  g_res_min->GetXaxis()->SetNdivisions(525, kTRUE);
  g_res_min->Draw("AP");

  TGraphErrors* g_res_max = new TGraphErrors(nmax_plot, residui_max_plot, n_err_plot, L_errmax_plot);
  g_res_max->SetMarkerStyle(20);
  g_res_max->SetMarkerColor(kRed);
  g_res_max->Draw("P");

  TF1* zero1 = new TF1("zero", "0", -100, 100);
  zero1->SetLineStyle(7);
  zero1->SetLineColor(kRed);

  cout << "-------------CHI QUADRO FIT LOGARITMICO---------------------\n";
  cout << "chi quadro min chi_log_min= " << fit_min->GetChisquare() << " con D.O.F= " << f_min->GetNDF() << endl;
  cout << "chi quadro max chi_log_max= " << fit_max->GetChisquare() << " con D.O.F= " << f_max->GetNDF() << endl;
  cout << "-------------------------------------------\n";
  output << "-------------CHI QUADRO FIT LOGARITMICO---------------------\n";
  output << "chi quadro min chi_log_min= " << fit_min->GetChisquare() << " con D.O.F= " << f_min->GetNDF() << endl;
  output << "chi quadro max chi_log_max= " << fit_max->GetChisquare() << " con D.O.F= " << f_max->GetNDF() << endl;
  output << "-------------------------------------------\n";

  zero1->Draw("SAME");
  c2->SaveAs("plot/fit_logaritmico.png");

  //////////////////////Da qui inizio a tirare fuori parametri///////////////////

  //calcolo di delta
  double delta_max = (-1) * (m_max / T_fin);
  double delta_errmax = sqrt((pow(m_max_err * (1 / T_fin), 2)) + (pow(T_err_fin * (m_max / pow(T_fin, 2)), 2)));
  double delta_min = (-1) * (m_min / T_fin);
  double delta_errmin = sqrt((pow(m_min_err * (1 / T_fin), 2)) + (pow(T_err_fin * (m_min / pow(T_fin, 2)), 2)));

  //media ponderata tra le delta dei fit logaritmici
  vector<double> delta_vec = { delta_max, delta_min };
  vector<double> delta_err_vec = { delta_errmax, delta_errmin };
  double delta_fin = media_ponderata(delta_vec, delta_err_vec);
  double delta_fin_err = media_ponderata_err(delta_vec, delta_err_vec);
  cout << "------------DELTA------------------------\n";
  cout << "delta_max=" << delta_max << "+-" << delta_errmax << endl;
  cout << "delta_min=" << delta_min << "+-" << delta_errmin << endl;
  cout << "delta finale=" << delta_fin << "+-" << delta_fin_err << endl;
  cout << "-----------------------------------------\n";
  output << "------------DELTA------------------------\n";
  output << "delta_max=" << delta_max << "+-" << delta_errmax << endl;
  output << "delta_min=" << delta_min << "+-" << delta_errmin << endl;
  output << "delta finale=" << delta_fin << "+-" << delta_fin_err << endl;
  output << "-----------------------------------------\n";

  double C = 3.46E-09;
  double err_C = 0.02E-09;

  double pulsazione = Puls(omegad, delta_fin);
  double err_puls = ErrP(omegad, omegad_err, delta_fin, delta_fin_err);
  cout << "------------PULSAZIONE omega_0-------------------\n";
  cout << "omega_0= " << pulsazione << "+-" << err_puls << endl;
  cout << "-------------------------------------------------\n";
  output << "------------PULSAZIONE omega_0-------------------\n";
  output << "omega_0= " << pulsazione << "+-" << err_puls << endl;
  output << "-------------------------------------------------\n";

  //////////SCRIVERE L E ALTRE COSE
  //
  //ostream output("NOME FILE OUTPUT");
  //output << "il coefficiente angolare è:  " << m << "+-" << fit->GetParameterError(0) <<endl;
  //output << "Il test del chi^2 è: " << fit->GetChisquare() << " con " << fit->GetNDF() << " gradi di libertà" << endl;
  //output << "l'intercetta è: " << q << "+-" << fit->GetParameterError(1) << endl;
  //output << "il delta è: " << delta << "+-" << delta err << endl;
}

double Puls(double w_d, double d)
{
  double pul = sqrt(pow(w_d, 2) + pow(d, 2));
  return pul;
}

double ErrP(double w_d, double errw_d, double d, double err_d)
{
  double errp = sqrt((pow((w_d / sqrt(pow(w_d, 2) + pow(d, 2))), 2) * pow(errw_d, 2)) + (pow((d / sqrt(pow(w_d, 2) + pow(d, 2))), 2) * pow(err_d, 2)));
  return errp;
}

double Indu(double pul, double C)
{
  double ind = 1 / (pow(pul, 2) * C);
  return ind;
}

double ErrL(double pul, double errp, double C, double errC)
{
  double errL = sqrt((pow(2 / (pow(pul, 3) * C), 2) * (pow(errp, 2))) + (pow(1 / (pow(pul, 2) * pow(C, 2)), 2) * (pow(errC, 2))));
  return errL;
}
