#include "../../programmi_utili/lib_analisi_dati.cxx"
#include "TGraphErrors.h"
#include "fstream"
#include "iostream"
#include "math.h"
#include "string"
#include "vector"
using namespace std;

void controllo_dati()
{
  //riempo i vettori con i dati
  vector<double> x0d;
  dati_in_txt("dati/dati_t-.txt", x0d, 3, 2);

  vector<double> x0c;
  dati_in_txt("dati/dati_t+.txt", x0c, 3, 2);

  vector<double> xmax;
  vector<double> ymax;
  dati_in_x_y_txt("dati/dati_t_max.txt", xmax, ymax, 5, 2, 4);

  vector<double> xmin;
  vector<double> ymin;
  dati_in_x_y_txt("dati/dati_t_min.txt", xmin, ymin, 5, 2, 4);

  vector<double> xosc;
  vector<double> yosc;
  dati_in_x_y_txt("dati/dati_oscilloscopio_001.txt", xosc, yosc, 2, 1, 2);

  TVectorD x0d_plot(x0d.size());
  TVectorD x0c_plot(x0c.size());
  TVectorD xmax_plot(xmax.size());
  TVectorD ymax_plot(xmax.size());
  TVectorD xmin_plot(xmin.size());
  TVectorD ymin_plot(xmin.size());
  TVectorD xosc_plot(xosc.size());
  TVectorD yosc_plot(xosc.size());
  TVectorD x_err_plot(x0d.size());

  TVector_convert(x0d, x0d_plot);
  TVector_convert(x0c, x0c_plot);
  TVector_convert(xmax, xmax_plot);
  TVector_convert(ymax, ymax_plot);
  TVector_convert(xmin, xmin_plot);
  TVector_convert(ymin, ymin_plot);
  TVector_convert(xosc, xosc_plot);
  TVector_convert(yosc, yosc_plot);
  TVector_fill_0(x_err_plot);

  TCanvas* c = new TCanvas("check", "check", 0, 400, 1000, 600);
  c->SetGrid();

  TGraph* g_0d = new TGraph(x0d_plot, x_err_plot);
  g_0d->SetMarkerStyle(24);
  g_0d->SetMarkerColor(kBlue);
  g_0d->SetTitle("Controllo dati ;t (s) ;V(t) (V)");
  g_0d->GetXaxis()->SetTitleSize(0.04);
  g_0d->GetXaxis()->SetLabelSize(0.047);
  g_0d->GetYaxis()->SetLabelSize(0.047);
  g_0d->SetMarkerSize(1.2);
  //g_0d->GetXaxis()->SetLimits(0, 8.5);
  g_0d->GetYaxis()->SetRangeUser(-10, 10);
  g_0d->GetXaxis()->SetTitleSize(0.052);
  g_0d->GetYaxis()->SetTitleSize(0.052);
  g_0d->GetYaxis()->SetTitleOffset(0.95);
  g_0d->GetXaxis()->SetMaxDigits(3);
  //g_0d->GetXaxis()->SetNdivisions(525, kTRUE);
  g_0d->Draw("AP");

  //faccio stampare tutti i punti
  TGraph* g_0c = new TGraph(x0c_plot, x_err_plot);
  g_0c->SetMarkerStyle(24);
  g_0c->SetMarkerColor(kRed);
  g_0c->Draw("P SAME");

  TGraph* g_max = new TGraph(xmax_plot, ymax_plot);
  g_max->SetMarkerStyle(20);
  g_max->SetMarkerColor(kRed);

  g_max->Draw("P SAME");

  TGraph* g_min = new TGraph(xmin_plot, ymin_plot);
  g_min->SetMarkerStyle(20);
  g_min->SetMarkerColor(kBlue);
  g_min->Draw("P SAME");

  TGraph* g_osc = new TGraph("dati/dati_oscilloscopio_001.txt", "%lg %lg");
  g_osc->SetMarkerStyle(20);
  g_osc->SetMarkerSize(0.2);
  g_osc->SetMarkerColor(kBlack);
  g_osc->Draw("P SAME");

  TLegend* legend_int_lin_T = new TLegend(0.7347, 0.57234783, 0.8747, 0.893043, NULL, "brNDC");
  legend_int_lin_T->AddEntry(g_0d, "t^{-}", "pe");
  legend_int_lin_T->AddEntry(g_0c, "t^{+}", "pe");
  legend_int_lin_T->AddEntry(g_min, "min", "pe");
  legend_int_lin_T->AddEntry(g_max, "max", "pe");
  legend_int_lin_T->AddEntry(g_osc, "Dati digitali", "p");
  legend_int_lin_T->Draw("SAME");

  c->SaveAs("plot/controllo_dati.png");
}
