#include "TAxis.h"
#include "TCanvas.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TH1.h"
#include "TLegend.h"
#include "TMath.h"
#include "TPad.h"
#include "TStyle.h"
#include "fstream"
#include "iostream"
#include "math.h"
#include "vector"

using namespace std;

void TVector_convert(vector<double>&, TVectorD&);

void fit_globale_sfasamenti_C()
{
    ifstream input("dati/dati_phi_C.txt");

    vector<double> f;
    vector<double> phi;
    vector<double> phi_err;
    vector<double> all;
    double u;
    while (input >> u)
    {
        all.push_back(u);
    }
    for (int i=0; i<all.size(); i=i+3)
    {
        f.push_back(all[i]);
        phi.push_back(all[i+1]);
        phi_err.push_back(all[i+2]);
    }
    TVectorD f_plot(f.size());
    TVectorD phi_plot(phi.size());
    TVectorD phi_err_plot(phi_err.size());
    TVector_convert(f, f_plot);
    TVector_convert(phi, phi_plot);
    TVector_convert(phi_err, phi_err_plot);


    TVectorD f_err_plot(f.size());
    for(int i=0; i< f.size(); i++)
    {
        f_err_plot(i)=f[i]*10e-06;
    }


 TCanvas* c_fit = new TCanvas("c_fit", "c_fit", 200, 100, 1000, 600);
    c_fit->SetGrid();
    c_fit->SetBottomMargin(0.12591304);

 TGraphErrors* g_0 = new TGraphErrors(f_plot, phi_plot, f_err_plot, phi_err_plot);
    g_0->SetTitle("Sfasamenti Condensatore; f (Hz); #Delta#phi");
    g_0->SetMarkerStyle(20);
    g_0->SetMarkerSize(0.5);
    g_0->GetXaxis()->SetLabelSize(0.047);
    g_0->GetYaxis()->SetLabelSize(0.047);
    g_0->GetXaxis()->SetTitleSize(0.052);
    g_0->GetYaxis()->SetTitleSize(0.052);
    g_0->GetYaxis()->SetTitleOffset(1);
    g_0->GetXaxis()->SetMaxDigits(3);
    g_0->Draw("AP");

    Double_t pi = M_PI;

 TF1* fit = new TF1("fit", "atan2( [0],( (([1]*[1])/(4*pi*x)) - ((2*pi*x)/2)) )", 30000, 80000);
    fit->SetParNames("Delta", "Omega");
    fit->SetParameter(0, 10000);
    fit->SetParameter(1, 304000);
    fit->SetNpx(400);
    fit->SetLineWidth(1);
    

  TFitResultPtr result_par_TC = g_0->Fit(fit, "S");

  TMatrixD cov_par_TC = result_par_TC->GetCovarianceMatrix();
  cov_par_TC.Print();


 TLegend* legend_PhiC = new TLegend(0.1547, 0.76234783, 0.2747, 0.893043, NULL, "brNDC");
    legend_PhiC->AddEntry(g_0, "#Delta#phi", "e");
    legend_PhiC->AddEntry(fit, "#Delta#phi (Fit)", "l");
    legend_PhiC->Draw("SAME");

    double Delta = fit->GetParameter(0);
    double Omega = fit->GetParameter(1);
    double chisq = fit->GetChisquare();


    ofstream output("risultati_fit_globale_sfasamenti.txx");
    output << "---------------- Risultati Fit Globale Sfasamenti -----------------" <<endl;
    output << "Il test del Chi^2 è: " << chisq << " con " << fit->GetNDF() << " gradi di libertà;" << endl;
    output << "Delta: " << Delta << "+-" << fit->GetParError(0) << endl;
    output << "Omega: " << Omega << "+-" << fit->GetParError(1) << endl;
    
 c_fit->SaveAs("fit_globale_sfasamenti_C.png");

}

void TVector_convert(vector<double>& v, TVectorD& v_root)
{
    for(int i = 0; i < v.size(); i++)
    {
        v_root(i) = v[i];
    }
}
