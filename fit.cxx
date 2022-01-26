using namespace std;
using namespace RooFit;
R__LOAD_LIBRARY(/Users/liang/Documents/Anti-neutron/RooFit/RooMyPdf/RooStudent_cxx.so);
R__LOAD_LIBRARY(/Users/liang/Documents/Anti-neutron/RooFit/RooMyPdf/RooBifurStudent_cxx.so);
double efficiency[12][20];
void PlotFit(RooPlot *xframe, RooPlot *pullframe);
void PlotFit(const char *name, const char *title, RooPlot *xframe, RooPlot *pullframe, int N_obj);
double getval(double x, RooAbsPdf *pdf);
void fit()
{

		TFile *f1_mc = new TFile("MCW0912_708_refine.root", "read");
		TH1F *h2 = (TH1F*)f1_mc->Get("hmlam_vx_fes");
	//	TTree *t2 = (TTree*)f1_mc->Get("bkg");
		
		gSystem->Load("libRooFit");

		RooRealVar x("hmlam_vx_fes", "M(#Lambda) (GeV/#font[52]{c}^{2})",1.10,1.130);
	//	x.setBins(400, "cache");

		RooDataHist data1("data1", "data1", x, h2);
	//	RooDataSet data("data", "data", x, Import(*t2));

		RooRealVar mean1("mean1", "mean1", 1.1157, 0.9, 1.25);	
		RooRealVar sigma1("sigma1", "sigma1", 1.28779e-02, 0, 0.5);
		RooGaussian ga1("ga1", "ga1", x, mean1, sigma1);

		RooRealVar mean2("mean2", "mean2", 1.1157, 0.9, 1.25);	
		RooRealVar sigma2("sigma2", "sigma2",3.94886e-02, 0, 0.5);
		RooGaussian ga2("ga2", "ga2", x, mean2, sigma2);

		RooRealVar mean3("mean3", "mean3", 1.1157, 0.9, 1.22);	
		RooRealVar sigma3("sigma3", "sigma3",3.94886e-02, 0, 0.9);
		RooGaussian ga3("ga3", "ga3", x, mean3, sigma3);

		RooRealVar bifurmean1("bifurmean1", "bifurmean1", 1.1157, 0.9, 1.25);
		RooRealVar bifursigma11("bifursigma11", "bifursigma11", 6.06036e-03, -0.9, 0.9);
		RooRealVar bifursigma12("bifursigma12", "bifursigma12", -3.90685e-03, -0.9, 0.9);
		RooBifurGauss bifur1("bifur1", "bifur1", x, bifurmean1, bifursigma11, bifursigma12);

		RooRealVar bifurmean2("bifurmean2", "bifurmean2", 1.1157, 0.9, 1.25);
		RooRealVar bifursigma21("bifursigma21", "bifursigma21", 3.29338e-03, -0.5, 0.5);
		RooRealVar bifursigma22("bifursigma22", "bifursigma22", -6.44787e-03, -0.5, 0.5);
		RooBifurGauss bifur2("bifur2", "bifur2", x, bifurmean2, bifursigma21, bifursigma22);

		RooRealVar CBm1("CBm1", "CBm1", 1.197, 1.18, 1.22);
		RooRealVar CBsigma1("CBsigma1", "CBsigma1", 0.1, -0.5, 0.5);
		RooRealVar CBalpha1("CBalpha1", "CBalpha1", 0.1, -0.5, 0.5);
		RooRealVar CBn1("CBn1", "CBn1", 1, 0, 20);
		RooCBShape CB1("CBshape1", "CBshape1", x, CBm1, CBsigma1, CBalpha1, CBn1);

		RooRealVar STm1("STm1", "STm1", 1.1157, 0.9, 1.25);
		RooRealVar STg1("STg1", "STg1", 3.64949e-03, 0, 1.22);
		RooRealVar STN1("STN1", "STN1", 9.49169e-01, 0, 2.22);
		RooStudent ST1("ST1", "ST1", x, STm1, STg1, STN1);

		RooRealVar STm2("STm2", "STm2", 1.197, 1.18, 1.22);
		RooRealVar STg2("STg2", "STg2", 10.5, 0, 100.22);
		RooRealVar STN2("STN2", "STN2", 50, -100.22, 100.22);
		RooStudent ST2("ST2", "ST2", x, STm2, STg2, STN2);

		RooRealVar bifurSTm1("bifurSTm1", "bifurSTm1", 1.1157);
		RooRealVar bifurSTgh1("bifurSTgh1", "bifurSTgh1", 10.197, -100, 100);
		RooRealVar bifurSTgl1("bifurSTgl1", "bifurSTgl1", 1.197, -100, 100);
		RooRealVar bifurSTNh1("bifurSTNh1", "bifurSTNh1", 100, -1000, 1000);
		RooRealVar bifurSTNl1("bifurSTNl1", "bifurSTNl1", 100.197, -1000, 1000);
		RooBifurStudent bifurST1("bifurST1", "bifurST1", x, bifurSTm1, bifurSTgh1, bifurSTgl1, bifurSTNh1, bifurSTNl1);


		RooRealVar bifurmean3("bifurmean3", "bifurmean3", 1.1157, 0.9, 1.22);
		RooRealVar bifursigma31("bifursigma31", "bifursigma31", 0.0034, -0.05, 0.05);
		RooRealVar bifursigma32("bifursigma32", "bifursigma32", -0.0043, -0.05, 0.05);
		RooBifurGauss bifur3("bifur3", "bifur3", x, bifurmean3, bifursigma31, bifursigma32);




		RooRealVar fra1("fra1", "fra1", 2.00767e+06, 10, 20000000800);
		RooRealVar fra2("fra2", "fra2", 2.12591e+05, 10, 200008000);
		RooRealVar fra3("fra3", "fra3", 5.89562e+05, 1000, 200080000);
		RooRealVar fra4("fra4", "fra4", 7.49341e+05, 1000, 200800000);
		RooRealVar fra5("fra5", "fra5", 1.67102e+06, 10000, 200000000);
		RooAddPdf signal("signal", "ga1+g2", RooArgList(ga1, ga2, bifur1, bifur2), RooArgList(fra1, fra2, fra3, fra4));
	//	RooAddPdf signal("signal", "ga1+g2", RooArgList(ga1, ga2, bifur1, bifur2, ST1), RooArgList(fra1, fra2, fra3, fra4, fra5));
	//	RooAddPdf signal("signal", "ga1+g2+bifur1+bifur2", RooArgList(ga1, ga2, ga3, bifur2), RooArgList(fra1, fra2, fra3, fra4));
	//	RooAddPdf signal("signal", "ga1+g2+bifur1+bifur2", RooArgList(ga1, bifur1, bifur2, bifur3), RooArgList(fra1, fra2, fra3, fra4));
	//	RooAddPdf signal("signal", "ga1+g2+bifur1+bifur2", RooArgList(ga1,  bifur1, bifur2), RooArgList(fra1,  fra3, fra4));

		RooFitResult *res = signal.fitTo(data1, Extended(), Save(1));

	/*	
		RooWorkspace *w = new RooWorkspace("w");
		w->import(signal);
		w->Print();
		w->writeToFile("MCFit.root");
*/



		new TCanvas();
		RooPlot* xframe = x.frame();
		data1.plotOn(xframe,  Name("SignalMC"));
		signal.plotOn(xframe, LineColor(4), Name("SUM"));
		signal.plotOn(xframe, Components("ga1"), LineColor(2), LineStyle(2), Name("GAUSS1"), Title("Gaussian 1"));
		signal.plotOn(xframe, Components("ga2"), LineColor(3), LineStyle(2), Name("GAUSS1"), Title("Gaussian 1"));
		signal.plotOn(xframe, Components("bifur1"), LineColor(4), LineStyle(2), Name("BIFUR1"), Title("Gaussian 2"));
		signal.plotOn(xframe, Components("bifur2"), LineColor(6), LineStyle(2), Name("BIFUR2"), Title("Gaussian 2"));
	//	signal.plotOn(xframe, Components("ST1"), LineColor(7), LineStyle(2), Name("STUDENT1"), Title("Gaussian 2"));
		xframe->Draw();
		xframe->SetXTitle("M(#Lambda) (GeV/c^{2})");
		xframe->SetYTitle("Events/(1.0 MeV/c^{2})");

	//	RooHist *hresid = xframe->pullHist("SignalMC", "SUM");
	//	RooPlot *frame2 = x.frame(Title("Residual Distribution"));
	//	frame2->addPlotable(hresid);

	//	RooHist *hpull = xframe->pullHist("SignalMC", "SUM");
	//	RooPlot *frame3 = x.frame(Title("Pull Distribution"));
	//	frame3->addPlotable(hpull, "P");

	//	TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
	//	c2->cd();
	//	frame2->GetYaxis()->SetTitleOffset(1.6);
	//	frame2->Draw();

	//	cout << h2->GetNbinsX() << endl;

	//	TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);
	//	c3->cd();
	//	frame3->GetYaxis()->SetTitleOffset(1.6);
	//	frame3->Draw();
	//	PlotFit("MC", "MC", xframe, frame3, 7);

/*
		auto legend = new TLegend(0.2,0.45,0.5,0.85);
		legend->SetFillColor(0);
		legend->SetBorderSize(0);
		legend->AddEntry(xframe->getObject(0),"DATA","elp");
		legend->AddEntry(xframe->getObject(1),"SUM","l");
		legend->AddEntry(xframe->getObject(2),"Gaussian 1","l");
		legend->AddEntry(xframe->getObject(3),"Gaussian 2","l");
		legend->AddEntry(xframe->getObject(4),"Bifur Gaussian 1", "l");
		legend->AddEntry(xframe->getObject(5),"Bifur Gaussian 2", "l");
		legend->AddEntry((TObject*)0,"#chi^{2}/ndf = 3.0", "");
		legend->Draw();
	//	c1->SaveAs("range1.eps");


		TCanvas *c4 = new TCanvas("c4", "c4", 800, 600);
		c4->cd();
		TPad *pad1 = new TPad("pad1", "Fit part",0,0.3,1,1);
		pad1->SetBottomMargin(0);
		pad1->SetTopMargin(-1);
		pad1->Draw();
		pad1->cd();
		pad1->Range(0.05128203,0.6,2.102564,1.231579);
		xframe->Draw();
		xframe->GetYaxis()->CenterTitle(true);
		xframe->GetYaxis()->SetLabelFont(42);
	//	xframe->GetYaxis()->SetLabelOffset(0.01);
		xframe->GetYaxis()->SetLabelSize(0.07);
		xframe->GetYaxis()->SetTitleSize(0.1);
		xframe->GetYaxis()->SetTitleOffset(0.8);
		xframe->GetYaxis()->SetTitleFont(42);
		legend->Draw();
		c4->cd();
		TPad *pad2 = new TPad("pad2", "Comparison part",0,0.05,1,0.3);
		pad2->Draw();
		pad2->cd();
		pad2->SetGridy();
		pad2->Range(0.05128203,-4.333333,2.102564,2.333333);
		pad2->SetBottomMargin(0.35);
		frame3->Draw();
		frame3->GetXaxis()->CenterTitle(true);
		frame3->GetXaxis()->SetLabelFont(42);
		frame3->GetXaxis()->SetLabelOffset(0.01);
		frame3->GetXaxis()->SetLabelSize(0.15);
		frame3->GetXaxis()->SetTitleSize(0.25);
		frame3->GetXaxis()->SetTitleOffset(0.60);
		frame3->GetXaxis()->SetTitleFont(42);
		frame3->GetYaxis()->SetNdivisions(505);
		frame3->GetYaxis()->SetLabelSize(0.16);
		frame3->GetYaxis()->SetTitle("#chi");
		frame3->GetYaxis()->SetTitleSize(0.3);
		frame3->GetYaxis()->SetTitleOffset(0.26);
		frame3->GetYaxis()->CenterTitle(true);
		frame3->SetMaximum(5.5);
		frame3->SetMinimum(-5.5);
*/


}

void PlotFit(const char *name, const char *title, RooPlot *xframe, RooPlot *pullframe, int N_obj){

		auto legend = new TLegend(0.2,0.45,0.5,0.85);
		legend->SetFillColor(0);
		legend->SetFillStyle(0);
		legend->SetBorderSize(0);
		legend->AddEntry(xframe->getObject(0),xframe->getObject(0)->GetName(),"elp");
		for(int i = 1; i < N_obj; i++){
			legend->AddEntry(xframe->getObject(i),xframe->getObject(i)->GetName(),"l");
		}
		char chi[20];
		sprintf(chi , "#chi^{2}/ndf = %1.2f", xframe->chiSquare( xframe->getObject(1)->GetName(), xframe->getObject(0)->GetName() ,18));
		legend->AddEntry((TObject*)0, chi, "");
		legend->Draw();
	//	c1->SaveAs("range1.eps");


		TCanvas *c4 = new TCanvas(name, title, 800, 600);
		c4->cd();
		TPad *pad1 = new TPad("pad1", "Fit part",0,0.3,1,1);
		pad1->SetBottomMargin(0);
		pad1->SetTopMargin(-1);
		pad1->Draw();
		pad1->cd();
		pad1->Range(0.05128203,0.6,2.102564,1.231579);
		xframe->Draw();
		xframe->GetYaxis()->CenterTitle(true);
		xframe->GetYaxis()->SetLabelFont(42);
	//	xframe->GetYaxis()->SetLabelOffset(0.01);
		xframe->GetYaxis()->SetLabelSize(0.07);
		xframe->GetYaxis()->SetTitleSize(0.1);
		xframe->GetYaxis()->SetTitleOffset(0.8);
		xframe->GetYaxis()->SetTitleFont(42);
		xframe->SetMinimum(0.001);
		legend->Draw();
		c4->cd();


		TPad *pad2 = new TPad("pad2", "Comparison part",0,0.05,1,0.3);
		pad2->Draw();
		pad2->cd();
		pad2->SetGridy();
		pad2->Range(0.05128203,-4.333333,2.102564,2.333333);
		pad2->SetBottomMargin(0.35);
		pullframe->Draw();
		pullframe->GetXaxis()->CenterTitle(true);
		pullframe->GetXaxis()->SetLabelFont(42);
		pullframe->GetXaxis()->SetLabelOffset(0.01);
		pullframe->GetXaxis()->SetLabelSize(0.15);
		pullframe->GetXaxis()->SetTitleSize(0.25);
		pullframe->GetXaxis()->SetTitleOffset(0.60);
		pullframe->GetXaxis()->SetTitleFont(42);
		pullframe->GetYaxis()->SetNdivisions(505);
		pullframe->GetYaxis()->SetLabelSize(0.16);
		pullframe->GetYaxis()->SetTitle("#chi");
		pullframe->GetYaxis()->SetTitleSize(0.3);
		pullframe->GetYaxis()->SetTitleOffset(0.26);
		pullframe->GetYaxis()->CenterTitle(true);
		pullframe->SetMaximum(5.5);
		pullframe->SetMinimum(-5.5);
		char outfile[50];
		sprintf(outfile, "%s.eps", name);
		c4->SaveAs(outfile);
		char outfile1[50];
		sprintf(outfile1, "%s.png", name);
		c4->SaveAs(outfile1);

}


