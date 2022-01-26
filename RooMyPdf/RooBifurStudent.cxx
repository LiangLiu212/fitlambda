#include "Riostream.h"
#include "RooBifurStudent.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
//#include "math.h"
#include "TMath.h"
ClassImp(RooBifurStudent)
	RooBifurStudent::RooBifurStudent(const char *name, const char *title,
			RooAbsReal& _x,
			RooAbsReal& _m,
			RooAbsReal& _gh,
			RooAbsReal& _gl,
			RooAbsReal& _Nh, 
			RooAbsReal& _Nl) :
		RooAbsPdf(name,title),
		x("x","x",this,_x),
		m("m","m",this,_m),
		gh("gh","gh",this,_gh),
		gl("gl","gl",this,_gl),
		Nh("Nh","Nh",this,_Nh),
		Nl("Nl","Nl",this,_Nl)
{
}


RooBifurStudent::RooBifurStudent(const RooBifurStudent& other, const char* name) :
	RooAbsPdf(other,name),
	x("x",this,other.x),
	m("m",this,other.m),
	gh("gh",this,other.gh),
	gl("gl",this,other.gl),
	Nh("Nh",this,other.Nh),
	Nl("Nl",this,other.Nl)
{
}



Double_t RooBifurStudent::evaluate() const
{

	double Ph = TMath::Gamma((Nh+1.0)/2.0) / ( (gh+gl) * TMath::Gamma(Nh/2.0) * sqrt(Nh) );
	double Pl = TMath::Gamma((Nl+1.0)/2.0) / ( (gh-gl) * TMath::Gamma(Nl/2.0) * sqrt(Nl) );

	double fra = 2*Ph*Pl/( (Ph + Pl) * sqrt(TMath::Pi()) );
	double value = 0;

	if( x >= m ){
		value = fra / TMath::Power(( 1.0 + (1.0/Nh)*( ((x-m)/(gl+gh)) * ((x-m)/(gl+gh)) ) ), ((Nh + 1.0)/2.0));
	}
	else{
		value = fra / TMath::Power(( 1.0 + (1.0/Nl)*( ((x-m)/(gh-gl)) * ((x-m)/(gh-gl)) ) ), ((Nl + 1.0)/2.0));
	}

	return value;
}

