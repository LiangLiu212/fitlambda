#include "Riostream.h"
#include "RooStudent.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
//#include "math.h"
#include "TMath.h"
ClassImp(RooStudent)
	RooStudent::RooStudent(const char *name, const char *title,
			RooAbsReal& _x,
			RooAbsReal& _m,
			RooAbsReal& _g,
			RooAbsReal& _N) :
		RooAbsPdf(name,title),
		x("x","x",this,_x),
		m("m","m",this,_m),
		g("g","g",this,_g),
		N("N","N",this,_N)
{
}


RooStudent::RooStudent(const RooStudent& other, const char* name) :
	RooAbsPdf(other,name),
	x("x",this,other.x),
	m("m",this,other.m),
	g("g",this,other.g),
	N("N",this,other.N)
{
}



Double_t RooStudent::evaluate() const
{
	double P = TMath::Gamma((N+1.0)/2.0)/TMath::Gamma(N/2.0)/g/sqrt(N);
	double part1 = P/sqrt(TMath::Pi());
	double part2 = TMath::Power((1.0+(1.0/N)*((x-m)/g)*((x-m)/g)),((N+1.0)/2.0));

	double value = part1/part2;
	return value;
}

