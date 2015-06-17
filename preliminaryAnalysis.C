#include <iostream>
#include <vector>

using namespace std;

#include "TFile.h"
#include "TTree.h"

int LengthV(vector<double> v) {
	return (int) v.size();
}

void preliminaryAnalysis(){

	TFile * f = new TFile("output.root");
	TTree * T_NaI = (TTree *)f->Get("TNaI");

	T_NaI->Draw("Sum$(edep)","LengthV(edep)>0");

}

