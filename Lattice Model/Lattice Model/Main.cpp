#include "Polymer.h"
#define L_PEPTIDE 0.01
#define L_GLYCAN 0.01
using namespace std;

#include<iostream>
#include<fstream>
#include<ctime>
#include<iomanip>

int main(void)
{
	ofstream OUTPUT("Output.txt");

	if (!OUTPUT.is_open())
	{
		cout << "Error: output file cannot be opened." << endl;
		system("pause");
		return 1;
	}
	else
	{
		cout << "Output file opened successfully." << endl;
	}

	srand(time(NULL));

	Polymer StaphAureous;

	for(int dist = 0; dist < 1000; dist = dist + 1)
	{
		
		StaphAureous.Find_Energy_To_Stretch_Peptide(dist*L_PEPTIDE);
		StaphAureous.Find_Energy_To_Stretch_Glycan(dist*L_GLYCAN);
		OUTPUT << dist*L_GLYCAN << "\t" << StaphAureous.Return_Energy_To_Stretch_Glycan() << "\t" << dist*L_PEPTIDE << "\t" << StaphAureous.Return_Energy_To_Stretch_Peptide() << "\t" << dist*sqrt(L_PEPTIDE*L_PEPTIDE+L_GLYCAN*L_GLYCAN) << "\t" << sqrt(pow(StaphAureous.Return_Energy_To_Stretch_Peptide(),2) + pow(StaphAureous.Return_Energy_To_Stretch_Glycan(),2)) << endl;
	}
	system("pause");
	return 0;

}