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

	/*for(int dist = 0; dist < 1000; dist = dist + 1)
	{
		
		StaphAureous.Find_Energy_To_Stretch_Peptide(dist*L_PEPTIDE);
		StaphAureous.Find_Energy_To_Stretch_Glycan(dist*L_GLYCAN);
		OUTPUT << dist*L_GLYCAN << "\t" << StaphAureous.Return_Energy_To_Stretch_Glycan() << "\t" << dist*L_PEPTIDE << "\t" << StaphAureous.Return_Energy_To_Stretch_Peptide() << "\t" << dist*sqrt(L_PEPTIDE*L_PEPTIDE+L_GLYCAN*L_GLYCAN) << "\t" << sqrt(pow(StaphAureous.Return_Energy_To_Stretch_Peptide(),2) + pow(StaphAureous.Return_Energy_To_Stretch_Glycan(),2)) << endl;
	}*/
	
	int Frequency_At_Length_Peptide[100];
	int Frequency_At_Length_Glycan[100];

	for (int i = 0; i < 100; i++)
	{
		for (int j = 0; j < 100; j++)
		{
			Frequency_At_Length_Peptide[(int)(StaphAureous.Return_Peptide_Length(i, j) / StaphAureous.Return_Splitter_Peptide())]++;
			Frequency_At_Length_Glycan[(int)(StaphAureous.Return_Glycan_Length(i, j) / StaphAureous.Return_Splitter_Glycan())]++;
		}
		
	}

	for (int i = 0; i < 100; i++)
	{
		OUTPUT << i*StaphAureous.Return_Splitter_Peptide() << "\t" << Frequency_At_Length_Peptide[i] << "\t" << i*StaphAureous.Return_Splitter_Glycan() << "\t" << Frequency_At_Length_Glycan[i] << endl;
	}
	
	system("pause");
	return 0;

}