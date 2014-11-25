#ifndef MONOMER_HEADER
#define MONOMER_HEADER

#include "Peptide.h"
#include "Glycan.h"

class Monomer
{
private:
	Glycan Glycan_1;
	Glycan Glycan_2;
	Peptide Peptide_1;
	Peptide Peptide_2;

	int Glycan_1_Bonds;
	int Glycan_2_Bonds;
	int Peptide_1_Bonds;
	int Peptide_2_Bonds;
	
	double Glycan_1_Length;
	double Glycan_2_Length;
	double Peptide_1_Length;
	double Peptide_2_Length;

	void Set_Glycan_1_Bonds();
	void Set_Peptide_1_Bonds();
	
	
public:
	Monomer();
	
	void Find_Set_Glycan_2_Bonds(int bonds);
	void Find_Set_Peptide_2_Bonds(int bonds);

	int Return_Glycan_1_Bonds();
	int Return_Glycan_2_Bonds();
	int Return_Peptide_1_Bonds();
	int Return_Peptide_2_Bonds();

	double Set_Glycan_1_Length(double extension);
	double Set_Glycan_2_Length(double extension);
	double Set_Peptide_1_Length(double extension);
	double Set_Peptide_2_Length(double extension);

	

};

#endif