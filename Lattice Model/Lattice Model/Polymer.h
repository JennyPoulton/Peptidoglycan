
#ifndef POLYMER_HEADER
#define POLYMER_HEADER

#include "Monomer.h"

class Polymer
{
private:
	Monomer Polymer1[100][100];
	double Glycan_Force[100][100];
	double Peptide_Force[100][100];
	int Current_Monomer;

	double Spring_Constant_Peptide_Direction;
	double Spring_Constant_Glycan_Direction;

	void Find_Spring_Constant_Peptide_Direction();
	void Find_Spring_Constant_Glycan_Direction();

	double Energy_To_Stretch_Peptide;
	double Energy_To_Stretch_Glycan;

	int Bar_Above_Right(int i, int j);
	int Bar_Above_Left(int i, int j);

	int Bar_Below_Right(int i, int j);
	int Bar_Below_Left(int i, int j);


	void Set_Joins();

public:
	Polymer();

	double Return_Energy_To_Stretch_Peptide();
	double Return_Energy_To_Stretch_Glycan();
	void Find_Energy_To_Stretch_Peptide(double x);
	void Find_Energy_To_Stretch_Glycan(double y);

	void Find_Force_Peptides(double Total_Force);
	void Find_Force_Glycans(double Total_Force);
	
	double Return_Spring_Constant_Peptide_Direction();
	double Return_Spring_Constant_Glycan_Direction();		

};

#endif