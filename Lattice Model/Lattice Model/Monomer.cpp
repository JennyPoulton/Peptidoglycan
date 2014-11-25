#include <cmath>
#include <stdlib.h>
#include "Monomer.h"

using namespace std;

Monomer::Monomer()
{
	Set_Glycan_1_Bonds();
	Set_Peptide_1_Bonds();
	return;
}

void Monomer::Set_Glycan_1_Bonds()
{
	if (rand()/(double)RAND_MAX<0.833333333)
	{
		Glycan_1.Set_Bonds(1);
	}
	else
	{
		Glycan_1.Set_Bonds(0);
	}
	return;
}

void Monomer::Set_Peptide_1_Bonds()
{
	if (rand()/(double)RAND_MAX<0.5)
	{
		Peptide_1.Set_Bonds(1);
	}
	else
	{
		Peptide_1.Set_Bonds(0);
	}
	return;
}

void Monomer::Find_Set_Glycan_2_Bonds(int bonds)
{
	Glycan_2.Set_Bonds(bonds);
}

void Monomer::Find_Set_Peptide_2_Bonds(int bonds)
{
	Peptide_2.Set_Bonds(bonds);
}

int Monomer::Return_Glycan_1_Bonds()
{
	int bonds = Glycan_1.Return_Number_Bonds();
	return bonds;
}
int Monomer::Return_Glycan_2_Bonds()
{
	int bonds = Glycan_2.Return_Number_Bonds();
	return bonds;
}

int Monomer::Return_Peptide_1_Bonds()
{
	int bonds = Peptide_1.Return_Number_Bonds();
	return bonds;
}

int Monomer::Return_Peptide_2_Bonds()
{
	int bonds = Peptide_2.Return_Number_Bonds();
	return bonds;
}