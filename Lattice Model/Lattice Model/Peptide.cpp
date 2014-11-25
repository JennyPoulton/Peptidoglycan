#include <cmath>
#include "Peptide.h"
using namespace std;

Peptide::Peptide()
{
	number_bonds_created = 0;
	return;
}

int Peptide::Return_Number_Bonds()
{
	return number_bonds_created;
}

void Peptide::Add_Bond()
{
	number_bonds_created++;
	return;
}

void Peptide::Set_Bonds(int bonds)
{
	number_bonds_created=bonds;
	return;
}

double Peptide::Return_Spring_Constant()
{
	return spring_constant;
}