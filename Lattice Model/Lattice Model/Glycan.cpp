#include <cmath>
#include "Glycan.h"
using namespace std;

Glycan::Glycan()
{
	number_bonds_created = 0;
	return;
}

int Glycan::Return_Number_Bonds()
{
	return number_bonds_created;
}

void Glycan::Add_Bond()
{
	number_bonds_created++;
	return;
}

void Glycan::Set_Bonds(int bonds)
{
	number_bonds_created = bonds;
	return;
}

double Glycan::Return_Spring_Constant()
{
	return spring_constant;
}