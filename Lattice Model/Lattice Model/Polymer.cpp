#include <cmath>
#include "Polymer.h"
using namespace std;
#define PEPTIDE 3
#define GLYCAN 5
#define L_PEPTIDE 0.01
#define L_GLYCAN 0.01

Polymer::Polymer()
{
	Set_Joins();
	Spring_Constant_Peptide_Direction = 0;
	Spring_Constant_Glycan_Direction = 0;
	Find_Spring_Constant_Peptide_Direction();
	Find_Spring_Constant_Glycan_Direction();
	Energy_To_Stretch_Glycan = 0;
	Energy_To_Stretch_Peptide = 0;

	return;
}

void Polymer::Set_Joins()
{
	for (int i = 0; i < 99; i++)
	{
		for (int j = 0; j < 99; j++)
		{
			Polymer1[i][j].Find_Set_Glycan_2_Bonds(Polymer1[i + 1][j].Return_Glycan_1_Bonds());
			Polymer1[i][j].Find_Set_Peptide_2_Bonds(Polymer1[i][j + 1].Return_Peptide_1_Bonds());
		}
	}
	for (int i = 0; i < 99; i++)
	{
		Polymer1[i][99].Find_Set_Glycan_2_Bonds(Polymer1[i + 1][99].Return_Glycan_1_Bonds());
		Polymer1[i][99].Find_Set_Peptide_2_Bonds(Polymer1[i][0].Return_Peptide_1_Bonds());
	}
	for (int j = 0; j < 99; j++)
	{
		Polymer1[99][j].Find_Set_Glycan_2_Bonds(Polymer1[0][j].Return_Glycan_1_Bonds());
		Polymer1[99][j].Find_Set_Peptide_2_Bonds(Polymer1[99][j + 1].Return_Peptide_1_Bonds());
	}
}

double Polymer::Return_Energy_To_Stretch_Glycan()
{
	return Energy_To_Stretch_Glycan;
}

double Polymer::Return_Energy_To_Stretch_Peptide()
{
	return Energy_To_Stretch_Peptide;
}

void Polymer::Find_Energy_To_Stretch_Peptide(double x)
{
	Energy_To_Stretch_Peptide = Spring_Constant_Peptide_Direction*x*x;
	return;
}

void Polymer::Find_Energy_To_Stretch_Glycan(double y)
{
	Energy_To_Stretch_Glycan = Spring_Constant_Glycan_Direction*y*y;
	return;
}

double Polymer::Return_Spring_Constant_Peptide_Direction()
{
	return Spring_Constant_Peptide_Direction;
}

void Polymer::Find_Spring_Constant_Peptide_Direction()
{
	int Number[100];
	for (int i = 0; i < 100; i++)
	{
		Number[i] = 0;
		for (int j = 0; j < 100; j++)
		{
			Number[i] = Number[i] + PEPTIDE*Polymer1[j][i].Return_Peptide_1_Bonds();
		}
	}

	double Spring_Constant = 0;

	for (int i = 0; i < 100; i++)
	{
		Spring_Constant = Spring_Constant + 1 / (double)Number[i];
	}

	Spring_Constant_Peptide_Direction = Spring_Constant;
	return;
}

double Polymer::Return_Spring_Constant_Glycan_Direction()
{
	return Spring_Constant_Glycan_Direction;
}

void Polymer::Find_Spring_Constant_Glycan_Direction()
{
	int Number[100];
	for (int i = 0; i < 100; i++)
	{
		Number[i] = 0;
		for (int j = 0; j < 100; j++)
		{
			Number[i] = Number[i] + GLYCAN*Polymer1[i][j].Return_Glycan_1_Bonds();
		}
	}

	double Spring_Constant = 0;

	for (int i = 0; i < 100; i++)
	{
		Spring_Constant = Spring_Constant + 1 / (double)Number[i];
	}

	Spring_Constant_Glycan_Direction = Spring_Constant;
	return;
}

void Polymer::Find_Force_Peptides(double Total_Force)
{

	for (int p = 0; p < 100; p++)
	{
		for (int q = 0; q < 100; q++)
		{
			double Fraction1 = 1;
			double Fraction2 = 1;

			for (int r = q; r > 0; r--)
			{
				int Count_Right_Up = Bar_Above_Right(p, r);
				int Count_Right_Down = Bar_Below_Right(p, r);
				int Count_Left_Up = Bar_Above_Left(p, r);
				int Count_Left_Down = Bar_Below_Left(p, r);

				int Count_Right_Both = 0;
				int Count_Left_Both = 0;


				if (Count_Right_Up < Count_Right_Down)
				{
					Count_Right_Both = Count_Right_Up;
				}
				else
				{
					Count_Right_Both = Count_Right_Down;
				}

				if (Count_Left_Up < Count_Left_Down)
				{
					Count_Left_Both = Count_Left_Up;
				}
				else
				{
					Count_Left_Both = Count_Left_Down;
				}

				if (r == q)
				{
					Fraction1 = Fraction1 / (double)(Count_Left_Up + Count_Right_Up);
				}
				else if (r == 0)
				{
					Fraction1 = Fraction1* (double)(Count_Left_Both + Count_Right_Both) / 100;
				}
				else
				{
					Fraction1 = Fraction1* (double)(Count_Left_Both + Count_Right_Both) / (double)(Count_Left_Up + Count_Right_Up);
				}
			}

			for (int r = q; r < 100; r++)
			{
				int Count_Right_Up = Bar_Above_Right(p, r);
				int Count_Right_Down = Bar_Below_Right(p, r);
				int Count_Left_Up = Bar_Above_Left(p, r);
				int Count_Left_Down = Bar_Below_Left(p, r);

				int Count_Right_Both = 0;
				int Count_Left_Both = 0;


				if (Count_Right_Up < Count_Right_Down)
				{
					Count_Right_Both = Count_Right_Up;
				}
				else
				{
					Count_Right_Both = Count_Right_Down;
				}

				if (Count_Left_Up < Count_Left_Down)
				{
					Count_Left_Both = Count_Left_Up;
				}
				else
				{
					Count_Left_Both = Count_Left_Down;
				}

				if (r == q)
				{
					Fraction2 = Fraction2 / (double)(Count_Left_Down + Count_Right_Down);
				}
				else if (r == 100)
				{
					Fraction2 = Fraction2* (double)(Count_Left_Both + Count_Right_Both) / 100;
				}
				else
				{
					Fraction2 = Fraction2* (double)(Count_Left_Both + Count_Right_Both) / (double)(Count_Left_Down + Count_Right_Down);
				}

				Peptide_Force[p][q] = Fraction1 + Fraction2;
			}
			
			Peptide_Force[p][q] = Fraction1 + Fraction2;
		}

	}
}


int Polymer::Bar_Above_Left(int i, int j)
{
	int count_left = 0;

	for (int k = 1; k < 100; k++)
	{
		if (Polymer1[i - k][j - 1].Return_Glycan_1_Bonds == 1)
		{
			count_left++;
		}
		else if (Polymer1[i - k][j - 1].Return_Glycan_1_Bonds == 0)
		{
			return count_left;
		}
	}
}

int Polymer::Bar_Above_Right(int i, int j)
{
	int count_left = 0;

	for (int k = 1; k < 100; k++)
	{
		if (Polymer1[i + k][j - 1].Return_Glycan_2_Bonds == 1)
		{
			count_left++;
		}
		else if (Polymer1[i - k][j - 1].Return_Glycan_2_Bonds == 0)
		{
			return count_left;
		}
	}
}


int Polymer::Bar_Below_Left(int i, int j)
{
	int count_left = 0;

	for (int k = 1; k < 100; k++)
	{
		if (Polymer1[i - k][j + 1].Return_Glycan_1_Bonds == 1)
		{
			count_left++;
		}
		else if (Polymer1[i - k][j + 1].Return_Glycan_1_Bonds == 0)
		{
			return count_left;
		}
	}
}


int Polymer::Bar_Below_Right(int i, int j)
{
	int count_right = 0;

	for (int k = 1; k < 100; k++)
	{
		if (Polymer1[i + k][j + 1].Return_Glycan_2_Bonds == 1)
		{
			count_right++;
		}
		else if (Polymer1[i - k][j + 1].Return_Glycan_2_Bonds == 0)
		{
			return count_right;
		}
	}
}



