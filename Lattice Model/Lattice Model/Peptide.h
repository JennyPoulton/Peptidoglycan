#ifndef PEPTIDE_HEADER
#define PEPTIDE_HEADER

class Peptide
{
private:
	int number_bonds_created;
	double spring_constant;
	double length;
public:
	Peptide();
	int Return_Number_Bonds();
	void Add_Bond();
	void Set_Bonds(int bonds);
	void Set_Length(double extension);

	double Return_Spring_Constant();
};

#endif