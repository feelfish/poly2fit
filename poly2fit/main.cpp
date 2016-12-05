#include"polyfit.h"


int main(){
	vector<double> timeseq = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
	vector<double> distseq = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
	vector<double> input_coe_mat;
	vector<double> output_coe_mat;

	 init1DCoMat(timeseq, distseq, input_coe_mat);
	 solve1DCoMat(input_coe_mat, output_coe_mat);

	 cout << output_coe_mat[0] << "   " << output_coe_mat[1] <<  endl;

	getchar();
	return 0;
}