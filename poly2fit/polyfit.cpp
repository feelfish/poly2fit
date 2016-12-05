#include"polyfit.h"

//using namespace std;


void linTimeSeq(vector<double> &curTimeSeq, double addSeq){
	curTimeSeq.push_back(addSeq);
}



void quadTimeSeq(vector<double> &curTimeSeq, vector<double>&quadSeq){
	for (vector<double>::iterator iterTimeSeq = curTimeSeq.begin(); iterTimeSeq != curTimeSeq.end(); ++iterTimeSeq){
		quadSeq.push_back((*iterTimeSeq) * (*iterTimeSeq));
	}
}


void cubiTimeSeq(vector<double> &curTimeSeq, vector<double>&cuibSeq){
	for (vector<double>::iterator iterTimeSeq = curTimeSeq.begin(); iterTimeSeq != curTimeSeq.end(); ++iterTimeSeq){
		cuibSeq.push_back((*iterTimeSeq)*(*iterTimeSeq)*(*iterTimeSeq));
	}
}



void quarTimeSeq(vector<double> &curTimeSeq, vector<double> &quarSeq){
	for (vector<double>::iterator iterTimeSeq = curTimeSeq.begin(); iterTimeSeq != curTimeSeq.end(); ++iterTimeSeq){
		quarSeq.push_back((*iterTimeSeq)*(*iterTimeSeq)*(*iterTimeSeq)*(*iterTimeSeq));
	}
}



void sumSeq(vector<double> &InputSeq, double &sumSeq){
	for (vector<double>::iterator iterSeq = InputSeq.begin(); iterSeq != InputSeq.end(); ++iterSeq){
		sumSeq = sumSeq + *iterSeq;
	}
}



void addDistSeq(vector<double> &curDistSeq, double addDist){
	curDistSeq.push_back(addDist);
}


void DistTime1Seq(vector<double>&curTimeSeq, vector<double> &curDistSeq, vector<double> &curDiTiSeq){
	if (curTimeSeq.size() != curDistSeq.size()) {
		cerr << "Length do not match" << endl; return;
	}
	for (vector<double>::iterator iterTi = curTimeSeq.begin(), iterDi = curDistSeq.begin(); iterTi != curTimeSeq.end(); ++iterTi, ++iterDi){
		curDiTiSeq.push_back((*iterTi)*(*iterDi));
	}
}


void DistTime2Seq(vector<double> &curTimeSeq, vector<double> &curDistSeq, vector<double> &DiTi2Seq){
	if (curTimeSeq.size() != curDistSeq.size()){
		cerr << "Size do not match" << endl;
		return;
	}

	for (vector<double>::iterator iterTi = curTimeSeq.begin(), iterDi = curDistSeq.begin(); iterTi != curTimeSeq.end(); ++iterTi, ++iterDi){
		DiTi2Seq.push_back((*iterDi)*(*iterTi)*(*iterTi));
	}
}

void init1DCoMat(vector<double>& Ti1Seq, vector<double>&DistSeq, vector<double>&CoeMat){
	vector<double> Ti2Seq, Di1Ti1Seq;
	quadTimeSeq(Ti1Seq, Ti2Seq);
	DistTime1Seq(Ti1Seq, DistSeq, Di1Ti1Seq);

	double sumTi0(0), sumTi1(0),sumTi2(0), sumDi1Ti0(0), sumDi1Ti1(0);
	sumTi0 = Ti1Seq.size();
	sumSeq(Ti1Seq, sumTi1);
	sumSeq(Ti2Seq, sumTi2);
	sumSeq(DistSeq, sumDi1Ti0);
	sumSeq(Di1Ti1Seq, sumDi1Ti1);

	CoeMat.push_back(sumTi0);					CoeMat.push_back(sumTi1);
	CoeMat.push_back(sumTi1);					CoeMat.push_back(sumTi2);

	// The Right Side Coefficients
	CoeMat.push_back(sumDi1Ti0);
	CoeMat.push_back(sumDi1Ti1);


}

void solve1DCoMat(vector<double>& Input_CoMat, vector<double> &Output_CoMat){
	int MatSize = 2 * 2;
	vector<double> adjuMat;
	adjuMat.push_back(Input_CoMat[3]);
	adjuMat.push_back(-Input_CoMat[2]);
	adjuMat.push_back(-Input_CoMat[1]);
	adjuMat.push_back(Input_CoMat[0]);

	double detInput_Co = Input_CoMat[0]*Input_CoMat[3] - Input_CoMat[1]*Input_CoMat[2];

	Output_CoMat.push_back((adjuMat[0] * Input_CoMat[4] + adjuMat[2] * Input_CoMat[5])/detInput_Co);
	Output_CoMat.push_back((adjuMat[1] * Input_CoMat[4] + adjuMat[3] * Input_CoMat[5])/detInput_Co);
}


void init2DCoMat(vector<double>&Ti1Seq, vector<double> &DistSeq, vector<double>&CoeMat){

	//Obtain Ti_n, DiTi_n Sequences
	vector<double> Ti2Seq, Ti3Seq, Ti4Seq, DiTi1Seq, DiTi2Seq;
	//Initialize time sequence in degree of T2,T3,T4
	quadTimeSeq(Ti1Seq, Ti2Seq);
	cubiTimeSeq(Ti1Seq, Ti3Seq);
	quarTimeSeq(Ti1Seq, Ti4Seq);
	//Initialize Distance*Time sequence in degree of D1T1, D1T2, 
	DistTime1Seq(Ti1Seq, DistSeq, DiTi1Seq);
	DistTime2Seq(Ti1Seq, DistSeq, DiTi2Seq);

	//Sum the Sequences above,
	double sumTi0(0), sumTi1(0), sumTi2(0), sumTi3(0), sumTi4(0), \
		sumDi1Ti0(0), sumDi1Ti1(0), sumDi1Ti2(0);
	//Sum Time Sequence in degree of 0,1,2,3,4
	sumTi0 = Ti1Seq.size();
	sumSeq(Ti1Seq, sumTi1);
	sumSeq(Ti2Seq, sumTi2);
	sumSeq(Ti3Seq, sumTi3);
	sumSeq(Ti4Seq, sumTi4);
	//Sum Distance*TimeSequance in degree of D1T0, D1T1, D1T2
	sumSeq(DistSeq, sumDi1Ti0);
	sumSeq(DiTi1Seq, sumDi1Ti1);
	sumSeq(DiTi2Seq, sumDi1Ti2);

	//Initialize Input Coefficients
	//The Left Side Coefficients 
	CoeMat.push_back(sumTi0);					CoeMat.push_back(sumTi1);				CoeMat.push_back(sumTi2);
	CoeMat.push_back(sumTi1);					CoeMat.push_back(sumTi2);				CoeMat.push_back(sumTi3);
	CoeMat.push_back(sumTi2);					CoeMat.push_back(sumTi3);				CoeMat.push_back(sumTi4);
	// The Right Side Coefficients
	CoeMat.push_back(sumDi1Ti0);
	CoeMat.push_back(sumDi1Ti1);
	CoeMat.push_back(sumDi1Ti2);

};



void solve2DCoMat(vector<double>& Input_CoMat, vector<double> &Output_CoMat){

	int MatSize = 3 * 3;
	//Generate Adjugate Matrix 
	vector<double> adjuMat;
	int i = 0;
	double temp;
	temp = pow(-1, i++ % 2) * (Input_CoMat[4] * Input_CoMat[8] - Input_CoMat[5] * Input_CoMat[7]);
	adjuMat.push_back(temp);
	temp = pow(-1, i++ % 2) * (Input_CoMat[3] * Input_CoMat[8] - Input_CoMat[5] * Input_CoMat[6]);
	adjuMat.push_back(temp);
	temp = pow(-1, i++ % 2) * (Input_CoMat[3] * Input_CoMat[7] - Input_CoMat[4] * Input_CoMat[6]);
	adjuMat.push_back(temp);
	temp = pow(-1, i++ % 2) * (Input_CoMat[1] * Input_CoMat[8] - Input_CoMat[2] * Input_CoMat[7]);
	adjuMat.push_back(temp),
		temp = pow(-1, i++ % 2) * (Input_CoMat[0] * Input_CoMat[8] - Input_CoMat[2] * Input_CoMat[6]);
	adjuMat.push_back(temp);
	temp = pow(-1, i++ % 2) * (Input_CoMat[0] * Input_CoMat[7] - Input_CoMat[1] * Input_CoMat[6]);
	adjuMat.push_back(temp);
	temp = pow(-1, i++ % 2) * (Input_CoMat[1] * Input_CoMat[5] - Input_CoMat[2] * Input_CoMat[4]);
	adjuMat.push_back(temp);
	temp = pow(-1, i++ % 2) * (Input_CoMat[0] * Input_CoMat[5] - Input_CoMat[2] * Input_CoMat[3]);
	adjuMat.push_back(temp);
	temp = pow(-1, i++ % 2) * (Input_CoMat[0] * Input_CoMat[4] - Input_CoMat[1] * Input_CoMat[3]);
	adjuMat.push_back(temp);

	//Obtain Determinant of Input_Coefficient Matrix
	i = 0;
	double det_InputCo;
	det_InputCo = \
		Input_CoMat[0] * pow(-1, i++ % 2)* (Input_CoMat[4] * Input_CoMat[8] - Input_CoMat[5] * Input_CoMat[7]) + \
		Input_CoMat[1] * pow(-1, i++ % 2)* (Input_CoMat[3] * Input_CoMat[8] - Input_CoMat[5] * Input_CoMat[6]) + \
		Input_CoMat[2] * pow(-1, i++ % 2)* (Input_CoMat[3] * Input_CoMat[7] - Input_CoMat[4] * Input_CoMat[6]);

	//Obtain Reverse
	Output_CoMat.push_back((adjuMat[0] * Input_CoMat[9] + adjuMat[3] * Input_CoMat[10] + adjuMat[6] * Input_CoMat[11]) / det_InputCo);
	Output_CoMat.push_back((adjuMat[1] * Input_CoMat[9] + adjuMat[4] * Input_CoMat[10] + adjuMat[7] * Input_CoMat[11]) / det_InputCo);
	Output_CoMat.push_back((adjuMat[2] * Input_CoMat[9] + adjuMat[5] * Input_CoMat[10] + adjuMat[8] * Input_CoMat[11]) / det_InputCo);

}