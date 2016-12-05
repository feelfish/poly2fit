#include<iostream>
#include<vector>
#include<math.h>

using namespace std;

void linTimeSeq(vector<double> &, double );
void quadTimeSeq(vector<double> &, vector<double>&);
void cubiTimeSeq(vector<double> &, vector<double>&);
void quarTimeSeq(vector<double> &, vector<double>&);


void	sumSeq(vector<double> &, double &);


void addDistSeq(vector<double> &, double &);

void DistTime1Seq(vector<double>&, vector<double> &, vector<double> &);
void DistTime2Seq(vector<double> &, vector<double> &, vector<double> &);

void init2DCoMat(vector<double>&, vector<double> &, vector<double>&);
void solve2DCoMat(vector<double>& , vector<double> &);

void init1DCoMat(vector<double>&, vector<double> &, vector<double>&);
void solve1DCoMat(vector<double>&, vector<double> &);