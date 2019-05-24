/***************************Matrix.h**********************************************
*矩阵类
*
*Contract:
*		wuguanbin@shao.ac.cn
*声明：	
*		Matrix A(2,3);				生成一个2x3的0零矩阵
*		Matrix B(A);				生成一个与A相同的矩阵  
*		Matrix C;					声明C是一个矩阵，未初始化
*									
*功能：								
*		A[i]=k;						对矩阵A的第i(int)个元素赋值k(int or double)
*		A(i,j)=k;					对矩阵A的第i(int，从0开始计数)行第j(int，从0开始计数)列元素赋值k(int or double)
*		B=A;						用A矩阵向B矩阵赋值
*		B==A;						判断B与A是否相等，相等返回true，否则返回false
*		A+B;						A、B矩阵相加，可赋值给C矩阵，如：C=A+B;
*		A-B;						A、B矩阵相减，可赋值给C矩阵，如：C=A-B;
*		A*B;						A、B矩阵相乘，可赋值给C矩阵，如：C=A*B;
*		A*a;						a(int or double)与矩阵A数乘，可赋值给C矩阵，如：C=A*a;	
*		A.trans();					矩阵A的转置，可赋值给矩阵C，如：C=A.trans();
*		A.inv();					矩阵A的逆矩阵，可赋值给矩阵C，如：C=A.inv();
*		A.getRow();					获取矩阵A的行数
*		A.getCol();					获取矩阵A的列数
*		A.print();					将矩阵A输出到屏幕
*		erase=EraseZerosCol(A,B);	清除矩阵A中所有元素为0的列，并将清除后的结果赋给矩阵B，返回矩阵A中被清除的列号给erase(set<int>)
*
*Last Update:
*		2019 05 24:finished the whole class
*
**********************************************************************************/
#pragma once
#include<cmath>
#include<cstring>
#include<stdio.h>
#include<new>
#include<iostream>
#include<set>
#define POS(i,j,k) ((i)*j+k)
class Matrix{
public:
	Matrix(int n,int m);
	Matrix();
	Matrix(Matrix &a);
	double operator=(double a);
	friend Matrix operator*(Matrix &a,Matrix &b);
	Matrix operator*(double a);
	friend Matrix operator*(double a, Matrix&b);
	Matrix operator+(Matrix &a);
	Matrix operator-(Matrix &a);
	Matrix trans();
	bool operator==(Matrix &a);
	double& operator[](int a);
	Matrix &operator=(Matrix& other);
	double& operator()(int i, int j);
	double getRow();
	double getCol();
	Matrix inv();
	friend std::set<int> EraseZerosCol(Matrix &input,Matrix &output);
	void print();
	~Matrix();

	
private:
	int row,col;
	double *elements;
	void initMatrix();
	int matinv(double *A, int n);
	int *imat(int n, int m);
	double *mat(int n, int m);
	void matcpy(double *A, const double *B, int n, int m);
	int ludcmp(double *A, int n, int *indx, double *d);
	void lubksb(const double *A, int n, const int *indx, double *b);
};
