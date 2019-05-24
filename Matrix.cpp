//***************************Matrix.cpp***********************************
#include"Matrix.h"
Matrix::Matrix(int n,int m){
	//构造函数：输入行列数，构造 n x m 0矩阵
	if(n<=0||m<=0)
		throw("输入行列输有误");
	elements=new double[n*m];
	this->row=n;
	this->col=m;
	initMatrix();
}
Matrix::Matrix():row(0),col(0),elements(nullptr){
	//构造函数：0x0 矩阵
}
Matrix::Matrix(Matrix&a)
{	//构造函数：矩阵初始化
	if (col + row > 0)
		delete[] elements;
	row = a.row;
	col = a.col;
	elements = new double[row*col];
	for (int i = 0;i < row*col;++i)
		elements[i] = a.elements[i];

}
Matrix::~Matrix(){
	//析构函数
	if (row + col > 0)
	{
		delete[]elements;
		row = 0;
		col = 0;
	}
	else
	{
		delete elements;
		row = 0;
		col = 0;
	}
}
void Matrix::print(){
	//输出矩阵到屏幕
	printf("\n");
	for(int i=0;i<row;++i)
		for(int j=0;j<col;++j)
		{
			if(j==col-1)
				printf(" %10.3f\n",elements[POS(i,col,j)]);
			else
				printf(" %10.3f ",elements[POS(i,col,j)]);
		}
	printf("\n");
}
Matrix Matrix::trans(){
	//矩阵转置函数
	Matrix C(col,row);
	for(int i=0;i<row;++i)
		for(int j=0;j<col;++j)
		{
			C.elements[POS(j,row,i)]=elements[POS(i,col,j)];
		}
	return C;
}
double Matrix::getRow(){
	//获得矩阵行数
	return row;
}
double Matrix::getCol(){
	//获得矩阵列数
	return col;
}
double  Matrix::operator=(double a){
	//重载=运算符，用于对矩阵元素赋值
	return a;
}
Matrix &Matrix::operator=(Matrix &other) {
	//重载=运算符，用于矩阵对矩阵赋值
	if (row != 0 || col != 0)
	{
		if (row != other.row || col != other.col)
			throw("矩阵维度不一致不能赋值！");
	}
	if (row == 0 && col == 0) {
		row = other.row;
		col = other.col;
		elements = new double[row*col];
	}
	elements = new double[row*col];
	for (int i = 0;i < row*col;++i)
		elements[i] = other.elements[i];
	return *this;
}
Matrix Matrix::operator*(double a){
	//重载*运算符，用于数乘矩阵
	Matrix C(*this);
	for(int i=0;i<row*col;++i)
		C.elements[i]=elements[i]*a;
	return C;
}
Matrix operator*(double a, Matrix &b) {
	//重载*运算符，用于数乘矩阵
	Matrix C(b);
	for (int i = 0;i<C.row*C.col;++i)
		C.elements[i] = b.elements[i] * a;
	return C;
}
Matrix Matrix::operator+(Matrix &a){
	//重载+运算符，用于矩阵相加
	if(row!=a.row||col!=a.col)
		throw("不同维度矩阵不能相加");
	Matrix C(a);
	for(int i=0;i<row*col;++i)
		C.elements[i]=elements[i]+a[i];
	return C;
}
Matrix Matrix::operator-(Matrix &a){
	//重载+运算符，用于矩阵相减
	if(row!=a.row||col!=a.col)
		throw("不同维度矩阵不能相减");
	Matrix C(a);
	for(int i=0;i<row*col;++i)
		C.elements[i]=elements[i]-a[i];
	return C;
}
bool Matrix::operator==(Matrix &a){
	//重载==运算符，判断矩阵是否相等
	if(row!=a.row||col!=a.col)
		return false;

	for(int i=0;i<row*col;i++)
	{
		if(fabs(elements[i]-a.elements[i])>1e-12)
			return false;
	}
	return true;
}
double& Matrix::operator()(int i, int j) {
	//重载()运算符，用于获取第i行第j列元素,i j均从0开始计数
	if (i > row || j > col || i < 0 || j < 0)
		throw("输入的行列号下标不合法！");
	return elements[POS(i, col, j)];
}
Matrix operator*(Matrix &a,Matrix &b){
	//重载*运算符，用于矩阵相乘
	if(a.col!=b.row)
		throw("矩阵相乘维度必须一致");
	else if(a.row<=0||a.col<=0||b.row<=0||b.col<=0)
	{
		throw("矩阵维度必须大于0");
	}
	Matrix C(a.row,b.col);
	C.elements=new double[a.row*b.row];
	double tmp=0;
	int see;
	int k=b.col;
	for(int i=0;i<a.row;i++)
	{
		for(int j=0;j<k;++j)
		{
			for(int l=0;l<a.col;++l)
			{
				tmp+=a.elements[POS(i,a.col,l)]*b.elements[POS(l,k,j)];
			}
			see=POS(i,k,j);
			C.elements[POS(i,k,j)]=tmp;
			tmp=0;
		}
	}
	return C;
}
double& Matrix::operator[](int i){
	//重载[]运算符，用于获取矩阵第i个元素
	if(i<0||i>=row*col)
		throw("下标输入有误");
	return elements[i];
}
Matrix Matrix::inv(){
	//矩阵求逆
	if(row!=col||row<=0)
		throw("求逆矩阵必须为方阵");
	Matrix tmp1(row,row),tmp2(row,row),tmp3(row,row);
	int flag;
	tmp1=trans();
	if(matinv(tmp1.elements,row))
		throw("求逆矩阵出错，请检查矩阵满秩");
	tmp3.elements=new double[row*row];
	tmp2=tmp1.trans();
	for(int i=0;i<row*row;++i)
		tmp3.elements[i]=tmp2.elements[i];
	return tmp3;
}
std::set<int> EraseZerosCol(Matrix& input,Matrix& output) {
	//清楚元素全为零的列，返回被清除 set<int> 列号
	if (input.row <= 0 || input.col <= 0)
		throw("输入矩阵维度不能少于等于0！");
	if (output.row + output.col>0)
	{
		delete[]output.elements;
		output.row = 0;
		output.col = 0;
	}
	std::set<int> erase;

	for (int i = 0;i < input.col;++i) {
		int flag = 0;
		for (int j = 0;j < input.row;++j) {
			if (fabs(input(j, i))> 1e-12)
			{
				flag = 1;
				break;
			}
		}
		if (flag == 0)
			erase.insert(i);
	}
	int row = input.row;
	int col = input.col - (int)erase.size();
	if (col == 0) {
		return erase;
	}
	output.row = row;
	output.col = col;
	output.elements = new double[row*col];
	std::set<int>::iterator iter;
	int k = 0;
	int r, c;
	for (int j = 0;j < input.col;++j)
	{ 
		iter = erase.find(j);
		if (iter != erase.end())
			continue;
		for (int i = 0;i < input.row;++i) {
			c = POS(i, col, k);
			output(i, k) = input(i, j);
		}
		k++;
	}
	return erase;

}
//****************************************  Private Function  ********************************************
int Matrix::matinv(double *A, int n)
{

    double d,*B;
    int i,j,*indx;
    
    indx=imat(n,1); B=mat(n,n); matcpy(B,A,n,n);
    if (ludcmp(B,n,indx,&d)) {free(indx); free(B); return -1;}
    for (j=0;j<n;j++) {
        for (i=0;i<n;i++) A[i+j*n]=0.0; A[j+j*n]=1.0;
        lubksb(B,n,indx,A+j*n);
    }
    free(indx); free(B);
    return 0;
}
int* Matrix::imat(int n, int m)
{
    int *p;
    
    if (n<=0||m<=0) return NULL;
    if (!(p=(int *)malloc(sizeof(int)*n*m))) {
        throw("integer matrix memory allocation error: n=%d,m=%d\n",n,m);
    }
    return p;
}
double* Matrix::mat(int n, int m)
{
    double *p;
    
    if (n<=0||m<=0) return NULL;
    if (!(p=(double *)malloc(sizeof(double)*n*m))) {
        throw("matrix memory allocation error: n=%d,m=%d\n",n,m);
    }
    return p;
}
void Matrix::matcpy(double *A, const double *B, int n, int m)
{
    memcpy(A,B,sizeof(double)*n*m);
}
int Matrix::ludcmp(double *A, int n, int *indx, double *d)
{
    double big,s,tmp,*vv=mat(n,1);
    int i,imax=0,j,k;
    
    *d=1.0;
    for (i=0;i<n;i++) {
        big=0.0; for (j=0;j<n;j++) if ((tmp=fabs(A[i+j*n]))>big) big=tmp;
        if (big>0.0) vv[i]=1.0/big; else {free(vv); return -1;}
    }
    for (j=0;j<n;j++) {
        for (i=0;i<j;i++) {
            s=A[i+j*n]; for (k=0;k<i;k++) s-=A[i+k*n]*A[k+j*n]; A[i+j*n]=s;
        }
        big=0.0;
        for (i=j;i<n;i++) {
            s=A[i+j*n]; for (k=0;k<j;k++) s-=A[i+k*n]*A[k+j*n]; A[i+j*n]=s;
            if ((tmp=vv[i]*fabs(s))>=big) {big=tmp; imax=i;}
        }
        if (j!=imax) {
            for (k=0;k<n;k++) {
                tmp=A[imax+k*n]; A[imax+k*n]=A[j+k*n]; A[j+k*n]=tmp;
            }
            *d=-(*d); vv[imax]=vv[j];
        }
        indx[j]=imax;
        if (A[j+j*n]==0.0) {free(vv); return -1;}
        if (j!=n-1) {
            tmp=1.0/A[j+j*n]; for (i=j+1;i<n;i++) A[i+j*n]*=tmp;
        }
    }
    free(vv);
    return 0;
}
void Matrix::lubksb(const double *A, int n, const int *indx, double *b)
{
    double s;
    int i,ii=-1,ip,j;
    
    for (i=0;i<n;i++) {
        ip=indx[i]; s=b[ip]; b[ip]=b[i];
        if (ii>=0) for (j=ii;j<i;j++) s-=A[i+j*n]*b[j]; else if (s) ii=i;
        b[i]=s;
    }
    for (i=n-1;i>=0;i--) {
        s=b[i]; for (j=i+1;j<n;j++) s-=A[i+j*n]*b[j]; b[i]=s/A[i+i*n];
    }
}
void Matrix::initMatrix() {

	memset(elements, 0, sizeof(double)*row*col);
}
