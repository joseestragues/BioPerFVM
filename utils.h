#ifndef __utils_h__
#define __utils_h__
#include<vector>

//Vector operations
std::vector<double> operator-( const std::vector<double>& v1 , const std::vector<double>& v2 );

std::vector<double> operator+( const std::vector<double>& v1 , const std::vector<double>& v2 );


std::vector<double> operator*( const std::vector<double>& v1 , const std::vector<double>& v2 );

std::vector<double> operator/( const std::vector<double>& v1 , const std::vector<double>& v2 );

std::vector<double> operator*( double d , const std::vector<double>& v1 );

std::vector<double> operator+( double d , const std::vector<double>& v1 );

std::vector<double> operator+( const std::vector<double>& v1 , double d );

std::vector<double> operator-( double d , const std::vector<double>& v1 );

std::vector<double> operator-( const std::vector<double>& v1 , double d  );

void operator+=( std::vector<double>& v1, const std::vector<double>& v2 );

void operator-=( std::vector<double>& v1, const std::vector<double>& v2 );

void operator/=( std::vector<double>& v1, const std::vector<double>& v2 );

void operator*=( std::vector<double>& v1, const double& a );

void operator*=( std::vector<double>& v1, const std::vector<double>& v2 );

void operator/=( std::vector<double>& v1, const double& a );

// y = y + a*x 
void axpy( std::vector<double>* y, double& a , std::vector<double>& x );
// y = y + a.*x
void axpy( std::vector<double>* y, std::vector<double>& a , std::vector<double>& x ); 

// y = y - a*x 
void naxpy( std::vector<double>* y, double& a , std::vector<double>& x );
// y = y - a.*x
void naxpy( std::vector<double>* y, std::vector<double>& a , std::vector<double>& x ); 
#endif