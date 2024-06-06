#include "./utils.h"
//Vector operations
std::vector<double> operator-( const std::vector<double>& v1 , const std::vector<double>& v2 )
{
 std::vector<double> v = v1;
 for( unsigned int i=0; i < v1.size() ; i++ )
 { v[i] -= v2[i]; }
 return v; 
}

std::vector<double> operator+( const std::vector<double>& v1 , const std::vector<double>& v2 )
{
 std::vector<double> v = v1;
 for( unsigned int i=0; i < v1.size() ; i++ )
 { v[i] += v2[i]; }
 return v; 
}

std::vector<double> operator*( const std::vector<double>& v1 , const std::vector<double>& v2 )
{
 std::vector<double> v = v1;
 for( unsigned int i=0; i < v1.size() ; i++ )
 { v[i] *= v2[i]; }
 return v; 
}

std::vector<double> operator/( const std::vector<double>& v1 , const std::vector<double>& v2 )
{
 std::vector<double> v = v1;
 for( unsigned int i=0; i < v1.size() ; i++ )
 { v[i] /= v2[i]; }
 return v; 
}

std::vector<double> operator*( double d , const std::vector<double>& v1 )
{
 std::vector<double> v = v1;
 for( unsigned int i=0; i < v1.size() ; i++ )
 { v[i] *= d; }
 return v; 
}

std::vector<double> operator+( double d , const std::vector<double>& v1 )
{
 std::vector<double> v = v1;
 for( unsigned int i=0; i < v1.size() ; i++ )
 { v[i] += d; }
 return v; 
}

std::vector<double> operator+( const std::vector<double>& v1 , double d )
{
 std::vector<double> v = v1;
 for( unsigned int i=0; i < v1.size() ; i++ )
 { v[i] += d; }
 return v; 
}

std::vector<double> operator-( double d , const std::vector<double>& v1 )
{
 std::vector<double> v = v1;
 for( unsigned int i=0; i < v1.size() ; i++ )
 { v[i] = d - v1[i]; }
 return v; 
}

std::vector<double> operator-( const std::vector<double>& v1 , double d  )
{
 std::vector<double> v = v1;
 for( unsigned int i=0; i < v1.size() ; i++ )
 { v[i] -= d; }
 return v; 
}

void operator+=( std::vector<double>& v1, const std::vector<double>& v2 )
{
 for( unsigned int i=0; i < v1.size() ; i++ )
 { v1[i] += v2[i]; }
 return; 
}

void operator-=( std::vector<double>& v1, const std::vector<double>& v2 )
{
 for( unsigned int i=0; i < v1.size() ; i++ )
 { v1[i] -= v2[i]; }
 return; 
}

void operator/=( std::vector<double>& v1, const std::vector<double>& v2 )
{
 for( unsigned int i=0; i < v1.size() ; i++ )
 { v1[i] /= v2[i]; }
 return;  
} 

void operator*=( std::vector<double>& v1, const double& a )
{
 for( unsigned int i=0; i < v1.size() ; i++ )
 { v1[i] *= a; }
 return; 
}

void operator*=( std::vector<double>& v1, const std::vector<double>& v2 )
{
 for( unsigned int i=0; i < v1.size() ; i++ )
 { v1[i] *= v2[i]; }
 return;  
}

void operator/=( std::vector<double>& v1, const double& a )
{
 for( unsigned int i=0; i < v1.size() ; i++ )
 { v1[i] /= a; }
 return;  
}

void axpy( std::vector<double>* y, double& a , std::vector<double>& x )
{
 for( unsigned int i=0; i < (*y).size() ; i++ )
 {
  (*y)[i] += a * x[i] ; 
 }
 return ; 
}

void axpy( std::vector<double>* y, std::vector<double>& a , std::vector<double>& x )
{
 for( unsigned int i=0; i < (*y).size() ; i++ )
 {
  (*y)[i] += a[i] * x[i] ; 
 }
 return; 
}

void naxpy( std::vector<double>* y, double& a , std::vector<double>& x )
{
 for( unsigned int i=0; i < (*y).size() ; i++ )
 {
  (*y)[i] -= a * x[i] ; 
 }
 return ; 
}

void naxpy( std::vector<double>* y, std::vector<double>& a , std::vector<double>& x )
{
 for( unsigned int i=0; i < (*y).size() ; i++ )
 {
  (*y)[i] -= a[i] * x[i] ; 
 }
 return; 
}
