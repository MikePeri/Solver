#pragma once
#include <complex>
using namespace std;
using std::complex; 

namespace solver
{
    class RealVariable{
        public:
            double _a,_b,_c;
        //Constructor:
        RealVariable()
        {
            _a=0;
            _b=1;
            _c=0;
        }//RealVariable
        RealVariable(double a,double b,double c)
        {
            _a=a;
            _b=b;
            _c=c;
        }//RealVariable
        //----------------------------------------
        // friend global binary operators
        //----------------------------------------
        friend const RealVariable operator+ (const RealVariable& c1, const RealVariable& c2);
        friend const RealVariable operator+ (const RealVariable& c1, const double& c2);
        friend const RealVariable operator+ (const double& c1, const RealVariable& c2);
        //----------------------------------------------------------------------------------------------
        friend const RealVariable operator- (const RealVariable& c1, const RealVariable& c2);
        friend const RealVariable operator- (const RealVariable& c1, const double& c2);
        friend const RealVariable operator- (const double& c1, const RealVariable& c2);
        friend const RealVariable operator- (const RealVariable& c1);
        //----------------------------------------------------------------------------------------------
        friend const RealVariable operator* (const RealVariable& c1, const RealVariable& c2);
        friend const RealVariable operator* (const RealVariable& c1, const double& c2);
        friend const RealVariable operator* (const double& c1, const RealVariable& c2);
        //----------------------------------------------------------------------------------------------
        friend const RealVariable operator/ (const RealVariable& c1, const double& c2);
        //----------------------------------------------------------------------------------------------
        friend const RealVariable operator^(const RealVariable& c1, const double& c2);
        //----------------------------------------------------------------------------------------------
        friend const RealVariable operator==(const RealVariable& c1, const RealVariable& c2);
        friend const RealVariable operator==(const RealVariable& c1, const double& c2);
        friend const RealVariable operator==(const double& c1, const RealVariable& c2);
        //----------------------------------
        // friend global IO operators
        //----------------------------------
        friend ostream& operator<< (ostream& os, const RealVariable& c);
        //friend istream& operator>> (istream& is, RealVariable& c);
        //-------------------------------------
    };//RealVariable

    class ComplexVariable
    {
        public:
            complex<double> _a,_b,_comp;
        ComplexVariable()
        {
            _a=complex<double>(0,0);
            _b=complex<double>(1,0);
            _comp=complex<double>(0,0);
        }//Empty constructor
         ComplexVariable(complex<double> a,complex<double> b,complex<double> comp)
        {
            _a=a;
            _b=b;
            _comp=comp;
        }//Build constructor
        
            //----------------------------------------
            // friend global binary operators
            //----------------------------------------
            friend const ComplexVariable operator+ (const ComplexVariable& c1, const ComplexVariable& c2);
            friend const ComplexVariable operator+ (const ComplexVariable& c1, const complex<double>& c2);
            friend const ComplexVariable operator+ (const complex<double>& c1, const ComplexVariable& c2);
            friend const ComplexVariable operator+ (const ComplexVariable& c1, const double& c2);
            friend const ComplexVariable operator+ (const double& c1, const ComplexVariable& c2);
            //----------------------------------------------------------------------------------------------
            friend const ComplexVariable operator- (const ComplexVariable& c1, const ComplexVariable& c2);
            friend const ComplexVariable operator- (const ComplexVariable& c1, const complex<double>& c2);
            friend const ComplexVariable operator- (const complex<double>& c1, const ComplexVariable& c2);
            friend const ComplexVariable operator- (const ComplexVariable& c1, const double& c2);
            friend const ComplexVariable operator- (const double& c1, const ComplexVariable& c2);
            friend const ComplexVariable operator- (const ComplexVariable& c1);
             //----------------------------------------------------------------------------------------------
            friend const ComplexVariable operator* (const ComplexVariable& c1, const ComplexVariable& c2);
            friend const ComplexVariable operator* (const ComplexVariable& c1, const complex<double>& c2);
            friend const ComplexVariable operator* (const complex<double>& c1, const ComplexVariable& c2);
            friend const ComplexVariable operator* (const ComplexVariable& c1, const double& c2);
            friend const ComplexVariable operator* (const double& c1, const ComplexVariable& c2);
             //----------------------------------------------------------------------------------------------
            //friend const ComplexVariable operator/ (const ComplexVariable& c1, const ComplexVariable& c2);
            //friend const ComplexVariable operator/ (const ComplexVariable& c1, const complex<double>& c2);
            friend const ComplexVariable operator/ (const ComplexVariable& c1, const double& c2);
            //friend const ComplexVariable operator/ (const double& c1, const ComplexVariable& c2);
             //----------------------------------------------------------------------------------------------
            friend const ComplexVariable operator^ (const ComplexVariable& c1, const double& c2);
             //----------------------------------------------------------------------------------------------
            friend const ComplexVariable operator== (const ComplexVariable& c1, const ComplexVariable& c2);
            friend const ComplexVariable operator== (const ComplexVariable& c1, const complex<double>& c2);
            friend const ComplexVariable operator== (const complex<double>& c1, const ComplexVariable& c2);
            friend const ComplexVariable operator== (const ComplexVariable& c1, const double& c2);
            friend const ComplexVariable operator== (const double& c1, const ComplexVariable& c2);
             //----------------------------------------------------------------------------------------------
            // friend global IO operators
            //----------------------------------
            
            friend ostream& operator<< (ostream& os, const ComplexVariable& c);
            //friend istream& operator>> (istream& is, ComplexVariable& c);
            //-------------------------------------
    };//ComplexVariable
    string toString(complex<double> comp);
    double solve(RealVariable equation);
    std::complex<double> solve(ComplexVariable equation);
}//solver

namespace complexExtend{
    const complex<double> operator+(double r,const complex<double>& c2);
    const complex<double> operator+(const complex<double>&c1,double r);
    const complex<double> operator-(double r,const complex<double>& c2);
    const complex<double> operator-(const complex<double>&c1,double r);
    const complex<double> operator*(double r,const complex<double>& c2);
    const complex<double> operator*(const complex<double>&c1,double r);
    const complex<double> operator/(const complex<double>&c1,double r);

}//complexExtend
