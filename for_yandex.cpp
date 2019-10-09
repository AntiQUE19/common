/* Калькулятор рациональных чисел( +, -, *, /)
пример ввода: 1/2 + 1/3 
вывод: 5/6*/
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <exception>

using namespace std;

class Rational {
public:
  Rational(){
		num = 0;
		denom = 1;
	};
	
  Rational(int numerator, int denominator){
	  
	if( denominator == 0 ){
		throw runtime_error("Invalid argument");
	}
	
	if( numerator == 0){
			num = 0;
			denom = 1;
			return;
	}
		
		int sign = MakeSign( numerator, denominator );
		int d = FindNOD( numerator, denominator);
		
		num = sign * numerator / d;
		denom = denominator / d;
	};

  int Numerator() const{
	return num;
	};
		
  int Denominator() const{
	return denom;
	};
	/*
	bool operator<( const Rational& a) const{
		return num * a.Denominator() - denom * a.Numerator() < 0;
	}	
	*/
	
private:
	int MakeSign( int &numerator, int &denominator ){		
		int sign;
		if( (numerator > 0 && denominator > 0) || (numerator < 0 && denominator < 0) ){
			sign = 1;
		}else{
			sign = -1;
		}
		numerator = abs(numerator);
		denominator = abs(denominator);
		return sign;
	};
	
	int FindNOD( int numerator, int denominator ){
		int result = 0;
		int temp;
		
		while( result == 0 ){		
			if( numerator > denominator )
			{
				temp = numerator;
				numerator = denominator;
				denominator = temp;
			}
	
			if( denominator % numerator == 0)
			{
				result = numerator;
			}
			denominator %= numerator;
		
		}
		return result;
	};
	



	int num;
	int denom;
};

Rational operator+( const Rational& a, const Rational& b){
	
	int num, denom;
	
	num = a.Numerator() * b.Denominator() + b.Numerator() * a.Denominator();
	denom = a.Denominator() * b.Denominator();
	
	Rational result(num, denom);
	
	return result;
}

Rational operator-( const Rational& a, const Rational& b){
	
	int num, denom;
	
	num = a.Numerator() * b.Denominator() - b.Numerator() * a.Denominator();
	denom = a.Denominator() * b.Denominator();
	
	Rational result(num, denom);
	
	return result;
}

bool operator==( const Rational& a, const Rational& b){
	return a.Numerator() == b.Numerator() && a.Denominator() == b.Denominator();
}

Rational operator*( const Rational& a, const Rational& b){
	
	int num, denom;
	
	num = a.Numerator() * b.Numerator();
	denom = a.Denominator() * b.Denominator();
	
	Rational result(num, denom);
	
	return result;
}

Rational operator/( const Rational& a, const Rational& b){
	
	if( b.Numerator() == 0 ){
		throw domain_error(" ");
	}
	
	int num, denom;
	
	num = a.Numerator() * b.Denominator();
	denom = a.Denominator() * b.Numerator();
	
	Rational result(num, denom);
	
	return result;
}

istream& operator>>( istream& in, Rational& r ){
	
	string value;
	int numerator, denominator;
	/*
	getline(in, value, '/');
	if(value.size() == 0){
		return in;
	}
	numerator = stoi(value);
	getline(in, value, ' ');
	denominator = stoi(value);
	*/
	in >> numerator;
	in.ignore(1);
	in >> denominator;
	
	r = Rational(numerator, denominator);
	return in;
}

ostream& operator<<( ostream& out, const Rational& r ){
	string value;
	value = to_string(r.Numerator()) + "/" + to_string(r.Denominator());
	out << value;
	
	return out;
}

Rational Calculate(const Rational& a, const Rational& b, const char& sign ){
		
	if( sign == '+' ){
		return a + b;
	}else if( sign == '-' ){
		return a - b;
	}else if( sign == '*' ){
		return a * b;
	}else if( sign == '/' ){
		if( b.Numerator() == 0 ){
			throw runtime_error("Division by zero");
		}
		return a / b;
	}else{
		throw runtime_error("");
	}
}

int main(){
	
	Rational a, b;
	char sign;
	
	try{
		cin >> a >> sign >> b; 
		cout << Calculate( a, b, sign ) << endl;
	} catch( runtime_error& e){
		cout << e.what() << endl;
	} 
	return 0;
}
