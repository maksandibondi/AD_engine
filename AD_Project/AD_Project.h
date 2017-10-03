#pragma once
#ifdef ADUTILITIES_EXPORTS  
#define ADUTILITIES_API __declspec(dllexport)   
#else  
#define ADUTILITIES_API __declspec(dllimport)   
#endif  

#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <numeric>

using namespace std;

namespace ADUTILITIES {
	// Class definitions
	class NodeOper;

	class Derivable {
	protected:
		NodeOper* _nodeOper = nullptr; // pointer to the operator node associated with this derivable
	public:
		double val;
		double deriv = 0;
		Derivable* symb_deriv = 0;
		bool flag = 0;
		//Constructors
		Derivable();
		Derivable(double value);
		Derivable(double value, double deriv);
		virtual ~Derivable() {}
		//Setters/getters
		void setNodeOper(NodeOper *_nodeOper);
		void configureNodeOper(char Type, Derivable f1, Derivable f2);
		void configureNodeOper(char Type, Derivable f1);
		NodeOper* getNodeAddress();
		void pushToTape(vector<NodeOper*> *tape);
		double getAdjoint();

		// friend function that can access protected members
		friend Derivable operator+(Derivable f1, Derivable f2);
		friend Derivable operator+(double f1, Derivable f2);
		friend Derivable operator+(Derivable f1, double f2);
		friend vector<Derivable> operator+(vector<Derivable> f1, vector<Derivable> f2);
		friend vector<Derivable> operator+(double f1, vector<Derivable> f2);
		friend vector<Derivable> operator+(vector<Derivable> f1, double f2);
		friend Derivable operator-(Derivable f1, Derivable f2);
		friend Derivable operator-(double f1, Derivable f2);
		friend Derivable operator-(Derivable f1, double f2);
		friend Derivable operator*(Derivable f1, Derivable f2);
		friend Derivable operator*(double f1, Derivable f2);
		friend Derivable operator*(Derivable f1, double f2);
		friend Derivable operator/(Derivable f1, Derivable f2);
		friend Derivable operator/(double f1, Derivable f2);
		friend Derivable operator/(Derivable f1, double f2);
		friend Derivable operator^(Derivable f1, Derivable f2);
		friend Derivable operator^(double f1, Derivable f2);
		friend Derivable operator^(Derivable f1, double f2);
		friend Derivable exponent(Derivable f);
		friend Derivable logarithm(Derivable f);
		friend Derivable Normal(Derivable f);
		friend Derivable abs(Derivable f);
	};

	class NodeOper {
	protected:
		NodeOper* node_f1 = nullptr;
		NodeOper* node_f2 = nullptr;
	public:
		char Type; // operator type
		double val;
		double adjoint = 0;

		//Constructors
		NodeOper();
		NodeOper(double value);
		virtual ~NodeOper() {}
		//Setters/getters
		NodeOper* get_f1_Adress();
		NodeOper* get_f2_Adress();
		void set_f1_adress(NodeOper* node_f1);
		void set_f2_adress(NodeOper* node_f2);
	};

	// singletone class containing global vars
	class GlobalVars {
	private:
		static GlobalVars* current_instance;
		GlobalVars();
	public:
		vector<NodeOper*>* globalTape = nullptr;
		char mode = NULL;
		static GlobalVars* instance();
	};
	

	// Algorithms we can put in separate library and create a good API

	void AAD(Derivable(*func)(vector<Derivable>), vector<Derivable> vars);

	Derivable FMAD(Derivable(*f)(vector<Derivable>, vector<double>), vector<Derivable> vec, vector<double> flag);

	vector<Derivable> vector_FMAD(vector<Derivable>(*f)(vector<vector<Derivable>>, vector<double>), vector<vector<Derivable>> vec, vector<double> flag);

	double delta_FDM(const double S, const double sigma, const double r, const double K, const double T, const double t, const double delta_S);
	double vega_FDM(const double S, const double sigma, const double r, const double K, const double T, const double t, const double delta_sigma);
	double rho_FDM(const double S, const double sigma, const double r, const double K, const double T, const double t, const double delta_r);
	double theta_FDM(const double S, const double sigma, const double r, const double K, const double T, const double t, const double delta_T);
	double thetaa_FDM(const double S, const double sigma, const double r, const double K, const double T, const double t, const double delta_t);
	double kaa_FDM(const double S, const double sigma, const double r, const double K, const double T, const double t, const double delta_k);
	// Functions to derive


	double BS(double S, double sigma, double r, double K, double T, double t);

	Derivable test(vector<Derivable> vec, vector<double> flag);

	Derivable test(vector<Derivable> vec);

	Derivable BS(vector<Derivable> vec, vector<double> flag);

	Derivable BS(vector<Derivable> vec);

	Derivable VarReduction(vector<Derivable> vec, vector<double> flag);

	Derivable VarReduction(vector<Derivable> vec);

	vector<Derivable> vectorBS(vector<vector<Derivable>> vec, vector<double> flag);


}