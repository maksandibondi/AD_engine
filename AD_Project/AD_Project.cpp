#include "AD_Project.h"


namespace ADUTILITIES {
	
	// Global Vars
	GlobalVars::GlobalVars(){
		globalTape = new vector<NodeOper*>();
		mode = char('fw');
	}
	GlobalVars *GlobalVars::current_instance = 0;
	GlobalVars* GlobalVars::instance(){
		if (!current_instance) {
			current_instance = new GlobalVars();
		}
		return current_instance;
	}

	// Derivable implementation
	Derivable::Derivable() {}
	Derivable::Derivable(double value) {
		this->_nodeOper = new NodeOper(value);
		this->val = value;
	}
	Derivable::Derivable(double value, double deriv) {
		this->val = value;
		this->deriv = deriv;
	}
	void Derivable::setNodeOper(NodeOper *_nodeOper) {
		this->_nodeOper = _nodeOper;
	}
	void Derivable::configureNodeOper(char Type, Derivable f1, Derivable f2) {
		this->_nodeOper->Type = Type;
		this->_nodeOper->set_f1_adress(f1.getNodeAddress());
		this->_nodeOper->set_f2_adress(f2.getNodeAddress());
	}
	void Derivable::configureNodeOper(char Type, Derivable f1) {
		this->_nodeOper->Type = Type;
		this->_nodeOper->set_f1_adress(f1.getNodeAddress());
		this->_nodeOper->set_f2_adress(nullptr);
	}
	NodeOper* Derivable::getNodeAddress() {
		return this->_nodeOper;
	}
	void Derivable::pushToTape(vector<NodeOper*> *tape) {
		tape->push_back(this->_nodeOper);
	}
	double Derivable::getAdjoint() {
		return (this->getNodeAddress()->adjoint);
	}


	// Implementations of Node Oper
	NodeOper::NodeOper() {
	}
	NodeOper::NodeOper(double value) {
		this->val = value;
		this->adjoint = 0;
	}
	NodeOper* NodeOper::get_f1_Adress() {
		return node_f1;
	}
	NodeOper* NodeOper::get_f2_Adress() {
		return node_f2;
	}
	void NodeOper::set_f1_adress(NodeOper* node_f1) {
		this->node_f1 = node_f1;
	}
	void NodeOper::set_f2_adress(NodeOper* node_f2) {
		this->node_f2 = node_f2;
	}


	// Operators overloading (elementary)
	Derivable operator+(Derivable f1, Derivable f2) {
		Derivable* fnew = nullptr;
		switch ((GlobalVars::instance())->mode) {
		case char('fw') :
			fnew = new Derivable(f1.val + f2.val, f1.deriv + f2.deriv);
			//(*fnew->symb_deriv) = f1.flag + f2.flag;
			break;
		case char('adj') :
			double value = (f1._nodeOper->val) + (f2._nodeOper->val);
			fnew = new Derivable(value);
			(*fnew).configureNodeOper('Plus', f1, f2);
			(*fnew).pushToTape(GlobalVars::instance()->globalTape);
		}
		return (*fnew);
	}
	Derivable operator+(double f1, Derivable f2) {
		Derivable* fnew = nullptr;
		Derivable* f1_new = nullptr;
		switch ((GlobalVars::instance())->mode) {
		case char('fw') :
			fnew = new Derivable(f1 + f2.val, f2.deriv);
			//(*fnew->symb_deriv) = f2.flag;
			break;
		case char('adj') :
			double value = f1 + (f2._nodeOper->val);
			f1_new = new Derivable(f1);
			fnew = new Derivable(value);
			(*fnew).configureNodeOper('Plus', *f1_new, f2);
			(*fnew).pushToTape(GlobalVars::instance()->globalTape);
			break;
		}
		return (*fnew);
	}
	Derivable operator+(Derivable f1, double f2) {
		Derivable* fnew = nullptr;
		Derivable* f2_new = nullptr;
		switch ((GlobalVars::instance())->mode) {
		case char('fw') :
			fnew = new Derivable(f1.val + f2, f1.deriv);
			//(*fnew->symb_deriv) = f1.flag;
			break;
		case char('adj') :
			double value = (f1._nodeOper->val) + f2;
			f2_new = new Derivable(f2);
			fnew = new Derivable(value);
			(*fnew).configureNodeOper('Plus', f1, *f2_new);
			(*fnew).pushToTape(GlobalVars::instance()->globalTape);
			break;
		}
		return (*fnew);
	}
	vector<Derivable> operator+(vector<Derivable> f1, vector<Derivable> f2) {
		vector<Derivable>* fnew = new vector<Derivable>();
		int sz = f1.size();
		for (int k = 0; k < sz; k++) {
			fnew->push_back(*(new Derivable(f1[k].val + f2[k].val, f1[k].deriv + f2[k].deriv)));
		}
		return (*fnew);
	}
	vector<Derivable> operator+(double f1, vector<Derivable> f2) {
		vector<Derivable>* fnew = new vector<Derivable>();
		int sz = f2.size();
		for (int k = 0; k < sz; k++) {
			fnew->push_back(*(new Derivable(f1 + f2[k].val, f2[k].deriv)));
		}
		return (*fnew);
	}
	vector<Derivable> operator+(vector<Derivable> f1, double f2) {
		vector<Derivable>* fnew = new vector<Derivable>();
		int sz = f1.size();
		for (int k = 0; k < sz; k++) {
			fnew->push_back(*(new Derivable(f1[k].val + f2, f1[k].deriv)));
		}
		return (*fnew);
	}

	Derivable operator-(Derivable f1, Derivable f2) {
		Derivable* fnew = nullptr;
		switch ((GlobalVars::instance())->mode) {
		case char('fw') :
			fnew = new Derivable(f1.val - f2.val, f1.deriv - f2.deriv);
			//(*fnew->symb_deriv) = f1.flag - f2.flag;
			break;
		case char('adj') :
			double value = (f1._nodeOper->val) - (f2._nodeOper->val);
			fnew = new Derivable(value);
			(*fnew).configureNodeOper('Min', f1, f2);
			(*fnew).pushToTape(GlobalVars::instance()->globalTape);
		}
		return (*fnew);
	}
	Derivable operator-(double f1, Derivable f2) {
		Derivable* fnew = nullptr;
		Derivable* f1_new = nullptr;
		switch ((GlobalVars::instance())->mode) {
		case char('fw') :
			fnew = new Derivable(f1 - f2.val, -f2.deriv);
			//(*fnew->symb_deriv) = -f2.flag;
			break;
		case char('adj') :
			double value = f1 - (f2._nodeOper->val);
			f1_new = new Derivable(f1);
			fnew = new Derivable(value);
			(*fnew).configureNodeOper('Min', *f1_new, f2);
			(*fnew).pushToTape(GlobalVars::instance()->globalTape);
			break;
		}
		return (*fnew);
	}
	Derivable operator-(Derivable f1, double f2) {
		Derivable* fnew = nullptr;
		Derivable* f2_new = nullptr;
		switch ((GlobalVars::instance())->mode) {
		case char('fw') :
			fnew = new Derivable(f1.val - f2, f1.deriv);
			//(*fnew->symb_deriv) = f1.flag;
			break;
		case char('adj') :
			double value = (f1._nodeOper->val) - f2;
			f2_new = new Derivable(f2);
			fnew = new Derivable(value);
			(*fnew).configureNodeOper('Min', f1, *f2_new);
			(*fnew).pushToTape(GlobalVars::instance()->globalTape);
			break;
		}
		return (*fnew);
	}
	vector<Derivable> operator-(vector<Derivable> f1, vector<Derivable> f2) {
		vector<Derivable>* fnew = new vector<Derivable>();
		int sz = f1.size();
		for (int k = 0; k < sz; k++) {
			fnew->push_back(*(new Derivable(f1[k].val - f2[k].val, f1[k].deriv - f2[k].deriv)));
		}
		return (*fnew);
	}
	vector<Derivable> operator-(double f1, vector<Derivable> f2) {
		vector<Derivable>* fnew = new vector<Derivable>();
		int sz = f2.size();
		for (int k = 0; k < sz; k++) {
			fnew->push_back(*(new Derivable(f1 - f2[k].val, -f2[k].deriv)));
		}
		return (*fnew);
	}
	vector<Derivable> operator-(vector<Derivable> f1, double f2) {
		vector<Derivable>* fnew = new vector<Derivable>();
		int sz = f1.size();
		for (int k = 0; k < sz; k++) {
			fnew->push_back(*(new Derivable(f1[k].val - f2, f1[k].deriv)));
		}
		return (*fnew);
	}


	Derivable operator*(Derivable f1, Derivable f2) {
		Derivable* fnew = nullptr;
		switch ((GlobalVars::instance())->mode) {
		case char('fw') :
			fnew = new Derivable(f1.val * f2.val, f1.deriv * f2.val + f1.val * f2.deriv);
			//(*fnew->symb_deriv) = f1.flag * f2 + f1 * f2.flag;
			break;
		case char('adj') :
			double value = (f1._nodeOper->val)*(f2._nodeOper->val);
			fnew = new Derivable(value);
			(*fnew).configureNodeOper('Mult', f1, f2);
			(*fnew).pushToTape(GlobalVars::instance()->globalTape);
			break;
		}
		return (*fnew);
	}
	Derivable operator*(double f1, Derivable f2) {
		Derivable* fnew = nullptr;
		Derivable* f1_new = nullptr;
		switch ((GlobalVars::instance())->mode) {
		case char('fw') :
			fnew = new Derivable(f1 * f2.val, f1 * (f2.deriv));
			//(*fnew->symb_deriv) = f1 * f2.flag;
			break;
		case char('adj') :
			double value = f1*(f2._nodeOper->val);
			f1_new = new Derivable(f1);
			fnew = new Derivable(value);
			(*fnew).configureNodeOper('Mult', *f1_new, f2);
			(*fnew).pushToTape(GlobalVars::instance()->globalTape);
			break;
		}
		return (*fnew);
	}
	Derivable operator*(Derivable f1, double f2) {
		Derivable* fnew = nullptr;
		Derivable* f2_new = nullptr;
		switch ((GlobalVars::instance())->mode) {
		case char('fw') :
			fnew = new Derivable(f1.val * f2, f1.deriv * f2);
			//(*fnew->symb_deriv) = f2 * f1.flag;
			break;
		case char('adj') :
			double value = (f1._nodeOper->val)*f2;
			f2_new = new Derivable(f2);
			fnew = new Derivable(value);
			(*fnew).configureNodeOper('Mult', f1, *f2_new);
			(*fnew).pushToTape(GlobalVars::instance()->globalTape);
			break;
		}
		return (*fnew);
	}
	vector<Derivable> operator*(vector<Derivable> f1, vector<Derivable> f2) {
		vector<Derivable>* fnew = new vector<Derivable>();
		int sz = f1.size();
		for (int k = 0; k < sz; k++) {
			fnew->push_back(*(new Derivable(f1[k].val*f2[k].val, f1[k].deriv*f2[k].val + f2[k].deriv*f1[k].val)));
		}
		return (*fnew);
	}
	vector<Derivable> operator*(double f1, vector<Derivable> f2) {
		vector<Derivable>* fnew = new vector<Derivable>();
		int sz = f2.size();
		for (int k = 0; k < sz; k++) {
			fnew->push_back(*(new Derivable(f1*f2[k].val, f1*(f2[k].deriv))));
		}
		return (*fnew);
	}
	vector<Derivable> operator*(vector<Derivable> f1, double f2) {
		vector<Derivable>* fnew = new vector<Derivable>();
		int sz = f1.size();
		for (int k = 0; k < sz; k++) {
			fnew->push_back(*(new Derivable(f1[k].val*f2, f2*(f1[k].deriv))));
		}
		return (*fnew);
	}


	Derivable operator/(Derivable f1, Derivable f2) {
		Derivable* fnew = nullptr;
		switch ((GlobalVars::instance())->mode) {
		case char('fw') :
			fnew = new Derivable(f1.val / f2.val, (f1.deriv * f2.val - f1.val * f2.deriv) / pow(f2.val, 2));
			//(*fnew->symb_deriv) = (f1.flag * (f2) - (f1) * f2.flag) / (f2) / (f2);
			break;
		case char('adj') :
			double value = (f1._nodeOper->val) / (f2._nodeOper->val);
			fnew = new Derivable(value);
			(*fnew).configureNodeOper('Div', f1, f2);
			(*fnew).pushToTape(GlobalVars::instance()->globalTape);
			break;
		}
		return (*fnew);
	}
	Derivable operator/(double f1, Derivable f2) {
		Derivable* fnew = nullptr;
		Derivable* f1_new = nullptr;
		switch ((GlobalVars::instance())->mode) {
		case char('fw') :
			fnew = new Derivable(f1 / f2.val, ((-f1) * f2.deriv) / pow(f2.val, 2));
			//(*fnew->symb_deriv) = (0 - f1*f2.flag) / (f2) / (f2);
			break;
		case char('adj') :
			double value = f1 / (f2._nodeOper->val);
			f1_new = new Derivable(f1);
			fnew = new Derivable(value);
			(*fnew).configureNodeOper('Div', *f1_new, f2);
			(*fnew).pushToTape(GlobalVars::instance()->globalTape);
			break;
		}
		return (*fnew);
	}
	Derivable operator/(Derivable f1, double f2) {
		Derivable* fnew = nullptr;
		Derivable* f2_new = nullptr;
		switch ((GlobalVars::instance())->mode) {
		case char('fw') :
			fnew = new Derivable(f1.val / f2, ((f1.deriv)/f2));
			//(*fnew->symb_deriv) = (f1.deriv*f2) / (f2) / (f2);
			break;
		case char('adj') :
			double value = (f1._nodeOper->val) / f2;
			f2_new = new Derivable(f2);
			fnew = new Derivable(value);
			(*fnew).configureNodeOper('Div', f1, *f2_new);
			(*fnew).pushToTape(GlobalVars::instance()->globalTape);

			break;
		}
		return (*fnew);
	}
	vector<Derivable> operator/(vector<Derivable> f1, vector<Derivable> f2) {
		vector<Derivable>* fnew = new vector<Derivable>();
		int sz = f1.size();
		for (int k = 0; k < sz; k++) {
			fnew->push_back(*(new Derivable(f1[k].val / f2[k].val, (f1[k].deriv * f2[k].val - f1[k].val * f2[k].deriv) / pow(f2[k].val, 2))));
		}
		return (*fnew);
	}
	vector<Derivable> operator/(double f1, vector<Derivable> f2) {
		vector<Derivable>* fnew = new vector<Derivable>();
		int sz = f2.size();
		for (int k = 0; k < sz; k++) {
			fnew->push_back(*(new Derivable(f1 / f2[k].val, ((-f1) * f2[k].deriv) / pow(f2[k].val, 2))));
		}
		return (*fnew);
	}
	vector<Derivable> operator/(vector<Derivable> f1, double f2) {
		vector<Derivable>* fnew = new vector<Derivable>();
		int sz = f1.size();
		for (int k = 0; k < sz; k++) {
			fnew->push_back(*(new Derivable(f1[k].val / f2, (f1[k].deriv * f2) / pow(f2, 2))));
		}
		return (*fnew);
	}


	Derivable operator^(Derivable f1, Derivable f2) {
		Derivable* fnew = nullptr;
		switch ((GlobalVars::instance())->mode) {
		case char('fw') :
			fnew = new Derivable(pow(f1.val, f2.val), pow(f1.val, f2.val)*log(f1.val)*(f2.deriv) + pow(f1.val, f2.val - 1)*f2.val*(f1.deriv));
			//(*fnew->symb_deriv) = f2*(f1 ^ (f2 - 1))*f2.flag;
			break;
		case char('adj') :
			double value = pow((f1._nodeOper->val), (f2._nodeOper->val));
			fnew = new Derivable(value); // Creating a new derivable. New node = nullptr. Val = value
			(*fnew).configureNodeOper('Pow', f1, f2);
			(*fnew).pushToTape(GlobalVars::instance()->globalTape);
			break;
		}
		return (*fnew);
	}
	Derivable operator^(double f1, Derivable f2) {
		Derivable* fnew = nullptr;
		Derivable* f1_new = nullptr;
		switch ((GlobalVars::instance())->mode) {
		case char('fw') :
			fnew = new Derivable(pow(f1, f2.val), pow(f1, f2.val)*log(f1)*(f2.deriv));
			//(*fnew->symb_deriv) = f2*(f1 ^ (f2 - 1))*f2.flag;
			break;
		case char('adj') :
			double value = pow(f1, (f2._nodeOper->val));
			f1_new = new Derivable(f1);
			fnew = new Derivable(value); // Creating a new derivable. New node = nullptr. Val = value
			(*fnew).configureNodeOper('Pow', *f1_new, f2);
			(*fnew).pushToTape(GlobalVars::instance()->globalTape);
			break;
		}
		return (*fnew);
	}
	Derivable operator^(Derivable f1, double f2) {
		Derivable* fnew = nullptr;
		Derivable* f2_new = nullptr;
		switch ((GlobalVars::instance())->mode) {
		case char('fw') :
			fnew = new Derivable(pow(f1.val, f2), f2*pow(f1.val, f2 - 1)*f1.deriv);
			break;
		case char('adj') :
			double value = pow((f1._nodeOper->val), f2);
			f2_new = new Derivable(f2);
			fnew = new Derivable(value); // Creating a new derivable. New node = nullptr. Val = value
			(*fnew).configureNodeOper('Pow', f1, *f2_new);
			(*fnew).pushToTape(GlobalVars::instance()->globalTape);
			break;
		}
		return (*fnew);
	}
	vector<Derivable> operator^(vector<Derivable> f1, vector<Derivable> f2) {
		vector<Derivable>* fnew = new vector<Derivable>();
		int sz = f1.size();
		for (int k = 0; k < sz; k++) {
			fnew->push_back(*(new Derivable(pow(f1[k].val, f2[k].val), pow(f1[k].val, f2[k].val)*log(f1[k].val)*f2[k].deriv + pow(f1[k].val, f2[k].val - 1)*f2[k].val*f1[k].deriv)));
		}
		return (*fnew);
	}
	vector<Derivable> operator^(double f1, vector<Derivable> f2) {
		vector<Derivable>* fnew = new vector<Derivable>();
		int sz = f2.size();
		for (int k = 0; k < sz; k++) {
			fnew->push_back(*(new Derivable(pow(f1, f2[k].val), pow(f1, f2[k].val)*log(f1)*f2[k].deriv)));
		}
		return (*fnew);
	}
	vector<Derivable> operator^(vector<Derivable> f1, double f2) {
		vector<Derivable>* fnew = new vector<Derivable>();
		int sz = f1.size();
		for (int k = 0; k < sz; k++) {
			fnew->push_back(*(new Derivable(pow(f1[k].val, f2), f2*pow(f1[k].val, f2 - 1)*f1[k].deriv)));
		}
		return (*fnew);
	}

	// exp for simple and vector modes
	Derivable exponent(Derivable f) {
		Derivable* fnew = nullptr;
		switch ((GlobalVars::instance())->mode) {
		case char('fw') :
			fnew = new Derivable(exp(f.val), exp(f.val)*(f.deriv));
			//(*fnew->symb_deriv) = exponent(*f.symb_deriv)*f.flag;
			break;
		case char('adj') :
			double value = exp(f._nodeOper->val);
			fnew = new Derivable(value);
			(*fnew).configureNodeOper('Exp', f);
			(*fnew).pushToTape(GlobalVars::instance()->globalTape);
			break;
		}
		return (*fnew);
	}
	vector<Derivable> exponent(vector<Derivable> f) {
		vector<Derivable>* fnew = new vector<Derivable>();
		int sz = f.size();
		for (int k = 0; k < sz; k++) {
			fnew->push_back(*(new Derivable(exp(f[k].val), exp(f[k].val)*(f[k].deriv))));
		}
		return (*fnew);
	}
	// abs for simple and vector modes
	Derivable abs(Derivable f) {

		double val = f.val;
		double value;
		double deriv;
		if (val >= 0) {
			value = val;
			deriv = f.deriv;
		}
		else {
			value = -val;
			deriv = -f.deriv;
		}

		Derivable* fnew = nullptr;
		switch ((GlobalVars::instance())->mode) {
		case char('fw') :
			fnew = new Derivable(value, deriv);
			//(*fnew->symb_deriv) = f.flag;
			break;
		case char('adj') :
			fnew = new Derivable(value);
			(*fnew).configureNodeOper('Abso', f);
			(*fnew).pushToTape(GlobalVars::instance()->globalTape);
			break;
		}
		return (*fnew);
	}
	vector<Derivable> abs(vector<Derivable> f) {
		vector<Derivable>* fnew = new vector<Derivable>();
		int sz = f.size();
		for (int k = 0; k < sz; k++) {
			double val = f[k].val;
			double value;
			double deriv;
			if (val >= 0) {
				value = val;
				deriv = f[k].deriv;
			}
			else {
				value = -val;
				deriv = -f[k].deriv;
			}
			fnew->push_back(*(new Derivable(value, deriv)));
		}
		return (*fnew);
	}
	// log for simple and vector modes
	Derivable logarithm(Derivable f) {
		Derivable* fnew = nullptr;
		switch ((GlobalVars::instance())->mode) {
		case char('fw') :
			fnew = new Derivable(log(f.val), (1 / f.val)*(f.deriv));
			//(*fnew->symb_deriv) = (1 / (*f.symb_deriv))*f.flag;
			break;
		case char('adj') :
			double value = log(f._nodeOper->val);
			fnew = new Derivable(value);
			(*fnew).configureNodeOper('Log', f);
			(*fnew).pushToTape(GlobalVars::instance()->globalTape);
			break;
		}
		return (*fnew);
	}
	vector<Derivable> logarithm(vector<Derivable> f) {
		vector<Derivable>* fnew = new vector<Derivable>();
		int sz = f.size();
		for (int k = 0; k < sz; k++) {
			fnew->push_back(*(new Derivable(log(f[k].val), (1 / f[k].val)*f[k].deriv)));
		}
		return (*fnew);
	}

	// Operators overloading (complex)
	Derivable Normal_ElHadj(Derivable u) {
		Derivable y;
		Derivable zinterm;
		Derivable z;
		Derivable num;
		Derivable den;
		switch ((GlobalVars::instance())->mode) {
		case char('fw') :
			y = abs(u);
			if (y.val > 35.0) {

				if (u.val > 0) {
					return 1;
				}
				else {
					return 0;
				}
			}
			else if (y.val <= 0.662912607) {
				//  evaluate erf() for |u| <= sqrt(2)*0.46875
				double a0 = 1.161110663653770e-2;
				double a1 = 3.951404679838207e-1;
				double a2 = 2.846603853776254e+1;
				double a3 = 1.887426188426510e+2;
				double a4 = 3.209377589138469e+3;

				double b0 = 1.767766952966369e-1;
				double b1 = 8.344316438579620;
				double b2 = 1.725514762600375e+2;
				double b3 = 1.813893686502485e+3;
				double b4 = 8.044716608901563e+3;

				Derivable z = y * y;
				y = u * ((((a0 * z + a1) * z + a2) * z + a3) * z + a4);
				y = y / ((((b0 * z + b1) * z + b2) * z + b3) * z + b4);
				return 0.5 + y; // y.val
			}
			zinterm = 0.5 * exponent((0 - y) * y / 2);
			if (y.val <= 4.0) {
				double c0 = 2.15311535474403846e-8;
				double c1 = 5.64188496988670089e-1;
				double c2 = 8.88314979438837594;
				double c3 = 6.61191906371416295e+1;
				double c4 = 2.98635138197400131e+2;
				double c5 = 8.81952221241769090e+2;
				double c6 = 1.71204761263407058e+3;
				double c7 = 2.05107837782607147e+3;
				double c8 = 1.23033935479799725e+3;

				double d0 = 1.0;
				double d1 = 1.57449261107098347e+1;
				double d2 = 1.17693950891312499e+2;
				double d3 = 5.37181101862009858e+2;
				double d4 = 1.62138957456669019e+3;
				double d5 = 3.29079923573345963e+3;
				double d6 = 4.36261909014324716e+3;
				double d7 = 3.43936767414372164e+3;
				double d8 = 1.23033935480374942e+3;

				// evaluate erfc() for sqrt(2)*0.46875 <= |u| <= sqrt(2)*4.0
				y = y / 1.4142135623730950488;
				num = ((((((((c0 * y + c1) * y + c2) * y + c3) * y + c4) * y + c5) * y + c6) * y + c7) * y + c8);
				den = ((((((((d0 * y + d1) * y + d2) * y + d3) * y + d4) * y + d5) * y + d6) * y + d7) * y + d8);

				y = num / den;
				y = zinterm * y;
			}
			else {
				double p0 = 1.63153871373020978e-2;
				double p1 = 3.05326634961232344e-1;
				double p2 = 3.60344899949804439e-1;
				double p3 = 1.25781726111229246e-1;
				double p4 = 1.60837851487422766e-2;
				double p5 = 6.58749161529837803e-4;

				double q0 = 1.00000000000000000;
				double q1 = 2.56852019228982242;
				double q2 = 1.87295284992346047;
				double q3 = 5.27905102951428412e-1;
				double q4 = 6.05183413124413191e-2;
				double q5 = 2.33520497626869185e-3;
				// evaluate erfc() for |u| > sqrt(2)*4.0
				z = zinterm * 1.41421356237309504880 / y;
				y = 2 / (y * y);
				y = y * (((((p0 * y + p1) * y + p2) * y + p3) * y + p4) * y + p5) / (((((q0 * y + q1) * y + q2) * y + q3) * y + q4) * y + q5);
				y = z*(0.564189583547756287 - y);
			}

			if (u.val < 0) {
				return y;
			}
			else {
				return 1 - y;
			}
			break;
		case char('adj') :
			y = abs(u);
			if (y.getNodeAddress()->val > 35.0) {

				if (u.getNodeAddress()->val > 0) {
					return 1;
				}
				else {
					return 0;
				}
			}
			else if (y.getNodeAddress()->val <= 0.662912607) {
				//  evaluate erf() for |u| <= sqrt(2)*0.46875
				double a0 = 1.161110663653770e-2;
				double a1 = 3.951404679838207e-1;
				double a2 = 2.846603853776254e+1;
				double a3 = 1.887426188426510e+2;
				double a4 = 3.209377589138469e+3;

				double b0 = 1.767766952966369e-1;
				double b1 = 8.344316438579620;
				double b2 = 1.725514762600375e+2;
				double b3 = 1.813893686502485e+3;
				double b4 = 8.044716608901563e+3;

				z = y * y;
				y = u * ((((a0 * z + a1) * z + a2) * z + a3) * z + a4);
				y = y / ((((b0 * z + b1) * z + b2) * z + b3) * z + b4);
				return 0.5 + y;// .getNodeAddress()->val;
			}
			zinterm = 0.5 * exponent((0 - y) * y / 2);
			if (y.getNodeAddress()->val <= 4.0) {
				double c0 = 2.15311535474403846e-8;
				double c1 = 5.64188496988670089e-1;
				double c2 = 8.88314979438837594;
				double c3 = 6.61191906371416295e+1;
				double c4 = 2.98635138197400131e+2;
				double c5 = 8.81952221241769090e+2;
				double c6 = 1.71204761263407058e+3;
				double c7 = 2.05107837782607147e+3;
				double c8 = 1.23033935479799725e+3;

				double d0 = 1.0;
				double d1 = 1.57449261107098347e+1;
				double d2 = 1.17693950891312499e+2;
				double d3 = 5.37181101862009858e+2;
				double d4 = 1.62138957456669019e+3;
				double d5 = 3.29079923573345963e+3;
				double d6 = 4.36261909014324716e+3;
				double d7 = 3.43936767414372164e+3;
				double d8 = 1.23033935480374942e+3;

				// evaluate erfc() for sqrt(2)*0.46875 <= |u| <= sqrt(2)*4.0
				y = y / 1.4142135623730950488;
				num = ((((((((c0 * y + c1) * y + c2) * y + c3) * y + c4) * y + c5) * y + c6) * y + c7) * y + c8);
				den = ((((((((d0 * y + d1) * y + d2) * y + d3) * y + d4) * y + d5) * y + d6) * y + d7) * y + d8);

				y = num / den;
				y = zinterm * y;
			}
			else {
				double p0 = 1.63153871373020978e-2;
				double p1 = 3.05326634961232344e-1;
				double p2 = 3.60344899949804439e-1;
				double p3 = 1.25781726111229246e-1;
				double p4 = 1.60837851487422766e-2;
				double p5 = 6.58749161529837803e-4;

				double q0 = 1.00000000000000000;
				double q1 = 2.56852019228982242;
				double q2 = 1.87295284992346047;
				double q3 = 5.27905102951428412e-1;
				double q4 = 6.05183413124413191e-2;
				double q5 = 2.33520497626869185e-3;
				// evaluate erfc() for |u| > sqrt(2)*4.0
				z = zinterm * 1.41421356237309504880 / y;
				y = 2 / (y * y);
				y = y * (((((p0 * y + p1) * y + p2) * y + p3) * y + p4) * y + p5) / (((((q0 * y + q1) * y + q2) * y + q3) * y + q4) * y + q5);
				y = z*(0.564189583547756287 - y);
			}

			if (u.getNodeAddress()->val < 0) {
				return y;
			}
			else {
				return 1 - y;
			}
			break;
		}
	}
	vector<Derivable> vector_Normal_ElHadj(vector<Derivable> u) {
		int sz = u.size();
		vector<Derivable>* y = new vector<Derivable>(sz, 0);
		vector<Derivable>* zinterm = new vector<Derivable>(sz, 0);
		vector<Derivable>* z = new vector<Derivable>(sz, 0);
		vector<Derivable>* num = new vector<Derivable>(sz, 0);
		vector<Derivable>* den = new vector<Derivable>(sz, 0);
		vector<Derivable>* result = new vector<Derivable>(sz, 0);
		for (int k = 0; k < sz; k++) {
			*y = abs(u);
			if ((*y)[k].val > 35.0) {

				if (u[k].val > 0) {
					(*result)[k] = 1;
					continue;
				}
				else {
					(*result)[k] = 0;
					continue;
				}
			}
			else if ((*y)[k].val <= 0.662912607) {
				//  evaluate erf() for |u| <= sqrt(2)*0.46875
				double a0 = 1.161110663653770e-2;
				double a1 = 3.951404679838207e-1;
				double a2 = 2.846603853776254e+1;
				double a3 = 1.887426188426510e+2;
				double a4 = 3.209377589138469e+3;

				double b0 = 1.767766952966369e-1;
				double b1 = 8.344316438579620;
				double b2 = 1.725514762600375e+2;
				double b3 = 1.813893686502485e+3;
				double b4 = 8.044716608901563e+3;

				(*z)[k] = (*y)[k] * (*y)[k];
				(*y)[k] = u[k] * ((((a0 * (*z)[k] + a1) * (*z)[k] + a2) * (*z)[k] + a3) * (*z)[k] + a4);
				(*y)[k] = (*y)[k] / ((((b0 * (*z)[k] + b1) * (*z)[k] + b2) * (*z)[k] + b3) * (*z)[k] + b4);
				(*result)[k] = (0.5 + (*y)[k]);
				continue;
			}
			(*zinterm)[k] = 0.5 * exponent((0 - (*y)[k]) * (*y)[k] / 2);
			if ((*y)[k].val <= 4.0) {
				double c0 = 2.15311535474403846e-8;
				double c1 = 5.64188496988670089e-1;
				double c2 = 8.88314979438837594;
				double c3 = 6.61191906371416295e+1;
				double c4 = 2.98635138197400131e+2;
				double c5 = 8.81952221241769090e+2;
				double c6 = 1.71204761263407058e+3;
				double c7 = 2.05107837782607147e+3;
				double c8 = 1.23033935479799725e+3;

				double d0 = 1.0;
				double d1 = 1.57449261107098347e+1;
				double d2 = 1.17693950891312499e+2;
				double d3 = 5.37181101862009858e+2;
				double d4 = 1.62138957456669019e+3;
				double d5 = 3.29079923573345963e+3;
				double d6 = 4.36261909014324716e+3;
				double d7 = 3.43936767414372164e+3;
				double d8 = 1.23033935480374942e+3;

				// evaluate erfc() for sqrt(2)*0.46875 <= |u| <= sqrt(2)*4.0
				(*y)[k] = (*y)[k] / 1.4142135623730950488;
				(*num)[k] = ((((((((c0 * (*y)[k] + c1) * (*y)[k] + c2) * (*y)[k] + c3) * (*y)[k] + c4) * (*y)[k] + c5) * (*y)[k] + c6) * (*y)[k] + c7) * (*y)[k] + c8);
				(*den)[k] = ((((((((d0 * (*y)[k] + d1) * (*y)[k] + d2) * (*y)[k] + d3) * (*y)[k] + d4) * (*y)[k] + d5) * (*y)[k] + d6) * (*y)[k] + d7) * (*y)[k] + d8);

				(*y)[k] = (*num)[k] / (*den)[k];
				(*y)[k] = (*zinterm)[k] * (*y)[k];
			}
			else {
				double p0 = 1.63153871373020978e-2;
				double p1 = 3.05326634961232344e-1;
				double p2 = 3.60344899949804439e-1;
				double p3 = 1.25781726111229246e-1;
				double p4 = 1.60837851487422766e-2;
				double p5 = 6.58749161529837803e-4;

				double q0 = 1.00000000000000000;
				double q1 = 2.56852019228982242;
				double q2 = 1.87295284992346047;
				double q3 = 5.27905102951428412e-1;
				double q4 = 6.05183413124413191e-2;
				double q5 = 2.33520497626869185e-3;
				// evaluate erfc() for |u| > sqrt(2)*4.0
				(*z)[k] = (*zinterm)[k] * 1.41421356237309504880 / (*y)[k];
				(*y)[k] = 2 / ((*y)[k] * (*y)[k]);
				(*y)[k] = (*y)[k] * (((((p0 * (*y)[k] + p1) * (*y)[k] + p2) * (*y)[k] + p3) * (*y)[k] + p4) * (*y)[k] + p5) / (((((q0 * (*y)[k] + q1) * (*y)[k] + q2) * (*y)[k] + q3) * (*y)[k] + q4) * (*y)[k] + q5);
				(*y)[k] = (*z)[k] * (0.564189583547756287 - (*y)[k]);
			}

			if (u[k].val < 0) {
				(*result)[k] = (*y)[k];
				continue;
			}
			else {
				(*result)[k] = (1 - (*y)[k]);
				continue;
			}
		}
		return (*result);
	}

	Derivable sqrt(Derivable f) {
		Derivable f2 = f ^ 0.5;
		return f2;
	}
	vector<Derivable> sqrt(vector<Derivable> f) {
		vector<Derivable>* fnew = new vector<Derivable>();
		int sz = f.size();
		for (int k = 0; k < sz; k++) {
			Derivable f2 = (f[k]) ^ 0.5;
			fnew->push_back(f2);
		}
		return (*fnew);
	}

	Derivable sqr(Derivable f) {
		Derivable f2 = f ^ 2;
		return f2;
	}
	vector<Derivable> sqr(vector<Derivable> f) {
		vector<Derivable>* fnew = new vector<Derivable>();
		int sz = f.size();
		for (int k = 0; k < sz; k++) {
			Derivable f2 = f[k] ^ 2;
			fnew->push_back(f2);
		}
		return (*fnew);
	}

	double factorial(double a) {
		double ans = 1;
		if (a != 0) {
			for (int i = 1; i <= a; i++) {
				ans = ans*i;
			}
		}
		else {
			ans = 1;
		}
		return ans;
	}
	double Normal_ElHadj(double u) {
		double y = std::abs(u);
		if (y > 35.0) {
			if (u > 0)
				return 1;
			else
				return 0;
		}
		if (y <= 0.662912607) {
			//  evaluate erf() for |u| <= sqrt(2)*0.46875
			double a0 = 1.161110663653770e-2;
			double a1 = 3.951404679838207e-1;
			double a2 = 2.846603853776254e+1;
			double a3 = 1.887426188426510e+2;
			double a4 = 3.209377589138469e+3;

			double b0 = 1.767766952966369e-1;
			double b1 = 8.344316438579620;
			double b2 = 1.725514762600375e+2;
			double b3 = 1.813893686502485e+3;
			double b4 = 8.044716608901563e+3;

			double z = y * y;
			y = u * ((((a0 * z + a1) * z + a2) * z + a3) * z + a4);
			y /= ((((b0 * z + b1) * z + b2) * z + b3) * z + b4);
			return 0.5 + y;
		}
		double zinterm = 0.5 * exp(-y * y / 2);
		if (y <= 4.0) {
			double c0 = 2.15311535474403846e-8;
			double c1 = 5.64188496988670089e-1;
			double c2 = 8.88314979438837594;
			double c3 = 6.61191906371416295e+1;
			double c4 = 2.98635138197400131e+2;
			double c5 = 8.81952221241769090e+2;
			double c6 = 1.71204761263407058e+3;
			double c7 = 2.05107837782607147e+3;
			double c8 = 1.23033935479799725e+3;

			double d0 = 1.0;
			double d1 = 1.57449261107098347e+1;
			double d2 = 1.17693950891312499e+2;
			double d3 = 5.37181101862009858e+2;
			double d4 = 1.62138957456669019e+3;
			double d5 = 3.29079923573345963e+3;
			double d6 = 4.36261909014324716e+3;
			double d7 = 3.43936767414372164e+3;
			double d8 = 1.23033935480374942e+3;

			// evaluate erfc() for sqrt(2)*0.46875 <= |u| <= sqrt(2)*4.0
			y = y / 1.4142135623730950488;
			double num = ((((((((c0 * y + c1) * y + c2) * y + c3) * y + c4) * y + c5) * y + c6) * y + c7) * y + c8);
			double den = ((((((((d0 * y + d1) * y + d2) * y + d3) * y + d4) * y + d5) * y + d6) * y + d7) * y + d8);

			y = num / den;
			y = zinterm * y;
		}
		else {
			double p0 = 1.63153871373020978e-2;
			double p1 = 3.05326634961232344e-1;
			double p2 = 3.60344899949804439e-1;
			double p3 = 1.25781726111229246e-1;
			double p4 = 1.60837851487422766e-2;
			double p5 = 6.58749161529837803e-4;

			double q0 = 1.00000000000000000;
			double q1 = 2.56852019228982242;
			double q2 = 1.87295284992346047;
			double q3 = 5.27905102951428412e-1;
			double q4 = 6.05183413124413191e-2;
			double q5 = 2.33520497626869185e-3;
			// evaluate erfc() for |u| > sqrt(2)*4.0
			double z = zinterm * 1.41421356237309504880 / y;
			y = 2 / (y * y);
			y = y * (((((p0 * y + p1) * y + p2) * y + p3) * y + p4) * y + p5) / (((((q0 * y + q1) * y + q2) * y + q3) * y + q4) * y + q5);
			y = z*(0.564189583547756287 - y);
		}

		if (u < 0)
			return y;
		return 1 - y;
	}
















	void AAD(Derivable(*func)(vector<Derivable>), vector<Derivable> vars) {
		vector<NodeOper*>* tape = GlobalVars::instance()->globalTape;
		GlobalVars::instance()->mode = char('adj'); // sets adj global mode
		Derivable function = func(vars); // fills the tape in fw mode
		int sz = tape->size();

		if (sz != 0) {

			(*tape)[sz - 1]->adjoint = 1; // setting last adjoint to 1

			for (int i = sz - 1; i >= 0; i--) {
				char Type = ((*tape)[i])->Type;

				// Obtain pointers to rhs and lhs nodes, current operator node
				NodeOper* f1 = ((*tape)[i])->get_f1_Adress(); //pointer to f1
				NodeOper* f2 = ((*tape)[i])->get_f2_Adress(); // pointer to f2
				NodeOper* current = (*tape)[i];

				// Initialize empty containers for partial derivs and new tape that will be used in recursion
				double derivf1 = 0;
				double derivf2 = 0;

				// check the case of operator node
				switch (Type) {

				case char('Plus') :
					derivf1 = 1;
					derivf2 = 1;
					break;

				case char('Min') :
					derivf1 = 1;
					derivf2 = -1;
					break;

				case char('Mult') :
					derivf1 = (f2->val); // here it must be derivable f2 instead of just value or casted double
					derivf2 = (f1->val); // deriv1 becomes a derivable sequence of operations to -> it has to has its own tape (has it be assigned explicitly?)
					break;

				case char('Div') :
					// Problem is here . it will be always 1/f2 because it does not know that 1 is constant
					derivf1 = 1 / (f2->val);
					derivf2 = (-f1->val) / pow((f2->val), 2);
					break;

				case char('Pow') :
					derivf1 = (f2->val)*pow(f1->val, f2->val - 1);
					derivf2 = pow(f1->val, f2->val)*log(f1->val);
					break;

				case char('Exp') :
					derivf1 = exp(f1->val);
					break;

				case char('Log') :
					derivf1 = 1 / (f1->val);
					break;

				case char('Abso') :
					if (f1 != nullptr) {
						if (f1->val >= 0) {
							derivf1 = 1;
						}
						else {
							derivf1 = -1;
						}
					}
					break;
				}

				// Assigning the new adjoint1s for  f1, f2

				if (f1 != nullptr) {
					(f1->adjoint) = (f1->adjoint) + (current->adjoint) * derivf1; // numerical adjoint
				}
				if (f2 != nullptr) {
					(f2->adjoint) = (f2->adjoint) + (current->adjoint) * derivf2;
				}

			}

		}
		
		/*cout << "value = " << function.val << endl;
		sz = vars.size();
		for (int i = 0; i < sz; i++) {
			cout << "derivative = " << vars[i].getNodeAddress()->adjoint << endl;
		}
		*/
	}

	Derivable FMAD(Derivable(*f)(vector<Derivable>, vector<double>), vector<Derivable> vec, vector<double> flag) {
		GlobalVars::instance()->mode = char('fw');
		Derivable result = f(vec, flag);
		return result;
	}

	vector<Derivable> vector_FMAD(vector<Derivable>(*f)(vector<vector<Derivable>>, vector<double>), vector<vector<Derivable>> vec, vector<double> flag) {
		//GlobalVars::instance()->mode = char('vfw');
		vector<Derivable> result = f(vec, flag);
		/*
		int sz = result.size();
		for (int i = 0; i < sz; i++) {
			cout << result[i].val << endl;
		}
		for (int i = 0; i < sz; i++) {
			cout << result[i].deriv << endl;
		}*/
		return result;
	}

	double delta_FDM(const double S, const double sigma, const double r, const double K, const double T, const double t, const double delta_S) {
		return (BS(S + delta_S, sigma, r, K, T, t) - BS(S, sigma, r, K, T, t)) / delta_S;
	}

	double vega_FDM(const double S, const double sigma, const double r, const double K, const double T, const double t, const double delta_sigma) {
		return (BS(S, sigma+delta_sigma, r, K, T, t) - BS(S, sigma, r, K, T, t)) / delta_sigma;
	}

	double rho_FDM(const double S, const double sigma, const double r, const double K, const double T, const double t, const double delta_r) {
		return (BS(S, sigma, r+delta_r, K, T, t) - BS(S, sigma, r, K, T, t)) / delta_r;
	}

	double theta_FDM(const double S, const double sigma, const double r, const double K, const double T, const double t, const double delta_T) {
		return -(BS(S, sigma, r, K, T+delta_T, t) - BS(S, sigma, r, K, T, t)) / delta_T;
	}

	double thetaa_FDM(const double S, const double sigma, const double r, const double K, const double T, const double t, const double delta_t) {
		return -(BS(S, sigma, r, K, T, t+delta_t) - BS(S, sigma, r, K, T, t)) / delta_t;
	}

	double kaa_FDM(const double S, const double sigma, const double r, const double K, const double T, const double t, const double delta_k) {
		return -(BS(S, sigma, r, K+delta_k, T, t) - BS(S, sigma, r, K, T, t)) / delta_k;
	}
	
	double test_complex_FDM(vector<double> vec, int der_position, double delta) {
		//vector<double> vec1 = vec;
		vec[der_position] = vec[der_position] + delta;

		double ans = (test_complex(vec) - 1) / delta;
		return ans;
	}

	double BS(double S, double sigma, double r, double K, double T, double t) {
		double d1 = (1 / (sigma * pow(T - t,0.5)))*(log(S / K) + (r + sigma*sigma)*(T - t));
		//cout << "d1 = " << d1.getValue() << endl;

		double d2 = d1 - sigma * pow(T - t, 0.5);
		//cout << "d2 = " << d2.getValue() << endl;

		double func = (Normal_ElHadj(d1)*S) - (Normal_ElHadj(d2)*K * exp((-r)*(T - t)));

		return func;
	}

	double test_complex(vector<double> vec) {
		int sz = vec.size();
		double func = vec[0];
		for (int i = 0; i < sz-1; i++) {
			func = func * vec[i + 1];
		}
		//double func = exp(vec[1])*vec[0] + vec[2] * log(vec[3])*std::abs(vec[4])*vec[5] * (vec[6] * vec[6])*vec[7] / vec[8] + vec[8] * log(vec[6])*vec[9] * vec[10] / vec[7] + exp(vec[10])*vec[9] * vec[3] * vec[11] * vec[12] * vec[13] / vec[14] * vec[5] * exp(vec[11]) + vec[14];
		return func;
	}

	


	Derivable test(vector<Derivable> vec, vector<double> flag) {
		int sz = vec.size();
		for (int i = 0; i < sz; i++) {
			vec[i].flag = flag[i];
			if (flag[i] == 1) {
				vec[i].deriv = 1;
			}
		}

		//Derivable func = Normal_ElHadj(vec[1]);
		Derivable func = exponent(vec[1])*vec[0] + vec[1];
		return func;

	}

	Derivable test(vector<Derivable> vec) {

		Derivable func = exponent(vec[1])*vec[0] + vec[1];

		return func;
	}

	Derivable test_complex(vector<Derivable> vec, vector<double> flag) {
		int sz = vec.size();
		for (int i = 0; i < sz; i++) {
			vec[i].flag = flag[i];
			if (flag[i] == 1) {
				vec[i].deriv = 1;
			}
		}

		Derivable func = vec[0];
		for (int i = 0; i < sz-1; i++) {
			func = func * vec[i + 1];
		}
		//Derivable func = exponent(vec[1])*vec[0] + vec[2]*logarithm(vec[3])*abs(vec[4])*vec[5]*(vec[6]^2)*vec[7]/vec[8]+vec[8]*logarithm(vec[6])*vec[9]*vec[10]/vec[7]+exponent(vec[10])*vec[9]*vec[3]*vec[11]*vec[12]*vec[13]/vec[14]*vec[5]*exponent(vec[11])+vec[14];
		return func;

	}

	Derivable test_complex(vector<Derivable> vec) {
		int sz = vec.size();
		Derivable func = vec[0];
		for (int i = 0; i < sz-1; i++) {
			func = func * vec[i + 1];
			//Derivable func = exponent(vec[1])*vec[0] + vec[2] * logarithm(vec[3])*abs(vec[4])*vec[5] * (vec[6] ^ 2)*vec[7] / vec[8] + vec[8] * logarithm(vec[6])*vec[9] * vec[10] / vec[7] + exponent(vec[10])*vec[9] * vec[3] * vec[11] * vec[12] * vec[13] / vec[14] * vec[5] * exponent(vec[11]) + vec[14];
		}
		return func;
	}

	Derivable BS(vector<Derivable> vec, vector<double> flag) {

		int sz = vec.size();
		for (int i = 0; i < sz; i++) {
			vec[i].flag = flag[i];
			if (flag[i] == 1) {
				vec[i].deriv = 1;
			}
		}


		Derivable d1 = (1 / (vec[1] * sqrt(vec[4] - vec[5])))*(logarithm(vec[0] / vec[3]) + (vec[2] + (vec[1] * vec[1] / 2))*(vec[4] - vec[5]));
		//cout << "d1 = " << d1.getValue() << endl;

		Derivable d2 = d1 - vec[1] * sqrt(vec[4] - vec[5]);
		//cout << "d2 = " << d2.getValue() << endl;

		Derivable func = (Normal_ElHadj(d1)*vec[0]) - (Normal_ElHadj(d2)*vec[3]*exponent((0 - vec[2])*(vec[4] - vec[5])));

		return func;
	}

	Derivable BS(vector<Derivable> vec) {

		Derivable d1 = (1 / (vec[1] * sqrt(vec[4] - vec[5])))*(logarithm(vec[0] / vec[3]) + (vec[2] + (sqr(vec[1])/2))*(vec[4] - vec[5]));
		//cout << "d1 = " << d1.getValue() << endl;

		Derivable d2 = d1 - vec[1] * sqrt(vec[4] - vec[5]);
		//cout << "d2 = " << d2.getValue() << endl;

		Derivable func = (Normal_ElHadj(d1)*vec[0]) - (Normal_ElHadj(d2)*vec[3] * exponent((0 - vec[2])*(vec[4] - vec[5])));

		return func;
	}

	Derivable VarReduction(vector<Derivable> vec, vector<double> flag) {

		int sz = vec.size();
		for (int i = 0; i < sz; i++) {
			vec[i].flag = flag[i];
			if (flag[i] == 1) {
				vec[i].deriv = 1;
			}
		}

		Derivable P = Normal_ElHadj((logarithm(vec[3] / vec[0]) - (vec[2] + (vec[1] ^ 2) / 2)*(vec[4] - vec[5])) / (vec[1] * sqrt(vec[4] - vec[5])));

		Derivable Q = Normal_ElHadj((logarithm(vec[3] / vec[0]) - (vec[2] - (vec[1] ^ 2) / 2)*(vec[4] - vec[5])) / (vec[1] * sqrt(vec[4] - vec[5])));

		Derivable func = vec[0] * (1 - 2 * P)*(exponent((vec[4] - vec[5])*(vec[1] ^ 2)) - 1) - 2 * vec[3] * (Q - P);


		return func;
	}

	Derivable VarReduction(vector<Derivable> vec) {

		Derivable P = Normal_ElHadj((logarithm(vec[3] / vec[0]) - (vec[2] + (vec[1] ^ 2) / 2)*(vec[4] - vec[5])) / (vec[1] * sqrt(vec[4] - vec[5])));

		Derivable Q = Normal_ElHadj((logarithm(vec[3] / vec[0]) - (vec[2] - (vec[1] ^ 2) / 2)*(vec[4] - vec[5])) / (vec[1] * sqrt(vec[4] - vec[5])));

		Derivable func = vec[0] * (1 - 2 * P)*(exponent((vec[4] - vec[5])*vec[1] ^ 2) - 1) - 2 * vec[3] * (Q - P);


		return func;
	}

	vector<Derivable> vectorBS(vector<vector<Derivable>> vec, vector<double> flag) {

		int sz = vec.size();
		for (int i = 0; i < sz; i++) {
			int sz2 = vec[i].size();
			for (int k = 0; k<sz2; k++) {
				vec[i][k].flag = flag[i];
				if (flag[i] == 1) {
					vec[i][k].deriv = 1;
				}
			}
		}


		vector<Derivable> d1 = (1 / (vec[1] * sqrt(vec[4] - vec[5])))*(logarithm(vec[0] / vec[3]) + (vec[2] + sqr(vec[1]) / 2)*(vec[4] - vec[5]));
		//cout << "d1 = " << d1.getValue() << endl;

		vector<Derivable> d2 = d1 - vec[1] * sqrt(vec[4] - vec[5]);
		//cout << "d2 = " << d2.getValue() << endl;

		vector<Derivable> func = (vector_Normal_ElHadj(d1)*vec[0]) - (vector_Normal_ElHadj(d2)*vec[3] * exponent(0 - vec[2])*(vec[4] - vec[5]));

		return func;
	}



}