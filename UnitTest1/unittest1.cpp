#include "stdafx.h"
#include "CppUnitTest.h"
#include "AD_Project.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;
using namespace ADUTILITIES;

namespace UnitTest1
{		
	TEST_CLASS(UnitTest1)
	{
	public:
		// test of simple function
		TEST_METHOD(TestEvalSimple) {
			double func = exp(2)*1 + 2;
		}

		TEST_METHOD(FDMSimple) {
			double da = ((exp(2)*(1 + 0.01) + 2) - (exp(2) * 1 + 2)) / 0.01;
			double db = ((exp(2+0.01)*1 + 2+0.01) - (exp(2) * 1 + 2)) / 0.01;
		}

		TEST_METHOD(TestFMADSimple) {
			Derivable a(1);
			Derivable b(2);
			Derivable result1 = FMAD(test, { a,b }, { 0 , 1});
			Derivable result2 = FMAD(test, { a,b }, { 1 , 0 });
			//double v = result.val;
			//double der = result.deriv;

			//system("pause");
		}

		TEST_METHOD(TestAADSimple) {
			Derivable a(1);
			Derivable b(2);
			AAD(test, { a,b });

			//double v = result.val;
			//double der = result.deriv;

			//system("pause");
		}

		// Test of functions initial
		TEST_METHOD(TestFMAD)
		{
			// TODO: Your test code here
			Derivable S = Derivable(17);
			Derivable sigma = Derivable(0.4);
			Derivable r = Derivable(0.04);
			Derivable K = Derivable(18);
			Derivable T = Derivable(1);
			Derivable t = Derivable(0);

			Derivable result = FMAD(BS, { S,sigma,r,K,T,t }, { 0,1,0,0,0,0 });
			double v = result.val;
			double der = result.deriv;
			//AAD(VarReduction, { S, sigma, r, K, T, t });


			//system("pause");


		}

		TEST_METHOD(TestVectorFMADBS)
		{
			// TODO: Your test code here

			vector<Derivable> _S = { Derivable(16), Derivable(22) };
			vector<Derivable> _sigma = { Derivable(0.1), Derivable(0.1) };
			vector<Derivable> _r = { Derivable(0.04),  Derivable(0.04) };
			vector<Derivable> _K = { Derivable(18),  Derivable(18) };
			vector<Derivable> _T = { Derivable(1), Derivable(1) };
			vector<Derivable> _t = { Derivable(0), Derivable(0) };

			// permits to obtain a greek graph in all the points simultaneously for BS
			vector<Derivable> result = vector_FMAD(vectorBS, { _S,_sigma,_r,_K,_T,_t }, {1,0,0,0,0,0});
	
		//	system("pause");


		}

		// Testing for a POC
		TEST_METHOD(TestEvalBS) {
			double S = 17;
			double sigma = 0.4;
			double r = 0.04;
			double K = 18;
			double T = 1;
			double t = 0;
			double res = BS(S, sigma, r, K, T, t);
		}

		TEST_METHOD(TestFMADForAllGreeks)
		{
			// TODO: Your test code here
			Derivable S = Derivable(17);
			Derivable sigma = Derivable(0.4);
			Derivable r = Derivable(0.04);
			Derivable K = Derivable(18);
			Derivable T = Derivable(1);
			Derivable t = Derivable(0);

			//vector<double> der(6, 0);

			for (int i = 0; i < 6; i++) {
				vector<double> flag(6, 0);
				flag[i] = 1;

				Derivable result = FMAD(BS, { S,sigma,r,K,T,t }, flag);
				//der[i] = result.deriv;
			}
		
			//system("pause");


		}

		TEST_METHOD(TestAAD)
		{
			// TODO: Your test code here
			Derivable S = Derivable(17);
			Derivable sigma = Derivable(0.4);
			Derivable r = Derivable(0.04);
			Derivable K = Derivable(18);
			Derivable T = Derivable(1);
			Derivable t = Derivable(0);

			AAD(BS, { S, sigma, r, K, T, t });
			//double d1 = S.getNodeAddress()->adjoint;
			//double d2 = sigma.getNodeAddress()->adjoint;
			//double d3 = r.getNodeAddress()->adjoint;
			//double d4 = K.getNodeAddress()->adjoint;
			//double d5 = T.getNodeAddress()->adjoint;
			//double d6 = t.getNodeAddress()->adjoint;
			//system("pause");


		}

		TEST_METHOD(TestFDMBSAllGreeks) {
			double S = 17;
			double sigma = 0.4;
			double r = 0.04;
			double K = 18;
			double T = 1;
			double t = 0;


			double delta = delta_FDM(S, sigma, r, K, T, t, S*0.01);
			double vega = vega_FDM(S, sigma, r, K, T, t, 0.01);
			double rho = rho_FDM(S, sigma, r, K, T, t, 0.01);
			double kaa = kaa_FDM(S, sigma, r, K, T, t, K*0.01);
			double theta = theta_FDM(S, sigma, r, K, T, t, T*0.01);
			double thetaa = thetaa_FDM(S, sigma, r, K, T, t, 0.01);

			//system("pause");
		}


	};
}