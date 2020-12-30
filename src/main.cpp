#include <stdio.h>
#include <array>
#include <complex>
#include <vector>
#include <cassert>
#include <limits>

#include "math.hpp"
#include "constants.hpp"
#include "opacity.hpp"

#define NDIM 3

const double h = 0.7;
const double omega_b = 0.05;
const double omega_c = 0.25;
const double Theta = 1.00;
const double omega_m = omega_b + omega_c;
//const double zeq = 2.5e+4 * std::pow(Theta, -4) * omega_m * h * h;
//const double aeq = 1.0 / (zeq + 1.0);
const double omega_gam = 5e-5;
const double omega_nu = 3.5e-5;
const double omega_r = omega_nu + omega_gam;
const double omega_lambda = 1.0 - omega_m - omega_r;
const double H0 = 100.0 * 1e5 / constants::clight;
const double ns = 1;
const double sigma_T = 2.3e-5;
const double hefrac = 0.24;
const double Yp = (1 - hefrac / 2.0);

double Hubble(double a) {
	return H0 * h * std::sqrt(omega_r / (a * a * a * a) + omega_m / (a * a * a) + omega_lambda);
}

#define NF 9

using state = std::array<double,NF>;

#define Theta0i 0
#define Theta1i 1
#define deltaci 2
#define ubi 3
#define deltabi 4
#define uci 5
#define Phii 6
#define N0i 7
#define N1i 8

double compute_Nele(double eta, double beta);
double compute_Npos(double eta, double beta);

double rho_b(double a) {
	using namespace constants;
	return (3.0 * omega_b * H0cgs * H0cgs * h * h) / (8.0 * M_PI * G * a * a * a);
}
double chemical_potential(double ne, double T) {
	using namespace constants;
	const auto beta = kb * T / (me * clight * clight);
	double eta = 0.0;
	double deta = 1.0e-7;
	double deta0;
	int iter = 0;
	double w = 1;
	do {
		const auto eta0 = eta - 0.5 * deta;
		const auto eta1 = eta + 0.5 * deta;
		const double f0 = ne - (compute_Nele(eta0, beta) - compute_Npos(eta0, beta));
		const double f1 = ne - (compute_Nele(eta1, beta) - compute_Npos(eta1, beta));
		const double dfdeta = (f1 - f0) / deta;
		const double f = 0.5 * (f0 + f1);
		deta0 = -f / dfdeta;
		eta += w * deta0;
		iter++;
		if (iter > 500) {
			printf("chem_pot : %e %e %e %e %e \n", ne, T, eta, deta0, w);
		}
		if (iter == 1000) {
			printf("maximum iterations exceeded in chemical_potential\n");
			abort();
		}
		w *= 0.98;
	} while (std::abs(deta0) > 1.0e-6);
	return eta;
}

void compute_chemistry(double rho, double T, double &H, double &Hp, double &He, double &Hep, double &Hepp, double &ne) {
	using namespace constants;
	constexpr auto evtoerg = 1.6021772e-12;
	constexpr auto eps_1_H = 13.59844 * evtoerg;
	constexpr auto eps_1_He = 24.58738 * evtoerg;
	constexpr auto eps_2_He = 54.41776 * evtoerg;
	constexpr auto g0_H = 2.0;
	constexpr auto g1_H = 1.0;
	constexpr auto g0_He = 1.0;
	constexpr auto g1_He = 2.0;
	constexpr auto g2_He = 1.0;
	const auto lambda3 = std::pow(hplanck * hplanck / (2 * me * kb * T), 1.5);
	const auto A0 = 2.0 / lambda3 * std::exp(-eps_1_H / (kb * T)) * g1_H / g0_H;
	const auto B0 = 2.0 / lambda3 * std::exp(-eps_1_He / (kb * T)) * g1_He / g0_He;
	const auto C0 = 2.0 / lambda3 * std::exp(-(eps_2_He - eps_1_He) / (kb * T)) * g2_He / g1_He;
	const double he_num_frac = hefrac / 4 / (hefrac / 4 + (1 - hefrac) / 1);

	const double n_nuc = rho / mh;

	double H0 = n_nuc * (1 - hefrac);
	double He0 = n_nuc * hefrac / 4;
	double err;
	ne = 0.5 * (H0 + 2 * He0);
	double Y;
#define Power std::pow
	int iters = 0;
//	printf("p&e begin\n");
	do {
//		printf( "%e %e\n", T, ne);
//		printf("%e\n", err);
		const auto ne0 = ne;
		const auto dne = -((ne - (A0 * H0) / (A0 + ne) - (B0 * He0 * (2 * C0 + ne)) / (Power(ne, 2) + B0 * (C0 + ne)))
				/ (1 + (A0 * H0) / Power(A0 + ne, 2) + (B0 * He0 * (B0 * C0 + ne * (4 * C0 + ne))) / Power(Power(ne, 2) + B0 * (C0 + ne), 2)));
		ne = std::min(std::max(0.5 * ne, ne + dne), 2 * ne);
		if (ne < 1e-100) {
			break;
		}
		err = std::abs(std::log(ne / ne0));
		iters++;
		if (iters > 1000) {
			printf("Max iters exceed in compute_electron_fraction\n");
			abort();
		}
	} while (err > 1.0e-6);
	Hp = A0 * H0 / (A0 + ne);
	H = H0 - Hp;
	Hep = B0 * ne * He0 / (B0 * C0 + B0 * ne + ne * ne);
	Hepp = B0 * C0 * He0 / (B0 * C0 + B0 * ne + ne * ne);
	He = He0 - Hep - Hepp;
//	printf("p&e end\n");
}
//
//void cooling_rate(double H, double Hp, double He, double Hep, double Hepp, double z, double T, double &C1, double &C2, double &C3, double &C4, double &C5,
//		double &C6, double &C7, double &C8, double &C9, double &C10, double &C11, double &C12, double &C13) {
//	const auto ne = Hp + Hep + 2 * Hepp;
//	const auto T3 = T / 1e3;
//	const auto T5 = T / 1e5;
//	const auto T6 = T / 1e6;
//	const auto tev = T * 8.61732814974056E-05;
//	const auto logtev = std::log(T);
//	const auto tiny = std::numeric_limits<double>::min();
//
//	const auto k1 = std::exp(
//			-32.71396786375 + 13.53655609057 * logtev - 5.739328757388 * std::pow(logtev, 2) + 1.563154982022 * std::pow(logtev, 3)
//					- 0.2877056004391 * std::pow(logtev, 4) + 0.03482559773736999 * std::pow(logtev, 5) - 0.00263197617559 * std::pow(logtev, 6)
//					+ 0.0001119543953861 * std::pow(logtev, 7) - 2.039149852002e-6 * std::pow(logtev, 8));
//
//	const auto k3 = std::exp(
//			-44.09864886561001 + 23.91596563469 * logtev - 10.75323019821 * std::pow(logtev, 2) + 3.058038757198 * std::pow(logtev, 3)
//					- 0.5685118909884001 * std::pow(logtev, 4) + 0.06795391233790001 * std::pow(logtev, 5) - 0.005009056101857001 * std::pow(logtev, 6)
//					+ 0.0002067236157507 * std::pow(logtev, 7) - 3.649161410833e-6 * std::pow(logtev, 8));
//
//	const auto k5 = std::exp(
//			-68.71040990212001 + 43.93347632635 * logtev - 18.48066993568 * std::pow(logtev, 2) + 4.701626486759002 * std::pow(logtev, 3)
//					- 0.7692466334492 * std::pow(logtev, 4) + 0.08113042097303 * std::pow(logtev, 5) - 0.005324020628287001 * std::pow(logtev, 6)
//					+ 0.0001975705312221 * std::pow(logtev, 7) - 3.165581065665e-6 * std::pow(logtev, 8));
//
//	if (T > 250) {
//		C1 = 7.50e-19 * std::exp(-118348 / T) / (1.0 + std::sqrt(T / 100000));
//	} else {
//		C1 = 0.0;
//	}
//	if (T > 25) {
//		C2 = 9.10e-27 * std::pow(T, -.1687) * std::exp(-13179 / T) / (1.0 + std::sqrt(T / 100000));
//	} else {
//		C2 = 0.0;
//	}
//	if (T > 1000) {
//		C3 = 5.54e-17 * std::pow(T, -.397) * std::exp(-473638 / T) / (1.0 + std::sqrt(T / 100000));
//	} else {
//		C3 = 0.0;
//	}
//	C4 = 2.18e-11 * k1;
//	C5 = 3.94e-11 * k3;
//	C6 = 8.72e-11 * k5;
//	if (T > 100) {
//		C7 = 5.01e-27 * std::pow(T, -.1687) * std::exp(-55338 / T) / (1.0 + std::sqrt(T / 100000));
//	} else {
//		C7 = 0.0;
//	}
//	C8 = 8.70e-27 * std::sqrt(T) * std::pow(T3, -0.2) / (1.0 + std::pow(T6, 0.7));
//	C9 = 1.55e-26 * std::pow(T, 0.3647);
//	if (T > 1000) {
//		C10 = 1.24e-13 * std::pow(T, -1.5) * (1 + 0.3 * std::exp(-94000 / T)) * std::exp(-470000 / T);
//	} else {
//		C10 = 0.0;
//	}
//	C11 = 3.48e-26 * std::sqrt(T) * std::pow(T / 1000, -0.2) / (1.0 + std::pow(T / 1000000, 0.7));
//	C12 = 1.43e-27 * std::sqrt(T) * (1.1 + 0.34 * std::exp(-std::pow(5.50 - std::log10(T), 2) / 3.0));
//	C13 = 5.64e-36 * std::pow(1 + z, 4) * (T);
//
//	C1 *= ne * H;
//	C2 *= ne * ne * He;
//	C3 *= ne * Hep;
//	C4 *= ne * H;
//	C5 *= ne * He;
//	C6 *= ne * Hep;
//	C7 *= ne * ne * Hep;
//	C8 *= ne * Hp;
//	C9 *= ne * Hep;
//	C10 *= ne * Hep;
//	C11 *= ne * Hepp;
//	C12 *= ne * (Hp + Hep + Hepp);
//	C13 *= ne;
//
//}

double fermi_dirac(double k, double eta, double beta) {
	const auto func = [k, eta, beta](double x) {
		if (x != 0.0) {
			double part1 = std::pow(x, k) * std::sqrt(1.0 + 0.5 * beta * x) / (std::exp(x - eta) + 1.0);
			double part2 = std::pow(x, -k - 2) * std::sqrt(1.0 + 0.5 * beta / x) / (std::exp(1.0 / x - eta) + 1.0);
			return part1 + part2;
		} else {
			return 0.0;
		}
	};
	return integrate(func, 0.0, 1.0);
}

double compute_Nele(double eta, double beta) {
	using namespace constants;
	constexpr double c0 = 8.0 * M_PI * std::sqrt(2) * std::pow(me * clight / hplanck, 3);
	return c0 * std::pow(beta, 1.5) * (fermi_dirac(0.5, eta, beta) + beta * fermi_dirac(1.5, eta, beta));
}

double compute_Npos(double eta, double beta) {
	using namespace constants;
	constexpr double c0 = 8.0 * M_PI * std::sqrt(2) * std::pow(me * clight / hplanck, 3);
	return c0 * std::pow(beta, 1.5) * (fermi_dirac(0.5, -eta - 2.0 / beta, beta) + beta * fermi_dirac(1.5, -eta - 2.0 / beta, beta));
}

void electron_state(double ne, double T, double &P, double &E) {
	using namespace constants;
	const auto beta = kb * T / (me * clight * clight);
	const auto eta = chemical_potential(ne, T);
	constexpr double c0 = 8.0 * M_PI * std::sqrt(2) * std::pow(me * clight / hplanck, 3) * me * clight * clight;
	const double F32ele = fermi_dirac(1.5, eta, beta);
	const double F52ele = fermi_dirac(2.5, eta, beta);
	const double F32pos = fermi_dirac(1.5, -eta - 2.0 / beta, beta);
	const double F52pos = fermi_dirac(2.5, -eta - 2.0 / beta, beta);
	const double Npos = compute_Npos(eta, beta);
	const double Pele = (2.0 / 3.0) * c0 * std::pow(beta, 2.5) * (F32ele + 0.5 * beta * F52ele);
	const double Ppos = (2.0 / 3.0) * c0 * std::pow(beta, 2.5) * (F32pos + 0.5 * beta * F52pos);
	const double Eele = c0 * std::pow(beta, 2.5) * (F32ele + beta * F52ele);
	const double Epos = c0 * std::pow(beta, 2.5) * (F32pos + beta * F52pos) + 2.0 * me * clight * clight * Npos;
	P = Ppos + Pele;
	E = Epos + Eele;
//	printf("%e %e\n", ne / compute_Nele(eta, beta), P/(kb*ne*T));

}

void pressure_and_energy(double rho, double T, double &P, double &eps) {
	using namespace constants;
	double H, Hp, He, Hep, Hepp, ne;
	compute_chemistry(rho, T, H, Hp, He, Hep, Hepp, ne);
//	printf("%e %e\n", ne, Ye);
	const double N = H + Hp + He + Hep + Hepp;
	double Pele, Eele, Pnuc, Enuc;
	if (ne / (H + Hp + He + Hep + Hepp) > 1.0e-6) {
		electron_state(ne, T, Pele, Eele);
	} else {
		Pele = Eele = 0.0;
	}
	Pnuc = kb * N * T;
	Enuc = 1.5 * Pnuc;
	eps = (Eele + Enuc) / rho;
	P = Pele + Pnuc;
}

void equation_of_state(double rho, double T, double &P, double &cs, double &eps) {
	using namespace constants;
	double Pp, Pn, epsp, epsn;
	double dT = T * 1.0e-4;
	double drho = rho * 1.0e-4;

	pressure_and_energy(rho, T + 0.5 * dT, Pp, epsp);
	pressure_and_energy(rho, T - 0.5 * dT, Pn, epsn);
	pressure_and_energy(rho, T, P, eps);
	const double dPdT = (Pp - Pn) / dT;
	const double Cv = (epsp - epsn) / dT;
	pressure_and_energy(rho + 0.5 * drho, T, Pp, epsp);
	pressure_and_energy(rho - 0.5 * drho, T, Pn, epsn);
	const double dPdrho = (Pp - Pn) / drho;
	const double cs2 = (P / (Cv * rho * rho) * dPdT + dPdrho) / (P / rho + eps + clight * clight);
	if (cs2 > 1) {
		printf("%e\n", std::sqrt(cs2));
	}
	assert(cs2 >= 0.0);
//	assert(cs2 < 1.0);
	cs = clight * std::sqrt(cs2);

}

void zero_order_universe(double amin, double amax) {
	using namespace constants;
	double minloga = std::log(amin);
	double maxloga = std::log(amax);
	double loga = minloga;
	constexpr int N = 100;
	double dloga = (maxloga - minloga) / N;
	double Trad = Theta * 2.73 / amin;
	double Tgas = Trad;
	for (int i = 0; i <= N; i++) {
		double loga = minloga + i * dloga;
		double a = std::exp(loga);
		double rho = rho_b(a) * 1000;
		double H, Hp, He, Hep, Hepp, ne;
		double C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, C11, C12, C13;
		Trad = Theta * 2.73 / a;
		compute_chemistry(rho, Tgas, H, Hp, He, Hep, Hepp, ne);
		const auto n = H + Hp + He + Hep + Hepp;
		const auto Y = ne / (H + Hp + 2 * (He + Hep + Hepp));
		const auto Cv = 1.5 * kb * (n + ne);
		const auto P = kb * (n + ne) * Tgas;
		const auto Pele = kb * ne * Tgas;
		const auto Pion = kb * n * Tgas;
		double Pfermi, csfermi, eps;
		equation_of_state(rho, Tgas, Pfermi, csfermi, eps);
		const auto cs = std::sqrt((5.0 / 3.0) * P / rho / (5.0 / 2.0 * P / rho + clight * clight));
		const auto Comp = 5.64e-36 * std::pow(1 / (a * Trad), 4) * ne;
		const auto hubble = Hubble(a) * H0cgs / H0;
		const auto c0 = Comp / (Cv * hubble);
		printf("%10.2e %10.2e %10.2e %10.2e %10.2e %10.2e  %10.2e  %10.2e %10.2e %10.2e  %10.2e \n", a, 1 / a - 1, rho, n, Trad, Tgas, Y, cs, csfermi / clight,
				P, Pfermi / P);
		Tgas = (Tgas + dloga * c0 * std::pow(Trad, 5)) / (1 + dloga * (2 + c0 * std::pow(Trad, 4)));

	}

}

//
auto build_electron_fraction_interpolation_function(double amin, double amax) {
//	double minloga = std::log(amin);
//	double maxloga = std::log(amax);
//	double dloga = 1.0 / 256.0;
//	int N = (maxloga - minloga) / dloga + 1;
//	dloga = (maxloga - minloga) / N;
//	N += 2;
////	printf( "electron fraction table size = %i\n", N);
//	maxloga += dloga;
//	minloga -= dloga;
//	dloga = (maxloga - minloga) / N;
//	std::vector<double> values(N + 1);
//	for (int i = 0; i <= N; i++) {
//		double loga = i * dloga + minloga;
//		double a = std::exp(loga);
//		values[i] = compute_electron_fraction(rho_b(a), 2.73 * Theta / a);
//	}
//
//	const auto func = [=](double loga) {
//		if (loga - dloga < minloga || loga + dloga > maxloga) {
//			printf("Error in build_electron_fraction_interpolation_function\n");
//			abort();
//		}
//		int i1 = std::min(std::max(1, int((loga - minloga) / (dloga))), N - 2);
//		int i0 = i1 - 1;
//		int i2 = i1 + 1;
//		int i3 = i2 + 1;
//		const double c0 = values[i1];
//		const double c1 = -values[i0] / 3.0 - 0.5 * values[i1] + values[i2] - values[i3] / 6.0;
//		const double c2 = 0.5 * values[i0] - values[i1] + 0.5 * values[i2];
//		const double c3 = -values[i0] / 6.0 + 0.5 * values[i1] - 0.5 * values[i2] + values[i3] / 6.0;
//		double x = (loga - i1 * dloga - minloga) / dloga;
////		return values[i2] * x + values[i1] * (1 - x);
//		return std::max(c0 + c1 * x + c2 * x * x + c3 * x * x * x, 0.0);
//	};
//	return func;
}

void advance(state &u, double k, double a0, double a1) {
//	auto Ye = build_electron_fraction_interpolation_function(1e-10, 1);
//	const double logamin = std::log(a0);
//	const double logamax = std::log(a1);
//	double loga = logamin;
//	while (loga < logamax) {
//		double a, H, Om, Or, eps, fb, fc, fgam, sigma, R, Onu, Ogam, fnu;
//		const auto compute_parameters = [&]() {
//			a = std::exp(loga);
//			H = Hubble(a);
//			Om = omega_m / (a * a * a) * std::pow(H0 * h / H, 2);
//			Or = omega_r / (a * a * a * a) * std::pow(H0 * h / H, 2);
//			eps = k / (a * H);
//			fb = omega_b / omega_m;
//			fgam = omega_gam / omega_r;
//			fnu = 1 - fgam;
//			fc = 1 - fb;
//			R = (3.0 * Om) / (4.0 * Or * fgam);
//		};
//		compute_parameters();
//		const auto lam_max = eps + 4.5 * Om + 12 * Or + eps * eps + 3.0;
//		const auto dloga = std::min(1.0 / lam_max, logamax - loga);
//		const auto mass = fc * u[deltaci] + fb * u[deltabi];
//		const auto rad = fnu * u[N0i] + fgam * u[Theta0i];
//		const auto mmom = fc * u[uci] + fb * u[ubi];
//		const auto rmom = fnu * u[N1i] + fgam * u[Theta1i];
//		const auto Phi = 1.5 / (eps * eps) * (4 * Or * (rad + rmom / eps) + Om * (mass + mmom / eps));
////		printf( "%e %e\n", u[Phii], Phi);
//		const auto compute_dudt_exp = [&](state u) {
//			state dudt;
//			const auto mass = fc * u[deltaci] + fb * u[deltabi];
//			const auto rad = fnu * u[N0i] + fgam * u[Theta0i];
//			dudt[Phii] = (0.5 * Om * mass + 2.0 * Or * rad - ((1.0 / 3.0) * eps * eps + 1.0) * u[Phii]);
//			dudt[Theta0i] = -eps * u[Theta1i] - dudt[Phii];
//			dudt[Theta1i] = eps / 3.0 * (u[Theta0i] - u[Phii]);
//			dudt[N0i] = -eps * u[N1i] - dudt[Phii];
//			dudt[N1i] = eps / 3.0 * (u[N0i] - u[Phii]);
//			dudt[deltaci] = -eps * u[uci] - 3 * dudt[Phii];
//			dudt[uci] = -u[uci] - eps * u[Phii];
//			dudt[deltabi] = -eps * u[ubi] - 3 * dudt[Phii];
//			dudt[ubi] = -u[ubi] - eps * u[Phii];
//			return dudt;
//		};
//		state u0 = u;
//		auto dudt = compute_dudt_exp(u0);
//		for (int i = 0; i < NF; i++) {
//			u[i] = u0[i] + dudt[i] * dloga;
//		}
//
//		sigma = sigma_T * omega_b * h * h * Ye(loga + dloga) / (a * a * a * H);
//		u[Theta1i] = (dloga * dloga * (3 * dudt[Theta1i] + dudt[ubi] * R) * sigma + 3 * R * u0[Theta1i]
//				+ dloga * (3 * dudt[Theta1i] * R + 3 * sigma * u0[Theta1i] + R * sigma * u0[ubi])) / (3. * (R + dloga * sigma + dloga * fgam * R * sigma));
//		u[ubi] = (dloga * dudt[ubi] * (R + dloga * fgam * R * sigma) + 3 * dloga * fgam * sigma * (dloga * dudt[Theta1i] + u0[Theta1i])
//				+ R * (1 + dloga * fgam * sigma) * u0[ubi]) / (R + dloga * sigma + dloga * fgam * R * sigma);
//		loga += dloga;
//	}
//
}

void initial_conditions(state &u, double k, double a) {
	const auto eps = k / (a * Hubble(a));
	u[Theta0i] = 0.5 * u[Phii];
	u[Theta1i] = -eps / 6.0 * u[Phii];
	u[N0i] = u[Theta0i];
	u[N1i] = u[Theta1i];
	u[deltaci] = 3.0 * u[Theta0i];
	u[uci] = 3.0 * u[Theta1i];
	u[deltabi] = 3.0 * u[Theta0i];
	u[ubi] = 3.0 * u[Theta1i];
}

int main() {
//	using namespace constants;
//	double Tgas = 1e6;
//	double eps;
//	for (double rho = 1e-14; rho < 1e+20; rho *= 2) {
//		double Pfermi, csfermi;
//		double P = rho * Tgas * kb / mh * (1 - hefrac + hefrac / 4);
//		equation_of_state(rho, Tgas, Pfermi, csfermi, eps);
//		printf("%e %e\n", rho, Pfermi/P);
//	}
////	abort();
////
	zero_order_universe(1.0e-10, 1);
//	for (double z = 0.001; z < 1.5; z += 0.1) {
//		printf("%e %e \n", z, hyper_geometric(0.5, 0.5, 1.5, z).real());
//	}
//	for (double x = 1; x < 4.0; x += 0.1) {
//		printf("%e %e\n", x, gamma_cmplx(cmplx(x, 0)).real());
//	}
//	double P, cs;
//	using namespace constants;
//	auto func = build_electron_fraction_interpolation_function(1e-12, 1);
//	for (double a = 1e-12; a <= 1.0; a *= 1.5) {
//		const double T = 2.73 * Theta / a;
//		equation_of_state(rho_b(a), T, P, cs);
//		printf("%e %e %e %e %e %e\n", a, rho_b(a), T, P, cs / clight, func(std::log(a)));
//	}
	//	auto func = build_electron_fraction_interpolation_function(1e-10, 1);
//	double max_err = 0.0;
//	int n = 0;
//	for (double a = 1.0e-10; a <= 1.0; a *= 1.01) {
//		n++;
//		const auto y0 = compute_electron_fraction(a);
//		const auto y1 = func(std::log(a));
//		max_err = std::max(max_err, std::abs(y0 - y1));
//		printf("%e %e %e %e\n", 1 / a - 1, y0, 2.73 / a, y1);
//	}
//	printf("Error = %e\n", max_err);
//	compute_electron_fraction(1, 1e10);
}
