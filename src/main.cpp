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

#define Power std::pow

//const double omega_b_h2 = 0.0240;
//const double omega_c_h2 = 0.1146;
//const double h = 0.697;
//const double z_eq = 3273;
//const double Neff = 3.84;
//const double hefrac = 0.24;
//const double Theta = 1.00;

const double h = 0.7;
const double omega_b_h2 = 0.05 * h * h;
const double omega_c_h2 = 0.25 * h * h;
const double z_eq = 3273;
const double Neff = 3.086;
const double hefrac = 0.24;
const double Theta = 1.00;

const double omega_b = omega_b_h2 / (h * h);
const double omega_c = omega_c_h2 / (h * h);
const double omega_m = omega_b + omega_c;
const double omega_r = 1.0 / (z_eq + 1.0) * omega_m;
const double omega_nu = omega_r * Neff / (8.0 / 7.0 * std::pow(11.0 / 4.0, 4.0 / 3.0) + Neff);
const double omega_gam = omega_r - omega_nu;
const double omega_lambda = 1.0 - omega_m - omega_r;
const double H0 = 100.0 * 1e5 / constants::clight;
const double ns = 1;
const double sigma_T = 2.3e-5;
const double Yp = (1 - hefrac / 2.0);

double Hubble(double a) {
	return H0 * h * std::sqrt(omega_r / (a * a * a * a) + omega_m / (a * a * a) + omega_lambda);
}
#define LMAX 4

#define Phii 0
#define deltai 1
#define ui 2
#define deltabi 3
#define ubi 4
#define Ni 5
#define Thetai (5+LMAX)
#define ThetaPi (5+LMAX*2)
#define NF (5+LMAX*3)
using state = std::array<double,NF>;

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
		w *= 0.99;
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
}

double fermi_dirac(double k, double eta, double beta) {
	if (eta < -5) {
		return std::exp(eta) * gamma_real(1 + k);
	} else {
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
}

double compute_Nele(double eta, double beta) {
	using namespace constants;
	constexpr double c0 = 8.0 * M_PI * std::sqrt(2) * std::pow(me * clight / hplanck, 3);
	double y = c0 * std::pow(beta, 1.5) * (fermi_dirac(0.5, eta, beta) + beta * fermi_dirac(1.5, eta, beta));
	return y;
}

double compute_Npos(double eta, double beta) {
	using namespace constants;
	constexpr double c0 = 8.0 * M_PI * std::sqrt(2) * std::pow(me * clight / hplanck, 3);
	double y = c0 * std::pow(beta, 1.5) * (fermi_dirac(0.5, -eta - 2.0 / beta, beta) + beta * fermi_dirac(1.5, -eta - 2.0 / beta, beta));
	return y;
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

}

void pressure_and_energy(double rho, double T, double &P, double &eps) {
	using namespace constants;
	double H, Hp, He, Hep, Hepp, ne;
	compute_chemistry(rho, T, H, Hp, He, Hep, Hepp, ne);
//	printf("%e %e\n", ne, Ye);
	const double N = H + Hp + He + Hep + Hepp;
	double Pele, Eele, Pnuc, Enuc;
	if (kb * T > 1.0e-2 * me * clight * clight && ne / (H + Hp + He + Hep + Hepp) > 1.0e-6) {
		electron_state(ne, T, Pele, Eele);
	} else {
//		printf("Shortcut\n");
		Pele = kb * ne * T;
		Eele = 1.5 * Pele;
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

void chemical_rates(double &k1, double &k2, double &k3, double &k4, double &k5, double &k6, double T) {

	const auto tev = T * 8.61732814974056E-05;
	const auto logtev = std::log(tev);
	k1 = exp(
			-32.71396786375 + 13.53655609057 * logtev - 5.739328757388 * std::pow(logtev, 2) + 1.563154982022 * std::pow(logtev, 3)
					- 0.2877056004391 * std::pow(logtev, 4) + 0.03482559773736999 * std::pow(logtev, 5) - 0.00263197617559 * std::pow(logtev, 6)
					+ 0.0001119543953861 * std::pow(logtev, 7) - 2.039149852002e-6 * std::pow(logtev, 8));

	k3 = exp(
			-44.09864886561001 + 23.91596563469 * logtev - 10.75323019821 * std::pow(logtev, 2) + 3.058038757198 * std::pow(logtev, 3)
					- 0.5685118909884001 * std::pow(logtev, 4) + 0.06795391233790001 * std::pow(logtev, 5) - 0.005009056101857001 * std::pow(logtev, 6)
					+ 0.0002067236157507 * std::pow(logtev, 7) - 3.649161410833e-6 * std::pow(logtev, 8));

	k4 = 3.92e-13 / std::pow(tev, 0.6353);

	if (tev > 0.1) {
		k4 += 1.54e-9 * (1. + 0.3 / std::exp(8.099328789667 / tev)) / (exp(40.49664394833662 / tev) * std::pow(tev, 1.5));
	}

	k5 = exp(
			-68.71040990212001 + 43.93347632635 * logtev - 18.48066993568 * std::pow(logtev, 2) + 4.701626486759002 * std::pow(logtev, 3)
					- 0.7692466334492 * std::pow(logtev, 4) + 0.08113042097303 * std::pow(logtev, 5) - 0.005324020628287001 * std::pow(logtev, 6)
					+ 0.0001975705312221 * std::pow(logtev, 7) - 3.165581065665e-6 * std::pow(logtev, 8));

	k2 = exp(
			-28.61303380689232 - 0.7241125657826851 * logtev - 0.02026044731984691 * std::pow(logtev, 2) - 0.002380861877349834 * std::pow(logtev, 3)
					- 0.0003212605213188796 * std::pow(logtev, 4) - 0.00001421502914054107 * std::pow(logtev, 5) + 4.989108920299513e-6 * std::pow(logtev, 6)
					+ 5.755614137575758e-7 * std::pow(logtev, 7) - 1.856767039775261e-8 * std::pow(logtev, 8) - 3.071135243196595e-9 * std::pow(logtev, 9));

	k6 = 3.36e-10 / sqrt(T) / std::pow((T / 1.e3), 0.2) / (1. + std::pow((T / 1.e6), 0.7));
}

void chemistry_update(double &H, double &Hp, double &He, double &Hep, double &Hepp, double &ne, double T, double Trad, double hubble, double dt) {
	const double H0 = H;
	const double Hp0 = Hp;
	const double He0 = He;
	const double Hep0 = Hep;
	const double Hepp0 = Hepp;
	double K1, K2, K3, K4, K5, K6;
	double J20, J21, J22;
	chemical_rates(K1, K2, K3, K4, K5, K6, T);
	photo_rates(J20, J21, J22, Trad);
//	J20 = J21 = J22 = 0.0;
	double nemax = (H0 + Hp0) + 2 * (He0 + Hep0 + Hepp0);
	double nemin = 0.0;
	double nemid;
	double err = 0.0;
	const auto compute_ne = [&](double ne) {
		using namespace constants;
		H = (H0 + dt * (H0 + Hp0) * K2 * ne) / (1 + dt * (J20 + (K1 + K2) * ne));
		Hp = H0 + Hp0 - H;
		const auto den = (-1 + Power(dt, 3) * J22 * (J21 + K3 * ne) * (J22 + (K5 - K6) * ne) - dt * (J21 + (K3 + K4 + K5 + K6) * ne)
				+ Power(dt, 2) * (Power(J22, 2) + J22 * (K4 + K5 - K6) * ne - ne * (J21 * (K5 + K6) + K4 * K6 * ne + K3 * (K5 + K6) * ne)));
		He = (He0 + dt * K4 * ne * (Hep0 + dt * (Hep0 + Hepp0) * K6 * ne) + dt * He0 * (J22 + ne * (K4 + K5 + K6 + dt * K4 * K6 * ne)))
				/ (1
						+ dt
								* (J21 + J22 + dt * J22 * K3 * ne + dt * J21 * (J22 + (K5 + K6) * ne)
										+ ne * (K3 + K4 + K5 + K6 + dt * (K3 * K5 + (K3 + K4) * K6) * ne)));
		Hep = (Hep0 + dt * He0 * J21 + dt * Hep0 * J21 + dt * ((He0 + Hep0) * K3 + (Hep0 + Hepp0 + dt * (He0 + Hep0 + Hepp0) * J21) * K6) * ne +
		Power(dt, 2) * (He0 + Hep0 + Hepp0) * K3 * K6 * Power(ne, 2))
				/ (1
						+ dt
								* (J21 + J22 + dt * J22 * K3 * ne + dt * J21 * (J22 + (K5 + K6) * ne)
										+ ne * (K3 + K4 + K5 + K6 + dt * (K3 * K5 + (K3 + K4) * K6) * ne)));
		Hepp = He0 + Hep0 + Hepp0 - He - Hep;
		return Hp + Hep + 2 * Hepp;
	};
	do {
		nemid = (nemax + nemin) / 2.0;
		double t1 = nemid - compute_ne(nemid);
		double t2 = nemax - compute_ne(nemax);
		if (t1 * t2 < 0.0) {
			nemin = nemid;
		} else {
			nemax = nemid;
		}
		err = 1.0 - nemin / nemax;
	} while (err > 1.0e-10 && nemid > 1e-100);
	ne = nemid;
	using namespace constants;
	constexpr auto evtoerg = 1.6021772e-12;
	constexpr auto eps_1_H = 13.59844 * evtoerg;
	constexpr auto g0_H = 2.0;
	constexpr auto g1_H = 1.0;
	const auto lambda3 = std::pow(hplanck * hplanck / (2 * me * kb * T), 1.5);
	const auto A0 = 2.0 / lambda3 * std::exp(-eps_1_H / (kb * T)) * g1_H / g0_H;
//	printf("%e %e %e %e\n", T, A0, J20 / K2, A0 / (J20/K2));
}

auto build_interpolation_function(const std::vector<double> values, double amin, double amax) {
	double minloga = std::log(amin);
	double maxloga = std::log(amax);
	int N = values.size() - 1;
	double dloga = (maxloga - minloga) / N;
	const auto func = [=](double loga) {
		if (loga < minloga || loga > maxloga) {
			printf("Error in build_electron_fraction_interpolation_function\n");
			abort();
		}
		int i1 = std::min(std::max(1, int((loga - minloga) / (dloga))), N - 2);
		int i0 = i1 - 1;
		int i2 = i1 + 1;
		int i3 = i2 + 1;
		const double c0 = values[i1];
		const double c1 = -values[i0] / 3.0 - 0.5 * values[i1] + values[i2] - values[i3] / 6.0;
		const double c2 = 0.5 * values[i0] - values[i1] + 0.5 * values[i2];
		const double c3 = -values[i0] / 6.0 + 0.5 * values[i1] - 0.5 * values[i2] + values[i3] / 6.0;
		double x = (loga - i1 * dloga - minloga) / dloga;
		return std::max(c0 + c1 * x + c2 * x * x + c3 * x * x * x, 0.0);
	};
	return func;
}

void zero_order_universe(double amin, double amax, std::function<double(double)> &csfunc, std::function<double(double)> &thomsonfunc) {
//	printf("omega_b = %e\n", omega_b);
//	printf("omega_c = %e\n", omega_c);
//	printf("omega_gam = %e\n", omega_gam);
//	printf("omega_nu = %e\n", omega_nu);
//	printf("omega_lam= %e\n", omega_lambda);

	using namespace constants;
	double minloga = std::log(amin);
	double maxloga = std::log(amax);
	double loga = minloga;
	constexpr int N = 1024;
	std::vector<double> sound_speed(N + 1), thomson(N + 1);
	double dloga = (maxloga - minloga) / N;
	double Trad = Theta * 2.73 / amin;
	double Tgas = Trad;
	double t = 0.0;
	double rho0 = rho_b(amin);
	double H, Hp, He, Hep, Hepp, ne;
	H = 0.0;
	Hp = (1.0 - hefrac) * rho0 / mh;
	He = 0.0;
	Hep = 0.0;
	Hepp = hefrac * rho0 / mh / 4;
	ne = Hp + 2 * Hepp;
	double last_a = amin;
	for (int i = 0; i <= N; i++) {
		double loga = minloga + i * dloga;
		double a = std::exp(loga);
		const auto hubble = Hubble(a) * H0cgs / H0;
		double fac = std::pow(last_a / a, 3);
		H *= fac;
		Hp *= fac;
		He *= fac;
		Hep *= fac;
		Hepp *= fac;
		ne *= fac;
		double rho = rho_b(a);
		Trad = Theta * 2.73 / a;
		const auto a3 = a * a * a;
		const auto dt = dloga / hubble;
		chemistry_update(H, Hp, He, Hep, Hepp, ne, Tgas, Trad, hubble, dt);
		const auto n = H + Hp + He + Hep + Hepp;
		const auto Y = ne / (H + Hp + 2 * (He + Hep + Hepp));
		const auto Cv = 1.5 * kb * (n + ne);
		const auto P = kb * (n + ne) * Tgas;
		const auto Pele = kb * ne * Tgas;
		const auto Pion = kb * n * Tgas;
		double Pfermi, cs, eps;
		equation_of_state(rho, Tgas, Pfermi, cs, eps);
		const auto Comp = 5.64e-36 * std::pow(1 / (a * Trad), 4) * ne;
		const auto sigma_T = 6.65e-25;
		const auto c0 = Comp / (Cv * hubble);
//		printf("%10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e  %10.2e  %10.2e  \n", t, a, 1 / a - 1, ne / n, H / n, Hp / n, He / n, Hep / n, Hepp / n);
		if (i > 0) {
			Tgas = (Tgas + dloga * c0 * std::pow(Trad, 5)) / (1 + dloga * (2 + c0 * std::pow(Trad, 4)));
		}
		sound_speed[i] = cs / clight;
		thomson[i] = clight * sigma_T * ne / hubble;
//		printf("%10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e  \n", t, a, 1 / a - 1, rho, n, Trad, Tgas, cs / clight, ne / n,
//				clight * sigma_T * ne / hubble);
		t += dt;
		//	t += dloga / hubble;
		last_a = a;
	}
	csfunc = build_interpolation_function(sound_speed, amin, amax);
	thomsonfunc = build_interpolation_function(thomson, amin, amax);
}

void advance(state &U, double k, double a0, double a1, std::function<double(double)> cs, std::function<double(double)> thomson) {
	const double logamin = std::log(a0);
	const double logamax = std::log(a1);
	double loga = logamin;
	double eta = 1.0 / (a0 * Hubble(a0));
	double cs2, etaaH, sigma;
	while (loga < logamax) {
		double a, H, Om, Or, eps, fb, fc, fgam, sigma, R, Onu, Ogam, fnu;
		const auto compute_parameters = [&]() {
			a = std::exp(loga);
			H = Hubble(a);
			Om = omega_m / (a * a * a) * std::pow(H0 * h / H, 2);
			Or = omega_r / (a * a * a * a) * std::pow(H0 * h / H, 2);
			eps = k / (a * H);
			fb = omega_b / omega_m;
			fgam = omega_gam / omega_r;
			fnu = 1 - fgam;
			fc = 1 - fb;
			R = (3.0 * Om) / (4.0 * Or * fgam);
			sigma = thomson(loga);
			cs2 = std::pow(cs(loga), 2);
			etaaH = eta * a * H;
		};
		compute_parameters();
		const auto lambda_max = std::max(std::max(std::max(1.0 + cs2 * eps, eps + 4.0 / etaaH), 0.5 * Om + 2.0 * Or), 1.0 + eps * eps);
		const auto dloga = std::min(0.5 / lambda_max, logamax - loga);

		auto &Phi = U[Phii];
		auto &delta = U[deltai];
		auto &u = U[ui];
		auto &deltab = U[deltabi];
		auto &ub = U[ubi];
		auto &Theta0 = U[Thetai + 0];
		auto &Theta1 = U[Thetai + 1];
		auto &Theta2 = U[Thetai + 2];
		auto &Theta3 = U[Thetai + 3];
		auto &ThetaP0 = U[ThetaPi + 0];
		auto &ThetaP1 = U[ThetaPi + 1];
		auto &ThetaP2 = U[ThetaPi + 2];
		auto &ThetaP3 = U[ThetaPi + 3];
		auto &N0 = U[Ni + 0];
		auto &N1 = U[Ni + 1];
		auto &N2 = U[Ni + 2];
		auto &N3 = U[Ni + 3];

		const auto Phi0 = U[Phii];
		const auto delta0 = U[deltai];
		const auto u0 = U[ui];
		const auto deltab0 = U[deltabi];
		const auto ub0 = U[ubi];
		const auto Theta00 = U[Thetai + 0];
		const auto Theta10 = U[Thetai + 1];
		const auto Theta20 = U[Thetai + 2];
		const auto Theta30 = U[Thetai + 3];
		const auto ThetaP00 = U[ThetaPi + 0];
		const auto ThetaP10 = U[ThetaPi + 1];
		const auto ThetaP20 = U[ThetaPi + 2];
		const auto ThetaP30 = U[ThetaPi + 3];
		const auto N00 = U[Ni + 0];
		const auto N10 = U[Ni + 1];
		const auto N20 = U[Ni + 2];
		const auto N30 = U[Ni + 3];
		const auto Psi = -12.0 * Or * (fgam * Theta2 + fnu * N2) / (eps * eps) - Phi;
//		printf("%10.2e ", a);
//		printf("%10.2e ", Psi);
//		for (int n = 0; n < NF; n++) {
//			printf("%10.2e ", U[n]);
//		}
//		printf("\n");

#define List(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16) a0; a1; a2; a3; a4; a5; a6; a7; a8; a9; a10; a11; a12; a13; a14; a15; a16
#define Rule(a,b) a = b
		List(
				Rule(Phi,(15*deltab0*dloga*Power(eps,2)*fb*Om*(5 + 6*dloga*sigma + Power(dloga,2)*Power(sigma,2)) + 15*delta0*dloga*Power(eps,2)*fc*Om*(5 + 6*dloga*sigma + Power(dloga,2)*Power(sigma,2)) - 2*(-75*Power(eps,2)*Phi0 + 36*Power(dloga,4)*eps*fnu*(2*N10 - 3*N30)*Or*Power(sigma,2) + 5*dloga*(5*Power(eps,4)*Phi0 - 3*Power(eps,2)*(10*fnu*N00*Or - 5*Phi0 + 6*Phi0*sigma + 10*fgam*Or*Theta00) + 180*Or*(fnu*N20 + fgam*Theta20)) + 15*Power(dloga,2)*(2*Power(eps,4)*Phi0*sigma - Power(eps,2)*sigma*(12*fnu*N00*Or + Phi0*(-6 + sigma) + 12*fgam*Or*Theta00) + 12*eps*Or*(2*fnu*N10 - 3*fnu*N30 + 2*fgam*Theta10 - 3*fgam*Theta30) + 6*Or*sigma*(12*fnu*N20 + fgam*(3*Theta20 + ThetaP00 + ThetaP20))) + Power(dloga,3)*sigma* (180*fnu*N20*Or*sigma + 5*Power(eps,4)*Phi0*sigma - 15*Power(eps,2)*sigma*(2*fnu*N00*Or - Phi0 + 2*fgam*Or*Theta00) + 54*eps*Or*(4*fnu*(2*N10 - 3*N30) + fgam*(2*Theta10 - 3*Theta30 - ThetaP10 - ThetaP30)))))/(30.*Power(eps,2)*(1 + dloga*sigma)*(5 + dloga*sigma))),
				Rule(delta,(-15*deltab0*dloga*Power(eps,2)*fb*Om*(5 + 6*dloga*sigma + Power(dloga,2)*Power(sigma,2)) - 5*delta0*Power(eps,2)*(-2 + 3*dloga*fc*Om)*(5 + 6*dloga*sigma + Power(dloga,2)*Power(sigma,2)) + 2*dloga*(5*Power(eps,4)*Phi0*(5 + 6*dloga*sigma + Power(dloga,2)*Power(sigma,2)) - 15*Power(eps,2)*(5 + 6*dloga*sigma + Power(dloga,2)*Power(sigma,2))*(2*fnu*N00*Or - Phi0 + 2*fgam*Or*Theta00) + 90*Or*(2*fnu*N20*(5 + 6*dloga*sigma + Power(dloga,2)*Power(sigma,2)) + fgam*(10 + 3*dloga*sigma)*Theta20 + dloga*fgam*sigma*(ThetaP00 + ThetaP20)) + 18*dloga*eps*Or*(2*fnu*(2*N10 - 3*N30)*(5 + 6*dloga*sigma + Power(dloga,2)*Power(sigma,2)) + fgam*(20 + 6*dloga*sigma)*Theta10 - 3*fgam*((10 + 3*dloga*sigma)*Theta30 + dloga*sigma*(ThetaP10 + ThetaP30))) - 5*Power(eps,3)*(5 + 6*dloga*sigma + Power(dloga,2)*Power(sigma,2))*u0))/ (10.*Power(eps,2)*(1 + dloga*sigma)*(5 + dloga*sigma))),
				Rule(u,(-15*deltab0*Power(dloga,2)*Power(eps,2)*fb*Om*(5 + 6*dloga*sigma + Power(dloga,2)*Power(sigma,2)) - 15*delta0*Power(dloga,2)*Power(eps,2)*fc*Om*(5 + 6*dloga*sigma + Power(dloga,2)*Power(sigma,2)) + 2*(36*Power(dloga,5)*eps*fnu*(2*N10 - 3*N30)*Or*Power(sigma,2) + Power(dloga,4)*sigma* (180*fnu*N20*Or*sigma + 5*Power(eps,4)*Phi0*sigma - 15*Power(eps,2)*sigma*(2*fnu*N00*Or - Phi0 + 2*fgam*Or*Theta00) - 18*eps*Or*(2*fnu*(2*N10 - 3*N30)*(-6 + sigma) + 3*fgam*(-2*Theta10 + 3*Theta30 + ThetaP10 + ThetaP30))) + 75*eps*u0 - 15*dloga*(60*fnu*N20*Or + 5*Power(eps,2)*Phi0 + 60*fgam*Or*Theta20 + eps*(5 - 6*sigma)*u0) + 5*Power(dloga,2)*(5*Power(eps,4)*Phi0 - 3*Power(eps,2)*(10*fnu*N00*Or + Phi0*(-5 + 6*sigma) + 10*fgam*Or*Theta00) - 18*Or*(2*fnu*N20*(-5 + 6*sigma) + fgam*(-10 + 3*sigma)*Theta20 + fgam*sigma*(ThetaP00 + ThetaP20)) - 3*eps*(12*fnu*(2*N10 - 3*N30)*Or + 12*fgam*Or*(2*Theta10 - 3*Theta30) - (-6 + sigma)*sigma*u0)) + 3*Power(dloga,3)*(10*Power(eps,4)*Phi0*sigma - 5*Power(eps,2)*sigma*(12*fnu*N00*Or + Phi0*(-6 + sigma) + 12*fgam*Or*Theta00) + 30*Or*sigma*(-2*fnu*N20*(-6 + sigma) + fgam*(3*Theta20 + ThetaP00 + ThetaP20)) + eps*(-12*fnu*(2*N10 - 3*N30)*Or*(-5 + 6*sigma) - 6*fgam*Or*((-20 + 6*sigma)*Theta10 - 3*((-10 + 3*sigma)*Theta30 + sigma*(ThetaP10 + ThetaP30))) - 5*Power(sigma,2)*u0))))/ (30.*eps*(1 + dloga*sigma)*(5 + dloga*sigma))),
				Rule(deltab,(-15*delta0*dloga*Power(eps,2)*fc*Om*(5 + 6*dloga*sigma + Power(dloga,2)*Power(sigma,2)) - 5*deltab0*Power(eps,2)*(-2 + 3*dloga*fb*Om)*(5 + 6*dloga*sigma + Power(dloga,2)*Power(sigma,2)) + 2*dloga*(5*Power(eps,4)*Phi0*(5 + 6*dloga*sigma + Power(dloga,2)*Power(sigma,2)) - 15*Power(eps,2)*(5 + 6*dloga*sigma + Power(dloga,2)*Power(sigma,2))*(2*fnu*N00*Or - Phi0 + 2*fgam*Or*Theta00) + 90*Or*(2*fnu*N20*(5 + 6*dloga*sigma + Power(dloga,2)*Power(sigma,2)) + fgam*(10 + 3*dloga*sigma)*Theta20 + dloga*fgam*sigma*(ThetaP00 + ThetaP20)) + 18*dloga*eps*Or*(2*fnu*(2*N10 - 3*N30)*(5 + 6*dloga*sigma + Power(dloga,2)*Power(sigma,2)) + fgam*(20 + 6*dloga*sigma)*Theta10 - 3*fgam*((10 + 3*dloga*sigma)*Theta30 + dloga*sigma*(ThetaP10 + ThetaP30))) - 5*Power(eps,3)*(5 + 6*dloga*sigma + Power(dloga,2)*Power(sigma,2))*ub0))/ (10.*Power(eps,2)*(1 + dloga*sigma)*(5 + dloga*sigma))),
				Rule(ub,-(75*delta0*Power(dloga,2)*Power(eps,2)*fc*Om*R + 300*Power(dloga,2)*Power(eps,2)*fnu*N00*Or*R + 720*Power(dloga,2)*eps*fnu*N10*Or*R - 720*Power(dloga,3)*eps*fnu*N10*Or*R + 1800*dloga*fnu*N20*Or*R - 1800*Power(dloga,2)*fnu*N20*Or*R - 1080*Power(dloga,2)*eps*fnu*N30*Or*R + 1080*Power(dloga,3)*eps*fnu*N30*Or*R + 150*dloga*Power(eps,2)*Phi0*R - 150*Power(dloga,2)*Power(eps,2)*Phi0*R - 50*Power(dloga,2)*Power(eps,4)*Phi0*R + 75*delta0*Power(dloga,3)*Power(eps,2)*fc*Om*sigma + 300*Power(dloga,3)*Power(eps,2)*fnu*N00*Or*sigma + 720*Power(dloga,3)*eps*fnu*N10*Or*sigma - 720*Power(dloga,4)*eps*fnu*N10*Or*sigma + 1800*Power(dloga,2)*fnu*N20*Or*sigma - 1800*Power(dloga,3)*fnu*N20*Or*sigma - 1080*Power(dloga,3)*eps*fnu*N30*Or*sigma + 1080*Power(dloga,4)*eps*fnu*N30*Or*sigma + 150*Power(dloga,2)*Power(eps,2)*Phi0*sigma - 150*Power(dloga,3)*Power(eps,2)*Phi0*sigma - 50*Power(dloga,3)*Power(eps,4)*Phi0*sigma + 165*delta0*Power(dloga,3)*Power(eps,2)*fc*Om*R*sigma + 660*Power(dloga,3)*Power(eps,2)*fnu*N00*Or*R*sigma + 1584*Power(dloga,3)*eps*fnu*N10*Or*R*sigma - 1584*Power(dloga,4)*eps*fnu*N10*Or*R*sigma + 3960*Power(dloga,2)*fnu*N20*Or*R*sigma - 3960*Power(dloga,3)*fnu*N20*Or*R*sigma - 2376*Power(dloga,3)*eps*fnu*N30*Or*R*sigma + 2376*Power(dloga,4)*eps*fnu*N30*Or*R*sigma + 330*Power(dloga,2)*Power(eps,2)*Phi0*R*sigma - 330*Power(dloga,3)*Power(eps,2)*Phi0*R*sigma - 110*Power(dloga,3)*Power(eps,4)*Phi0*R*sigma + 90*delta0*Power(dloga,4)*Power(eps,2)*fc*Om*Power(sigma,2) + 360*Power(dloga,4)*Power(eps,2)*fnu*N00*Or*Power(sigma,2) + 864*Power(dloga,4)*eps*fnu*N10*Or*Power(sigma,2) - 864*Power(dloga,5)*eps*fnu*N10*Or*Power(sigma,2) + 2160*Power(dloga,3)*fnu*N20*Or*Power(sigma,2) - 2160*Power(dloga,4)*fnu*N20*Or*Power(sigma,2) - 1296*Power(dloga,4)*eps*fnu*N30*Or*Power(sigma,2) + 1296*Power(dloga,5)*eps*fnu*N30*Or*Power(sigma,2) + 180*Power(dloga,3)*Power(eps,2)*Phi0*Power(sigma,2) - 180*Power(dloga,4)*Power(eps,2)*Phi0*Power(sigma,2) - 60*Power(dloga,4)*Power(eps,4)*Phi0*Power(sigma,2) + 105*delta0*Power(dloga,4)*Power(eps,2)*fc*Om*R*Power(sigma,2) + 420*Power(dloga,4)*Power(eps,2)*fnu*N00*Or*R*Power(sigma,2) + 1008*Power(dloga,4)*eps*fnu*N10*Or*R*Power(sigma,2) - 1008*Power(dloga,5)*eps*fnu*N10*Or*R*Power(sigma,2) + 2520*Power(dloga,3)*fnu*N20*Or*R*Power(sigma,2) - 2520*Power(dloga,4)*fnu*N20*Or*R*Power(sigma,2) - 1512*Power(dloga,4)*eps*fnu*N30*Or*R*Power(sigma,2) + 1512*Power(dloga,5)*eps*fnu*N30*Or*R*Power(sigma,2) + 210*Power(dloga,3)*Power(eps,2)*Phi0*R*Power(sigma,2) - 210*Power(dloga,4)*Power(eps,2)*Phi0*R*Power(sigma,2) - 70*Power(dloga,4)*Power(eps,4)*Phi0*R*Power(sigma,2) + 15*delta0*Power(dloga,5)*Power(eps,2)*fc*Om*Power(sigma,3) + 60*Power(dloga,5)*Power(eps,2)*fnu*N00*Or*Power(sigma,3) + 144*Power(dloga,5)*eps*fnu*N10*Or*Power(sigma,3) - 144*Power(dloga,6)*eps*fnu*N10*Or*Power(sigma,3) + 360*Power(dloga,4)*fnu*N20*Or*Power(sigma,3) - 360*Power(dloga,5)*fnu*N20*Or*Power(sigma,3) - 216*Power(dloga,5)*eps*fnu*N30*Or*Power(sigma,3) + 216*Power(dloga,6)*eps*fnu*N30*Or*Power(sigma,3) + 30*Power(dloga,4)*Power(eps,2)*Phi0*Power(sigma,3) - 30*Power(dloga,5)*Power(eps,2)*Phi0*Power(sigma,3) - 10*Power(dloga,5)*Power(eps,4)*Phi0*Power(sigma,3) + 15*delta0*Power(dloga,5)*Power(eps,2)*fc*Om*R*Power(sigma,3) + 60*Power(dloga,5)*Power(eps,2)*fnu*N00*Or*R*Power(sigma,3) + 144*Power(dloga,5)*eps*fnu*N10*Or*R*Power(sigma,3) - 144*Power(dloga,6)*eps*fnu*N10*Or*R*Power(sigma,3) + 360*Power(dloga,4)*fnu*N20*Or*R*Power(sigma,3) - 360*Power(dloga,5)*fnu*N20*Or*R*Power(sigma,3) - 216*Power(dloga,5)*eps*fnu*N30*Or*R*Power(sigma,3) + 216*Power(dloga,6)*eps*fnu*N30*Or*R*Power(sigma,3) + 30*Power(dloga,4)*Power(eps,2)*Phi0*R*Power(sigma,3) - 30*Power(dloga,5)*Power(eps,2)*Phi0*R*Power(sigma,3) - 10*Power(dloga,5)*Power(eps,4)*Phi0*R*Power(sigma,3) - 30*cs2*deltab0*dloga*R*(5 + dloga*sigma)*Power(eps + dloga*eps*sigma,2) + 15*deltab0*Power(dloga,2)*Power(eps,2)*fb*Om*(R + dloga*sigma + dloga*R*sigma)*(5 + 6*dloga*sigma + Power(dloga,2)*Power(sigma,2)) + 300*Power(dloga,2)*Power(eps,2)*fgam*Or*R*Theta00 - 150*Power(dloga,2)*Power(eps,2)*sigma*Theta00 + 300*Power(dloga,3)*Power(eps,2)*fgam*Or*sigma*Theta00 + 660*Power(dloga,3)*Power(eps,2)*fgam*Or*R*sigma*Theta00 - 180*Power(dloga,3)*Power(eps,2)*Power(sigma,2)*Theta00 + 360*Power(dloga,4)*Power(eps,2)*fgam*Or*Power(sigma,2)*Theta00 + 420*Power(dloga,4)*Power(eps,2)*fgam*Or*R*Power(sigma,2)*Theta00 - 30*Power(dloga,4)*Power(eps,2)*Power(sigma,3)*Theta00 + 60*Power(dloga,5)*Power(eps,2)*fgam*Or*Power(sigma,3)*Theta00 + 60*Power(dloga,5)*Power(eps,2)*fgam*Or*R*Power(sigma,3)*Theta00 + 720*Power(dloga,2)*eps*fgam*Or*R*Theta10 - 720*Power(dloga,3)*eps*fgam*Or*R*Theta10 - 450*dloga*eps*sigma*Theta10 + 720*Power(dloga,3)*eps*fgam*Or*sigma*Theta10 - 720*Power(dloga,4)*eps*fgam*Or*sigma*Theta10 + 936*Power(dloga,3)*eps*fgam*Or*R*sigma*Theta10 - 936*Power(dloga,4)*eps*fgam*Or*R*sigma*Theta10 - 540*Power(dloga,2)*eps*Power(sigma,2)*Theta10 + 216*Power(dloga,4)*eps*fgam*Or*Power(sigma,2)*Theta10 - 216*Power(dloga,5)*eps*fgam*Or*Power(sigma,2)*Theta10 + 216*Power(dloga,4)*eps*fgam*Or*R*Power(sigma,2)*Theta10 - 216*Power(dloga,5)*eps*fgam*Or*R*Power(sigma,2)*Theta10 - 90*Power(dloga,3)*eps*Power(sigma,3)*Theta10 + 1800*dloga*fgam*Or*R*Theta20 - 1800*Power(dloga,2)*fgam*Or*R*Theta20 + 300*Power(dloga,2)*Power(eps,2)*sigma*Theta20 + 1800*Power(dloga,2)*fgam*Or*sigma*Theta20 - 1800*Power(dloga,3)*fgam*Or*sigma*Theta20 + 2340*Power(dloga,2)*fgam*Or*R*sigma*Theta20 - 2340*Power(dloga,3)*fgam*Or*R*sigma*Theta20 + 360*Power(dloga,3)*Power(eps,2)*Power(sigma,2)*Theta20 + 540*Power(dloga,3)*fgam*Or*Power(sigma,2)*Theta20 - 540*Power(dloga,4)*fgam*Or*Power(sigma,2)*Theta20 + 540*Power(dloga,3)*fgam*Or*R*Power(sigma,2)*Theta20 - 540*Power(dloga,4)*fgam*Or*R*Power(sigma,2)*Theta20 + 60*Power(dloga,4)*Power(eps,2)*Power(sigma,3)*Theta20 - 1080*Power(dloga,2)*eps*fgam*Or*R*Theta30 + 1080*Power(dloga,3)*eps*fgam*Or*R*Theta30 - 1080*Power(dloga,3)*eps*fgam*Or*sigma*Theta30 + 1080*Power(dloga,4)*eps*fgam*Or*sigma*Theta30 - 1404*Power(dloga,3)*eps*fgam*Or*R*sigma*Theta30 + 1404*Power(dloga,4)*eps*fgam*Or*R*sigma*Theta30 - 324*Power(dloga,4)*eps*fgam*Or*Power(sigma,2)*Theta30 + 324*Power(dloga,5)*eps*fgam*Or*Power(sigma,2)*Theta30 - 324*Power(dloga,4)*eps*fgam*Or*R*Power(sigma,2)*Theta30 + 324*Power(dloga,5)*eps*fgam*Or*R*Power(sigma,2)*Theta30 + 180*Power(dloga,2)*fgam*Or*R*sigma*ThetaP00 - 180*Power(dloga,3)*fgam*Or*R*sigma*ThetaP00 + 180*Power(dloga,3)*fgam*Or*Power(sigma,2)*ThetaP00 - 180*Power(dloga,4)*fgam*Or*Power(sigma,2)*ThetaP00 + 180*Power(dloga,3)*fgam*Or*R*Power(sigma,2)*ThetaP00 - 180*Power(dloga,4)*fgam*Or*R*Power(sigma,2)*ThetaP00 - 108*Power(dloga,3)*eps*fgam*Or*R*sigma*ThetaP10 + 108*Power(dloga,4)*eps*fgam*Or*R*sigma*ThetaP10 - 108*Power(dloga,4)*eps*fgam*Or*Power(sigma,2)*ThetaP10 + 108*Power(dloga,5)*eps*fgam*Or*Power(sigma,2)*ThetaP10 - 108*Power(dloga,4)*eps*fgam*Or*R*Power(sigma,2)*ThetaP10 + 108*Power(dloga,5)*eps*fgam*Or*R*Power(sigma,2)*ThetaP10 + 180*Power(dloga,2)*fgam*Or*R*sigma*ThetaP20 - 180*Power(dloga,3)*fgam*Or*R*sigma*ThetaP20 + 180*Power(dloga,3)*fgam*Or*Power(sigma,2)*ThetaP20 - 180*Power(dloga,4)*fgam*Or*Power(sigma,2)*ThetaP20 + 180*Power(dloga,3)*fgam*Or*R*Power(sigma,2)*ThetaP20 - 180*Power(dloga,4)*fgam*Or*R*Power(sigma,2)*ThetaP20 - 108*Power(dloga,3)*eps*fgam*Or*R*sigma*ThetaP30 + 108*Power(dloga,4)*eps*fgam*Or*R*sigma*ThetaP30 - 108*Power(dloga,4)*eps*fgam*Or*Power(sigma,2)*ThetaP30 + 108*Power(dloga,5)*eps*fgam*Or*Power(sigma,2)*ThetaP30 - 108*Power(dloga,4)*eps*fgam*Or*R*Power(sigma,2)*ThetaP30 + 108*Power(dloga,5)*eps*fgam*Or*R*Power(sigma,2)*ThetaP30 - 150*eps*R*ub0 + 150*dloga*eps*R*ub0 - 330*dloga*eps*R*sigma*ub0 + 330*Power(dloga,2)*eps*R*sigma*ub0 - 210*Power(dloga,2)*eps*R*Power(sigma,2)*ub0 + 210*Power(dloga,3)*eps*R*Power(sigma,2)*ub0 - 30*Power(dloga,3)*eps*R*Power(sigma,3)*ub0 + 30*Power(dloga,4)*eps*R*Power(sigma,3)*ub0)/ (30.*eps*(1 + dloga*sigma)*(5 + dloga*sigma)*(R + dloga*sigma + dloga*R*sigma))),
				Rule(Theta0,(-15*deltab0*dloga*Power(eps,2)*fb*Om*(5 + 6*dloga*sigma + Power(dloga,2)*Power(sigma,2)) - 15*delta0*dloga*Power(eps,2)*fc*Om*(5 + 6*dloga*sigma + Power(dloga,2)*Power(sigma,2)) + 2*(36*Power(dloga,4)*eps*fnu*(2*N10 - 3*N30)*Or*Power(sigma,2) + 75*Power(eps,2)*Theta00 + 5*dloga*(5*Power(eps,4)*Phi0 - 3*Power(eps,2)*(10*fnu*N00*Or - 5*Phi0 + 10*fgam*Or*Theta00 - 6*sigma*Theta00) - 15*Power(eps,3)*Theta10 + 180*Or*(fnu*N20 + fgam*Theta20)) + 15*Power(dloga,2)*(2*Power(eps,4)*Phi0*sigma + Power(eps,2)*sigma*(-12*fnu*N00*Or + 6*Phi0 - 12*fgam*Or*Theta00 + sigma*Theta00) - 6*Power(eps,3)*sigma*Theta10 + 12*eps*Or*(2*fnu*N10 - 3*fnu*N30 + 2*fgam*Theta10 - 3*fgam*Theta30) + 6*Or*sigma*(12*fnu*N20 + fgam*(3*Theta20 + ThetaP00 + ThetaP20))) + Power(dloga,3)*sigma*(180*fnu*N20*Or*sigma + 5*Power(eps,4)*Phi0*sigma - 15*Power(eps,2)*sigma*(2*fnu*N00*Or - Phi0 + 2*fgam*Or*Theta00) - 15*Power(eps,3)*sigma*Theta10 + 54*eps*Or*(8*fnu*N10 - 12*fnu*N30 + 2*fgam*Theta10 - 3*fgam*Theta30 - fgam*ThetaP10 - fgam*ThetaP30))))/(30.*Power(eps,2)*(1 + dloga*sigma)*(5 + dloga*sigma))),
				Rule(Theta1,(-15*delta0*Power(dloga,2)*Power(eps,2)*fc*Om*(R + dloga*sigma + dloga*R*sigma)*(5 + 6*dloga*sigma + Power(dloga,2)*Power(sigma,2)) - 15*deltab0*Power(dloga,2)*Power(eps,2)*(5 + 6*dloga*sigma + Power(dloga,2)*Power(sigma,2))*(-2*cs2*R*sigma + fb*Om*(R + dloga*sigma + dloga*R*sigma)) + 2*(36*Power(dloga,6)*eps*fnu*(2*N10 - 3*N30)*Or*(1 + R)*Power(sigma,3) + 225*eps*R*Theta10 + Power(dloga,5)*Power(sigma,2)*(180*fnu*N20*Or*(1 + R)*sigma + 5*Power(eps,4)*Phi0*(1 + R)*sigma - 36*eps*fnu*(2*N10 - 3*N30)*Or*(-6 + R*(-7 + sigma) + sigma) - 15*Power(eps,2)*(1 + R)*sigma*(2*fnu*N00*Or - Phi0 + 2*fgam*Or*Theta00) + 54*eps*fgam*Or*(1 + R)*(2*Theta10 - 3*Theta30 - ThetaP10 - ThetaP30)) - 15*dloga*(60*fnu*N20*Or*R + 60*fgam*Or*R*Theta20 + 5*Power(eps,2)*R*(Phi0 - Theta00 + 2*Theta20) - eps*sigma*(3*(5 + 6*R)*Theta10 + 5*R*ub0)) + Power(dloga,4)*sigma*(5*Power(eps,4)*Phi0*(6 + 7*R)*sigma - 15*Power(eps,2)*sigma* (2*fnu*N00*Or*(6 + 7*R) + Phi0*(-6 + R*(-7 + sigma) + sigma) + 12*fgam*Or*Theta00 + 14*fgam*Or*R*Theta00 - sigma*Theta00 + 2*sigma*Theta20) + 90*Or*sigma*(-2*fnu*N20*(-6 + R*(-7 + sigma) + sigma) + fgam*(1 + R)*(3*Theta20 + ThetaP00 + ThetaP20)) - 3*eps*(12*fnu*(2*N10 - 3*N30)*Or*(-5 + 6*sigma + R*(-11 + 7*sigma)) + 6*fgam*Or*((-20 + 6*sigma + R*(-26 + 6*sigma))*Theta10 - 3*((-10 + 3*sigma + R*(-13 + 3*sigma))*Theta30 + (R*(-1 + sigma) + sigma)*(ThetaP10 + ThetaP30))) + 5*R*Power(sigma,2)*ub0)) + Power(dloga,3)*(5*Power(eps,4)*Phi0*(5 + 11*R)*sigma - 15*Power(eps,2)*sigma* (2*fnu*N00*Or*(5 + 11*R) + Phi0*(-5 - 11*R + 6*sigma + 7*R*sigma) + 10*fgam*Or*Theta00 + 22*fgam*Or*R*Theta00 - 6*sigma*Theta00 - R*sigma*Theta00 + 12*sigma*Theta20 + 2*R*sigma*Theta20) - 90*Or*sigma*(2*fnu*N20*(-5 + 6*sigma + R*(-11 + 7*sigma)) + fgam*((-10 + 3*sigma + R*(-13 + 3*sigma))*Theta20 + (R*(-1 + sigma) + sigma)*(ThetaP00 + ThetaP20))) + eps*(-36*fnu*(2*N10 - 3*N30)*Or*(5*sigma + R*(-5 + 11*sigma)) - 18*fgam*Or* (R*(-20 + 26*sigma)*Theta10 + 10*sigma*(2*Theta10 - 3*Theta30) - 3*R*((-10 + 13*sigma)*Theta30 + sigma*(ThetaP10 + ThetaP30))) + 15*Power(sigma,2)*(3*sigma*Theta10 - 6*R*ub0 + R*sigma*ub0) )) + 5*Power(dloga,2)*(5*Power(eps,4)*Phi0*R - 3*Power(eps,2)*(10*fnu*N00*Or*R + Phi0*(-5*R + 5*sigma + 11*R*sigma) + 10*fgam*Or*R*Theta00 - 5*sigma*Theta00 - 6*R*sigma*Theta00 + 10*sigma*Theta20 + 12*R*sigma*Theta20) - 18*Or*(2*fnu*N20*(5*sigma + R*(-5 + 11*sigma)) + fgam*(10*sigma*Theta20 + R*(-10 + 13*sigma)*Theta20 + R*sigma*(ThetaP00 + ThetaP20))) - 3*eps*(12*fnu*(2*N10 - 3*N30)*Or*R + 12*fgam*Or*R*(2*Theta10 - 3*Theta30) - sigma*(3*(6 + R)*sigma*Theta10 - 5*R*ub0 + 6*R*sigma*ub0)))))/ (90.*eps*(1 + dloga*sigma)*(5 + dloga*sigma)*(R + dloga*sigma + dloga*R*sigma))),
				Rule(Theta2,(50*Theta20 + 5*dloga*(eps*(4*Theta10 - 6*Theta30) + sigma*(3*Theta20 + ThetaP00 + ThetaP20)) + 3*Power(dloga,2)*eps*sigma*(2*Theta10 - 3*Theta30 - ThetaP10 - ThetaP30))/ (10.*(1 + dloga*sigma)*(5 + dloga*sigma))),
				Rule(Theta3,(dloga*eps*etaaH*Theta20 - 4*dloga*Theta30 + etaaH*Theta30)/(etaaH + dloga*etaaH*sigma)),
				Rule(ThetaP0,(10*ThetaP00 - 10*dloga*eps*ThetaP10 + dloga*sigma*(5*Theta20 + 7*ThetaP00 + 5*ThetaP20) + Power(dloga,2)*eps*sigma*(2*Theta10 - 3*Theta30 - 5*ThetaP10 - 3*ThetaP30))/ (2.*(1 + dloga*sigma)*(5 + dloga*sigma))),
				Rule(ThetaP1,(3*ThetaP10 + dloga*eps*(ThetaP00 - 2*ThetaP20))/(3 + 3*dloga*sigma)),
				Rule(ThetaP2,(25*ThetaP20 + 5*dloga*(sigma*(Theta20 + ThetaP00 + 2*ThetaP20) + eps*(2*ThetaP10 - 3*ThetaP30)) + Power(dloga,2)*eps*sigma*(2*Theta10 - 3*Theta30 - ThetaP10 - 6*ThetaP30))/ (5.*(1 + dloga*sigma)*(5 + dloga*sigma))),
				Rule(ThetaP3,(dloga*eps*etaaH*ThetaP20 - 4*dloga*ThetaP30 + etaaH*ThetaP30)/(etaaH + dloga*etaaH*sigma)),
				Rule(N0,(-30*dloga*Power(eps,3)*N10*(5 + 6*dloga*sigma + Power(dloga,2)*Power(sigma,2)) + 10*dloga*Power(eps,4)*Phi0*(5 + 6*dloga*sigma + Power(dloga,2)*Power(sigma,2)) - 15*Power(eps,2)*(5 + 6*dloga*sigma + Power(dloga,2)*Power(sigma,2))*(N00*(-2 + 4*dloga*fnu*Or) + dloga*(deltab0*fb*Om + delta0*fc*Om - 2*Phi0 + 4*fgam*Or*Theta00)) + 180*dloga*Or*(2*fnu*N20*(5 + 6*dloga*sigma + Power(dloga,2)*Power(sigma,2)) + fgam*(10 + 3*dloga*sigma)*Theta20 + dloga*fgam*sigma*(ThetaP00 + ThetaP20)) + 36*Power(dloga,2)*eps*Or*(2*fnu*(2*N10 - 3*N30)*(5 + 6*dloga*sigma + Power(dloga,2)*Power(sigma,2)) + fgam*(20 + 6*dloga*sigma)*Theta10 - 3*fgam*((10 + 3*dloga*sigma)*Theta30 + dloga*sigma*(ThetaP10 + ThetaP30))))/(30.*Power(eps,2)*(1 + dloga*sigma)*(5 + dloga*sigma))),
				Rule(N1,(450*eps*N10 + 72*Power(dloga,5)*eps*fnu*(2*N10 - 3*N30)*Or*Power(sigma,2) + 30*dloga*(5*Power(eps,2)*(N00 - 2*N20 - Phi0) + 18*eps*N10*sigma - 60*Or*(fnu*N20 + fgam*Theta20)) - 5*Power(dloga,2)*(15*deltab0*Power(eps,2)*fb*Om + 15*delta0*Power(eps,2)*fc*Om - 10*Power(eps,4)*Phi0 + 6*Power(eps,2)*(10*fnu*N00*Or - 5*Phi0 - 6*N00*sigma + 12*N20*sigma + 6*Phi0*sigma + 10*fgam*Or*Theta00) + 18*eps*(4*fnu*(2*N10 - 3*N30)*Or - N10*Power(sigma,2) + 8*fgam*Or*Theta10 - 12*fgam*Or*Theta30) + 36*Or*(2*fnu*N20*(-5 + 6*sigma) + fgam*(-10 + 3*sigma)*Theta20 + fgam*sigma*(ThetaP00 + ThetaP20))) + Power(dloga,4)*sigma*(360*fnu*N20*Or*sigma + 10*Power(eps,4)*Phi0*sigma - 15*Power(eps,2)*sigma*(deltab0*fb*Om + delta0*fc*Om + 4*fnu*N00*Or - 2*Phi0 + 4*fgam*Or*Theta00) - 36*eps*Or*(2*fnu*(2*N10 - 3*N30)*(-6 + sigma) + 3*fgam*(-2*Theta10 + 3*Theta30 + ThetaP10 + ThetaP30))) + 6*Power(dloga,3)*(10*Power(eps,4)*Phi0*sigma - 5*Power(eps,2)*sigma*(3*deltab0*fb*Om + 3*delta0*fc*Om + 12*fnu*N00*Or - 6*Phi0 - N00*sigma + 2*N20*sigma + Phi0*sigma + 12*fgam*Or*Theta00) + 30*Or*sigma*(-2*fnu*N20*(-6 + sigma) + fgam*(3*Theta20 + ThetaP00 + ThetaP20)) - 6*eps*Or*(2*fnu*(2*N10 - 3*N30)*(-5 + 6*sigma) + fgam*(-20 + 6*sigma)*Theta10 - 3*fgam*((-10 + 3*sigma)*Theta30 + sigma*(ThetaP10 + ThetaP30)))))/(90.*eps*(1 + dloga*sigma)*(5 + dloga*sigma))),
				Rule(N2,N20 + (dloga*eps*(2*N10 - 3*N30))/5.), Rule(N3,dloga*eps*N20 + N30 - (4*dloga*N30)/etaaH));
		eta += dloga / (a * H);
		loga += dloga;
	}

}

void initial_conditions(state &u, double k, double a) {
	const auto eps = k / (a * Hubble(a));
	const auto Psii = 1.0;
	u[Thetai + 0] = -0.5 * Psii;
	u[Thetai + 1] = eps / 6.0 * Psii;
	u[Ni + 0] = u[Thetai + 0];
	u[Ni + 1] = u[Thetai + 1];
	u[Ni + 2] = eps * eps / 30.0 * Psii;
	u[deltai] = 3.0 * u[Thetai + 0];
	u[ui] = 3.0 * u[Thetai + 1];
	u[deltabi] = 3.0 * u[Thetai + 0];
	u[ubi] = 3.0 * u[Thetai + 1];
	u[Phii] = -12.0 * omega_nu / (a * omega_m + omega_r) * (u[Ni + 2]) / (eps * eps) - Psii;
	for (int l = 0; l < LMAX; l++) {
		u[ThetaPi + l] = 0.0;
	}
	for (int l = 2; l < LMAX; l++) {
		u[Thetai + l] = 0.0;
	}
	for (int l = 3; l < LMAX; l++) {
		u[Ni + l] = 0.0;
	}
}

int main() {
	std::array<std::array<double, 3>, 3> A = { { { -1, -2, -3 }, { 4, 5, 6 }, { 7, 8, 0 } } };
//	compute_eigenvalues<3>(A)
	//	printf("%e\n", matrix_determinant<3>(A));

//	const auto f = [](cmplx x) {
//		return (x + 1.0) * (x + 1.0) * x;
//	};
//	const auto dfdx = [](cmplx x) {
//		return (x + 1.0) * (x + 1.0) + 2.0 * (x + 1.0) * x;
//	};
//	auto roots = find_roots(f, 3);
//	printf("%e %e %e %e %e %e \n", roots[0].real(), roots[0].imag(), roots[1].real(), roots[1].imag(), roots[2].real(), roots[2].imag());

	double amin = 1e-8;
	double amax = 1.0;
	std::function<double(double)> cs, thomson;
	zero_order_universe(amin, amax, cs, thomson);
	for (double k = 1e-4; k <= 1; k *= 1.1) {
		state u;
		initial_conditions(u, k, amin);
		advance(u, k, amin, amax, cs, thomson);
	//	break;
		printf("%e %e\n", k, k * std::pow(u[Phii], 2));
	}
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
//	std::function<double(double)> cs, thomson;
//	zero_order_universe(1.0e-8, 1, cs, thomson);
//	for( double a = 1e-8; a <= 1.0; a *= 1.01) {
//		printf( "%e %e %e\n", 1/a-1, cs(std::log(a)), thomson(std::log(a)));
//	}
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
