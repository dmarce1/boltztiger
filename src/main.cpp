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
//#define LMAX 8
//
//#define Phii 0
//#define deltai 1
//#define ui 2
//#define deltabi 3
//#define ubi 4
//#define Ni 5
//#define Thetai (5+LMAX)
//#define ThetaPi (5+LMAX*2)
//#define NF (5+LMAX*3)
//using state = std::array<double,NF>;
//
//double compute_Nele(double eta, double beta);
//double compute_Npos(double eta, double beta);

double rho_b(double a) {
	using namespace constants;
	return (3.0 * omega_b * H0cgs * H0cgs * h * h) / (8.0 * M_PI * G * a * a * a);
}

void saha(double rho, double T, double &H, double &Hp, double &He, double &Hep, double &Hepp, double &ne) {
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

#define Sqrt(a) std::sqrt(a)

void chemistry_update(double &H, double &Hp, double &He, double &Hep, double &Hepp, double &ne, double T, double a, double dt) {
	bool use_saha;
	double H1 = H;
	double Hp1 = Hp;
	double He1 = He;
	double Hep1 = Hep;
	double Hepp1 = Hepp;
	double ne1 = ne;
	double H0 = H;
	double Hp0 = Hp;
	double He0 = He;
	double Hep0 = Hep;
	double Hepp0 = Hepp;
	double ne0 = ne;
	using namespace constants;
	double rho = ((H + Hp) + (He + Hep + Hepp) * 4) * mh;
	if (ne > (H + Hp)) {
		saha(rho, T, H1, Hp1, He1, Hep1, Hepp1, ne1);
		if (ne1 > (H1 + Hp1)) {
			use_saha = true;
		} else {
			use_saha = false;
		}
		use_saha = true;
	} else {
		use_saha = false;
	}
	if (use_saha) {
		H = H1;
		Hp = Hp1;
		He = He1;
		Hep = Hep1;
		Hepp = Hepp1;
		ne = ne1;
	} else {
		double nH = H + Hp;
		double x0 = Hp / nH;
		const auto dxdt = [=](double x0, double dt) {
			double hubble = Hubble(a) * H0cgs / H0;
			using namespace constants;
			const auto B1 = 13.6 * evtoerg;
			const auto phi2 = std::max(0.448 * std::log(B1 / (kb * T)), 0.0);
			const auto alpha2 = 64.0 * M_PI / std::sqrt(27.0 * M_PI) * B1 * 2.0 * std::pow(hbar, 2) / std::pow(me * clight, 3) / std::sqrt(kb * T / B1) * phi2;
			const auto beta = std::pow((me * kb * T) / (2 * M_PI * hbar * hbar), 1.5) * std::exp(-B1 / (kb * T)) * alpha2;
			const auto lambda_a = 8.0 * M_PI * hbar * clight / (3.0 * B1);
			const auto num = hplanck * clight / lambda_a / kb / T;
			const auto beta2 = beta * std::exp(std::min(num, 80.0));
			const auto La = 8.0 * M_PI * hubble / (a * std::pow(lambda_a, 3) * nH);
			const auto L2s = 8.227;
			const auto func = [=](double x) {
				return (x - (dt * (L2s + La / (1 - x)) * (beta * (1 - x) - alpha2 * nH * Power(x, 2))) / (beta2 + L2s + La / (1 - x)) - x0) * (1 - x);
			};
			double x = find_root(func);
			return (x - x0) / dt;
		};
		double gam = 1.0 - 1.0 / std::sqrt(2);
		double dx1 = dxdt(x0, gam * dt);
		double dx2 = dxdt(x0 + (1 - 2 * gam) * dx1 * dt, gam * dt);
		double x = (x0 + 0.5 * (dx1 * dt + dx2 * dt));
		He = He0 + Hep0 + Hepp0;
		Hep = 0.0;
		Hepp = 0.0;
		H = (1.0 - x) * (H0 + Hp0);
		Hp = x * (H0 + Hp0);
		ne = Hp;
	}
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
	constexpr int N = 10000;
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
	double P0, rho, P, rho1, P1;
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
		rho0 = rho1;
		P0 = P1;
		rho1 = rho;
		P1 = P;
		rho = rho_b(a);
		Trad = Theta * 2.73 / a;
		const auto dt = dloga / hubble;
		double cs;
		const auto sigma_T = 6.65e-25;
		double sigmaC;
		if (i > 0) {
			constexpr auto gamma = 1.0 - 1.0 / std::sqrt(2);
			chemistry_update(H, Hp, He, Hep, Hepp, ne, Tgas, a, 0.5 * dt);
			const auto mu = (H + Hp + 4 * He + 4 * Hep + 4 * Hepp) * mh / (H + Hp + He + Hep + Hepp + ne);
			sigmaC = mu / me * clight * (8.0 / 3.0) * omega_gam / (a * omega_m) * sigma_T * ne / hubble;
			const auto dTgasdT1 = ((Tgas + gamma * dloga * sigmaC * Trad) / (1 + gamma * dloga * (2 + sigmaC)) - Tgas) / (gamma * dloga);
			const auto T1 = Tgas + (1 - 2 * gamma) * dTgasdT1 * dloga;
			const auto dTgasdT2 = ((T1 + gamma * dloga * sigmaC * Trad) / (1 + gamma * dloga * (2 + sigmaC)) - T1) / (gamma * dloga);
			Tgas += 0.5 * (dTgasdT1 + dTgasdT2) * dloga;
			chemistry_update(H, Hp, He, Hep, Hepp, ne, Tgas, a, 0.5 * dt);
		}
		const auto n = H + Hp + He + Hep + Hepp;
		P = kb * (n + ne) * Tgas;
		if (i > 1 && i <= N) {
			sound_speed[i - 1] = std::sqrt((P - P0) / (rho - rho0)) / clight;
		} else if (i == 1) {
			sound_speed[0] = std::sqrt((P - P1) / (rho - rho1)) / clight;
		}
		if (i == N) {
			sound_speed[N - 1] = std::sqrt((P - P1) / (rho - rho1)) / clight;
		}
		sound_speed[i] = cs / clight;
		thomson[i] = clight * sigma_T * ne / hubble;
//		printf("%e %e %e %e %e %e\n", 1 / a - 1, Tgas, Trad, (Hp + Hep + 2 * Hepp) / (H + Hp + 2 * He + 2 * Hep + 2 * Hepp), thomson[i], sigmaC);
		t += dt;
		last_a = a;
	}
	csfunc = build_interpolation_function(sound_speed, amin, amax);
	thomsonfunc = build_interpolation_function(thomson, amin, amax);
//	abort();
}

//void advance(state &U, double k, double a0, double a1, std::function<double(double)> cs, std::function<double(double)> thomson) {
//	const double logamin = std::log(a0);
//	const double logamax = std::log(a1);
//	double loga = logamin;
//	double eta = 1.0 / (a0 * Hubble(a0));
//	double cs2, etaaH, sigma;
//	while (loga < logamax) {
//		double a, H, Om, Or, eps, fb, fc, fgam, sigma, R, Onu, Ogam, fnu;
//		const auto compute_parameters = [&](double loga) {
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
//			sigma = thomson(loga);
//			cs2 = std::pow(cs(loga), 2);
//			etaaH = eta * a * H;
//		};
//		compute_parameters(loga);
//		const auto lambda_max = std::max(std::max(std::max(1.0 + cs2 * eps, eps + 4.0 / etaaH), 0.5 * Om + 2.0 * Or), (1 + eps * eps / 3.0));
//		const auto dloga = std::min(4 / lambda_max, logamax - loga);
//
//		const auto dudt_exp = [&](state U) {
//			state dudt;
//			dudt[Phii] = 0.0;
//			dudt[deltai] = -eps * U[ui];
//			dudt[ui] = -U[ui];
//			dudt[deltabi] = -eps * U[ubi];
//			dudt[ubi] = -U[ubi] + eps * cs2 * U[deltabi];
//			dudt[Thetai + 0] = -eps * U[Thetai + 1];
//			dudt[ThetaPi + 0] = -eps * U[ThetaPi + 1];
//			dudt[Ni + 0] = -eps * U[Ni + 1];
//			for (int l = 1; l < LMAX - 1; l++) {
//				dudt[Thetai + l] = (eps / (2 * l + 1)) * (l * U[Thetai + l - 1] - (l + 1) * U[Thetai + l + 1]);
//				dudt[ThetaPi + l] = (eps / (2 * l + 1)) * (l * U[ThetaPi + l - 1] - (l + 1) * U[ThetaPi + l + 1]);
//				dudt[Ni + l] = (eps / (2 * l + 1)) * (l * U[Ni + l - 1] - (l + 1) * U[Ni + l + 1]);
//			}
//			dudt[Thetai + LMAX - 1] = eps * U[Thetai + LMAX - 2];
//			dudt[ThetaPi + LMAX - 1] = eps * U[ThetaPi + LMAX - 2];
//			dudt[Ni + LMAX - 1] = eps * U[Ni + LMAX - 2];
//			return dudt;
//		};
//
//		const auto dudt_imp =
//				[&](state U, double dloga) {
//					double Phi, delta, u, deltab, ub, Theta0, Theta1, Theta2, ThetaP0, ThetaP1, ThetaP2, N0, N1, N2;
//					const auto Phi0 = U[Phii];
//					const auto delta0 = U[deltai];
//					const auto u0 = U[ui];
//					const auto deltab0 = U[deltabi];
//					const auto ub0 = U[ubi];
//					const auto Theta00 = U[Thetai + 0];
//					const auto Theta10 = U[Thetai + 1];
//					const auto Theta20 = U[Thetai + 2];
//					const auto ThetaP00 = U[ThetaPi + 0];
//					const auto ThetaP10 = U[ThetaPi + 1];
//					const auto ThetaP20 = U[ThetaPi + 2];
//					const auto N00 = U[Ni + 0];
//					const auto N10 = U[Ni + 1];
//					const auto N20 = U[Ni + 2];
//					state dudt;
//#define List(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13) a0; a1; a2; a3; a4; a5; a6; a7; a8; a9; a10; a11; a12; a13
//#define Rule(a,b) a = b
//					List(
//							Rule(Phi,(3*(40*dloga*Power(eps,2)*fnu*N00*Or - 240*dloga*fnu*N20*Or + 20*Power(eps,2)*Phi0 + 30*dloga*Power(eps,2)*fb*Om*Phi0 + 30*dloga*Power(eps,2)*fc*Om*Phi0 + 40*dloga*Power(eps,2)*fgam*Or*Phi0 + 40*dloga*Power(eps,2)*fnu*Or*Phi0 + 52*Power(dloga,2)*Power(eps,2)*fnu*N00*Or*sigma - 312*Power(dloga,2)*fnu*N20*Or*sigma + 26*dloga*Power(eps,2)*Phi0*sigma + 39*Power(dloga,2)*Power(eps,2)*fb*Om*Phi0*sigma + 39*Power(dloga,2)*Power(eps,2)*fc*Om*Phi0*sigma + 52*Power(dloga,2)*Power(eps,2)*fgam*Or*Phi0*sigma + 52*Power(dloga,2)*Power(eps,2)*fnu*Or*Phi0*sigma + 12*Power(dloga,3)*Power(eps,2)*fnu*N00*Or*Power(sigma,2) - 72*Power(dloga,3)*fnu*N20*Or*Power(sigma,2) + 6*Power(dloga,2)*Power(eps,2)*Phi0*Power(sigma,2) + 9*Power(dloga,3)*Power(eps,2)*fb*Om*Phi0*Power(sigma,2) + 9*Power(dloga,3)*Power(eps,2)*fc*Om*Phi0*Power(sigma,2) + 12*Power(dloga,3)*Power(eps,2)*fgam*Or*Phi0*Power(sigma,2) + 12*Power(dloga,3)*Power(eps,2)*fnu*Or*Phi0*Power(sigma,2) + deltab0*dloga*Power(eps,2)*fb*Om*(10 + 13*dloga*sigma + 3*Power(dloga,2)*Power(sigma,2)) + delta0*dloga*Power(eps,2)*fc*Om*(10 + 13*dloga*sigma + 3*Power(dloga,2)*Power(sigma,2)) + 40*dloga*Power(eps,2)*fgam*Or*Theta00 + 52*Power(dloga,2)*Power(eps,2)*fgam*Or*sigma*Theta00 + 12*Power(dloga,3)*Power(eps,2)*fgam*Or*Power(sigma,2)*Theta00 - 240*dloga*fgam*Or*Theta20 - 96*Power(dloga,2)*fgam*Or*sigma*Theta20 - 24*Power(dloga,2)*fgam*Or*sigma*ThetaP00 - 24*Power(dloga,2)*fgam*Or*sigma*ThetaP20))/ (Power(eps,2)*(6 + dloga*(6 + 2*Power(eps,2) + 9*fb*Om + 9*fc*Om + 12*fgam*Or + 12*fnu*Or))*(1 + dloga*sigma)*(10 + 3*dloga*sigma))),
//							Rule(delta,(delta0*Power(eps,2)*(6 + dloga*(6 + 2*Power(eps,2) + 9*fb*Om + 12*fgam*Or + 12*fnu*Or))*(10 + 13*dloga*sigma + 3*Power(dloga,2)*Power(sigma,2)) - 3*dloga*(3*deltab0*Power(eps,2)*fb*Om*(10 + 13*dloga*sigma + 3*Power(dloga,2)*Power(sigma,2)) - 2*Power(eps,4)*Phi0*(10 + 13*dloga*sigma + 3*Power(dloga,2)*Power(sigma,2)) + 6*Power(eps,2)*(10 + 13*dloga*sigma + 3*Power(dloga,2)*Power(sigma,2))*(2*fnu*N00*Or - Phi0 + 2*fgam*Or*Theta00) - 72*Or*(fnu*N20*(10 + 13*dloga*sigma + 3*Power(dloga,2)*Power(sigma,2)) + 2*fgam*(5 + 2*dloga*sigma)*Theta20 + dloga*fgam*sigma*(ThetaP00 + ThetaP20))))/ (Power(eps,2)*(6 + dloga*(6 + 2*Power(eps,2) + 9*fb*Om + 9*fc*Om + 12*fgam*Or + 12*fnu*Or))*(1 + dloga*sigma)*(10 + 3*dloga*sigma))),
//							Rule(u,-((30*delta0*Power(dloga,2)*Power(eps,2)*fc*Om + 120*Power(dloga,2)*Power(eps,2)*fnu*N00*Or + 720*dloga*fnu*N20*Or + 240*Power(dloga,2)*Power(eps,2)*fnu*N20*Or + 1080*Power(dloga,2)*fb*fnu*N20*Om*Or + 1080*Power(dloga,2)*fc*fnu*N20*Om*Or + 1440*Power(dloga,2)*fgam*fnu*N20*Power(Or,2) + 1440*Power(dloga,2)*Power(fnu,2)*N20*Power(Or,2) + 60*dloga*Power(eps,2)*Phi0 + 90*Power(dloga,2)*Power(eps,2)*fb*Om*Phi0 + 90*Power(dloga,2)*Power(eps,2)*fc*Om*Phi0 + 120*Power(dloga,2)*Power(eps,2)*fgam*Or*Phi0 + 120*Power(dloga,2)*Power(eps,2)*fnu*Or*Phi0 + 39*delta0*Power(dloga,3)*Power(eps,2)*fc*Om*sigma + 156*Power(dloga,3)*Power(eps,2)*fnu*N00*Or*sigma + 936*Power(dloga,2)*fnu*N20*Or*sigma + 312*Power(dloga,3)*Power(eps,2)*fnu*N20*Or*sigma + 1404*Power(dloga,3)*fb*fnu*N20*Om*Or*sigma + 1404*Power(dloga,3)*fc*fnu*N20*Om*Or*sigma + 1872*Power(dloga,3)*fgam*fnu*N20*Power(Or,2)*sigma + 1872*Power(dloga,3)*Power(fnu,2)*N20*Power(Or,2)*sigma + 78*Power(dloga,2)*Power(eps,2)*Phi0*sigma + 117*Power(dloga,3)*Power(eps,2)*fb*Om*Phi0*sigma + 117*Power(dloga,3)*Power(eps,2)*fc*Om*Phi0*sigma + 156*Power(dloga,3)*Power(eps,2)*fgam*Or*Phi0*sigma + 156*Power(dloga,3)*Power(eps,2)*fnu*Or*Phi0*sigma + 9*delta0*Power(dloga,4)*Power(eps,2)*fc*Om*Power(sigma,2) + 36*Power(dloga,4)*Power(eps,2)*fnu*N00*Or*Power(sigma,2) + 216*Power(dloga,3)*fnu*N20*Or*Power(sigma,2) + 72*Power(dloga,4)*Power(eps,2)*fnu*N20*Or*Power(sigma,2) + 324*Power(dloga,4)*fb*fnu*N20*Om*Or*Power(sigma,2) + 324*Power(dloga,4)*fc*fnu*N20*Om*Or*Power(sigma,2) + 432*Power(dloga,4)*fgam*fnu*N20*Power(Or,2)*Power(sigma,2) + 432*Power(dloga,4)*Power(fnu,2)*N20*Power(Or,2)*Power(sigma,2) + 18*Power(dloga,3)*Power(eps,2)*Phi0*Power(sigma,2) + 27*Power(dloga,4)*Power(eps,2)*fb*Om*Phi0*Power(sigma,2) + 27*Power(dloga,4)*Power(eps,2)*fc*Om*Phi0*Power(sigma,2) + 36*Power(dloga,4)*Power(eps,2)*fgam*Or*Phi0*Power(sigma,2) + 36*Power(dloga,4)*Power(eps,2)*fnu*Or*Phi0*Power(sigma,2) + 3*deltab0*Power(dloga,2)*Power(eps,2)*fb*Om*(10 + 13*dloga*sigma + 3*Power(dloga,2)*Power(sigma,2)) + 120*Power(dloga,2)*Power(eps,2)*fgam*Or*Theta00 + 156*Power(dloga,3)*Power(eps,2)*fgam*Or*sigma*Theta00 + 36*Power(dloga,4)*Power(eps,2)*fgam*Or*Power(sigma,2)*Theta00 + 720*dloga*fgam*Or*Theta20 + 240*Power(dloga,2)*Power(eps,2)*fgam*Or*Theta20 + 1080*Power(dloga,2)*fb*fgam*Om*Or*Theta20 + 1080*Power(dloga,2)*fc*fgam*Om*Or*Theta20 + 1440*Power(dloga,2)*Power(fgam,2)*Power(Or,2)*Theta20 + 1440*Power(dloga,2)*fgam*fnu*Power(Or,2)*Theta20 + 288*Power(dloga,2)*fgam*Or*sigma*Theta20 + 96*Power(dloga,3)*Power(eps,2)*fgam*Or*sigma*Theta20 + 432*Power(dloga,3)*fb*fgam*Om*Or*sigma*Theta20 + 432*Power(dloga,3)*fc*fgam*Om*Or*sigma*Theta20 + 576*Power(dloga,3)*Power(fgam,2)*Power(Or,2)*sigma*Theta20 + 576*Power(dloga,3)*fgam*fnu*Power(Or,2)*sigma*Theta20 + 72*Power(dloga,2)*fgam*Or*sigma*ThetaP00 + 24*Power(dloga,3)*Power(eps,2)*fgam*Or*sigma*ThetaP00 + 108*Power(dloga,3)*fb*fgam*Om*Or*sigma*ThetaP00 + 108*Power(dloga,3)*fc*fgam*Om*Or*sigma*ThetaP00 + 144*Power(dloga,3)*Power(fgam,2)*Power(Or,2)*sigma*ThetaP00 + 144*Power(dloga,3)*fgam*fnu*Power(Or,2)*sigma*ThetaP00 + 72*Power(dloga,2)*fgam*Or*sigma*ThetaP20 + 24*Power(dloga,3)*Power(eps,2)*fgam*Or*sigma*ThetaP20 + 108*Power(dloga,3)*fb*fgam*Om*Or*sigma*ThetaP20 + 108*Power(dloga,3)*fc*fgam*Om*Or*sigma*ThetaP20 + 144*Power(dloga,3)*Power(fgam,2)*Power(Or,2)*sigma*ThetaP20 + 144*Power(dloga,3)*fgam*fnu*Power(Or,2)*sigma*ThetaP20 - cs2*dloga*(deltab0*Power(eps,2)*(6 + dloga*(6 + 2*Power(eps,2) + 9*fc*Om + 12*fgam*Or + 12*fnu*Or))*(10 + 13*dloga*sigma + 3*Power(dloga,2)*Power(sigma,2)) - 3*dloga*(3*delta0*Power(eps,2)*fc*Om*(10 + 13*dloga*sigma + 3*Power(dloga,2)*Power(sigma,2)) - 2*Power(eps,4)*Phi0*(10 + 13*dloga*sigma + 3*Power(dloga,2)*Power(sigma,2)) + 6*Power(eps,2)*(10 + 13*dloga*sigma + 3*Power(dloga,2)*Power(sigma,2))*(2*fnu*N00*Or - Phi0 + 2*fgam*Or*Theta00) - 72*Or*(fnu*N20*(10 + 13*dloga*sigma + 3*Power(dloga,2)*Power(sigma,2)) + 2*fgam*(5 + 2*dloga*sigma)*Theta20 + dloga*fgam*sigma*(ThetaP00 + ThetaP20)))) - 60*eps*u0 - 60*dloga*eps*u0 - 20*dloga*Power(eps,3)*u0 - 90*dloga*eps*fb*Om*u0 - 90*dloga*eps*fc*Om*u0 - 120*dloga*eps*fgam*Or*u0 - 120*dloga*eps*fnu*Or*u0 - 78*dloga*eps*sigma*u0 - 78*Power(dloga,2)*eps*sigma*u0 - 26*Power(dloga,2)*Power(eps,3)*sigma*u0 - 117*Power(dloga,2)*eps*fb*Om*sigma*u0 - 117*Power(dloga,2)*eps*fc*Om*sigma*u0 - 156*Power(dloga,2)*eps*fgam*Or*sigma*u0 - 156*Power(dloga,2)*eps*fnu*Or*sigma*u0 - 18*Power(dloga,2)*eps*Power(sigma,2)*u0 - 18*Power(dloga,3)*eps*Power(sigma,2)*u0 - 6*Power(dloga,3)*Power(eps,3)*Power(sigma,2)*u0 - 27*Power(dloga,3)*eps*fb*Om*Power(sigma,2)*u0 - 27*Power(dloga,3)*eps*fc*Om*Power(sigma,2)*u0 - 36*Power(dloga,3)*eps*fgam*Or*Power(sigma,2)*u0 - 36*Power(dloga,3)*eps*fnu*Or*Power(sigma,2)*u0)/ (eps*(6 + dloga*(6 + 2*Power(eps,2) + 9*fb*Om + 9*fc*Om + 12*fgam*Or + 12*fnu*Or))*(1 + dloga*sigma)*(10 + 3*dloga*sigma)))),
//							Rule(deltab,(deltab0*Power(eps,2)*(6 + dloga*(6 + 2*Power(eps,2) + 9*fc*Om + 12*fgam*Or + 12*fnu*Or))*(10 + 13*dloga*sigma + 3*Power(dloga,2)*Power(sigma,2)) - 3*dloga*(3*delta0*Power(eps,2)*fc*Om*(10 + 13*dloga*sigma + 3*Power(dloga,2)*Power(sigma,2)) - 2*Power(eps,4)*Phi0*(10 + 13*dloga*sigma + 3*Power(dloga,2)*Power(sigma,2)) + 6*Power(eps,2)*(10 + 13*dloga*sigma + 3*Power(dloga,2)*Power(sigma,2))*(2*fnu*N00*Or - Phi0 + 2*fgam*Or*Theta00) - 72*Or*(fnu*N20*(10 + 13*dloga*sigma + 3*Power(dloga,2)*Power(sigma,2)) + 2*fgam*(5 + 2*dloga*sigma)*Theta20 + dloga*fgam*sigma*(ThetaP00 + ThetaP20))))/ (Power(eps,2)*(6 + dloga*(6 + 2*Power(eps,2) + 9*fb*Om + 9*fc*Om + 12*fgam*Or + 12*fnu*Or))*(1 + dloga*sigma)*(10 + 3*dloga*sigma))),
//							Rule(ub,-((120*Power(dloga,2)*Power(eps,2)*fnu*N00*Or*R + 720*dloga*fnu*N20*Or*R + 240*Power(dloga,2)*Power(eps,2)*fnu*N20*Or*R + 1080*Power(dloga,2)*fb*fnu*N20*Om*Or*R + 1080*Power(dloga,2)*fc*fnu*N20*Om*Or*R + 1440*Power(dloga,2)*fgam*fnu*N20*Power(Or,2)*R + 1440*Power(dloga,2)*Power(fnu,2)*N20*Power(Or,2)*R + 60*dloga*Power(eps,2)*Phi0*R + 90*Power(dloga,2)*Power(eps,2)*fb*Om*Phi0*R + 90*Power(dloga,2)*Power(eps,2)*fc*Om*Phi0*R + 120*Power(dloga,2)*Power(eps,2)*fgam*Or*Phi0*R + 120*Power(dloga,2)*Power(eps,2)*fnu*Or*Phi0*R + 120*Power(dloga,3)*Power(eps,2)*fnu*N00*Or*sigma + 720*Power(dloga,2)*fnu*N20*Or*sigma + 240*Power(dloga,3)*Power(eps,2)*fnu*N20*Or*sigma + 1080*Power(dloga,3)*fb*fnu*N20*Om*Or*sigma + 1080*Power(dloga,3)*fc*fnu*N20*Om*Or*sigma + 1440*Power(dloga,3)*fgam*fnu*N20*Power(Or,2)*sigma + 1440*Power(dloga,3)*Power(fnu,2)*N20*Power(Or,2)*sigma + 60*Power(dloga,2)*Power(eps,2)*Phi0*sigma + 90*Power(dloga,3)*Power(eps,2)*fb*Om*Phi0*sigma + 90*Power(dloga,3)*Power(eps,2)*fc*Om*Phi0*sigma + 120*Power(dloga,3)*Power(eps,2)*fgam*Or*Phi0*sigma + 120*Power(dloga,3)*Power(eps,2)*fnu*Or*Phi0*sigma + 276*Power(dloga,3)*Power(eps,2)*fnu*N00*Or*R*sigma + 1656*Power(dloga,2)*fnu*N20*Or*R*sigma + 552*Power(dloga,3)*Power(eps,2)*fnu*N20*Or*R*sigma + 2484*Power(dloga,3)*fb*fnu*N20*Om*Or*R*sigma + 2484*Power(dloga,3)*fc*fnu*N20*Om*Or*R*sigma + 3312*Power(dloga,3)*fgam*fnu*N20*Power(Or,2)*R*sigma + 3312*Power(dloga,3)*Power(fnu,2)*N20*Power(Or,2)*R*sigma + 138*Power(dloga,2)*Power(eps,2)*Phi0*R*sigma + 207*Power(dloga,3)*Power(eps,2)*fb*Om*Phi0*R*sigma + 207*Power(dloga,3)*Power(eps,2)*fc*Om*Phi0*R*sigma + 276*Power(dloga,3)*Power(eps,2)*fgam*Or*Phi0*R*sigma + 276*Power(dloga,3)*Power(eps,2)*fnu*Or*Phi0*R*sigma + 156*Power(dloga,4)*Power(eps,2)*fnu*N00*Or*Power(sigma,2) + 936*Power(dloga,3)*fnu*N20*Or*Power(sigma,2) + 312*Power(dloga,4)*Power(eps,2)*fnu*N20*Or*Power(sigma,2) + 1404*Power(dloga,4)*fb*fnu*N20*Om*Or*Power(sigma,2) + 1404*Power(dloga,4)*fc*fnu*N20*Om*Or*Power(sigma,2) + 1872*Power(dloga,4)*fgam*fnu*N20*Power(Or,2)*Power(sigma,2) + 1872*Power(dloga,4)*Power(fnu,2)*N20*Power(Or,2)*Power(sigma,2) + 78*Power(dloga,3)*Power(eps,2)*Phi0*Power(sigma,2) + 117*Power(dloga,4)*Power(eps,2)*fb*Om*Phi0*Power(sigma,2) + 117*Power(dloga,4)*Power(eps,2)*fc*Om*Phi0*Power(sigma,2) + 156*Power(dloga,4)*Power(eps,2)*fgam*Or*Phi0*Power(sigma,2) + 156*Power(dloga,4)*Power(eps,2)*fnu*Or*Phi0*Power(sigma,2) + 192*Power(dloga,4)*Power(eps,2)*fnu*N00*Or*R*Power(sigma,2) + 1152*Power(dloga,3)*fnu*N20*Or*R*Power(sigma,2) + 384*Power(dloga,4)*Power(eps,2)*fnu*N20*Or*R*Power(sigma,2) + 1728*Power(dloga,4)*fb*fnu*N20*Om*Or*R*Power(sigma,2) + 1728*Power(dloga,4)*fc*fnu*N20*Om*Or*R*Power(sigma,2) + 2304*Power(dloga,4)*fgam*fnu*N20*Power(Or,2)*R*Power(sigma,2) + 2304*Power(dloga,4)*Power(fnu,2)*N20*Power(Or,2)*R*Power(sigma,2) + 96*Power(dloga,3)*Power(eps,2)*Phi0*R*Power(sigma,2) + 144*Power(dloga,4)*Power(eps,2)*fb*Om*Phi0*R*Power(sigma,2) + 144*Power(dloga,4)*Power(eps,2)*fc*Om*Phi0*R*Power(sigma,2) + 192*Power(dloga,4)*Power(eps,2)*fgam*Or*Phi0*R*Power(sigma,2) + 192*Power(dloga,4)*Power(eps,2)*fnu*Or*Phi0*R*Power(sigma,2) + 36*Power(dloga,5)*Power(eps,2)*fnu*N00*Or*Power(sigma,3) + 216*Power(dloga,4)*fnu*N20*Or*Power(sigma,3) + 72*Power(dloga,5)*Power(eps,2)*fnu*N20*Or*Power(sigma,3) + 324*Power(dloga,5)*fb*fnu*N20*Om*Or*Power(sigma,3) + 324*Power(dloga,5)*fc*fnu*N20*Om*Or*Power(sigma,3) + 432*Power(dloga,5)*fgam*fnu*N20*Power(Or,2)*Power(sigma,3) + 432*Power(dloga,5)*Power(fnu,2)*N20*Power(Or,2)*Power(sigma,3) + 18*Power(dloga,4)*Power(eps,2)*Phi0*Power(sigma,3) + 27*Power(dloga,5)*Power(eps,2)*fb*Om*Phi0*Power(sigma,3) + 27*Power(dloga,5)*Power(eps,2)*fc*Om*Phi0*Power(sigma,3) + 36*Power(dloga,5)*Power(eps,2)*fgam*Or*Phi0*Power(sigma,3) + 36*Power(dloga,5)*Power(eps,2)*fnu*Or*Phi0*Power(sigma,3) + 36*Power(dloga,5)*Power(eps,2)*fnu*N00*Or*R*Power(sigma,3) + 216*Power(dloga,4)*fnu*N20*Or*R*Power(sigma,3) + 72*Power(dloga,5)*Power(eps,2)*fnu*N20*Or*R*Power(sigma,3) + 324*Power(dloga,5)*fb*fnu*N20*Om*Or*R*Power(sigma,3) + 324*Power(dloga,5)*fc*fnu*N20*Om*Or*R*Power(sigma,3) + 432*Power(dloga,5)*fgam*fnu*N20*Power(Or,2)*R*Power(sigma,3) + 432*Power(dloga,5)*Power(fnu,2)*N20*Power(Or,2)*R*Power(sigma,3) + 18*Power(dloga,4)*Power(eps,2)*Phi0*R*Power(sigma,3) + 27*Power(dloga,5)*Power(eps,2)*fb*Om*Phi0*R*Power(sigma,3) + 27*Power(dloga,5)*Power(eps,2)*fc*Om*Phi0*R*Power(sigma,3) + 36*Power(dloga,5)*Power(eps,2)*fgam*Or*Phi0*R*Power(sigma,3) + 36*Power(dloga,5)*Power(eps,2)*fnu*Or*Phi0*R*Power(sigma,3) + 3*deltab0*Power(dloga,2)*Power(eps,2)*fb*Om*(R + dloga*sigma + dloga*R*sigma)* (10 + 13*dloga*sigma + 3*Power(dloga,2)*Power(sigma,2)) + 3*delta0*Power(dloga,2)*Power(eps,2)*fc*Om*(R + dloga*sigma + dloga*R*sigma)*(10 + 13*dloga*sigma + 3*Power(dloga,2)*Power(sigma,2)) + 120*Power(dloga,2)*Power(eps,2)*fgam*Or*R*Theta00 + 120*Power(dloga,3)*Power(eps,2)*fgam*Or*sigma*Theta00 + 276*Power(dloga,3)*Power(eps,2)*fgam*Or*R*sigma*Theta00 + 156*Power(dloga,4)*Power(eps,2)*fgam*Or*Power(sigma,2)*Theta00 + 192*Power(dloga,4)*Power(eps,2)*fgam*Or*R*Power(sigma,2)*Theta00 + 36*Power(dloga,5)*Power(eps,2)*fgam*Or*Power(sigma,3)*Theta00 + 36*Power(dloga,5)*Power(eps,2)*fgam*Or*R*Power(sigma,3)*Theta00 - 180*dloga*eps*sigma*Theta10 - 180*Power(dloga,2)*eps*sigma*Theta10 - 60*Power(dloga,2)*Power(eps,3)*sigma*Theta10 - 270*Power(dloga,2)*eps*fb*Om*sigma*Theta10 - 270*Power(dloga,2)*eps*fc*Om*sigma*Theta10 - 360*Power(dloga,2)*eps*fgam*Or*sigma*Theta10 - 360*Power(dloga,2)*eps*fnu*Or*sigma*Theta10 - 234*Power(dloga,2)*eps*Power(sigma,2)*Theta10 - 234*Power(dloga,3)*eps*Power(sigma,2)*Theta10 - 78*Power(dloga,3)*Power(eps,3)*Power(sigma,2)*Theta10 - 351*Power(dloga,3)*eps*fb*Om*Power(sigma,2)*Theta10 - 351*Power(dloga,3)*eps*fc*Om*Power(sigma,2)*Theta10 - 468*Power(dloga,3)*eps*fgam*Or*Power(sigma,2)*Theta10 - 468*Power(dloga,3)*eps*fnu*Or*Power(sigma,2)*Theta10 - 54*Power(dloga,3)*eps*Power(sigma,3)*Theta10 - 54*Power(dloga,4)*eps*Power(sigma,3)*Theta10 - 18*Power(dloga,4)*Power(eps,3)*Power(sigma,3)*Theta10 - 81*Power(dloga,4)*eps*fb*Om*Power(sigma,3)*Theta10 - 81*Power(dloga,4)*eps*fc*Om*Power(sigma,3)*Theta10 - 108*Power(dloga,4)*eps*fgam*Or*Power(sigma,3)*Theta10 - 108*Power(dloga,4)*eps*fnu*Or*Power(sigma,3)*Theta10 + 720*dloga*fgam*Or*R*Theta20 + 240*Power(dloga,2)*Power(eps,2)*fgam*Or*R*Theta20 + 1080*Power(dloga,2)*fb*fgam*Om*Or*R*Theta20 + 1080*Power(dloga,2)*fc*fgam*Om*Or*R*Theta20 + 1440*Power(dloga,2)*Power(fgam,2)*Power(Or,2)*R*Theta20 + 1440*Power(dloga,2)*fgam*fnu*Power(Or,2)*R*Theta20 + 720*Power(dloga,2)*fgam*Or*sigma*Theta20 + 240*Power(dloga,3)*Power(eps,2)*fgam*Or*sigma*Theta20 + 1080*Power(dloga,3)*fb*fgam*Om*Or*sigma*Theta20 + 1080*Power(dloga,3)*fc*fgam*Om*Or*sigma*Theta20 + 1440*Power(dloga,3)*Power(fgam,2)*Power(Or,2)*sigma*Theta20 + 1440*Power(dloga,3)*fgam*fnu*Power(Or,2)*sigma*Theta20 + 1008*Power(dloga,2)*fgam*Or*R*sigma*Theta20 + 336*Power(dloga,3)*Power(eps,2)*fgam*Or*R*sigma*Theta20 + 1512*Power(dloga,3)*fb*fgam*Om*Or*R*sigma*Theta20 + 1512*Power(dloga,3)*fc*fgam*Om*Or*R*sigma*Theta20 + 2016*Power(dloga,3)*Power(fgam,2)*Power(Or,2)*R*sigma*Theta20 + 2016*Power(dloga,3)*fgam*fnu*Power(Or,2)*R*sigma*Theta20 + 288*Power(dloga,3)*fgam*Or*Power(sigma,2)*Theta20 + 96*Power(dloga,4)*Power(eps,2)*fgam*Or*Power(sigma,2)*Theta20 + 432*Power(dloga,4)*fb*fgam*Om*Or*Power(sigma,2)*Theta20 + 432*Power(dloga,4)*fc*fgam*Om*Or*Power(sigma,2)*Theta20 + 576*Power(dloga,4)*Power(fgam,2)*Power(Or,2)*Power(sigma,2)*Theta20 + 576*Power(dloga,4)*fgam*fnu*Power(Or,2)*Power(sigma,2)*Theta20 + 288*Power(dloga,3)*fgam*Or*R*Power(sigma,2)*Theta20 + 96*Power(dloga,4)*Power(eps,2)*fgam*Or*R*Power(sigma,2)*Theta20 + 432*Power(dloga,4)*fb*fgam*Om*Or*R*Power(sigma,2)*Theta20 + 432*Power(dloga,4)*fc*fgam*Om*Or*R*Power(sigma,2)*Theta20 + 576*Power(dloga,4)*Power(fgam,2)*Power(Or,2)*R*Power(sigma,2)*Theta20 + 576*Power(dloga,4)*fgam*fnu*Power(Or,2)*R*Power(sigma,2)*Theta20 + 72*Power(dloga,2)*fgam*Or*R*sigma*ThetaP00 + 24*Power(dloga,3)*Power(eps,2)*fgam*Or*R*sigma*ThetaP00 + 108*Power(dloga,3)*fb*fgam*Om*Or*R*sigma*ThetaP00 + 108*Power(dloga,3)*fc*fgam*Om*Or*R*sigma*ThetaP00 + 144*Power(dloga,3)*Power(fgam,2)*Power(Or,2)*R*sigma*ThetaP00 + 144*Power(dloga,3)*fgam*fnu*Power(Or,2)*R*sigma*ThetaP00 + 72*Power(dloga,3)*fgam*Or*Power(sigma,2)*ThetaP00 + 24*Power(dloga,4)*Power(eps,2)*fgam*Or*Power(sigma,2)*ThetaP00 + 108*Power(dloga,4)*fb*fgam*Om*Or*Power(sigma,2)*ThetaP00 + 108*Power(dloga,4)*fc*fgam*Om*Or*Power(sigma,2)*ThetaP00 + 144*Power(dloga,4)*Power(fgam,2)*Power(Or,2)*Power(sigma,2)*ThetaP00 + 144*Power(dloga,4)*fgam*fnu*Power(Or,2)*Power(sigma,2)*ThetaP00 + 72*Power(dloga,3)*fgam*Or*R*Power(sigma,2)*ThetaP00 + 24*Power(dloga,4)*Power(eps,2)*fgam*Or*R*Power(sigma,2)*ThetaP00 + 108*Power(dloga,4)*fb*fgam*Om*Or*R*Power(sigma,2)*ThetaP00 + 108*Power(dloga,4)*fc*fgam*Om*Or*R*Power(sigma,2)*ThetaP00 + 144*Power(dloga,4)*Power(fgam,2)*Power(Or,2)*R*Power(sigma,2)*ThetaP00 + 144*Power(dloga,4)*fgam*fnu*Power(Or,2)*R*Power(sigma,2)*ThetaP00 + 72*Power(dloga,2)*fgam*Or*R*sigma*ThetaP20 + 24*Power(dloga,3)*Power(eps,2)*fgam*Or*R*sigma*ThetaP20 + 108*Power(dloga,3)*fb*fgam*Om*Or*R*sigma*ThetaP20 + 108*Power(dloga,3)*fc*fgam*Om*Or*R*sigma*ThetaP20 + 144*Power(dloga,3)*Power(fgam,2)*Power(Or,2)*R*sigma*ThetaP20 + 144*Power(dloga,3)*fgam*fnu*Power(Or,2)*R*sigma*ThetaP20 + 72*Power(dloga,3)*fgam*Or*Power(sigma,2)*ThetaP20 + 24*Power(dloga,4)*Power(eps,2)*fgam*Or*Power(sigma,2)*ThetaP20 + 108*Power(dloga,4)*fb*fgam*Om*Or*Power(sigma,2)*ThetaP20 + 108*Power(dloga,4)*fc*fgam*Om*Or*Power(sigma,2)*ThetaP20 + 144*Power(dloga,4)*Power(fgam,2)*Power(Or,2)*Power(sigma,2)*ThetaP20 + 144*Power(dloga,4)*fgam*fnu*Power(Or,2)*Power(sigma,2)*ThetaP20 + 72*Power(dloga,3)*fgam*Or*R*Power(sigma,2)*ThetaP20 + 24*Power(dloga,4)*Power(eps,2)*fgam*Or*R*Power(sigma,2)*ThetaP20 + 108*Power(dloga,4)*fb*fgam*Om*Or*R*Power(sigma,2)*ThetaP20 + 108*Power(dloga,4)*fc*fgam*Om*Or*R*Power(sigma,2)*ThetaP20 + 144*Power(dloga,4)*Power(fgam,2)*Power(Or,2)*R*Power(sigma,2)*ThetaP20 + 144*Power(dloga,4)*fgam*fnu*Power(Or,2)*R*Power(sigma,2)*ThetaP20 - 60*eps*R*ub0 - 60*dloga*eps*R*ub0 - 20*dloga*Power(eps,3)*R*ub0 - 90*dloga*eps*fb*Om*R*ub0 - 90*dloga*eps*fc*Om*R*ub0 - 120*dloga*eps*fgam*Or*R*ub0 - 120*dloga*eps*fnu*Or*R*ub0 - 138*dloga*eps*R*sigma*ub0 - 138*Power(dloga,2)*eps*R*sigma*ub0 - 46*Power(dloga,2)*Power(eps,3)*R*sigma*ub0 - 207*Power(dloga,2)*eps*fb*Om*R*sigma*ub0 - 207*Power(dloga,2)*eps*fc*Om*R*sigma*ub0 - 276*Power(dloga,2)*eps*fgam*Or*R*sigma*ub0 - 276*Power(dloga,2)*eps*fnu*Or*R*sigma*ub0 - 96*Power(dloga,2)*eps*R*Power(sigma,2)*ub0 - 96*Power(dloga,3)*eps*R*Power(sigma,2)*ub0 - 32*Power(dloga,3)*Power(eps,3)*R*Power(sigma,2)*ub0 - 144*Power(dloga,3)*eps*fb*Om*R*Power(sigma,2)*ub0 - 144*Power(dloga,3)*eps*fc*Om*R*Power(sigma,2)*ub0 - 192*Power(dloga,3)*eps*fgam*Or*R*Power(sigma,2)*ub0 - 192*Power(dloga,3)*eps*fnu*Or*R*Power(sigma,2)*ub0 - 18*Power(dloga,3)*eps*R*Power(sigma,3)*ub0 - 18*Power(dloga,4)*eps*R*Power(sigma,3)*ub0 - 6*Power(dloga,4)*Power(eps,3)*R*Power(sigma,3)*ub0 - 27*Power(dloga,4)*eps*fb*Om*R*Power(sigma,3)*ub0 - 27*Power(dloga,4)*eps*fc*Om*R*Power(sigma,3)*ub0 - 36*Power(dloga,4)*eps*fgam*Or*R*Power(sigma,3)*ub0 - 36*Power(dloga,4)*eps*fnu*Or*R*Power(sigma,3)*ub0)/ (eps*(6 + dloga*(6 + 2*Power(eps,2) + 9*fb*Om + 9*fc*Om + 12*fgam*Or + 12*fnu*Or))*(1 + dloga*sigma)*(10 + 3*dloga*sigma)*(R + dloga*sigma + dloga*R*sigma)))),
//							Rule(Theta0,(-120*dloga*Power(eps,2)*fnu*N00*Or + 720*dloga*fnu*N20*Or + 60*dloga*Power(eps,2)*Phi0 + 20*dloga*Power(eps,4)*Phi0 - 156*Power(dloga,2)*Power(eps,2)*fnu*N00*Or*sigma + 936*Power(dloga,2)*fnu*N20*Or*sigma + 78*Power(dloga,2)*Power(eps,2)*Phi0*sigma + 26*Power(dloga,2)*Power(eps,4)*Phi0*sigma - 36*Power(dloga,3)*Power(eps,2)*fnu*N00*Or*Power(sigma,2) + 216*Power(dloga,3)*fnu*N20*Or*Power(sigma,2) + 18*Power(dloga,3)*Power(eps,2)*Phi0*Power(sigma,2) + 6*Power(dloga,3)*Power(eps,4)*Phi0*Power(sigma,2) - 3*deltab0*dloga*Power(eps,2)*fb*Om*(10 + 13*dloga*sigma + 3*Power(dloga,2)*Power(sigma,2)) - 3*delta0*dloga*Power(eps,2)*fc*Om*(10 + 13*dloga*sigma + 3*Power(dloga,2)*Power(sigma,2)) + 60*Power(eps,2)*Theta00 + 60*dloga*Power(eps,2)*Theta00 + 20*dloga*Power(eps,4)*Theta00 + 90*dloga*Power(eps,2)*fb*Om*Theta00 + 90*dloga*Power(eps,2)*fc*Om*Theta00 + 120*dloga*Power(eps,2)*fnu*Or*Theta00 + 78*dloga*Power(eps,2)*sigma*Theta00 + 78*Power(dloga,2)*Power(eps,2)*sigma*Theta00 + 26*Power(dloga,2)*Power(eps,4)*sigma*Theta00 + 117*Power(dloga,2)*Power(eps,2)*fb*Om*sigma*Theta00 + 117*Power(dloga,2)*Power(eps,2)*fc*Om*sigma*Theta00 + 156*Power(dloga,2)*Power(eps,2)*fnu*Or*sigma*Theta00 + 18*Power(dloga,2)*Power(eps,2)*Power(sigma,2)*Theta00 + 18*Power(dloga,3)*Power(eps,2)*Power(sigma,2)*Theta00 + 6*Power(dloga,3)*Power(eps,4)*Power(sigma,2)*Theta00 + 27*Power(dloga,3)*Power(eps,2)*fb*Om*Power(sigma,2)*Theta00 + 27*Power(dloga,3)*Power(eps,2)*fc*Om*Power(sigma,2)*Theta00 + 36*Power(dloga,3)*Power(eps,2)*fnu*Or*Power(sigma,2)*Theta00 + 720*dloga*fgam*Or*Theta20 + 288*Power(dloga,2)*fgam*Or*sigma*Theta20 + 72*Power(dloga,2)*fgam*Or*sigma*ThetaP00 + 72*Power(dloga,2)*fgam*Or*sigma*ThetaP20)/ (Power(eps,2)*(6 + dloga*(6 + 2*Power(eps,2) + 9*fb*Om + 9*fc*Om + 12*fgam*Or + 12*fnu*Or))*(1 + dloga*sigma)*(10 + 3*dloga*sigma))),
//							Rule(Theta1,-(120*Power(dloga,2)*Power(eps,2)*fnu*N00*Or*R + 720*dloga*fnu*N20*Or*R + 240*Power(dloga,2)*Power(eps,2)*fnu*N20*Or*R + 1080*Power(dloga,2)*fb*fnu*N20*Om*Or*R + 1080*Power(dloga,2)*fc*fnu*N20*Om*Or*R + 1440*Power(dloga,2)*fgam*fnu*N20*Power(Or,2)*R + 1440*Power(dloga,2)*Power(fnu,2)*N20*Power(Or,2)*R + 60*dloga*Power(eps,2)*Phi0*R + 90*Power(dloga,2)*Power(eps,2)*fb*Om*Phi0*R + 90*Power(dloga,2)*Power(eps,2)*fc*Om*Phi0*R + 120*Power(dloga,2)*Power(eps,2)*fgam*Or*Phi0*R + 120*Power(dloga,2)*Power(eps,2)*fnu*Or*Phi0*R + 120*Power(dloga,3)*Power(eps,2)*fnu*N00*Or*sigma + 720*Power(dloga,2)*fnu*N20*Or*sigma + 240*Power(dloga,3)*Power(eps,2)*fnu*N20*Or*sigma + 1080*Power(dloga,3)*fb*fnu*N20*Om*Or*sigma + 1080*Power(dloga,3)*fc*fnu*N20*Om*Or*sigma + 1440*Power(dloga,3)*fgam*fnu*N20*Power(Or,2)*sigma + 1440*Power(dloga,3)*Power(fnu,2)*N20*Power(Or,2)*sigma + 60*Power(dloga,2)*Power(eps,2)*Phi0*sigma + 90*Power(dloga,3)*Power(eps,2)*fb*Om*Phi0*sigma + 90*Power(dloga,3)*Power(eps,2)*fc*Om*Phi0*sigma + 120*Power(dloga,3)*Power(eps,2)*fgam*Or*Phi0*sigma + 120*Power(dloga,3)*Power(eps,2)*fnu*Or*Phi0*sigma + 276*Power(dloga,3)*Power(eps,2)*fnu*N00*Or*R*sigma + 1656*Power(dloga,2)*fnu*N20*Or*R*sigma + 552*Power(dloga,3)*Power(eps,2)*fnu*N20*Or*R*sigma + 2484*Power(dloga,3)*fb*fnu*N20*Om*Or*R*sigma + 2484*Power(dloga,3)*fc*fnu*N20*Om*Or*R*sigma + 3312*Power(dloga,3)*fgam*fnu*N20*Power(Or,2)*R*sigma + 3312*Power(dloga,3)*Power(fnu,2)*N20*Power(Or,2)*R*sigma + 138*Power(dloga,2)*Power(eps,2)*Phi0*R*sigma + 207*Power(dloga,3)*Power(eps,2)*fb*Om*Phi0*R*sigma + 207*Power(dloga,3)*Power(eps,2)*fc*Om*Phi0*R*sigma + 276*Power(dloga,3)*Power(eps,2)*fgam*Or*Phi0*R*sigma + 276*Power(dloga,3)*Power(eps,2)*fnu*Or*Phi0*R*sigma + 156*Power(dloga,4)*Power(eps,2)*fnu*N00*Or*Power(sigma,2) + 936*Power(dloga,3)*fnu*N20*Or*Power(sigma,2) + 312*Power(dloga,4)*Power(eps,2)*fnu*N20*Or*Power(sigma,2) + 1404*Power(dloga,4)*fb*fnu*N20*Om*Or*Power(sigma,2) + 1404*Power(dloga,4)*fc*fnu*N20*Om*Or*Power(sigma,2) + 1872*Power(dloga,4)*fgam*fnu*N20*Power(Or,2)*Power(sigma,2) + 1872*Power(dloga,4)*Power(fnu,2)*N20*Power(Or,2)*Power(sigma,2) + 78*Power(dloga,3)*Power(eps,2)*Phi0*Power(sigma,2) + 117*Power(dloga,4)*Power(eps,2)*fb*Om*Phi0*Power(sigma,2) + 117*Power(dloga,4)*Power(eps,2)*fc*Om*Phi0*Power(sigma,2) + 156*Power(dloga,4)*Power(eps,2)*fgam*Or*Phi0*Power(sigma,2) + 156*Power(dloga,4)*Power(eps,2)*fnu*Or*Phi0*Power(sigma,2) + 192*Power(dloga,4)*Power(eps,2)*fnu*N00*Or*R*Power(sigma,2) + 1152*Power(dloga,3)*fnu*N20*Or*R*Power(sigma,2) + 384*Power(dloga,4)*Power(eps,2)*fnu*N20*Or*R*Power(sigma,2) + 1728*Power(dloga,4)*fb*fnu*N20*Om*Or*R*Power(sigma,2) + 1728*Power(dloga,4)*fc*fnu*N20*Om*Or*R*Power(sigma,2) + 2304*Power(dloga,4)*fgam*fnu*N20*Power(Or,2)*R*Power(sigma,2) + 2304*Power(dloga,4)*Power(fnu,2)*N20*Power(Or,2)*R*Power(sigma,2) + 96*Power(dloga,3)*Power(eps,2)*Phi0*R*Power(sigma,2) + 144*Power(dloga,4)*Power(eps,2)*fb*Om*Phi0*R*Power(sigma,2) + 144*Power(dloga,4)*Power(eps,2)*fc*Om*Phi0*R*Power(sigma,2) + 192*Power(dloga,4)*Power(eps,2)*fgam*Or*Phi0*R*Power(sigma,2) + 192*Power(dloga,4)*Power(eps,2)*fnu*Or*Phi0*R*Power(sigma,2) + 36*Power(dloga,5)*Power(eps,2)*fnu*N00*Or*Power(sigma,3) + 216*Power(dloga,4)*fnu*N20*Or*Power(sigma,3) + 72*Power(dloga,5)*Power(eps,2)*fnu*N20*Or*Power(sigma,3) + 324*Power(dloga,5)*fb*fnu*N20*Om*Or*Power(sigma,3) + 324*Power(dloga,5)*fc*fnu*N20*Om*Or*Power(sigma,3) + 432*Power(dloga,5)*fgam*fnu*N20*Power(Or,2)*Power(sigma,3) + 432*Power(dloga,5)*Power(fnu,2)*N20*Power(Or,2)*Power(sigma,3) + 18*Power(dloga,4)*Power(eps,2)*Phi0*Power(sigma,3) + 27*Power(dloga,5)*Power(eps,2)*fb*Om*Phi0*Power(sigma,3) + 27*Power(dloga,5)*Power(eps,2)*fc*Om*Phi0*Power(sigma,3) + 36*Power(dloga,5)*Power(eps,2)*fgam*Or*Phi0*Power(sigma,3) + 36*Power(dloga,5)*Power(eps,2)*fnu*Or*Phi0*Power(sigma,3) + 36*Power(dloga,5)*Power(eps,2)*fnu*N00*Or*R*Power(sigma,3) + 216*Power(dloga,4)*fnu*N20*Or*R*Power(sigma,3) + 72*Power(dloga,5)*Power(eps,2)*fnu*N20*Or*R*Power(sigma,3) + 324*Power(dloga,5)*fb*fnu*N20*Om*Or*R*Power(sigma,3) + 324*Power(dloga,5)*fc*fnu*N20*Om*Or*R*Power(sigma,3) + 432*Power(dloga,5)*fgam*fnu*N20*Power(Or,2)*R*Power(sigma,3) + 432*Power(dloga,5)*Power(fnu,2)*N20*Power(Or,2)*R*Power(sigma,3) + 18*Power(dloga,4)*Power(eps,2)*Phi0*R*Power(sigma,3) + 27*Power(dloga,5)*Power(eps,2)*fb*Om*Phi0*R*Power(sigma,3) + 27*Power(dloga,5)*Power(eps,2)*fc*Om*Phi0*R*Power(sigma,3) + 36*Power(dloga,5)*Power(eps,2)*fgam*Or*Phi0*R*Power(sigma,3) + 36*Power(dloga,5)*Power(eps,2)*fnu*Or*Phi0*R*Power(sigma,3) + 3*deltab0*Power(dloga,2)*Power(eps,2)*fb*Om*(R + dloga*sigma + dloga*R*sigma)* (10 + 13*dloga*sigma + 3*Power(dloga,2)*Power(sigma,2)) + 3*delta0*Power(dloga,2)*Power(eps,2)*fc*Om*(R + dloga*sigma + dloga*R*sigma)*(10 + 13*dloga*sigma + 3*Power(dloga,2)*Power(sigma,2)) + 120*Power(dloga,2)*Power(eps,2)*fgam*Or*R*Theta00 + 120*Power(dloga,3)*Power(eps,2)*fgam*Or*sigma*Theta00 + 276*Power(dloga,3)*Power(eps,2)*fgam*Or*R*sigma*Theta00 + 156*Power(dloga,4)*Power(eps,2)*fgam*Or*Power(sigma,2)*Theta00 + 192*Power(dloga,4)*Power(eps,2)*fgam*Or*R*Power(sigma,2)*Theta00 + 36*Power(dloga,5)*Power(eps,2)*fgam*Or*Power(sigma,3)*Theta00 + 36*Power(dloga,5)*Power(eps,2)*fgam*Or*R*Power(sigma,3)*Theta00 - 180*eps*R*Theta10 - 180*dloga*eps*R*Theta10 - 60*dloga*Power(eps,3)*R*Theta10 - 270*dloga*eps*fb*Om*R*Theta10 - 270*dloga*eps*fc*Om*R*Theta10 - 360*dloga*eps*fgam*Or*R*Theta10 - 360*dloga*eps*fnu*Or*R*Theta10 - 180*dloga*eps*sigma*Theta10 - 180*Power(dloga,2)*eps*sigma*Theta10 - 60*Power(dloga,2)*Power(eps,3)*sigma*Theta10 - 270*Power(dloga,2)*eps*fb*Om*sigma*Theta10 - 270*Power(dloga,2)*eps*fc*Om*sigma*Theta10 - 360*Power(dloga,2)*eps*fgam*Or*sigma*Theta10 - 360*Power(dloga,2)*eps*fnu*Or*sigma*Theta10 - 234*dloga*eps*R*sigma*Theta10 - 234*Power(dloga,2)*eps*R*sigma*Theta10 - 78*Power(dloga,2)*Power(eps,3)*R*sigma*Theta10 - 351*Power(dloga,2)*eps*fb*Om*R*sigma*Theta10 - 351*Power(dloga,2)*eps*fc*Om*R*sigma*Theta10 - 468*Power(dloga,2)*eps*fgam*Or*R*sigma*Theta10 - 468*Power(dloga,2)*eps*fnu*Or*R*sigma*Theta10 - 234*Power(dloga,2)*eps*Power(sigma,2)*Theta10 - 234*Power(dloga,3)*eps*Power(sigma,2)*Theta10 - 78*Power(dloga,3)*Power(eps,3)*Power(sigma,2)*Theta10 - 351*Power(dloga,3)*eps*fb*Om*Power(sigma,2)*Theta10 - 351*Power(dloga,3)*eps*fc*Om*Power(sigma,2)*Theta10 - 468*Power(dloga,3)*eps*fgam*Or*Power(sigma,2)*Theta10 - 468*Power(dloga,3)*eps*fnu*Or*Power(sigma,2)*Theta10 - 54*Power(dloga,2)*eps*R*Power(sigma,2)*Theta10 - 54*Power(dloga,3)*eps*R*Power(sigma,2)*Theta10 - 18*Power(dloga,3)*Power(eps,3)*R*Power(sigma,2)*Theta10 - 81*Power(dloga,3)*eps*fb*Om*R*Power(sigma,2)*Theta10 - 81*Power(dloga,3)*eps*fc*Om*R*Power(sigma,2)*Theta10 - 108*Power(dloga,3)*eps*fgam*Or*R*Power(sigma,2)*Theta10 - 108*Power(dloga,3)*eps*fnu*Or*R*Power(sigma,2)*Theta10 - 54*Power(dloga,3)*eps*Power(sigma,3)*Theta10 - 54*Power(dloga,4)*eps*Power(sigma,3)*Theta10 - 18*Power(dloga,4)*Power(eps,3)*Power(sigma,3)*Theta10 - 81*Power(dloga,4)*eps*fb*Om*Power(sigma,3)*Theta10 - 81*Power(dloga,4)*eps*fc*Om*Power(sigma,3)*Theta10 - 108*Power(dloga,4)*eps*fgam*Or*Power(sigma,3)*Theta10 - 108*Power(dloga,4)*eps*fnu*Or*Power(sigma,3)*Theta10 + 720*dloga*fgam*Or*R*Theta20 + 240*Power(dloga,2)*Power(eps,2)*fgam*Or*R*Theta20 + 1080*Power(dloga,2)*fb*fgam*Om*Or*R*Theta20 + 1080*Power(dloga,2)*fc*fgam*Om*Or*R*Theta20 + 1440*Power(dloga,2)*Power(fgam,2)*Power(Or,2)*R*Theta20 + 1440*Power(dloga,2)*fgam*fnu*Power(Or,2)*R*Theta20 + 720*Power(dloga,2)*fgam*Or*sigma*Theta20 + 240*Power(dloga,3)*Power(eps,2)*fgam*Or*sigma*Theta20 + 1080*Power(dloga,3)*fb*fgam*Om*Or*sigma*Theta20 + 1080*Power(dloga,3)*fc*fgam*Om*Or*sigma*Theta20 + 1440*Power(dloga,3)*Power(fgam,2)*Power(Or,2)*sigma*Theta20 + 1440*Power(dloga,3)*fgam*fnu*Power(Or,2)*sigma*Theta20 + 1008*Power(dloga,2)*fgam*Or*R*sigma*Theta20 + 336*Power(dloga,3)*Power(eps,2)*fgam*Or*R*sigma*Theta20 + 1512*Power(dloga,3)*fb*fgam*Om*Or*R*sigma*Theta20 + 1512*Power(dloga,3)*fc*fgam*Om*Or*R*sigma*Theta20 + 2016*Power(dloga,3)*Power(fgam,2)*Power(Or,2)*R*sigma*Theta20 + 2016*Power(dloga,3)*fgam*fnu*Power(Or,2)*R*sigma*Theta20 + 288*Power(dloga,3)*fgam*Or*Power(sigma,2)*Theta20 + 96*Power(dloga,4)*Power(eps,2)*fgam*Or*Power(sigma,2)*Theta20 + 432*Power(dloga,4)*fb*fgam*Om*Or*Power(sigma,2)*Theta20 + 432*Power(dloga,4)*fc*fgam*Om*Or*Power(sigma,2)*Theta20 + 576*Power(dloga,4)*Power(fgam,2)*Power(Or,2)*Power(sigma,2)*Theta20 + 576*Power(dloga,4)*fgam*fnu*Power(Or,2)*Power(sigma,2)*Theta20 + 288*Power(dloga,3)*fgam*Or*R*Power(sigma,2)*Theta20 + 96*Power(dloga,4)*Power(eps,2)*fgam*Or*R*Power(sigma,2)*Theta20 + 432*Power(dloga,4)*fb*fgam*Om*Or*R*Power(sigma,2)*Theta20 + 432*Power(dloga,4)*fc*fgam*Om*Or*R*Power(sigma,2)*Theta20 + 576*Power(dloga,4)*Power(fgam,2)*Power(Or,2)*R*Power(sigma,2)*Theta20 + 576*Power(dloga,4)*fgam*fnu*Power(Or,2)*R*Power(sigma,2)*Theta20 + 72*Power(dloga,2)*fgam*Or*R*sigma*ThetaP00 + 24*Power(dloga,3)*Power(eps,2)*fgam*Or*R*sigma*ThetaP00 + 108*Power(dloga,3)*fb*fgam*Om*Or*R*sigma*ThetaP00 + 108*Power(dloga,3)*fc*fgam*Om*Or*R*sigma*ThetaP00 + 144*Power(dloga,3)*Power(fgam,2)*Power(Or,2)*R*sigma*ThetaP00 + 144*Power(dloga,3)*fgam*fnu*Power(Or,2)*R*sigma*ThetaP00 + 72*Power(dloga,3)*fgam*Or*Power(sigma,2)*ThetaP00 + 24*Power(dloga,4)*Power(eps,2)*fgam*Or*Power(sigma,2)*ThetaP00 + 108*Power(dloga,4)*fb*fgam*Om*Or*Power(sigma,2)*ThetaP00 + 108*Power(dloga,4)*fc*fgam*Om*Or*Power(sigma,2)*ThetaP00 + 144*Power(dloga,4)*Power(fgam,2)*Power(Or,2)*Power(sigma,2)*ThetaP00 + 144*Power(dloga,4)*fgam*fnu*Power(Or,2)*Power(sigma,2)*ThetaP00 + 72*Power(dloga,3)*fgam*Or*R*Power(sigma,2)*ThetaP00 + 24*Power(dloga,4)*Power(eps,2)*fgam*Or*R*Power(sigma,2)*ThetaP00 + 108*Power(dloga,4)*fb*fgam*Om*Or*R*Power(sigma,2)*ThetaP00 + 108*Power(dloga,4)*fc*fgam*Om*Or*R*Power(sigma,2)*ThetaP00 + 144*Power(dloga,4)*Power(fgam,2)*Power(Or,2)*R*Power(sigma,2)*ThetaP00 + 144*Power(dloga,4)*fgam*fnu*Power(Or,2)*R*Power(sigma,2)*ThetaP00 + 72*Power(dloga,2)*fgam*Or*R*sigma*ThetaP20 + 24*Power(dloga,3)*Power(eps,2)*fgam*Or*R*sigma*ThetaP20 + 108*Power(dloga,3)*fb*fgam*Om*Or*R*sigma*ThetaP20 + 108*Power(dloga,3)*fc*fgam*Om*Or*R*sigma*ThetaP20 + 144*Power(dloga,3)*Power(fgam,2)*Power(Or,2)*R*sigma*ThetaP20 + 144*Power(dloga,3)*fgam*fnu*Power(Or,2)*R*sigma*ThetaP20 + 72*Power(dloga,3)*fgam*Or*Power(sigma,2)*ThetaP20 + 24*Power(dloga,4)*Power(eps,2)*fgam*Or*Power(sigma,2)*ThetaP20 + 108*Power(dloga,4)*fb*fgam*Om*Or*Power(sigma,2)*ThetaP20 + 108*Power(dloga,4)*fc*fgam*Om*Or*Power(sigma,2)*ThetaP20 + 144*Power(dloga,4)*Power(fgam,2)*Power(Or,2)*Power(sigma,2)*ThetaP20 + 144*Power(dloga,4)*fgam*fnu*Power(Or,2)*Power(sigma,2)*ThetaP20 + 72*Power(dloga,3)*fgam*Or*R*Power(sigma,2)*ThetaP20 + 24*Power(dloga,4)*Power(eps,2)*fgam*Or*R*Power(sigma,2)*ThetaP20 + 108*Power(dloga,4)*fb*fgam*Om*Or*R*Power(sigma,2)*ThetaP20 + 108*Power(dloga,4)*fc*fgam*Om*Or*R*Power(sigma,2)*ThetaP20 + 144*Power(dloga,4)*Power(fgam,2)*Power(Or,2)*R*Power(sigma,2)*ThetaP20 + 144*Power(dloga,4)*fgam*fnu*Power(Or,2)*R*Power(sigma,2)*ThetaP20 - 60*dloga*eps*R*sigma*ub0 - 60*Power(dloga,2)*eps*R*sigma*ub0 - 20*Power(dloga,2)*Power(eps,3)*R*sigma*ub0 - 90*Power(dloga,2)*eps*fb*Om*R*sigma*ub0 - 90*Power(dloga,2)*eps*fc*Om*R*sigma*ub0 - 120*Power(dloga,2)*eps*fgam*Or*R*sigma*ub0 - 120*Power(dloga,2)*eps*fnu*Or*R*sigma*ub0 - 78*Power(dloga,2)*eps*R*Power(sigma,2)*ub0 - 78*Power(dloga,3)*eps*R*Power(sigma,2)*ub0 - 26*Power(dloga,3)*Power(eps,3)*R*Power(sigma,2)*ub0 - 117*Power(dloga,3)*eps*fb*Om*R*Power(sigma,2)*ub0 - 117*Power(dloga,3)*eps*fc*Om*R*Power(sigma,2)*ub0 - 156*Power(dloga,3)*eps*fgam*Or*R*Power(sigma,2)*ub0 - 156*Power(dloga,3)*eps*fnu*Or*R*Power(sigma,2)*ub0 - 18*Power(dloga,3)*eps*R*Power(sigma,3)*ub0 - 18*Power(dloga,4)*eps*R*Power(sigma,3)*ub0 - 6*Power(dloga,4)*Power(eps,3)*R*Power(sigma,3)*ub0 - 27*Power(dloga,4)*eps*fb*Om*R*Power(sigma,3)*ub0 - 27*Power(dloga,4)*eps*fc*Om*R*Power(sigma,3)*ub0 - 36*Power(dloga,4)*eps*fgam*Or*R*Power(sigma,3)*ub0 - 36*Power(dloga,4)*eps*fnu*Or*R*Power(sigma,3)*ub0)/ (3.*eps*(6 + dloga*(6 + 2*Power(eps,2) + 9*fb*Om + 9*fc*Om + 12*fgam*Or + 12*fnu*Or))*(1 + dloga*sigma)*(10 + 3*dloga*sigma)*(R + dloga*sigma + dloga*R*sigma))),
//							Rule(Theta2,(2*(5 + 2*dloga*sigma)*Theta20 + dloga*sigma*(ThetaP00 + ThetaP20))/((1 + dloga*sigma)*(10 + 3*dloga*sigma))),
//							Rule(ThetaP0,(10*ThetaP00 + dloga*sigma*(5*Theta20 + 8*ThetaP00 + 5*ThetaP20))/((1 + dloga*sigma)*(10 + 3*dloga*sigma))),
//							Rule(ThetaP1,ThetaP10/(1 + dloga*sigma)),
//							Rule(ThetaP2,(10*ThetaP20 + dloga*sigma*(Theta20 + ThetaP00 + 4*ThetaP20))/((1 + dloga*sigma)*(10 + 3*dloga*sigma))),
//							Rule(N0,(2*dloga*Power(eps,4)*(N00 + Phi0)*(10 + 13*dloga*sigma + 3*Power(dloga,2)*Power(sigma,2)) + 3*Power(eps,2)*(10 + 13*dloga*sigma + 3*Power(dloga,2)*Power(sigma,2))*(N00*(2 + dloga*(2 + 3*fb*Om + 3*fc*Om + 4*fgam*Or)) - dloga*(deltab0*fb*Om + delta0*fc*Om - 2*Phi0 + 4*fgam*Or*Theta00)) + 72*dloga*Or*(fnu*N20*(10 + 13*dloga*sigma + 3*Power(dloga,2)*Power(sigma,2)) + 2*fgam*(5 + 2*dloga*sigma)*Theta20 + dloga*fgam*sigma*(ThetaP00 + ThetaP20)))/ (Power(eps,2)*(6 + dloga*(6 + 2*Power(eps,2) + 9*fb*Om + 9*fc*Om + 12*fgam*Or + 12*fnu*Or))*(1 + dloga*sigma)*(10 + 3*dloga*sigma))),
//							Rule(N1,(2*dloga*Power(eps,3)*N10*(10 + 13*dloga*sigma + 3*Power(dloga,2)*Power(sigma,2)) + 3*eps*N10*(2 + dloga*(2 + 3*fb*Om + 3*fc*Om + 4*fgam*Or + 4*fnu*Or))*(10 + 13*dloga*sigma + 3*Power(dloga,2)*Power(sigma,2)) - dloga*Power(eps,2)*(40*dloga*fnu*N00*Or + 80*dloga*fnu*N20*Or + 20*Phi0 + 30*dloga*fb*Om*Phi0 + 30*dloga*fc*Om*Phi0 + 40*dloga*fgam*Or*Phi0 + 40*dloga*fnu*Or*Phi0 + 52*Power(dloga,2)*fnu*N00*Or*sigma + 104*Power(dloga,2)*fnu*N20*Or*sigma + 26*dloga*Phi0*sigma + 39*Power(dloga,2)*fb*Om*Phi0*sigma + 39*Power(dloga,2)*fc*Om*Phi0*sigma + 52*Power(dloga,2)*fgam*Or*Phi0*sigma + 52*Power(dloga,2)*fnu*Or*Phi0*sigma + 12*Power(dloga,3)*fnu*N00*Or*Power(sigma,2) + 24*Power(dloga,3)*fnu*N20*Or*Power(sigma,2) + 6*Power(dloga,2)*Phi0*Power(sigma,2) + 9*Power(dloga,3)*fb*Om*Phi0*Power(sigma,2) + 9*Power(dloga,3)*fc*Om*Phi0*Power(sigma,2) + 12*Power(dloga,3)*fgam*Or*Phi0*Power(sigma,2) + 12*Power(dloga,3)*fnu*Or*Phi0*Power(sigma,2) + deltab0*dloga*fb*Om*(10 + 13*dloga*sigma + 3*Power(dloga,2)*Power(sigma,2)) + delta0*dloga*fc*Om*(10 + 13*dloga*sigma + 3*Power(dloga,2)*Power(sigma,2)) + 40*dloga*fgam*Or*Theta00 + 52*Power(dloga,2)*fgam*Or*sigma*Theta00 + 12*Power(dloga,3)*fgam*Or*Power(sigma,2)*Theta00 + 80*dloga*fgam*Or*Theta20 + 32*Power(dloga,2)*fgam*Or*sigma*Theta20 + 8*Power(dloga,2)*fgam*Or*sigma*ThetaP00 + 8*Power(dloga,2)*fgam*Or*sigma*ThetaP20) - 12*dloga*Or*(2 + dloga*(3*fb*Om + 3*fc*Om + 4*(fgam + fnu)*Or))* (fnu*N20*(10 + 13*dloga*sigma + 3*Power(dloga,2)*Power(sigma,2)) + 2*fgam*(5 + 2*dloga*sigma)*Theta20 + dloga*fgam*sigma*(ThetaP00 + ThetaP20)))/ (eps*(6 + dloga*(6 + 2*Power(eps,2) + 9*fb*Om + 9*fc*Om + 12*fgam*Or + 12*fnu*Or))*(1 + dloga*sigma)*(10 + 3*dloga*sigma))),
//							Rule(N2,N20));
//					dudt[Phii] = (Phi - Phi0) / dloga;
//					dudt[deltai] = (delta - delta0) / dloga;
//					dudt[ui] = (u - u0) / dloga;
//					dudt[deltabi] = (deltab - deltab) / dloga;
//					dudt[ubi] = (ub - ub0) / dloga;
//					dudt[Thetai + 0] = (Theta0 - Theta00) / dloga;
//					dudt[Thetai + 1] = (Theta1 - Theta10) / dloga;
//					dudt[Thetai + 2] = (Theta2 - Theta20) / dloga;
//					dudt[ThetaPi + 0] = (ThetaP0 - ThetaP00) / dloga;
//					dudt[ThetaPi + 1] = (ThetaP1 - ThetaP10) / dloga;
//					dudt[ThetaPi + 2] = (ThetaP2 - ThetaP20) / dloga;
//					dudt[Ni + 0] = (N0 - N00) / dloga;
//					dudt[Ni + 1] = (N1 - N10) / dloga;
//					dudt[Ni + 2] = (N2 - N20) / dloga;
//					for (int l = 3; l < LMAX - 1; l++) {
//						dudt[Thetai + l] = -sigma / (1 + sigma * dloga) * U[Thetai + l];
//						dudt[ThetaPi + l] = -sigma / (1 + sigma * dloga) * U[ThetaPi + l];
//						dudt[Ni + l] = 0.0;
//					}
//					dudt[Thetai + LMAX - 1] = -sigma / (1 + sigma * dloga + LMAX / etaaH * dloga) * U[Thetai + LMAX - 1];
//					dudt[ThetaPi + LMAX - 1] = -sigma / (1 + sigma * dloga + LMAX / etaaH * dloga) * U[ThetaPi + LMAX - 1];
//					dudt[Ni + LMAX - 1] = -1 / (1 + LMAX / etaaH * dloga) * U[Ni + LMAX - 1];
//					return dudt;
//				};
//
//		constexpr int Nstage = 3;
//		constexpr double gamma = (3.0 + std::sqrt(3)) / 6.0;
//		std::array<std::array<double, Nstage>, Nstage> a_imp = { { { 0, 0, 0 }, { 0, gamma, 0 }, { 0, 1 - 2 * gamma, gamma } } };
//		std::array<std::array<double, Nstage>, Nstage> a_exp = { { { 0, 0, 0 }, { gamma, 0, 0 }, { gamma - 1.0, 2 * (1 - gamma), 0 } } };
//		std::array<double, Nstage> c_imp = { 0, gamma, 1 - gamma };
//		std::array<double, Nstage> c_exp = { 0, gamma, 1 - gamma };
//		std::array<double, Nstage> w_imp = { 0, 0.5, 0.5 };
//		std::array<double, Nstage> w_exp = { 0, 0.5, 0.5 };
//		std::array<state, Nstage> du_exp;
//		std::array<state, Nstage> du_imp;
//		state U0 = U;
//		for (int s = 0; s < Nstage; s++) {
//			state U1 = U0;
//			for (int m = 0; m < s; m++) {
//				for (int n = 0; n < NF; n++) {
//					U1[n] += a_exp[s][m] * du_exp[m][n] * dloga;
//					U1[n] += a_imp[s][m] * du_imp[m][n] * dloga;
//				}
//			}
//			if (a_imp[s][s] == 0.0) {
//				for (int n = 0; n < NF; n++) {
//					du_imp[s][n] = 0.0;
//				}
//			} else {
//				compute_parameters(loga + c_imp[s] * dloga);
//				du_imp[s] = dudt_imp(U1, dloga * a_imp[s][s]);
//			}
//			for (int n = 0; n < NF; n++) {
//				U1[n] += a_imp[s][s] * du_imp[s][n] * dloga;
//			}
//			compute_parameters(loga + c_exp[s] * dloga);
//			du_exp[s] = dudt_exp(U1);
//
//		}
//		U = U0;
//		for (int s = 0; s < Nstage; s++) {
//			for (int n = 0; n < NF; n++) {
//				U[n] += w_exp[s] * du_exp[s][n] * dloga;
//				U[n] += w_imp[s] * du_imp[s][n] * dloga;
//			}
//		}
//
//		eta += dloga / (a * H);
//		loga += dloga;
//	}
//
//}
//
//void initial_conditions(state &u, double k, double a) {
//	const auto eps = k / (a * Hubble(a));
//	const auto Psii = 1.0;
//	u[Thetai + 0] = -0.5 * Psii;
//	u[Thetai + 1] = eps / 6.0 * Psii;
//	u[Ni + 0] = u[Thetai + 0];
//	u[Ni + 1] = u[Thetai + 1];
//	u[Ni + 2] = eps * eps / 30.0 * Psii;
//	u[deltai] = 3.0 * u[Thetai + 0];
//	u[ui] = 3.0 * u[Thetai + 1];
//	u[deltabi] = 3.0 * u[Thetai + 0];
//	u[ubi] = 3.0 * u[Thetai + 1];
//	u[Phii] = -12.0 * omega_nu / (a * omega_m + omega_r) * (u[Ni + 2]) / (eps * eps) - Psii;
//	for (int l = 0; l < LMAX; l++) {
//		u[ThetaPi + l] = 0.0;
//	}
//	for (int l = 2; l < LMAX; l++) {
//		u[Thetai + l] = 0.0;
//	}
//	for (int l = 3; l < LMAX; l++) {
//		u[Ni + l] = 0.0;
//	}
//}

#define LMAX 50

#define hdoti 0
#define deltaci 1
#define deltabi 2
#define thetabi 3
#define FLi 4
#define GLi (4+LMAX)
#define NLi (4+2*LMAX)

#define deltagami (FLi+0)
#define thetagami (FLi+1)
#define F2i (FLi+2)
#define deltanui (NLi+0)
#define thetanui (NLi+1)
#define N2i (NLi+2)
#define G0i (GLi+0)
#define G1i (GLi+1)
#define G2i (GLi+2)

#define NFIELD (4+(3*LMAX))

using cos_state = std::array<double,NFIELD>;

double advance2(cos_state &U, double k, double amin, double amax, std::function<double(double)> &csfunc, std::function<double(double)> &thomsonfunc) {
	cos_state U0;
	auto loga = std::log(amin);
	const auto logamax = std::log(amax);
	const double eps = k / (amin * Hubble(amin));
	const double C = 1.0 * std::pow(eps, -1.5);
	double a = std::exp(loga);
	double hubble = Hubble(a);
	double Or = omega_r / (omega_r + a * omega_m + (a * a * a * a) * (1.0 - omega_m - omega_r));
	double Om = omega_m / (omega_r / a + omega_m + (a * a * a) * (1.0 - omega_m - omega_r));
	double Ogam = omega_gam * Or / omega_r;
	double Onu = omega_nu * Or / omega_r;
	double Ob = omega_b * Om / omega_m;
	double Oc = omega_c * Om / omega_m;
	double tau = 1.0 / (amin * hubble);
	double Rnu = Onu / Or;
	double eta = 2.0 * C - C * (5 + 4 * Rnu) / (6 * (15 + 4 * Rnu)) * eps * eps;
	U[deltanui] = U[deltagami] = -2.0 / 3.0 * C * eps * eps;
	U[deltabi] = U[deltaci] = 3.0 / 4.0 * U[deltagami];
	U[thetabi] = U[thetagami] = -C / 18.0 * eps * eps * eps;
	U[thetanui] = (23 + 4 * Rnu) / (15 + 4 * Rnu) * U[thetagami];
	U[N2i] = 0.5 * (4.0 * C) / (3.0 * (15 + 4 * Rnu)) * eps * eps;
	U[hdoti] = 2.0 * C * eps * eps;
	U[G0i] = U[G1i] = U[G2i] = U[F2i] = 0.0;
	for (int l = 3; l < LMAX; l++) {
		U[FLi + l] = U[NLi + l] = U[GLi + l] = 0.0;
	}
	eta = (0.5 * U[hdoti] - (1.5 * (Oc * U[deltaci] + Ob * U[deltabi]) + 1.5 * (Ogam * U[deltagami] + Onu * U[deltanui]))) / (eps * eps);
	while (loga < logamax) {
		double a = std::exp(loga);
		double hubble = Hubble(a);
		double eps = k / (a * hubble);
		Or = omega_r / (omega_r + a * omega_m + (a * a * a * a) * (1.0 - omega_m - omega_r));
		Om = omega_m / (omega_r / a + omega_m + (a * a * a) * (1.0 - omega_m - omega_r));
		Ogam = omega_gam * Or / omega_r;
		Onu = omega_nu * Or / omega_r;
		Ob = omega_b * Om / omega_m;
		Oc = omega_c * Om / omega_m;
		double cs2 = std::pow(csfunc(std::log(a)), 2);
		double lambda_max = 0.0;
		lambda_max = std::max(lambda_max, eps * std::max(1.0, std::sqrt(cs2)));
		lambda_max = std::max(lambda_max, std::sqrt(3 * eps * eps + 8 * Or) / sqrt(5));
		double dloga = std::min(std::min(1.0e-2, 1.0 / lambda_max), logamax - loga);
		U0 = U;
		cos_state dudt;
		double beta[3] = { 1, 0.25, (2.0 / 3.0) };
		double tm[3] = { 0, 1, 0.5 };
		double loga0 = loga;
		double tau0 = tau;
		double eta0 = eta;
//		printf("%10.2e ", a);
//		for (int f = 0; f < NFIELD; f++) {
//			printf("%10.2e ", U[f]);
//		}
//		printf("\n");
		double test1 = eps * eps * eta;
		double test2 = -0.5 * U[hdoti];
		double test3 = (1.5 * (Oc * U[deltaci] + Ob * U[deltabi]) + 1.5 * (Ogam * U[deltagami] + Onu * U[deltanui]));
//		printf("%e %e\n", tau,(test1 + test2 + test3)/test1);
		for (int i = 0; i < 3; i++) {
			loga = loga0 + tm[i] * dloga;
			a = std::exp(loga);
			hubble = Hubble(a);
			eps = k / (a * hubble);
			Or = omega_r / (omega_r + a * omega_m + (a * a * a * a) * (1.0 - omega_m - omega_r));
			Om = omega_m / (omega_r / a + omega_m + (a * a * a) * (1.0 - omega_m - omega_r));
			Ogam = omega_gam * Or / omega_r;
			Onu = omega_nu * Or / omega_r;
			Ob = omega_b * Om / omega_m;
			Oc = omega_c * Om / omega_m;
			cs2 = std::pow(csfunc(std::log(a)), 2);
			double dtau = 1.0 / (a * hubble);
			double etadot = (1.5 * ((Ob * U[thetabi]) + (4.0 / 3.0) * (Ogam * U[thetagami] + Onu * U[thetanui])) / eps);
			double factor = ((a * omega_m) + 4 * a * a * a * a * (1 - omega_m - omega_r))
					/ (2 * a * omega_m + 2 * omega_r + 2 * a * a * a * a * (1 - omega_m - omega_r));
			dudt[hdoti] = (-factor * U[hdoti] - (3.0 * (Oc * U[deltaci] + Ob * U[deltabi]) + 6.0 * (Ogam * U[deltagami] + Onu * U[deltanui])));

			dudt[deltaci] = -0.5 * U[hdoti];

			dudt[deltabi] = -eps * U[thetabi] - 0.5 * U[hdoti];

			dudt[deltagami] = -4.0 / 3.0 * eps * U[thetagami] - (2.0 / 3.0) * U[hdoti];

			dudt[deltanui] = -4.0 / 3.0 * eps * U[thetanui] - (2.0 / 3.0) * U[hdoti];

			dudt[thetabi] = -U[thetabi] + cs2 * eps * U[deltabi];

			dudt[thetagami] = eps * (0.25 * U[deltagami] - 0.5 * U[F2i]);

			dudt[thetanui] = eps * (0.25 * U[deltanui] - 0.5 * U[N2i]);

			dudt[F2i] = (8.0 / 15.0) * eps * U[thetagami] + (4.0 / 15.0) * U[hdoti] + (8.0 / 5.0) * etadot - (3.0 / 5.0) * eps * U[FLi + 3];

			dudt[N2i] = (8.0 / 15.0) * eps * U[thetanui] + (4.0 / 15.0) * U[hdoti] + (8.0 / 5.0) * etadot - (3.0 / 5.0) * eps * U[NLi + 3];

			dudt[GLi + 0] = -eps * U[GLi + 1];
			dudt[GLi + 1] = eps / (3) * (U[GLi + 0] - 2 * U[GLi + 2]);
			dudt[GLi + 2] = eps / (5) * (2 * U[GLi + 1] - 3 * U[GLi + 3]);
			for (int l = 3; l < LMAX - 1; l++) {
				dudt[FLi + l] = eps / (2 * l + 1) * (l * U[FLi - 1 + l] - (l + 1) * U[FLi + 1 + l]);
				dudt[NLi + l] = eps / (2 * l + 1) * (l * U[NLi - 1 + l] - (l + 1) * U[NLi + 1 + l]);
				dudt[GLi + l] = eps / (2 * l + 1) * (l * U[GLi - 1 + l] - (l + 1) * U[GLi + 1 + l]);
			}
			dudt[FLi + LMAX - 1] = (eps * U[FLi + LMAX - 2] - LMAX / (tau * a * hubble) * U[FLi + LMAX - 1]) / (2 * LMAX - 1);
			dudt[NLi + LMAX - 1] = (eps * U[NLi + LMAX - 2] - LMAX / (tau * a * hubble) * U[NLi + LMAX - 1]) / (2 * LMAX - 1);
			dudt[GLi + LMAX - 1] = (eps * U[GLi + LMAX - 2] - LMAX / (tau * a * hubble) * U[GLi + LMAX - 1]) / (2 * LMAX - 1);
			for (int f = 0; f < NFIELD; f++) {
				U[f] = (1 - beta[i]) * U0[f] + beta[i] * (U[f] + dudt[f] * dloga);
			}
			tau = (1 - beta[i]) * tau0 + beta[i] * (tau + dtau * dloga);
			eta = (1 - beta[i]) * eta0 + beta[i] * (eta + etadot * dloga);
		}
		loga = loga0 + dloga;
		a = std::exp(loga);
		double &thetab = U[thetabi];
		double &thetagam = U[thetagami];
		double &F2 = U[F2i];
		double &G0 = U[G0i];
		double &G1 = U[G1i];
		double &G2 = U[G2i];
		double thetab0 = thetab;
		double thetagam0 = thetagam;
		double F20 = F2;
		double G00 = G0;
		double G10 = G1;
		double G20 = G2;
		double sigma = thomsonfunc(std::log(a));
#define Rule(a, b) a = b
//		Rule(thetab,
//				-((-3 * Ob * thetab0 - 3 * dloga * Ob * sigma * thetab0 - 4 * dloga * Ogam * sigma * thetagam0)
//						/ (3 * Ob + 3 * dloga * Ob * sigma + 4 * dloga * Ogam * sigma)));
//		Rule(thetagam,
//				-((-3 * dloga * Ob * sigma * thetab0 - 3 * Ob * thetagam0 - 4 * dloga * Ogam * sigma * thetagam0)
//						/ (3 * Ob + 3 * dloga * Ob * sigma + 4 * dloga * Ogam * sigma)));
//		Rule(F2, -((-10 * F20 - 4 * dloga * F20 * sigma - dloga * G00 * sigma - dloga * G20 * sigma) / ((1 + dloga * sigma) * (10 + 3 * dloga * sigma))));
//		Rule(G0,
//				-((-10 * G00 - 5 * dloga * F20 * sigma - 8 * dloga * G00 * sigma - 5 * dloga * G20 * sigma) / ((1 + dloga * sigma) * (10 + 3 * dloga * sigma))));
//		Rule(G1, G10 / (1 + dloga * sigma));
//		Rule(G2, -((-10 * G20 - dloga * F20 * sigma - dloga * G00 * sigma - 4 * dloga * G20 * sigma) / ((1 + dloga * sigma) * (10 + 3 * dloga * sigma))));
//		for (int l = 3; l < LMAX; l++) {
//			U[GLi + l] /= (1 + dloga * sigma);
//			U[FLi + l] /= (1 + dloga * sigma);
//		}
	}
	return 0;
}
//
//double advance3(double &hdot, double &delta, double &Theta0, double &Theta1, double k, double amin, double amax) {
//	auto loga = std::log(amin);
//	const auto logamax = std::log(amax);
//	const double eps = k / (amin * Hubble(amin));
//	const double C = 1/(eps*eps);
//	double h = C * eps * eps;
//	Theta0 = -1.0 / 6.0 * C * eps * eps;
//	delta = 0;
//	Theta1 = 0;
//	double eta = 2.0 * C;
//	const double a = std::exp(loga);
//	double Or = omega_r / (omega_r + a * omega_m + (a * a * a * a) * (1.0 - omega_m - omega_r));
//	double Om = omega_m / (omega_r / a + omega_m + (a * a * a) * (1.0 - omega_m - omega_r));
//	hdot = 2.0 * eps * eps * eta + 3.0 * (Om * delta + 4.0 * Or * Theta0);
//	double Theta00 = Theta0;
//	double delta0 = delta;
//	double Theta10 = Theta1;
//	double hdot0 = hdot;
//	while (loga < logamax) {
//		double a = std::exp(loga);
//		double hubble = Hubble(a);
//		double eps = k / (a * hubble);
//		double Or = omega_r / (omega_r + a * omega_m + (a * a * a * a) * (1.0 - omega_m - omega_r));
//		double Om = omega_m / (omega_r / a + omega_m + (a * a * a) * (1.0 - omega_m - omega_r));
//		double lambda_max = 0.0;
//		lambda_max = std::max(lambda_max, eps + 1.0 / 6.0 * (6.0 * Om + 32.0 * Or));
//		lambda_max = std::max(lambda_max, 1.0 / 2.0 * (6.0 * Om + 32.0 * Or));
//		lambda_max = std::max(lambda_max, 1.0 + (6.0 * Om + 32.0 * Or));
//		const double dloga = 0.0001 / lambda_max;
//		const double dhdot = -hdot - 3.0 * (2.0 * Om * delta + 32.0 / 3.0 * Or * Theta0);
//		const double ddelta = -0.5 * hdot;
//		const double dTheta0 = -eps * Theta1 - 1.0 / 6.0 * hdot;
//		const double dTheta1 = 1.0 / 3.0 * eps * Theta0;
//		hdot += dhdot * dloga;
//		delta += ddelta * dloga;
//		Theta0 += dTheta0 * dloga;
//		Theta1 += dTheta1 * dloga;
//		loga += dloga;
//		printf( "%e %e %e %e\n", a, eta, hdot, delta);
//		a = std::exp(loga);
//		hubble = Hubble(a);
//		eps = k / (a * hubble);
//		Or = omega_r / (omega_r + a * omega_m + (a * a * a * a) * (1.0 - omega_m - omega_r));
//		Om = omega_m / (omega_r / a + omega_m + (a * a * a) * (1.0 - omega_m - omega_r));
//		const double dhdot2 = -hdot - 3.0 * (2.0 * Om * delta + 32.0 / 3.0 * Or * Theta0);
//		const double ddelta2 = -0.5 * hdot;
//		const double dTheta02 = -eps * Theta1 - 1.0 / 6.0 * hdot;
//		const double dTheta12 = 1.0 / 3.0 * eps * Theta0;
//		hdot += 0.5 * (dhdot2 - dhdot) * dloga;
//		delta += 0.5 * (ddelta2 - ddelta) * dloga;
//		Theta0 += 0.5 * (dTheta02 - dTheta0) * dloga;
//		Theta1 += 0.5 * (dTheta12 - dTheta1) * dloga;
//		eta = (-1.5 * (Om * delta + 4.0 * Or * Theta0) + 0.5 * hdot)/(eps*eps);
//	}
//	Theta0 /= Theta00;
//	Theta1 /= Theta10;
//	delta /= delta0;
//	hdot /= hdot0;
//	abort();
//	return eta;
//}

int main() {
	double hdot, deltac, deltar, thetar;

//	advance2(h, eta, deltac, deltar, thetar, 1.0e-1, 1.0e-8, 1.0);
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
	zero_order_universe(amin / 1.1, amax * 1.1, cs, thomson);
	for (double k = 1e-4; k <= 2.5; k *= 1.1) {
		cos_state u;
		const auto phi = advance2(u, k, amin, amax, cs, thomson);
		printf("%e %e\n", k, u[deltaci] * u[deltaci]);
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
