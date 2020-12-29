#include <stdio.h>
#include <array>
#include <complex>
#include <vector>
#include <cassert>

#define NDIM 3

const double h = 0.7;
const double omega_b = 0.05;
const double omega_c = 0.25;
const double Theta = 1.01;
const double omega_m = omega_b + omega_c;
//const double zeq = 2.5e+4 * std::pow(Theta, -4) * omega_m * h * h;
//const double aeq = 1.0 / (zeq + 1.0);
const double omega_gam = 5e-5;
const double omega_nu = 3.5e-5;
const double omega_r = omega_nu + omega_gam;
const double omega_lambda = 1.0 - omega_m - omega_r;
constexpr double clight = 2.9979e10;
const double H0 = 100.0 * 1e5 / clight;
const double ns = 1;
constexpr auto hplanck = 6.6260755e-27;
const double sigma_T = 2.3e-5;
const double hefrac = 0.24;
constexpr double H0cgs = 3.24e-18;
const double Yp = (1 - hefrac / 2.0);
constexpr auto kb = 1.380658e-16;
constexpr auto me = 9.1093897e-28;
constexpr double G = 6.67e-8;
constexpr double mh = 1.67e-24;

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

double chemical_potential(double ne, double T) {
	const auto beta = kb * T / (me * clight * clight);
	double eta = 0.0;
	double deta = 1.0e-6;
	double deta0;
	int iter = 0;
	do {
		const auto eta0 = eta - 0.5 * deta;
		const auto eta1 = eta + 0.5 * deta;
		const double f0 = ne - (compute_Nele(eta0, beta) - compute_Npos(eta0, beta));
		const double f1 = ne - (compute_Nele(eta1, beta) - compute_Npos(eta1, beta));
		const double dfdeta = (f1 - f0) / deta;
		const double f = 0.5 * (f0 + f1);
		deta0 = -f / dfdeta;
		eta += deta0;
		iter++;
		if (iter > 500) {
			printf("chem_pot : %e %e %e %e \n", ne, T, eta, deta0);
		}
		if (iter == 1000) {
			printf("maximum iterations exceeded in chemical_potential\n");
			abort();
		}
	} while (std::abs(deta0) > 1.0e-6);
	return eta;
}

double compute_electron_fraction(double rho, double T) {
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

	double H = n_nuc * (1 - hefrac);
	double He = n_nuc * hefrac / 4;
	double err;
	double ne = 0.5 * (H + 2 * He);
	double Y;
#define Power std::pow
	int iters = 0;
//	printf("p&e begin\n");
	do {
//		printf("%e\n", err);
		const auto ne0 = ne;
		const auto dne = -((ne - (A0 * H) / (A0 + ne) - (B0 * He * (2 * C0 + ne)) / (Power(ne, 2) + B0 * (C0 + ne)))
				/ (1 + (A0 * H) / Power(A0 + ne, 2) + (B0 * He * (B0 * C0 + ne * (4 * C0 + ne))) / Power(Power(ne, 2) + B0 * (C0 + ne), 2)));
		ne = std::min(std::max(0.5 * ne, ne + dne), 2 * ne);
		err = std::abs(std::log(ne / ne0));
		Y = ne / (H + 4 * He);
		iters++;
		if (Y < 1e-100) {
			break;
		}
		if (iters > 500) {
			printf("efrac: %e %e %e %e \n", err, ne, Y, T);
		}
		if (iters > 1000) {
			printf("Max iters exceed in compute_electron_fraction\n");
			abort();
		}
	} while (err > 1.0e-6);
//	printf("p&e end\n");
	return Y;
}

double fermi_dirac(double k, double eta, double beta) {
	double sum = 0.0;
	int N = 256;
	double dx = 1.0 / N;
	for (int i = 0; i < N; i++) {
		double x = (i + 0.5) * dx;
		double sum1 = std::pow(x, k) * std::sqrt(1.0 + beta * x) / (std::exp(x - eta) + 1.0);
		double sum2 = std::pow(x, -k - 2) * std::sqrt(1.0 + beta / x) / (std::exp(1.0 / x - eta) + 1.0);
		sum += (sum1 + sum2) * dx;
	}
	return sum;
}

double compute_Nele(double eta, double beta) {
	constexpr double c0 = 8.0 * M_PI * std::sqrt(2) * std::pow(me * clight / hplanck, 3);
	return c0 * std::pow(beta, 1.5) * (fermi_dirac(0.5, eta, beta) + fermi_dirac(1.5, eta, beta));
}

double compute_Npos(double eta, double beta) {
	constexpr double c0 = 8.0 * M_PI * std::sqrt(2) * std::pow(me * clight / hplanck, 3);
	return c0 * std::pow(beta, 1.5) * (fermi_dirac(0.5, -eta - 2.0 / beta, beta) + beta * fermi_dirac(1.5, -eta - 2.0 / beta, beta));
}

void electron_state(double ne, double T, double &P, double &E) {
	const auto beta = kb * T / (me * clight * clight);
	const auto eta = chemical_potential(ne, T);
	constexpr double c0 = 8.0 * M_PI * std::sqrt(2) * std::pow(me * clight / hplanck, 3) * me * clight * clight;
	const double F32ele = fermi_dirac(1.5, eta, beta);
	const double F52ele = fermi_dirac(2.5, eta, beta);
	const double F32pos = fermi_dirac(1.5, -eta - 2.0 / beta, beta);
	const double F52pos = fermi_dirac(2.5, -eta - 2.0 / beta, beta);
	const double Npos = compute_Npos(eta, beta);
	const double Pele = 2.0 * c0 * std::pow(beta, 2.5) * (F32ele + 0.5 * beta * F52ele);
	const double Ppos = 2.0 * c0 * std::pow(beta, 2.5) * (F32pos + 0.5 * beta * F52pos);
	const double Eele = c0 * std::pow(beta, 2.5) * (F32ele + beta * F52ele);
	const double Epos = c0 * std::pow(beta, 2.5) * (F32pos + beta * F52pos) + 2.0 * me * clight * clight * Npos;
	P = Ppos + Pele;
	E = Epos + Eele;

}

void pressure_and_energy(double rho, double T, double &P, double &eps) {
	const double Ye = compute_electron_fraction(rho, T);
	const double ne = rho / mh * Ye;
//	printf("%e %e\n", ne, Ye);
	const double N = (rho / mh) * (1 - 3.0 * hefrac / 4.0);
	double Pele, Eele, Pnuc, Enuc;
	electron_state(ne, T, Pele, Eele);
	Pnuc = kb * N * T;
	Enuc = 1.5 * Pnuc;
	eps = (Eele + Enuc) / rho;
	P = Pele + Pnuc;
}

void equation_of_state(double rho, double T, double &P, double &cs) {
	double Pp, Pn, epsp, epsn, eps;
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

double rho_b(double a) {
	return (3.0 * omega_b * H0cgs * H0cgs * h * h) / (8.0 * M_PI * G * a * a * a);
}

auto build_electron_fraction_interpolation_function(double amin, double amax) {
	double minloga = std::log(amin);
	double maxloga = std::log(amax);
	double dloga = 1.0 / 256.0;
	int N = (maxloga - minloga) / dloga + 1;
	dloga = (maxloga - minloga) / N;
	N += 2;
//	printf( "electron fraction table size = %i\n", N);
	maxloga += dloga;
	minloga -= dloga;
	dloga = (maxloga - minloga) / N;
	std::vector<double> values(N + 1);
	for (int i = 0; i <= N; i++) {
		double loga = i * dloga + minloga;
		double a = std::exp(loga);
		values[i] = compute_electron_fraction(rho_b(a), 2.73 * Theta / a);
	}

	const auto func = [=](double loga) {
		if (loga - dloga < minloga || loga + dloga > maxloga) {
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
//		return values[i2] * x + values[i1] * (1 - x);
		return std::max(c0 + c1 * x + c2 * x * x + c3 * x * x * x, 0.0);
	};
	return func;
}

void advance(state &u, double k, double a0, double a1) {
	auto Ye = build_electron_fraction_interpolation_function(1e-10, 1);
	const double logamin = std::log(a0);
	const double logamax = std::log(a1);
	double loga = logamin;
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
		};
		compute_parameters();
		const auto lam_max = eps + 4.5 * Om + 12 * Or + eps * eps + 3.0;
		const auto dloga = std::min(1.0 / lam_max, logamax - loga);
		const auto mass = fc * u[deltaci] + fb * u[deltabi];
		const auto rad = fnu * u[N0i] + fgam * u[Theta0i];
		const auto mmom = fc * u[uci] + fb * u[ubi];
		const auto rmom = fnu * u[N1i] + fgam * u[Theta1i];
		const auto Phi = 1.5 / (eps * eps) * (4 * Or * (rad + rmom / eps) + Om * (mass + mmom / eps));
//		printf( "%e %e\n", u[Phii], Phi);
		const auto compute_dudt_exp = [&](state u) {
			state dudt;
			const auto mass = fc * u[deltaci] + fb * u[deltabi];
			const auto rad = fnu * u[N0i] + fgam * u[Theta0i];
			dudt[Phii] = (0.5 * Om * mass + 2.0 * Or * rad - ((1.0 / 3.0) * eps * eps + 1.0) * u[Phii]);
			dudt[Theta0i] = -eps * u[Theta1i] - dudt[Phii];
			dudt[Theta1i] = eps / 3.0 * (u[Theta0i] - u[Phii]);
			dudt[N0i] = -eps * u[N1i] - dudt[Phii];
			dudt[N1i] = eps / 3.0 * (u[N0i] - u[Phii]);
			dudt[deltaci] = -eps * u[uci] - 3 * dudt[Phii];
			dudt[uci] = -u[uci] - eps * u[Phii];
			dudt[deltabi] = -eps * u[ubi] - 3 * dudt[Phii];
			dudt[ubi] = -u[ubi] - eps * u[Phii];
			return dudt;
		};
		state u0 = u;
		auto dudt = compute_dudt_exp(u0);
		for (int i = 0; i < NF; i++) {
			u[i] = u0[i] + dudt[i] * dloga;
		}

		sigma = sigma_T * omega_b * h * h * Ye(loga + dloga) / (a * a * a * H);
		u[Theta1i] = (dloga * dloga * (3 * dudt[Theta1i] + dudt[ubi] * R) * sigma + 3 * R * u0[Theta1i]
				+ dloga * (3 * dudt[Theta1i] * R + 3 * sigma * u0[Theta1i] + R * sigma * u0[ubi])) / (3. * (R + dloga * sigma + dloga * fgam * R * sigma));
		u[ubi] = (dloga * dudt[ubi] * (R + dloga * fgam * R * sigma) + 3 * dloga * fgam * sigma * (dloga * dudt[Theta1i] + u0[Theta1i])
				+ R * (1 + dloga * fgam * sigma) * u0[ubi]) / (R + dloga * sigma + dloga * fgam * R * sigma);
		loga += dloga;
	}

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
	double P, cs;
	auto func = build_electron_fraction_interpolation_function(1e-12, 1);
	for (double a = 1e-12; a <= 1.0; a *= 1.5) {
		const double T = 2.73 * Theta / a;
		equation_of_state(rho_b(a), T, P, cs);
		printf("%e %e %e %e %e %e\n", a, rho_b(a), T, P, cs / clight, func(std::log(a)));
	}
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
