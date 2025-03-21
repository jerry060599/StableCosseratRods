#pragma once
//-------------------------------------------------------------------------------

#include <cmath>
#include <intrin.h>
#include <glm/glm.hpp>

//-------------------------------------------------------------------------------
namespace cy {
	//-------------------------------------------------------------------------------

	template<typename T, typename S> inline T MultSign(T v, S sign) { return v * (sign < 0 ? T(-1) : T(1)); }	//!< Multiplies the given value with the given sign
	template <> inline float  MultSign<float, float >(float  v, float  sign) { return _mm_cvtss_f32(_mm_xor_ps(_mm_set_ss(v), _mm_and_ps(__m128 ({ -0.0f,-0.0f,-0.0f,-0.0f }), _mm_set_ss(sign)))); }
	template <> inline double MultSign<double, double>(double v, double sign) { return _mm_cvtsd_f64(_mm_xor_pd(_mm_set1_pd(v), _mm_and_pd(__m128d({ -0.0, -0.0 }), _mm_set1_pd(sign)))); }

	template <typename T> inline void Sort2Asc(T r[2], T const v[2]) { r[0] = std::min(v[0], v[1]); r[1] = std::max(v[0], v[1]); }

	//-------------------------------------------------------------------------------

	template <typename ftype>
	constexpr ftype PolynomialDefaultEpsilon() { return ftype(1e-4); }

	template <>
	constexpr double PolynomialDefaultEpsilon<double>() { return 1e-12; }

	//-------------------------------------------------------------------------------

	enum RootFinder {
		NEWTON,
		BISECTION,
		ITP
	};

	template <typename ftype> inline int LinearRoots(ftype& root, ftype const coef[2], ftype rangeMin = 0, ftype rangeMax = 1);
	template <typename ftype> inline int QuadraticRoots(ftype roots[2], ftype const coef[3], ftype rangeMin = 0, ftype rangeMax = 1);
	template <typename ftype> inline int CubicRoots(ftype roots[3], ftype const coef[4], ftype rangeMin = 0, ftype rangeMax = 1, ftype epsilon = PolynomialDefaultEpsilon<ftype>());

	template <unsigned int N, RootFinder rootFinder = NEWTON, typename ftype>
	inline int PolynomialRoots(ftype roots[N], ftype const coef[N + 1], ftype rangeMin = 0, ftype rangeMax = 1, ftype epsilon = PolynomialDefaultEpsilon<ftype>());

	//-------------------------------------------------------------------------------

	template <unsigned int N, typename ftype>
	inline ftype PolynomialEval(ftype const coef[N + 1], ftype x) {
		ftype r = coef[N];
		for (int i = N - 1; i >= 0; --i) {
			r = r * x + coef[i];
		}
		return r;
	}

	//-------------------------------------------------------------------------------

	template <unsigned int N, typename ftype>
	inline void PolynomialEvalWithDeriv(ftype& val, ftype& deriv, ftype const coef[N + 1], ftype x) {
		ftype  p = coef[N] * x + coef[N - 1];
		ftype dp = coef[N];
		for (int i = N - 2; i >= 0; --i) {
			dp = dp * x + p;
			p = p * x + coef[i];
		}
		val = p;
		deriv = dp;
	}

	//-------------------------------------------------------------------------------

	template <int N, typename ftype>
	struct Polynomial {
		ftype coef[N + 1];

		ftype const& operator [] (int i) const { return coef[i]; }
		ftype& operator [] (int i) { return coef[i]; }

		ftype Eval(ftype x) const { return PolynomialEval<N>(coef, x); }
		void  EvalWithDeriv(ftype& val, ftype& deriv, ftype x) const { PolynomialEvalWithDeriv<N>(val, deriv, coef, x); }

		Polynomial<N, ftype> operator + (Polynomial<N, ftype> const& p) const { Polynomial<N, ftype> r; for (int i = 0; i <= N; ++i) r[i] = coef[i] + p[i]; return r; }
		Polynomial<N, ftype> operator - (Polynomial<N, ftype> const& p) const { Polynomial<N, ftype> r; for (int i = 0; i <= N; ++i) r[i] = coef[i] - p[i]; return r; }
		void operator += (Polynomial<N, ftype> const& p) const { for (int i = 0; i <= N; ++i) coef[i] += p[i]; }
		void operator -= (Polynomial<N, ftype> const& p) const { for (int i = 0; i <= N; ++i) coef[i] -= p[i]; }

		template <int M> Polynomial<N + M, ftype> operator * (Polynomial<M, ftype> const& p) const {
			Polynomial<N + M, ftype> result;
			for (int i = 0; i <= N + M; ++i) result[i] = ftype(0);
			for (int i = 0; i <= N; ++i) {
				for (int j = 0; j <= M; ++j) {
					result[i + j] += coef[i] * p.coef[j];
				}
			}
			return result;
		}

		Polynomial<2 * N, ftype> Squared() const {
			Polynomial<2 * N, ftype> result;
			for (int i = 0; i <= 2 * N; ++i) result[i] = ftype(0);
			for (int i = 0; i <= N; ++i) {
				result[2 * i] += coef[i] * coef[i];
				for (int j = i + 1; j <= N; ++j) {
					result[i + j] += 2 * coef[i] * coef[j];
				}
			}
			return result;
		}

		Polynomial<N - 1, ftype> Deriv() const {
			Polynomial<N - 1, ftype> result;
			result[0] = coef[1];
#pragma warning(push)
#pragma warning( disable: 4244 )
			for (int i = 2; i <= N; ++i) result[i - 1] = i * coef[i];
#pragma warning(pop)
			return result;
		}

		template <RootFinder rootFinder = NEWTON>
		int Roots(ftype roots[N], ftype rangeMin = 0, ftype rangeMax = 1, ftype epsilon = PolynomialDefaultEpsilon<ftype>()) const {
			return Roots<N, rootFinder>(roots, rangeMin, rangeMax, epsilon);
		}

		//private:

		template <int N, RootFinder rootFinder>
		int Roots(ftype roots[N], ftype rangeMin, ftype rangeMax, ftype epsilon) const {
			if (coef[N] != 0) {
				Polynomial<N - 1, ftype> deriv;
				deriv[0] = coef[1];
				for (int i = 2; i <= N; ++i) deriv[i - 1] = i * coef[i];

				ftype derivRoots[N];
				int nd = deriv.Roots(derivRoots, rangeMin, rangeMax, epsilon);
				derivRoots[nd] = rangeMax;
				ftype x0 = rangeMin;
				ftype y0 = PolynomialEval<N>(coef, x0);
				int nr = 0;
				for (int i = 0; i <= nd; ++i) {
					ftype x1 = derivRoots[i];
					ftype y1 = PolynomialEval<N>(coef, x1);
					if (y0 != y1 && (y0 > 0) != (y1 > 0)) {
						bool increasing = y1 > 0;
						ftype xr;
						if (rootFinder == NEWTON) {
							xr = NewtonSafe<N>(deriv, x0, x1, y0, y1, increasing, epsilon);
						}
						if (rootFinder == BISECTION) {
							xr = Bisection<N>(x0, x1, y0, y1, increasing, epsilon);
						}
						if (rootFinder == ITP) {
							xr = ITPSteps<N>(x0, x1, y0, y1, increasing, epsilon);
						}
						roots[nr++] = xr;
					}
					x0 = x1;
					y0 = y1;
				}
				return nr;

			}
			else {
				return Roots<N - 1, rootFinder>(roots, rangeMin, rangeMax, epsilon);
			}
		}

		template <>
		int Roots<2, NEWTON>(ftype roots[2], ftype rangeMin, ftype rangeMax, ftype epsilon) const {
			return QuadraticRoots<ftype>(roots, coef, rangeMin, rangeMax);
		}
		template <>
		int Roots<2, BISECTION>(ftype roots[2], ftype rangeMin, ftype rangeMax, ftype epsilon) const {
			return QuadraticRoots<ftype>(roots, coef, rangeMin, rangeMax);
		}
		template <>
		int Roots<2, ITP>(ftype roots[2], ftype rangeMin, ftype rangeMax, ftype epsilon) const {
			return QuadraticRoots<ftype>(roots, coef, rangeMin, rangeMax);
		}

		template <>
		int Roots<3, NEWTON>(ftype roots[3], ftype rangeMin, ftype rangeMax, ftype epsilon) const {
			return CubicRoots<ftype>(roots, coef, rangeMin, rangeMax, epsilon);
		}

		template <unsigned int N>
		inline ftype NewtonSafe(Polynomial<N - 1, ftype> const& deriv, ftype x0, ftype x1, ftype y0, ftype y1, int increasing, ftype epsilon) const {
			ftype len = x1 - x0;
			ftype xr = (x0 + x1) * ftype(0.5);
			ftype ep2 = 2 * epsilon;
			if (len <= ep2) return xr;

			ftype xbounds[2] = { x0, x1 };
			ftype ybounds[2] = { y0, y1 };

			ftype yr = PolynomialEval<N>(coef, xr);
			int i = increasing ^ (yr <= 0);
			while (true) {

				xbounds[i] = xr;
				ybounds[i] = yr;
				if (xbounds[1] - xbounds[0] <= ep2) return (xbounds[0] + xbounds[1]) * ftype(0.5);

				ftype dy = deriv.Eval(xr);
				ftype dx = yr / dy;
				ftype xn = xr - dx;
				if (std::abs(dx) > epsilon) {
					if (xn > xbounds[0] && xn < xbounds[1]) {
						xr = xn;
					}
					else {	// if Newton step failed (also catches NANs)
						xr = (xbounds[0] + xbounds[1]) * ftype(0.5);
					}
					yr = PolynomialEval<N>(coef, xr);
					i = increasing ^ (yr <= 0);
				}
				else {
					xr = xr - MultSign(epsilon, i - 1);
					yr = PolynomialEval<N>(coef, xr);
					int j = increasing ^ (yr <= 0);
					if (i != j) return xn;
				}
			}
		}

		template <>
		inline ftype NewtonSafe<3>(Polynomial<2, ftype> const& deriv, ftype x0, ftype x1, ftype y0, ftype y1, int increasing, ftype epsilon) const {
			ftype len = x1 - x0;
			ftype xr = (x0 + x1) * ftype(0.5);
			ftype ep2 = 2 * epsilon;
			if (len <= ep2) return xr;

			ftype xbounds[2] = { x0, x1 };

			ftype yr = PolynomialEval<3>(coef, xr);
			int i = increasing ^ (yr <= 0);

			xbounds[i] = xr;
			ftype xb0 = xbounds[0] + epsilon;
			ftype xb1 = xbounds[1] - epsilon;
			xr -= yr / deriv.Eval(xr);
			glm::clamp(xr, xb0, xb1);
			xr -= PolynomialEval<3>(coef, xr) / deriv.Eval(xr);
			glm::clamp(xr, xb0, xb1);
			xr -= PolynomialEval<3>(coef, xr) / deriv.Eval(xr);
			glm::clamp(xr, xb0, xb1);
			ftype xn = xr - PolynomialEval<3>(coef, xr) / deriv.Eval(xr);
			glm::clamp(xn, xb0, xb1);
			while (std::abs(xr - xn) > epsilon) {
				xr = xn;
				ftype yr = PolynomialEval<3>(coef, xn);
				ftype dx = yr / deriv.Eval(xn);
				xn = xr - dx;
				glm::clamp(xn, xb0, xb1);
			}
			return xn;
		}

		template <unsigned int N>
		inline ftype Bisection(ftype x0, ftype x1, ftype y0, ftype y1, int increasing, ftype epsilon) const {
			ftype len = x1 - x0;
			ftype xr = (x0 + x1) * ftype(0.5);
			ftype ep2 = 2 * epsilon;
			if (len <= ep2) return xr;

			ftype xbounds[2] = { x0, x1 };
			ftype ybounds[2] = { y0, y1 };

			ftype yr = PolynomialEval<N>(coef, xr);
			int i = increasing ^ (yr <= 0);
			xbounds[i] = xr;
			ybounds[i] = yr;
			len = xbounds[1] - xbounds[0];

			while (len > ep2) {
				xr = (xbounds[0] + xbounds[1]) * ftype(0.5);
				yr = PolynomialEval<N>(coef, xr);
				i = increasing ^ (yr <= 0);
				xbounds[i] = xr;
				ybounds[i] = yr;
				len = xbounds[1] - xbounds[0];
			}

			return (xbounds[0] + xbounds[1]) * ftype(0.5);
		}

		static ftype RegulaFalsi(ftype x0, ftype x1, ftype y0, ftype y1) { return (y1 * x0 - y0 * x1) / (y1 - y0); }
		static ftype RegulaFalsi(ftype const x[2], ftype const y[2]) { return RegulaFalsi(x[0], x[1], y[0], y[1]); }

		template <unsigned int N>
		inline ftype ITPSteps(ftype x0, ftype x1, ftype y0, ftype y1, int increasing, ftype epsilon) const {
			const ftype k1 = ftype(1.0);
			const ftype k2 = ftype(2.0);
			const ftype n0 = ftype(0.0);

			int j = 0;
			ftype nmax = std::ceil(std::log2((x1 - x0) / (2 * epsilon))) + n0;

			ftype xbounds[2], ybounds[2];
			xbounds[0] = x0;
			xbounds[1] = x1;
			ybounds[0] = y0;
			ybounds[1] = y1;
			ftype xr = (xbounds[0] + xbounds[1]) * ftype(0.5);

			while (xbounds[1] - xbounds[0] > 2 * epsilon) {
				ftype xm = xr;
				ftype r = epsilon * std::pow(ftype(2), nmax - j) - (xbounds[1] - xbounds[0]) * ftype(0.5);
				ftype delta = k1 * std::pow(std::abs(xbounds[1] - xbounds[0]), k2);
				j++;

				// Interpolation
				xr = RegulaFalsi(xbounds, ybounds);

				// Truncation
				ftype sigma_val = xm - xr;
				if (delta <= std::abs(sigma_val)) {
					xr += MultSign(delta, sigma_val);
				}
				else {
					xr = xm;
				}

				// Projection
				if (std::abs(xr - xm) > r) {
					xr = xm - MultSign(r, sigma_val);
				}

				ftype yr = PolynomialEval<N>(coef, xr);
				int i = increasing ^ (yr <= 0);
				xbounds[i] = xr;
				ybounds[i] = yr;

				xr = (xbounds[0] + xbounds[1]) * ftype(0.5);
			}
			return xr;
		}

	};

	//-------------------------------------------------------------------------------

	template <typename ftype>
	inline int LinearRoots(ftype& root, ftype const coef[2], ftype rangeMin, ftype rangeMax) {
		if (coef[1] != 0) {
			ftype r = -coef[0] / coef[1];
			root = r;
			return (r >= rangeMin && r <= rangeMax);
		}
		else {
			root = (rangeMin + rangeMax) / 2;
			return coef[0] == 0;
		}
	}

	//-------------------------------------------------------------------------------

	template <typename ftype>
	inline int QuadraticRoots(ftype roots[2], ftype const coef[3], ftype rangeMin, ftype rangeMax) {
		ftype c = coef[0];
		ftype b = coef[1];
		ftype a = coef[2];
		ftype delta = b * b - 4 * a * c;
		if (delta >= 0) {
			ftype d = sqrt(delta);
			ftype rv[2];
			ftype t = ftype(-0.5) * (b + MultSign(d, b));
			rv[0] = t / a;
			rv[1] = c / t;
			ftype r[2];
			Sort2Asc(r, rv);
			int r0i = (r[0] >= rangeMin) & (r[0] <= rangeMax);
			int r1i = (r[1] >= rangeMin) & (r[1] <= rangeMax);
			roots[0] = r[0];
			roots[r0i] = r[1];
			return r0i + r1i;
		}
		return 0;
	}

	template <>
	inline int QuadraticRoots<float>(float roots[2], float const coef[3], float rangeMin, float rangeMax) {
		__m128 _0abc = _mm_set_ps(0.0f, coef[2], coef[1], coef[0]);
		__m128 _02a2b2c = _mm_add_ps(_0abc, _0abc);
		__m128 _2a2c_bb = _mm_shuffle_ps(_0abc, _02a2b2c, _MM_SHUFFLE(2, 0, 1, 1));
		__m128 _2c2a_bb = _mm_shuffle_ps(_0abc, _02a2b2c, _MM_SHUFFLE(0, 2, 1, 1));
		__m128 _4ac_b2 = _mm_mul_ps(_2a2c_bb, _2c2a_bb);
		__m128 _4ac = _mm_shuffle_ps(_4ac_b2, _4ac_b2, _MM_SHUFFLE(2, 2, 2, 2));
		if (_mm_comige_ss(_4ac_b2, _4ac)) {
			__m128 delta = _mm_sub_ps(_4ac_b2, _4ac);
			__m128 sqrtd = _mm_sqrt_ss(delta);
			__m128 signb = _mm_set_ps(-0.0f, -0.0f, -0.0f, -0.0f);
			__m128 db = _mm_xor_ps(sqrtd, _mm_and_ps(_2a2c_bb, signb));
			__m128 b_db = _mm_add_ss(_2a2c_bb, db);
			__m128 _2t = _mm_xor_ps(b_db, signb);
			__m128 _2c_2t = _mm_shuffle_ps(_2t, _02a2b2c, _MM_SHUFFLE(0, 0, 0, 0));
			__m128 _2t_2a = _mm_shuffle_ps(_02a2b2c, _2t, _MM_SHUFFLE(0, 0, 2, 2));
			__m128 rv = _mm_div_ps(_2c_2t, _2t_2a);
			__m128 r0 = _mm_min_ps(rv, _mm_shuffle_ps(rv, rv, _MM_SHUFFLE(3, 2, 1, 2)));
			__m128 r = _mm_max_ps(r0, _mm_shuffle_ps(r0, r0, _MM_SHUFFLE(3, 2, 2, 0)));
			__m128 range = _mm_set_ps(rangeMax, rangeMax, rangeMin, rangeMin);
			__m128 minT = _mm_cmpge_ps(r, range);
			__m128 maxT = _mm_cmple_ps(r, _mm_shuffle_ps(range, range, _MM_SHUFFLE(3, 2, 2, 2)));
			__m128 valid = _mm_and_ps(minT, _mm_shuffle_ps(maxT, _mm_setzero_ps(), _MM_SHUFFLE(3, 2, 1, 0)));
			__m128 rr = _mm_blendv_ps(_mm_shuffle_ps(r, r, _MM_SHUFFLE(3, 2, 0, 1)), r, valid);
			roots[0] = _mm_cvtss_f32(rr);
			roots[1] = _mm_cvtss_f32(_mm_shuffle_ps(rr, rr, _MM_SHUFFLE(3, 2, 0, 1)));
			return _mm_popcnt_u32(_mm_movemask_ps(valid));
		}
		return 0;
	}

	//-------------------------------------------------------------------------------

	template <typename ftype>
	inline int CubicRoots(ftype roots[3], ftype const coef[4], ftype rangeMin, ftype rangeMax, ftype epsilon) {
		class Newton {
			ftype const* coef;
			ftype epsilon;
		public:
			Newton(ftype const* c, ftype eps) : coef(c), epsilon(eps) {}
			inline ftype Iter(ftype const* deriv, ftype x0, ftype x1) {
				ftype len = x1 - x0;
				ftype xr = (x0 + x1) * ftype(0.5);
				ftype ep2 = 2 * epsilon;
				if (len <= ep2) return xr;
				xr -= PolynomialEval<3>(coef, xr) / PolynomialEval<2>(deriv, xr);
				glm::clamp(xr, x0, x1);
				xr -= PolynomialEval<3>(coef, xr) / PolynomialEval<2>(deriv, xr);
				glm::clamp(xr, x0, x1);
				xr -= PolynomialEval<3>(coef, xr) / PolynomialEval<2>(deriv, xr);
				glm::clamp(xr, x0, x1);
				ftype xn = xr - PolynomialEval<3>(coef, xr) / PolynomialEval<2>(deriv, xr);
				glm::clamp(xn, x0, x1);
				while (std::abs(xr - xn) > epsilon) {
					xr = xn;
					xn = xr - PolynomialEval<3>(coef, xn) / PolynomialEval<2>(deriv, xn);
					glm::clamp(xn, x0, x1);
				}
				return xn;
			}

			inline ftype IterSafe(ftype const* deriv, ftype x0, ftype x1, ftype y0) {
				ftype len = x1 - x0;
				ftype xr = (x0 + x1) * ftype(0.5);
				ftype ep2 = 2 * epsilon;
				if (len <= ep2) return xr;
				ftype yr = PolynomialEval<3>(coef, xr);
				ftype xbounds[2] = { x0, x1 };
				int i = (y0 <= 0) ^ (yr <= 0);
				xbounds[i] = xr;
				ftype xb0 = xbounds[0] + epsilon;
				ftype xb1 = xbounds[1] - epsilon;
				xr -= yr / PolynomialEval<2>(deriv, xr);
				glm::clamp(xr, xb0, xb1);
				xr -= PolynomialEval<3>(coef, xr) / PolynomialEval<2>(deriv, xr);
				glm::clamp(xr, xb0, xb1);
				xr -= PolynomialEval<3>(coef, xr) / PolynomialEval<2>(deriv, xr);
				glm::clamp(xr, xb0, xb1);
				ftype xn = xr - PolynomialEval<3>(coef, xr) / PolynomialEval<2>(deriv, xr);
				glm::clamp(xn, xb0, xb1);
				while (std::abs(xr - xn) > epsilon) {
					xr = xn;
					xn = xr - PolynomialEval<3>(coef, xn) / PolynomialEval<2>(deriv, xn);
					glm::clamp(xn, xb0, xb1);
				}
				return xn;
			}
		};

		Newton newton(coef, epsilon);

		ftype y0 = PolynomialEval<3>(coef, rangeMin);
		ftype y1 = PolynomialEval<3>(coef, rangeMax);

		ftype c = coef[1];
		ftype b = 2 * coef[2];
		ftype a = 3 * coef[3];
		ftype delta = b * b - 4 * a * c;
		if (delta >= 0) {
			ftype d = sqrt(delta);
			ftype rv[2];
			ftype t = ftype(-0.5) * (b + MultSign(d, b));
			rv[0] = t / a;
			rv[1] = c / t;
			ftype r[2];
			Sort2Asc(r, rv);
			int r00 = r[0] >= rangeMin;
			int r10 = r[1] >= rangeMin;
			int r01 = r[0] >= rangeMax;
			int r11 = r[1] >= rangeMax;
			int rcase = (r11 << 3) | (r01 << 2) | (r10 << 1) | r00;
			ftype deriv[4] = { c, b, a, 0 };
			ftype yr0, yr1;
			int numRoots;
			switch (rcase) {
				// (r[0] < rangeMin && r[1] < rangeMin) || (r[0] >= rangeMax && r[1] >= rangeMin)
			case 0:
			case 15:
				if ((y0 > 0) != (y1 > 0)) {
					roots[0] = newton.Iter(deriv, rangeMin, rangeMax);
					return 1;
				}
				return 0;
				// (r[0] < rangeMin && r[1] >= rangeMin && r[1] < rangeMax)
			case 2:
				yr1 = PolynomialEval<3>(coef, r[1]);
				numRoots = 0;
				if ((y0 > 0) != (yr1 > 0)) {
					roots[0] = newton.IterSafe(deriv, rangeMin, r[1], y0);
					numRoots++;
				}
				if ((yr1 > 0) != (y1 > 0)) {
					roots[numRoots++] = newton.Iter(deriv, r[1], rangeMax);
				}
				return numRoots;
				// (r[0] >= rangeMin && r[0] < rangeMax && r[1] >= rangeMin && r[1] < rangeMax)
			case 3:
				yr0 = PolynomialEval<3>(coef, r[0]);
				yr1 = PolynomialEval<3>(coef, r[1]);
				numRoots = 0;
				if ((y0 > 0) != (yr0 > 0)) {
					roots[0] = newton.Iter(deriv, rangeMin, r[0]);
					numRoots++;
				}
				if ((yr0 > 0) != (yr1 > 0)) {
					roots[numRoots++] = newton.IterSafe(deriv, r[0], r[1], yr0);
				}
				if ((yr1 > 0) != (y1 > 0)) {
					roots[numRoots++] = newton.Iter(deriv, r[1], rangeMax);
				}
				return numRoots;
				// (r[0] < rangeMin && r[1] >= rangeMax)
			case 10:
				if ((y0 > 0) != (y1 > 0)) {
					roots[0] = newton.IterSafe(deriv, rangeMin, rangeMax, y0);
					return 1;
				}
				return 0;
				// (r[0] >= rangeMin && r[0] < rangeMax && r[1] >= rangeMax)
			case 11:
				yr0 = PolynomialEval<3>(coef, r[0]);
				numRoots = 0;
				if ((y0 > 0) != (yr0 > 0)) {
					roots[0] = newton.Iter(deriv, rangeMin, r[0]);
					numRoots++;
				}
				if ((yr0 > 0) != (y1 > 0)) {
					roots[numRoots++] = newton.IterSafe(deriv, r[0], rangeMax, yr0);
				}
				return numRoots;
				//nodefault;
			}
		}
		else {
			if ((y0 > 0) != (y1 > 0)) {
				ftype deriv[4] = { c, b, a, 0 };
				roots[0] = newton.Iter(deriv, rangeMin, rangeMax);
				return 1;
			}
			return 0;
		}
	}


	/*
	template <>
	inline int CubicRoots<float>( float roots[3], float const coef[4], float rangeMin, float rangeMax, float epsilon )
	{
		class Newton
		{
			__m128 abcd;
			__m128 bdcb;
			__m128 ca_3a0;
			float epsilon;
		public:
			Newton( __m128 _abcd, __m128 _3a2bc0, float eps ) : abcd(_abcd), epsilon(eps)
			{
				bdcb = _mm_shuffle_ps( abcd, abcd, _MM_SHUFFLE(2,0,1,2) );
				ca_3a0 = _mm_shuffle_ps( _3a2bc0, abcd, _MM_SHUFFLE(1,3,3,0) );
			}

			inline __m128 PolyDeriv( __m128 _ca3ax )
			{
				__m128 _xxxx            = _mm_shuffle_ps( _ca3ax, _ca3ax, _MM_SHUFFLE(0,0,0,0) );
				__m128 _cx_ax_3ax_x2    = _mm_mul_ps( _ca3ax, _xxxx );
				__m128 _b_ax_3ax_b      = _mm_blend_ps( _cx_ax_3ax_x2, bdcb, 0b1001 );
				__m128 _x2_x2_x_x       = _mm_shuffle_ps( _xxxx, _cx_ax_3ax_x2, _MM_SHUFFLE(0,0,0,0) );
				__m128 _bx2_ax3_3ax2_bx = _mm_mul_ps( _b_ax_3ax_b, _x2_x2_x_x );
				__m128 _cx_ax3_3ax2_bx  = _mm_blend_ps( _bx2_ax3_3ax2_bx, _cx_ax_3ax_x2, 0b1000 );
				__m128 _bx2_d_c_bx      = _mm_blend_ps( _bx2_ax3_3ax2_bx, bdcb, 0b0110 );
				__m128 _cxbx2_ax3d_3ax2c_2bx = _mm_add_ps( _cx_ax3_3ax2_bx, _bx2_d_c_bx );
				__m128 _ax3d_cxbx2_2bx_3ax2c = _mm_shuffle_ps( _cxbx2_ax3d_3ax2c_2bx, _cxbx2_ax3d_3ax2c_2bx, _MM_SHUFFLE(2,3,0,1) );
				__m128 _p_p_d_d = _mm_add_ps( _cxbx2_ax3d_3ax2c_2bx, _ax3d_cxbx2_2bx_3ax2c );
				return _p_p_d_d;
			}
			inline __m128 TakeStep( __m128 _ca3ax, __m128 _p_p_d_d, __m128 _x0, __m128 _x1 )
			{
				__m128 _deltax  = _mm_div_ss( _mm_shuffle_ps(_p_p_d_d,_p_p_d_d,_MM_SHUFFLE(2,2,2,2)), _p_p_d_d );
				__m128 _ca3ax0 = _mm_sub_ss( _ca3ax, _deltax );
				__m128 _ca3ax1 = _mm_max_ss( _ca3ax0, _x0 );
				__m128 _ca3axn = _mm_min_ss( _ca3ax1, _x1 );
				return _ca3axn;
			}

			inline __m128 Step( __m128 _ca3ax, __m128 _x0, __m128 _x1 )
			{
				__m128 _p_p_d_d = PolyDeriv( _ca3ax );
				return TakeStep( _ca3ax, _p_p_d_d, _x0, _x1 );
			}

			inline float Iter( float x0, float x1 ) {
				float len = x1 - x0;
				float xr = (x0 + x1) * float(0.5);
				float ep2 = 2*epsilon;
				if ( len <= ep2 ) return xr;
				__m128 _xr = _mm_set_ss( xr );
				__m128 _x0 = _mm_set_ss( x0 );
				__m128 _x1 = _mm_set_ss( x1 );
				__m128 _eps  = _mm_set_ss(epsilon);
				__m128 _neps = _mm_set_ss(-epsilon);
				__m128 _ca3ax = _mm_blend_ps( ca_3a0, _xr, 0b0001 );
				__m128 _ca3axn = Step( _ca3ax, _x0, _x1 );
				__m128 dx = _mm_sub_ss( _ca3ax, _ca3axn );
				_ca3ax = _ca3axn;
				while ( _mm_comigt_ss( dx, _eps ) || _mm_comilt_ss( dx, _neps ) ) {
					_ca3axn = Step( _ca3ax, _x0, _x1 );
					dx = _mm_sub_ss( _ca3ax, _ca3axn );
					_ca3ax = _ca3axn;
				}
				return _mm_cvtss_f32(_ca3ax);
			}

			inline float IterSafe( float x0, float x1, float y0 )
			{
				//return Iter(x0,x1);
				float len = x1 - x0;
				float xr = (x0 + x1) * float(0.5);
				float ep2 = 2*epsilon;
				if ( len <= ep2 ) return xr;
				__m128 _xr = _mm_set_ss( xr );
				__m128 _eps  = _mm_set_ss(epsilon);
				__m128 _neps = _mm_set_ss(-epsilon);
				__m128 _ca3ax = _mm_blend_ps( ca_3a0, _xr, 0b0001 );
				__m128 _p_p_d_d = PolyDeriv( _ca3ax );
				float yr = _mm_cvtss_f32( _mm_shuffle_ps(_p_p_d_d,_p_p_d_d,_MM_SHUFFLE(2,2,2,2)) );
				float xbounds[2] = { x0, x1 };
				int i = (y0<=0) ^ (yr<=0);
				xbounds[i] = xr;
				float xb0 = xbounds[0] + epsilon;
				float xb1 = xbounds[1] - epsilon;
				__m128 _x0 = _mm_set_ss( xb0 );
				__m128 _x1 = _mm_set_ss( xb1 );
				__m128 _ca3axn = TakeStep( _ca3ax, _p_p_d_d, _x0, _x1 );
				__m128 dx = _mm_sub_ss( _ca3ax, _ca3axn );
				_ca3ax = _ca3axn;
				while ( _mm_comigt_ss( dx, _eps ) || _mm_comilt_ss( dx, _neps ) ) {
					_ca3axn = Step( _ca3ax, _x0, _x1 );
					dx = _mm_sub_ss( _ca3ax, _ca3axn );
					_ca3ax = _ca3axn;
				}
				return _mm_cvtss_f32(_ca3ax);
			}

		};

		__m128 _abcd = _mm_loadu_ps(coef);

		float y0 = PolynomialEval<3>( coef, rangeMin );
		float y1 = PolynomialEval<3>( coef, rangeMax );

		// coefficients of the derivative
		__m128 _3210     = _mm_set_ps( 3.0f, 2.0f, 1.0f, 0.0f );
		__m128 q_abc0    = _mm_mul_ps( _abcd, _3210 );
		Newton newton( _abcd, q_abc0, epsilon );

		__m128 q_2a2b2c0 = _mm_add_ps( q_abc0, q_abc0 );
		__m128 q_2a2c_bb = _mm_shuffle_ps( q_abc0, q_2a2b2c0, _MM_SHUFFLE(3,1,2,2) );
		__m128 q_2c2a_bb = _mm_shuffle_ps( q_abc0, q_2a2b2c0, _MM_SHUFFLE(1,3,2,2) );
		__m128 q_4ac_b2  = _mm_mul_ps( q_2a2c_bb, q_2c2a_bb );
		__m128 q_4ac     = _mm_shuffle_ps( q_4ac_b2, q_4ac_b2, _MM_SHUFFLE(2,2,2,2) );
		if ( _mm_comige_ss( q_4ac_b2, q_4ac ) ) {
			__m128 delta  = _mm_sub_ps( q_4ac_b2, q_4ac );
			__m128 sqrtd  = _mm_sqrt_ss(delta);
			__m128 signb  = _mm_set_ps(-0.0f,-0.0f,-0.0f,-0.0f);
			__m128 db     = _mm_xor_ps( sqrtd, _mm_and_ps( q_2a2c_bb, signb ) );
			__m128 b_db   = _mm_add_ss( q_2a2c_bb, db );
			__m128 _2t    = _mm_xor_ps( b_db, signb );
			__m128 _2c_2t = _mm_shuffle_ps( _2t,   q_2a2c_bb, _MM_SHUFFLE(2,2,0,0) );
			__m128 _2t_2a = _mm_shuffle_ps( q_2a2c_bb, _2t,   _MM_SHUFFLE(0,0,3,3) );
			__m128 _rv    = _mm_div_ps( _2c_2t, _2t_2a );
			__m128 _r0    = _mm_min_ps( _rv, _mm_shuffle_ps( _rv, _rv, _MM_SHUFFLE(3,0,1,2) ) );
			__m128 _r     = _mm_max_ps( _r0, _mm_shuffle_ps( _r0, _r0, _MM_SHUFFLE(1,2,3,0) ) );
			__m128 range  = _mm_set_ps( rangeMax, rangeMax, rangeMin, rangeMin );
			__m128 _rcmp  = _mm_cmpge_ps( _r, range );
			int rcase = _mm_movemask_ps(_rcmp);
			float yr0, yr1;
			int numRoots;
			float r0 = _mm_cvtss_f32(_r);
			float r1 = _mm_cvtss_f32(_mm_shuffle_ps(_r,_r,_MM_SHUFFLE(1,1,1,1)));
			switch ( rcase ) {
				// (r[0] < rangeMin && r[1] < rangeMin) || (r[0] >= rangeMax && r[1] >= rangeMin)
				case 0:
				case 15:
					if ( (y0>0) != (y1>0) ) {
						roots[0] = newton.Iter( rangeMin, rangeMax );
						return 1;
					}
					return 0;
				// (r[0] < rangeMin && r[1] >= rangeMin && r[1] < rangeMax)
				case 2:
					yr1 = PolynomialEval<3>( coef, r1 );
					numRoots = 0;
					if ( (y0>0) != (yr1>0) ) {
						roots[0] = newton.IterSafe( rangeMin, r1, y0 );
						numRoots++;
					}
					if ( (yr1>0) != (y1>0) ) {
						roots[numRoots++] = newton.Iter( r1, rangeMax );
					}
					return numRoots;
				// (r[0] >= rangeMin && r[0] < rangeMax && r[1] >= rangeMin && r[1] < rangeMax)
				case 3:
					yr0 = PolynomialEval<3>( coef, r0 );
					yr1 = PolynomialEval<3>( coef, r1 );
					numRoots = 0;
					if ( (y0>0) != (yr0>0) ) {
						roots[0] = newton.Iter( rangeMin, r0 );
						numRoots++;
					}
					if ( (yr0>0) != (yr1>0) ) {
						roots[numRoots++] = newton.IterSafe( r0, r1, yr0 );
					}
					if ( (yr1>0) != (y1>0) ) {
						roots[numRoots++] = newton.Iter( r1, rangeMax );
					}
					return numRoots;
				// (r[0] < rangeMin && r[1] >= rangeMax)
				case 10:
					if ( (y0>0) != (y1>0) ) {
						roots[0] = newton.IterSafe( rangeMin, rangeMax, y0 );
						return 1;
					}
					return 0;
				// (r[0] >= rangeMin && r[0] < rangeMax && r[1] >= rangeMax)
				case 11:
					yr0 = PolynomialEval<3>( coef, r0 );
					numRoots = 0;
					if ( (y0>0) != (yr0>0) ) {
						roots[0] = newton.Iter( rangeMin, r0 );
						numRoots++;
					}
					if ( (yr0>0) != (y1>0) ) {
						roots[numRoots++] = newton.IterSafe( r0, rangeMax, yr0 );
					}
					return numRoots;
				nodefault;
			}
		} else {
			if ( (y0>0) != (y1>0) ) {
				roots[0] = newton.Iter( rangeMin, rangeMax );
				return 1;
			}
			return 0;
		}
	}
	//*/

	//-------------------------------------------------------------------------------

	template <unsigned int N, RootFinder rootFinder, typename ftype>
	inline int PolynomialRoots(ftype roots[N], ftype const coef[N + 1], ftype rangeMin, ftype rangeMax, ftype epsilon) {
		Polynomial<N, ftype> poly;
		MemCopy(poly.coef, coef, N + 1);
		return poly.Roots<rootFinder>(roots, rangeMin, rangeMax, epsilon);
	}

	template <> inline int PolynomialRoots<1, NEWTON, float >(float  roots[1], float  const coef[2], float  rangeMin, float  rangeMax, float  epsilon) { return LinearRoots(*roots, coef, rangeMin, rangeMax); }
	template <> inline int PolynomialRoots<1, BISECTION, float >(float  roots[1], float  const coef[2], float  rangeMin, float  rangeMax, float  epsilon) { return LinearRoots(*roots, coef, rangeMin, rangeMax); }
	template <> inline int PolynomialRoots<1, ITP, float >(float  roots[1], float  const coef[2], float  rangeMin, float  rangeMax, float  epsilon) { return LinearRoots(*roots, coef, rangeMin, rangeMax); }
	template <> inline int PolynomialRoots<1, NEWTON, double>(double roots[1], double const coef[2], double rangeMin, double rangeMax, double epsilon) { return LinearRoots(*roots, coef, rangeMin, rangeMax); }
	template <> inline int PolynomialRoots<1, BISECTION, double>(double roots[1], double const coef[2], double rangeMin, double rangeMax, double epsilon) { return LinearRoots(*roots, coef, rangeMin, rangeMax); }
	template <> inline int PolynomialRoots<1, ITP, double>(double roots[1], double const coef[2], double rangeMin, double rangeMax, double epsilon) { return LinearRoots(*roots, coef, rangeMin, rangeMax); }


	template <> inline int PolynomialRoots<2, NEWTON, float >(float  roots[2], float  const coef[3], float  rangeMin, float  rangeMax, float  epsilon) { return QuadraticRoots(roots, coef, rangeMin, rangeMax); }
	template <> inline int PolynomialRoots<2, BISECTION, float >(float  roots[2], float  const coef[3], float  rangeMin, float  rangeMax, float  epsilon) { return QuadraticRoots(roots, coef, rangeMin, rangeMax); }
	template <> inline int PolynomialRoots<2, ITP, float >(float  roots[2], float  const coef[3], float  rangeMin, float  rangeMax, float  epsilon) { return QuadraticRoots(roots, coef, rangeMin, rangeMax); }
	template <> inline int PolynomialRoots<2, NEWTON, double>(double roots[2], double const coef[3], double rangeMin, double rangeMax, double epsilon) { return QuadraticRoots(roots, coef, rangeMin, rangeMax); }
	template <> inline int PolynomialRoots<2, BISECTION, double>(double roots[2], double const coef[3], double rangeMin, double rangeMax, double epsilon) { return QuadraticRoots(roots, coef, rangeMin, rangeMax); }
	template <> inline int PolynomialRoots<2, ITP, double>(double roots[2], double const coef[3], double rangeMin, double rangeMax, double epsilon) { return QuadraticRoots(roots, coef, rangeMin, rangeMax); }

	template <> inline int PolynomialRoots<3, NEWTON, float >(float  roots[3], float  const coef[4], float  rangeMin, float  rangeMax, float  epsilon) { return CubicRoots(roots, coef, rangeMin, rangeMax, epsilon); }
	template <> inline int PolynomialRoots<3, NEWTON, double>(double roots[3], double const coef[4], double rangeMin, double rangeMax, double epsilon) { return CubicRoots(roots, coef, rangeMin, rangeMax, epsilon); }

	//-------------------------------------------------------------------------------
} // namespace hf
//-------------------------------------------------------------------------------
