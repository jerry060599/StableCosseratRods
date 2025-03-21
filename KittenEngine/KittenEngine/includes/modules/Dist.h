#pragma once
#include <math.h>
#include <glm/glm.hpp>

// Jerry Hsu 2022

namespace Kitten {
	/// <summary>
	/// Represents a distribution tracked through its mean and variance.
	/// </summary>
	struct Dist {
	public:
		double min, max;
		double X = 0;
		double XX = 0;
		int num = 0;

		Dist(int n = 0) {
			X = XX = 0;
			num = n;
			min = INFINITY;
			max = -INFINITY;
		}

		void accu(double x) {
			++num;
			X += x;
			XX += x * x;
			min = fmin(min, x);
			max = fmax(max, x);
		}

		void accu(Dist other) {
			num += other.num;
			X += other.X;
			XX += other.XX;
		}

		double mean() {
			return X / num;
		}

		double var() {
			if (num < 2)
				return 0;
			return (XX - X * X / num) / (num - 1);
		}

		double sd() {
			return sqrt(var());
		}

		static Dist zero() {
			return Dist(0);
		}
	};

	// A light weight running fit of a small degree polynomial up to a cubic
	template <int N = 4, typename T = float>
	struct PolyFit {
		typedef glm::vec<N, T, glm::defaultp> v_type;
		typedef glm::mat<N, N, T, glm::defaultp> m_type;

		m_type AtA;
		v_type Atb;

		PolyFit() {
			AtA = m_type(0);
			Atb = v_type(0);
		}

		// Accumulate a new data point
		void accu(T x, T y) {
			v_type row;
			row[0] = 1;
			for (int i = 1; i < N; i++)
				row[i] = row[i - 1] * x;

			AtA += glm::outerProduct(row, row);
			Atb += row * y;
		}

		// Returns the coefficients of the polynomial in increasing order of degree
		v_type coeff() {
			return inverse(AtA) * Atb;
		}
	};
}