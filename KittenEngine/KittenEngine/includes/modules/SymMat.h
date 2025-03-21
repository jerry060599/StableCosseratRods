#pragma once
#include "Common.h"

namespace Kitten {
	// Simple 3x3 hessian.
	// There is a templated version of this but the compiler doesnt unroll it properly
	// So here is a simple version where whats needed is hand unrolled.
	struct hess3 {
		union {
			float dat[6];
			struct {
				vec3 diag;
				vec3 upperTriangle;
			};
		};

		KITTEN_FUNC_DECL hess3() {}
		KITTEN_FUNC_DECL hess3(float d) : dat{ d, d, d } {}
		KITTEN_FUNC_DECL hess3(const hess3& m) {
			for (int i = 0; i < 6; i++)
				dat[i] = m[i];
		}

		KITTEN_FUNC_DECL hess3& operator=(const hess3& m) {
			for (int i = 0; i < 6; i++)
				dat[i] = m[i];
			return *this;
		}

		KITTEN_FUNC_DECL float& operator[](const int i) {
			return dat[i];
		}

		KITTEN_FUNC_DECL const float& operator[](const int i) const {
			return dat[i];
		}

		KITTEN_FUNC_DECL hess3 operator+(const hess3& m) const {
			hess3 r;
			for (int i = 0; i < 6; i++)
				r[i] = dat[i] + m[i];
			return r;
		}

		KITTEN_FUNC_DECL hess3 operator-(const hess3& m) const {
			hess3 r;
			for (int i = 0; i < 6; i++)
				r[i] = dat[i] - m[i];
			return r;
		}

		KITTEN_FUNC_DECL hess3& operator+=(const hess3& m) {
			for (int i = 0; i < 6; i++)
				dat[i] += m[i];
			return *this;
		}

		KITTEN_FUNC_DECL hess3& operator-=(const hess3& m) {
			for (int i = 0; i < 6; i++)
				dat[i] -= m[i];
			return *this;
		}

		KITTEN_FUNC_DECL hess3 operator*(const float s) const {
			hess3 r;
			for (int i = 0; i < 6; i++)
				r[i] = dat[i] * s;
			return r;
		}

		KITTEN_FUNC_DECL hess3 operator/(const float s) const {
			float invS = 1.0f / s;
			hess3 r;
			for (int i = 0; i < 6; i++)
				r[i] = dat[i] * invS;
			return r;
		}

		KITTEN_FUNC_DECL hess3 operator-() const {
			hess3 r;
			for (int i = 0; i < 6; i++)
				r[i] = -dat[i];
			return r;
		}

		KITTEN_FUNC_DECL static hess3 outer(vec3 v) {
			hess3 m;
			m[0] = v[0] * v[0];
			m[1] = v[1] * v[1];
			m[2] = v[2] * v[2];
			m[3] = v[0] * v[1];
			m[4] = v[0] * v[2];
			m[5] = v[1] * v[2];
			return m;
		}

		KITTEN_FUNC_DECL explicit operator mat3() const {
			mat3 m;
			m[0][0] = dat[0];
			m[1][1] = dat[1];
			m[2][2] = dat[2];
			m[0][1] = m[1][0] = dat[3];
			m[0][2] = m[2][0] = dat[4];
			m[1][2] = m[2][1] = dat[5];
			return m;
		}
	};

	// Simple 4x4 hessian.
	// There is a templated version of this but the compiler doesnt unroll it properly
	// So here is a simple version where whats needed is hand unrolled.
	struct hess4 {
		union {
			float dat[10];
			struct {
				vec4 diag;
				float upperTriangle[6];
			};
		};

		KITTEN_FUNC_DECL hess4() {}
		KITTEN_FUNC_DECL hess4(float d) : dat{ d, d, d, d } {}
		KITTEN_FUNC_DECL hess4(const hess4& m) {
			for (int i = 0; i < 10; i++)
				dat[i] = m[i];
		}

		KITTEN_FUNC_DECL hess4& operator=(const hess4& m) {
			for (int i = 0; i < 10; i++)
				dat[i] = m[i];
			return *this;
		}

		KITTEN_FUNC_DECL float& operator[](const int i) {
			return dat[i];
		}

		KITTEN_FUNC_DECL const float& operator[](const int i) const {
			return dat[i];
		}

		KITTEN_FUNC_DECL hess4 operator+(const hess4& m) const {
			hess4 r;
			for (int i = 0; i < 10; i++)
				r[i] = dat[i] + m[i];
			return r;
		}

		KITTEN_FUNC_DECL hess4 operator-(const hess4& m) const {
			hess4 r;
			for (int i = 0; i < 10; i++)
				r[i] = dat[i] - m[i];
			return r;
		}

		KITTEN_FUNC_DECL hess4& operator+=(const hess4& m) {
			for (int i = 0; i < 10; i++)
				dat[i] += m[i];
			return *this;
		}

		KITTEN_FUNC_DECL hess4& operator-=(const hess4& m) {
			for (int i = 0; i < 10; i++)
				dat[i] -= m[i];
			return *this;
		}

		KITTEN_FUNC_DECL hess4 operator*(const float s) const {
			hess4 r;
			for (int i = 0; i < 10; i++)
				r[i] = dat[i] * s;
			return r;
		}

		KITTEN_FUNC_DECL hess4 operator/(const float s) const {
			float invS = 1.0f / s;
			hess4 r;
			for (int i = 0; i < 10; i++)
				r[i] = dat[i] * invS;
			return r;
		}

		KITTEN_FUNC_DECL hess4 operator-() const {
			hess4 r;
			for (int i = 0; i < 10; i++)
				r[i] = -dat[i];
			return r;
		}

		KITTEN_FUNC_DECL static hess4 outer(vec4 v) {
			hess4 m;
			m[0] = v[0] * v[0];
			m[1] = v[1] * v[1];
			m[2] = v[2] * v[2];
			m[3] = v[3] * v[3];
			m[4] = v[0] * v[1];
			m[5] = v[0] * v[2];
			m[6] = v[0] * v[3];
			m[7] = v[1] * v[2];
			m[8] = v[1] * v[3];
			m[9] = v[2] * v[3];
			return m;
		}

		KITTEN_FUNC_DECL explicit operator mat4() const {
			mat4 m;
			m[0][0] = dat[0];
			m[1][1] = dat[1];
			m[2][2] = dat[2];
			m[3][3] = dat[3];
			m[0][1] = m[1][0] = dat[4];
			m[0][2] = m[2][0] = dat[5];
			m[0][3] = m[3][0] = dat[6];
			m[1][2] = m[2][1] = dat[7];
			m[1][3] = m[3][1] = dat[8];
			m[2][3] = m[3][2] = dat[9];
			return m;
		}
	};

	// A more compact representation for symetric matrices
	template <int N, typename T>
	struct SymMat {
		typedef glm::mat<N, N, T, glm::defaultp> mat_type;
		typedef glm::vec<N, T, glm::defaultp> vec_type;
		static constexpr int DATA_LEN = (N * (N + 1)) / 2;

		// Always laid out as diagonal, then the upper triangle going column major
		union {
			T data[DATA_LEN];
			vec_type diag;		// The diagonal
		};

		// Constructors and conversions
		KITTEN_FUNC_DECL SymMat() {}

		constexpr KITTEN_FUNC_DECL SymMat(const T d) : data{} {
			for (int i = 0; i < N; i++)
				data[i] = d;
		}

		constexpr KITTEN_FUNC_DECL SymMat(const vec_type d) : data{} {
			for (int i = 0; i < N; i++)
				data[i] = d[i];
		}

		constexpr KITTEN_FUNC_DECL SymMat(const SymMat& other) : data{} {
			for (int i = 0; i < DATA_LEN; ++i)
				data[i] = other.data[i];
		}

		template <int M, typename U>
		constexpr KITTEN_FUNC_DECL explicit SymMat(const SymMat<M, U>& other) {
			constexpr int K = M < N ? M : N;
			// Copy diaganol
			for (int i = 0; i < K; i++)
				data[i] = other.data[i];
			// Copy upper triangle
			constexpr int KLim = (K * (K + 1)) / 2 - K;
			for (int i = 0; i < KLim; i++)
				data[N + i] = other.data[M + i];
			// Zero out the rest
			for (int i = K + KLim; i < DATA_LEN; i++)
				data[i] = T(0);
		}

		KITTEN_FUNC_DECL SymMat<N, T>& operator=(const SymMat<N, T>& other) {
			for (int i = 0; i < DATA_LEN; ++i)
				data[i] = other.data[i];
			return *this;
		}

		constexpr KITTEN_FUNC_DECL explicit operator mat_type() const {
			mat_type m;
			for (int i = 0; i < N; i++)
				m[i][i] = data[i];

			for (int i = 1; i < N; i++)
				for (int j = 0; j < i; j++) {
					T v = data[N + j + (i * (i - 1)) / 2];
					m[i][j] = v;
					m[j][i] = v;
				}
			return m;
		}

		constexpr KITTEN_FUNC_DECL const vec_type col(const int i) const {
			vec_type v;
			// Grab before diagonal
			for (int j = 0; j < i; j++)
				v[j] = data[j + N + (i * (i - 1)) / 2];
			// Grab diagonal
			v[i] = data[i];
			// Grab after diagonal
			for (int j = i + 1; j < N; j++)
				v[j] = data[i + N + (j * (j - 1)) / 2];

			return v;
		}

		// Linear algebra
		constexpr KITTEN_FUNC_DECL SymMat<N, T> operator+(const T scalar) const {
			SymMat<N, T> result(*this);
			for (int i = 0; i < N; i++)
				result.data[i] = data[i] + scalar;
			return result;
		}

		constexpr KITTEN_FUNC_DECL SymMat<N, T> operator-(const T scalar) const {
			SymMat<N, T> result(*this);
			for (int i = 0; i < N; i++)
				result.data[i] = data[i] - scalar;
			return result;
		}

		constexpr KITTEN_FUNC_DECL SymMat<N, T> operator+(const SymMat<N, T>& other) const {
			SymMat<N, T> result;
			for (int i = 0; i < DATA_LEN; i++)
				result.data[i] = data[i] + other.data[i];
			return result;
		}

		constexpr KITTEN_FUNC_DECL SymMat<N, T> operator-(const SymMat<N, T>& other) const {
			SymMat<N, T> result;
			for (int i = 0; i < DATA_LEN; i++)
				result.data[i] = data[i] - other.data[i];
			return result;
		}

		constexpr KITTEN_FUNC_DECL SymMat<N, T> operator*(const T scalar) const {
			SymMat<N, T> result;
			for (int i = 0; i < DATA_LEN; i++)
				result.data[i] = data[i] * scalar;
			return result;
		}

		constexpr KITTEN_FUNC_DECL SymMat<N, T> operator/(const T scalar) const {
			T invScalar = 1 / scalar;
			SymMat<N, T> result;
			for (int i = 0; i < DATA_LEN; i++)
				result.data[i] = data[i] * invScalar;
			return result;
		}

		constexpr KITTEN_FUNC_DECL SymMat<N, T> operator*(const SymMat<N, T>& other) const {
			SymMat<N, T> result;
			for (int i = 0; i < N; i++)
				result[i] = dot(col(i), other.col(i));
			for (int i = 1; i < N; i++)
				for (int j = 0; j < i; j++)
					result[N + j + (i * (i - 1)) / 2] = dot(col(i), other.col(j));
			return result;
		}

		constexpr KITTEN_FUNC_DECL vec_type operator*(const vec_type& v) const {
			vec_type result = v[0] * col(0);
			for (int i = 1; i < N; i++)
				result += v[i] * col(i);
			return result;
		}

		constexpr KITTEN_FUNC_DECL SymMat<N, T> operator-() const {
			SymMat<N, T> result;
			for (int i = 0; i < DATA_LEN; i++)
				result.data[i] = -data[i];
			return result;
		}

		constexpr KITTEN_FUNC_DECL SymMat<N, T>& operator+=(const T scaler) {
			for (int i = 0; i < N; i++)
				data[i] += scaler;
			return *this;
		}

		constexpr KITTEN_FUNC_DECL SymMat<N, T>& operator-=(const T scaler) {
			for (int i = 0; i < N; i++)
				data[i] -= scaler;
			return *this;
		}

		constexpr KITTEN_FUNC_DECL SymMat<N, T>& operator+=(const SymMat<N, T>& other) {
			for (int i = 0; i < DATA_LEN; i++)
				data[i] += other.data[i];
			return *this;
		}

		constexpr KITTEN_FUNC_DECL SymMat<N, T>& operator-=(const SymMat<N, T>& other) {
			for (int i = 0; i < DATA_LEN; i++)
				data[i] -= other.data[i];
			return *this;
		}

		constexpr KITTEN_FUNC_DECL SymMat<N, T>& operator*=(const T scalar) {
			for (int i = 0; i < DATA_LEN; i++)
				data[i] *= scalar;
			return *this;
		}

		constexpr KITTEN_FUNC_DECL SymMat<N, T>& operator/=(const T scalar) {
			T invScalar = 1 / scalar;
			for (int i = 0; i < DATA_LEN; i++)
				data[i] *= invScalar;
			return *this;
		}

		constexpr KITTEN_FUNC_DECL SymMat<N, T>& operator*=(const SymMat<N, T>& other) {
			*this = *this * other;
			return *this;
		}

		constexpr KITTEN_FUNC_DECL T& operator[](const int i) {
			return data[i];
		}

		constexpr KITTEN_FUNC_DECL const T& operator[](const int i) const {
			return data[i];
		}

		constexpr KITTEN_FUNC_DECL static const SymMat<N, T> outer(const vec_type v) {
			SymMat<N, T> result(0);
			for (int i = 0; i < N; i++)
				result[i] = v[i] * v[i];
			for (int i = 1; i < N; i++)
				for (int j = 0; j < i; j++)
					result[N + j + (i * (i - 1)) / 2] = v[i] * v[j];
			return result;
		}
	};

	typedef SymMat<2, float> symmat2;
	typedef SymMat<3, float> symmat3;
	typedef SymMat<4, float> symmat4;

	typedef SymMat<2, double> symdmat2;
	typedef SymMat<3, double> symdmat3;
	typedef SymMat<4, double> symdmat4;
}