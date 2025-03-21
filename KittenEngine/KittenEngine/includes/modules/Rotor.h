#pragma once

#include "Common.h"

namespace Kitten {
	template <typename T>
	struct RotorX {
		typedef glm::vec<3, T, glm::defaultp> q_type;
		typedef glm::vec<4, T, glm::defaultp> v_type;

		union {
			v_type v;
			struct {
				q_type q;	// The bivector part laid out in the { y^z, z^x, x^y } basis
				T w;		// The scaler part
			};
			struct {
				T x, y, z, s;
			};
		};

		// Create a quaternion from an angle (in radians) and rotation axis
		KITTEN_FUNC_DECL static RotorX<T> angleAxis(T rad, q_type axis) {
			rad *= 0.5;
			return RotorX<T>(sin(rad) * axis, cos(rad));
		}

		// Create a quaternion from an angle (in degrees) and rotation axis
		KITTEN_FUNC_DECL static RotorX<T> angleAxisDeg(T deg, q_type axis) {
			deg *= 0.00872664625997165;
			return RotorX<T>(sin(deg) * axis, cos(deg));
		}

		// Create a quaternion from euler angles in radians
		KITTEN_FUNC_DECL static RotorX<T> eulerAngles(q_type rad) {
			rad *= 0.5;
			q_type c = cos(rad);
			q_type s = sin(rad);
			return RotorX<T>(q_type(0, 0, s.z), c.z) * RotorX<T>(q_type(0, s.y, 0), c.y) * RotorX<T>(q_type(s.x, 0, 0), c.x);
		}

		// Create a quaternion from euler angles in radians
		KITTEN_FUNC_DECL static RotorX<T> eulerAngles(T x, T y, T z) {
			return eulerAngles(q_type(x, y, z));
		}

		// Create a quaternion from euler angles in degrees
		KITTEN_FUNC_DECL static RotorX<T> eulerAnglesDeg(q_type deg) {
			return eulerAngles(deg * 0.0174532925199432958f);
		}

		// Create a quaternion from euler angles in degrees
		KITTEN_FUNC_DECL static RotorX<T> eulerAnglesDeg(T x, T y, T z) {
			return eulerAnglesDeg(q_type(x, y, z));
		}

		KITTEN_FUNC_DECL static RotorX<T> fromTo(q_type from, q_type to) {
			q_type h = (from + to) * 0.5f;
			T l = length2(h);
			if (l > 0) h *= inversesqrt(l);
			else h = orthoBasisX((vec3)from)[1];

			return RotorX<T>(cross(from, h), dot(from, h));
		}

		// Returns the multiplicative identity rotor
		KITTEN_FUNC_DECL static RotorX<T> identity() {
			return RotorX<T>();
		}

		KITTEN_FUNC_DECL RotorX(T x, T y = 0, T z = 0, T w = 0) : v(x, y, z, w) {}
		KITTEN_FUNC_DECL RotorX(q_type q, T w = 0) : q(q), w(w) {}
		KITTEN_FUNC_DECL RotorX(v_type v) : v(v) {}
		KITTEN_FUNC_DECL RotorX() : v(0, 0, 0, 1) {}

		KITTEN_FUNC_DECL RotorX(const RotorX<T>& other) : v(other.v) {}

		template<typename U>
		KITTEN_FUNC_DECL explicit RotorX(const RotorX<U>& other) : v((v_type)other.v) {}

		KITTEN_FUNC_DECL RotorX<T>& operator=(const RotorX<T>& rhs) {
			v = rhs.v;
			return *this;
		}

		// Get the multiplicative inverse
		KITTEN_FUNC_DECL RotorX<T> inverse() const {
			return RotorX<T>(-q, w);
		}

		KITTEN_FUNC_DECL RotorX<T> operator- () const {
			return inverse();
		}

		// Rotate a vector by this rotor
		KITTEN_FUNC_DECL q_type rotate(q_type v) const {
			// Calculate v * ab
			q_type a = w * v + cross(q, v);	// The vector
			T c = dot(v, q);			// The trivector

			// Calculate (w - q) * (a + c). Ignoring the scaler-trivector parts
			return w * a		// The scaler-vector product
				+ cross(q, a)	// The bivector-vector product
				+ c * q;		// The bivector-trivector product
		}

		KITTEN_FUNC_DECL mat<3, 3, T, defaultp> matrix() {
			mat3 cm = crossMatrix(q);
			return abT(q, q) + mat<3, 3, T, defaultp>(w * w) + 2 * w * cm + cm * cm;
		}

		KITTEN_FUNC_DECL static RotorX<T> fromMatrix(mat<3, 3, T, defaultp> m) {
			RotorX<T> q;
			T t;
			if (m[2][2] < 0) {
				if (m[0][0] > m[1][1]) {
					t = 1 + m[0][0] - m[1][1] - m[2][2];
					q = RotorX<T>(t, m[1][0] + m[0][1], m[0][2] + m[2][0], m[2][1] - m[1][2]);
				}
				else {
					t = 1 - m[0][0] + m[1][1] - m[2][2];
					q = RotorX<T>(m[1][0] + m[0][1], t, m[2][1] + m[1][2], m[0][2] - m[2][0]);
				}
			}
			else {
				if (m[0][0] < -m[1][1]) {
					t = 1 - m[0][0] - m[1][1] + m[2][2];
					q = RotorX<T>(m[0][2] + m[2][0], m[2][1] + m[1][2], t, m[1][0] - m[0][1]);
				}
				else {
					t = 1 + m[0][0] + m[1][1] + m[2][2];
					q = RotorX<T>(m[2][1] - m[1][2], m[0][2] - m[2][0], m[1][0] - m[0][1], t);
				}
			}

			return RotorX<T>(((0.5f / glm::sqrt(t)) * q.v)).inverse();
		}

		// Get the euler angle in radians
		KITTEN_FUNC_DECL q_type euler() {
			return q_type(
				atan2(2 * (w * q.x + q.y * q.z), 1 - 2 * (q.x * q.x + q.y * q.y)),
				asin(2 * (w * q.y - q.x * q.z)),
				atan2(2 * (w * q.z + q.x * q.y), 1 - 2 * (q.y * q.y + q.z * q.z))
			);
		}

		// Get the euler angle in degrees
		KITTEN_FUNC_DECL q_type eulerDeg() {
			return euler() * (T)57.29577951308232;
		}

		// Returns both the axis and rotation angle in radians
		KITTEN_FUNC_DECL q_type axis(T& angle) {
			T l = length(q);
			if (l == 0) {
				angle = 0;
				return q_type(1, 0, 0);
			}

			angle = 2 * atan2(l, w);
			return q / l;
		}

		// Returns the axis of rotation
		KITTEN_FUNC_DECL q_type axis() { T a; return axis(a); }

		// Returns both the axis and rotation angle in degrees
		KITTEN_FUNC_DECL q_type axisDeg(T& angle) {
			q_type a = axis(angle);
			angle *= 57.29577951308232f;
			return a;
		}

		// Returns the angle of rotation in radians
		KITTEN_FUNC_DECL T angle() { T a; axis(a); return a; }

		// Returns the angle of rotation in degrees
		KITTEN_FUNC_DECL T angleDeg() { T a; axis(a); return a * 57.29577951308232; }

		KITTEN_FUNC_DECL friend q_type operator*(RotorX<T> lhs, const q_type& rhs) {
			return lhs.rotate(rhs);
		}

		KITTEN_FUNC_DECL friend RotorX<T> operator*(RotorX<T> lhs, const RotorX<T>& rhs) {
			return RotorX<T>(lhs.w * rhs.q + rhs.w * lhs.q + cross(lhs.q, rhs.q),
				lhs.w * rhs.w - dot(lhs.q, rhs.q));
		}

		KITTEN_FUNC_DECL RotorX<T>& operator+=(const RotorX<T>& rhs) {
			v += rhs.v;
			return *this;
		}

		KITTEN_FUNC_DECL friend RotorX<T> operator+(RotorX<T> lhs, const RotorX<T>& rhs) {
			return RotorX<T>(lhs.v + rhs.v);
		}

		KITTEN_FUNC_DECL friend RotorX<T> operator*(T lhs, const RotorX<T>& rhs) {
			if (rhs.w == 1) return RotorX<T>();
			T na = lhs * acos(rhs.w);	// New angle
			T nw = cos(na);				// New cosine
			T s = sqrt((1 - nw * nw) / length2(rhs.q));
			if (fract(na * 0.1591549430918954) > 0.5) s = -s;
			return RotorX<T>(rhs.q * s, nw);
		}

		// Gets the vec4 repersentation laid out in { y^z, z^x, x^y, scaler }
		KITTEN_FUNC_DECL explicit operator v_type() const {
			return v;
		}

		KITTEN_FUNC_DECL T& operator[](std::size_t idx) {
			return v[idx];
		}
	};

	// mix rotors a to b from t=[0, 1] (unclamped)
	template<typename T>
	KITTEN_FUNC_DECL inline RotorX<T> mix(RotorX<T> a, RotorX<T> b, T t) {
		return (t * (b * a.inverse())) * a;
	}

	template<typename T>
	KITTEN_FUNC_DECL inline T dot(RotorX<T> a, RotorX<T> b) {
		return glm::dot(a.v, b.v);
	}

	template<typename T>
	KITTEN_FUNC_DECL inline RotorX<T> normalize(RotorX<T> a) {
		return RotorX<T>(glm::normalize(a.v));
	}

	// Projects x onto the constraint q*e = d
	template<typename T>
	KITTEN_FUNC_DECL inline RotorX<T> projectRotor(RotorX<T> x, vec<3, T, glm::defaultp> e, vec<3, T, glm::defaultp> d) {
		auto q = RotorX<T>::fromTo(e, d);
		auto qp = q.inverse() * x;
		qp.y = qp.z = 0;
		return q * normalize(qp);
	}

	template <typename T>
	KITTEN_FUNC_DECL void print(RotorX<T> v, const char* format = "%.4f") {
		printf("{");
		for (int i = 0; i < 4; i++) {
			printf(format, v[i]);
			if (i != 3) printf(", ");
		}
		if (abs(length2(v.v) - 1) < 1e-3) {
			printf("}, euler: {");
			auto a = v.eulerDeg();
			for (int i = 0; i < 3; i++) {
				printf(format, a[i]);
				if (i != 2) printf(", ");
			}
			printf("} deg\n");
		}
		else
			printf("}\n");
	}

	using Rotor = RotorX<float>;
	using RotorD = RotorX<double>;
}