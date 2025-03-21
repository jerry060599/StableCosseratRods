#include "../Cosserat.h"
#include "KittenEngine/includes/modules/SymMat.h"
#include <Eigen/Eigen>

namespace Cosserat {
	void Sim::iterateVertVBD() {
		const float h = meta.h;
		const float damping = meta.damping / h;
		const bool proportionalDamping = meta.useProportionalDamping;

#pragma omp parallel for if(multithreaded) schedule(dynamic, COSSERAT_THREAD_BATCH)
		for (int tid = 0; tid < verts.size(); tid++) {
			// Linear change
			auto v0 = verts[tid];

			if (v0.invMass > 0) {
				// Hessian H but we only track the diagonal
				float H = 1 / (v0.invMass * h * h);
				// vel has been overwritten to contain y - pos
				vec3 dx = dxs[tid];
				vec3 f = (1 / (h * h * v0.invMass)) * (vels[tid] - dx);

				// Cosserat stretching energy
				const int start = vertSegIndices[tid];
				const int end = vertSegIndices[tid + 1];
				for (int s = start; s < end; s++) {
					int segId = vertSegIds[s];
					auto seg = segs[segId];

					int oid = (seg.i.x == tid) ? seg.i.y : seg.i.x;
					auto p1 = verts[oid].pos;
					auto p1dx = dxs[oid];

					float invl = 1 / seg.l;
					float order = (seg.i.x == tid) ? 1 : -1;
					vec3 c = ((p1 - v0.pos) + (p1dx - dx)) * (order * invl) - seg.q * vec3(1, 0, 0);

					float k = abs(seg.kss) * invl;
					float d = k * invl;
					float damp = damping * (proportionalDamping ? d : invl);
					f += order * k * c - damp * dx;
					H += d + damp;
				}

				// Local solve
				vec3 delta = f * (1 / H);
				dx += delta;

				dxs[tid] = dx;
			}
		}
	}

	// Return A s.t. A a = q a for some quaternion q 
	inline mat4 leftQMulMatrix(Rotor q) {
		return mat4(
			q.w, q.v.z, -q.v.y, -q.v.x,
			-q.v.z, q.w, q.v.x, -q.v.y,
			q.v.y, -q.v.x, q.w, -q.v.z,
			q.v.x, q.v.y, q.v.z, q.v.w
		);
	}

	// Return A s.t. A a = a q for some quaternion q 
	inline mat4 rightQMulMatrix(Rotor q) {
		return mat4(
			q.w, -q.v.z, q.v.y, -q.v.x,
			q.v.z, q.w, -q.v.x, -q.v.y,
			-q.v.y, q.v.x, q.w, -q.v.z,
			q.v.x, q.v.y, q.v.z, q.v.w
		);
	}

	template <int n>
	glm::vec<n, double, glm::defaultp> eigenValuesSolve(
		const glm::mat<n, n, double, glm::defaultp>& M,
		glm::mat<n, n, double, glm::defaultp>& eigenVectors) {
		// Map the glm matrix to an Eigen matrix
		Eigen::Map<const Eigen::Matrix<double, n, n, Eigen::RowMajor>> eigenMatrix(&M[0][0]);

		// Compute the eigenvalues and eigenvectors using Eigen
		Eigen::EigenSolver<Eigen::Matrix<double, n, n>> solver(eigenMatrix);
		if (solver.info() != Eigen::Success)
			throw std::runtime_error("Eigen decomposition failed.");

		// Extract the eigenvalues (real and imaginary parts)
		Eigen::Matrix<double, n, 1> realPart = solver.eigenvalues().real();
		Eigen::Matrix<double, n, 1> imagPart = solver.eigenvalues().imag();

		// Extract the eigenvectors
		Eigen::Matrix<double, n, n> eigenVecs = solver.eigenvectors().real();

		// Convert the eigenvectors from Eigen matrix to glm matrix
		for (int row = 0; row < n; ++row)
			for (int col = 0; col < n; ++col)
				eigenVectors[col][row] = eigenVecs(row, col); // Note glm is column-major

		// Convert the eigenvalues from Eigen vector to glm vector
		glm::vec<n, double, glm::defaultp> glmEigenValues;
		for (int i = 0; i < n; ++i) {
			glmEigenValues[i] = (imagPart[i] == 0.0) ? realPart[i] : std::numeric_limits<double>::quiet_NaN();
		}

		return glmEigenValues;
	}


	void Sim::iterateSegPSDPLCVBD() {
		// VBD but with no inertia in the local variational energy
		// Essentially local newton iterations modified to use the least norm solution
#pragma omp parallel for if(multithreaded) schedule(dynamic, COSSERAT_THREAD_BATCH)
		for (int tid = 0; tid < segs.size(); tid++) {
			Seg seg = segs[tid];

			if (seg.kss < 0) continue;

			const vec3 v0 = verts[seg.i.x].pos;
			const vec3 dx = dxs[seg.i.x];
			const vec3 p1 = verts[seg.i.y].pos;
			const vec3 p1dx = dxs[seg.i.y];

			// Torque and hessian for the stretching constraint
			vec3 cs = ((p1 - v0) + (p1dx - dx)) / seg.l - seg.q * vec3(1, 0, 0);
			vec4 r = (Kit::Rotor(2 * seg.kss * cs) * seg.q * Kit::Rotor(1)).v;

			mat4x3 H43 = rightQMulMatrix(Kit::Rotor(1) * seg.q.inverse());
			mat4 H = (4.f * seg.kss) * (transpose(H43) * H43)
				+ leftQMulMatrix(Kit::Rotor(2 * seg.kss * cs)) * rightQMulMatrix(Kit::Rotor(1));

			// Torque and hessian for the bending constraints
			const int start = segAngleIndices[tid];
			const int end = segAngleIndices[tid + 1];
			for (int s = start; s < end; s++) {
				int angId = segAngleIds[s];
				RestAngle ang = restAngles[angId];

				const int oid = (ang.i.x == tid) ? ang.i.y : ang.i.x;
				const auto oq = segs[oid].q;
				auto qq = oq.inverse() * seg.q;
				if (ang.i.x == tid) qq = qq.inverse();

				float phi = (dot(qq.v, ang.qRest.v) > 0) ? 1 : -1;

				if (ang.i.x == tid) ang.qRest = ang.qRest.inverse();
				r -= (phi * ang.kbt) * (oq * Kit::Rotor(ang.qRest)).v;
				H[0][0] += ang.kbt; H[1][1] += ang.kbt; H[2][2] += ang.kbt; H[3][3] += ang.kbt;
			}

			if (start == end)
				seg.q = Rotor::fromTo(vec3(1, 0, 0), normalize((p1 - v0) + (p1dx - dx)));
			else {
				dmat4 eigenVectors;
				dvec4 eigenValues = eigenValuesSolve((dmat4)H, eigenVectors);
				eigenValues = glm::max(eigenValues, dvec4(1e2));	// We use 1e2 here because any less will cause numerical instability
				H = eigenVectors * diag(eigenValues) * transpose(eigenVectors);

				r -= dot(seg.q.v, r) * seg.q.v;
				// This needs to be computed as double. Otherwise it will be unstable
				auto AtAInv = (mat4)inverse((dmat4)transpose(H) * (dmat4)H);
				float lambda = dot(seg.q.v, AtAInv * (-r * H)) / dot(seg.q.v, AtAInv * seg.q.v);
				vec4 d = AtAInv * (-r * H - seg.q.v * lambda);
				auto newQ = normalize(Kit::Rotor(seg.q.v + d));

				if (!glm::any(isnan(newQ.v))) seg.q = newQ;	// Last line of defense against numerical instabilities.
			}
			segs[tid].q = seg.q;
		}
	}

	void Sim::iterateSegLCVBD() {
		// VBD but with no inertia in the local variational energy
		// Essentially local newton iterations modified to use the least norm solution
#pragma omp parallel for if(multithreaded) schedule(dynamic, COSSERAT_THREAD_BATCH)
		for (int tid = 0; tid < segs.size(); tid++) {
			Seg seg = segs[tid];

			if (seg.kss < 0) continue;

			const vec3 v0 = verts[seg.i.x].pos;
			const vec3 dx = dxs[seg.i.x];
			const vec3 p1 = verts[seg.i.y].pos;
			const vec3 p1dx = dxs[seg.i.y];

			// Torque and hessian for the stretching constraint
			vec3 cs = ((p1 - v0) + (p1dx - dx)) / seg.l - seg.q * vec3(1, 0, 0);
			vec4 r = (Kit::Rotor(2 * seg.kss * cs) * seg.q * Kit::Rotor(1)).v;

			mat4x3 H43 = rightQMulMatrix(Kit::Rotor(1) * seg.q.inverse());
			mat4 H = (4.f * seg.kss) * (transpose(H43) * H43)
				+ leftQMulMatrix(Kit::Rotor(2 * seg.kss * cs)) * rightQMulMatrix(Kit::Rotor(1));

			// Torque and hessian for the bending constraints
			const int start = segAngleIndices[tid];
			const int end = segAngleIndices[tid + 1];
			for (int s = start; s < end; s++) {
				int angId = segAngleIds[s];
				RestAngle ang = restAngles[angId];

				const int oid = (ang.i.x == tid) ? ang.i.y : ang.i.x;
				const auto oq = segs[oid].q;
				auto qq = oq.inverse() * seg.q;
				if (ang.i.x == tid) qq = qq.inverse();

				float phi = (dot(qq.v, ang.qRest.v) > 0) ? 1 : -1;

				if (ang.i.x == tid) ang.qRest = ang.qRest.inverse();
				r -= (phi * ang.kbt) * (oq * Kit::Rotor(ang.qRest)).v;
				H[0][0] += ang.kbt; H[1][1] += ang.kbt; H[2][2] += ang.kbt; H[3][3] += ang.kbt;
			}

			if (start == end)
				seg.q = Rotor::fromTo(vec3(1, 0, 0), normalize((p1 - v0) + (p1dx - dx)));
			else {
				r -= dot(seg.q.v, r) * seg.q.v;
				// This needs to be computed as double. Otherwise it will be unstable
				auto AtAInv = (mat4)inverse((dmat4)transpose(H) * (dmat4)H);
				float lambda = dot(seg.q.v, AtAInv * (-r * H)) / dot(seg.q.v, AtAInv * seg.q.v);
				vec4 d = AtAInv * (-r * H - seg.q.v * lambda);
				auto newQ = normalize(Kit::Rotor(seg.q.v + d));

				if (!glm::any(isnan(newQ.v))) seg.q = newQ;	// Last line of defense against numerical instabilities.
			}
			segs[tid].q = seg.q;
		}
	}

	void Sim::iterateSegVBD() {
		// VBD but with no inertia in the local variational energy
		// Essentially local newton iterations
#pragma omp parallel for if(multithreaded) schedule(dynamic, COSSERAT_THREAD_BATCH)
		for (int tid = 0; tid < segs.size(); tid++) {
			Seg seg = segs[tid];

			if (seg.kss < 0) continue;

			const vec3 v0 = verts[seg.i.x].pos;
			const vec3 dx = dxs[seg.i.x];
			const vec3 p1 = verts[seg.i.y].pos;
			const vec3 p1dx = dxs[seg.i.y];

			// Torque and hessian for the stretching constraint
			vec3 cs = ((p1 - v0) + (p1dx - dx)) / seg.l - seg.q * vec3(1, 0, 0);
			vec4 r = (Kit::Rotor(2 * seg.kss * cs) * seg.q * Kit::Rotor(1)).v;

			mat4x3 H43 = rightQMulMatrix(Kit::Rotor(1) * seg.q.inverse());
			mat4 H = (4.f * seg.kss) * (transpose(H43) * H43)
				+ leftQMulMatrix(Kit::Rotor(2 * seg.kss * cs)) * rightQMulMatrix(Kit::Rotor(1));

			// Torque and hessian for the bending constraints
			const int start = segAngleIndices[tid];
			const int end = segAngleIndices[tid + 1];
			for (int s = start; s < end; s++) {
				int angId = segAngleIds[s];
				RestAngle ang = restAngles[angId];

				const int oid = (ang.i.x == tid) ? ang.i.y : ang.i.x;
				const auto oq = segs[oid].q;
				auto qq = oq.inverse() * seg.q;
				if (ang.i.x == tid) qq = qq.inverse();

				float phi = (dot(qq.v, ang.qRest.v) > 0) ? 1 : -1;

				if (ang.i.x == tid) ang.qRest = ang.qRest.inverse();
				r -= (phi * ang.kbt) * (oq * Kit::Rotor(ang.qRest)).v;
				H[0][0] += ang.kbt; H[1][1] += ang.kbt; H[2][2] += ang.kbt; H[3][3] += ang.kbt;
			}

			if (start == end)
				seg.q = Rotor::fromTo(vec3(1, 0, 0), normalize((p1 - v0) + (p1dx - dx)));
			else {
				vec4 d = -(inverse(H) * r);

				auto newQ = normalize(Kit::Rotor(seg.q.v + d));
				if (!glm::any(isnan(newQ.v))) seg.q = newQ;	// Last line of defense against numerical instabilities.
			}
			segs[tid].q = seg.q;
		}
	}
}