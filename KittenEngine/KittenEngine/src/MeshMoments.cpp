#include "../includes/modules/Mesh.h"

namespace Kitten {
	float Mesh::zerothMoment() {
		double vol = 0;

#pragma omp parallel
		{
			double lvol = 0;
#pragma omp for schedule(static, 2048)
			for (int i = 0; i < indices.size(); i += 3)
				lvol += determinant(mat3(vertices[indices[i + 0]].pos, vertices[indices[i + 1]].pos, vertices[indices[i + 2]].pos));
#pragma omp critical
			vol += lvol;
		}

		return (float)(vol / 6);
	}

	vec3 Mesh::firstMoment() {
		dvec3 f(0);

#pragma omp parallel
		{
			dvec3 lf(0);
#pragma omp for schedule(static, 2048)
			for (int i = 0; i < indices.size(); i += 3) {
				mat3 tet(vertices[indices[i + 0]].pos, vertices[indices[i + 1]].pos, vertices[indices[i + 2]].pos);
				lf += determinant(tet) * (tet * vec3(1));
			}
#pragma omp critical
			f += lf;
		}

		return vec3((1. / 24.) * f);
	}

	mat3 Mesh::secondMoment() {
		dmat3 m(0);

#pragma omp parallel
		{
			dmat3 lm(0);
#pragma omp for schedule(static, 2048)
			for (int i = 0; i < indices.size(); i += 3) {
				mat3 tet(vertices[indices[i + 0]].pos, vertices[indices[i + 1]].pos, vertices[indices[i + 2]].pos);
				mat3 a;

				for (int k = 0; k < 3; k++)
					a[k][k] = (pow2(tet[0][k]) + pow2(tet[1][k]) + pow2(tet[2][k]) +
						tet[0][k] * tet[1][k] + tet[0][k] * tet[2][k] + tet[1][k] * tet[2][k]) / 60;

				tet = transpose(tet);
				a[0][1] = a[1][0] = dot(tet[0], tet[1] + dot(tet[1], vec3(1))) / 120;
				a[0][2] = a[2][0] = dot(tet[0], tet[2] + dot(tet[2], vec3(1))) / 120;
				a[1][2] = a[2][1] = dot(tet[1], tet[2] + dot(tet[2], vec3(1))) / 120;

				lm += determinant(tet) * a;
			}
#pragma omp critical
			m += lm;
		}

		return mat3(m);
	}

	vec3 Mesh::centerOfMass() {
		return firstMoment() / zerothMoment();
	}
}