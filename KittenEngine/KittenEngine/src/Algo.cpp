#include "../includes/modules/Algo.h"

using namespace Eigen;

std::vector<glm::vec3> Kitten::bluenoiseSample(vector<glm::vec3>& samples, int N) {
	using namespace glm;
	using namespace std;
	vector<vec3> newSamples;
	newSamples.reserve(N);
	newSamples.push_back(samples.back());
	samples.pop_back();

	vector<float> closestDist(samples.size(), std::numeric_limits<float>::infinity());

	for (int i = 1; i < N; i++) {
		float maxDist = 0;
		int furthestIdx = -1;

		vec3 lastPicked = newSamples.back();
		for (int j = 0; j < samples.size(); j++) {
			vec3 p = samples[j];

			float nd = length2(p - lastPicked);
			float minDist = closestDist[j];
			if (nd < minDist)
				closestDist[j] = minDist = nd;

			if (minDist > maxDist) {
				maxDist = minDist;
				furthestIdx = j;
			}
		}

		newSamples.push_back(samples[furthestIdx]);
		std::swap(samples[furthestIdx], samples.back());
		samples.pop_back();
		std::swap(closestDist[furthestIdx], closestDist.back());
		closestDist.pop_back();
	}
	return newSamples;
}

std::vector<float> Kitten::polylineUniformSample(std::function<vec3(float)> f, float a, float b, const int numSamples, const int numItr, const float learningRate) {
	std::vector<float> ts(numSamples);
	for (size_t i = 0; i < numSamples; i++)
		ts[i] = glm::mix(a, b, i / (float)(numSamples - 1));

	vec3 root = f(a);
	for (size_t k = 0; k < numItr; k++) {
		float lt = a;
		vec3 left = root;
		float mt = ts[1];
		vec3 mid = f(mt);

		for (size_t i = 1; i < numSamples - 1; i++) {
			vec3 right = f(ts[i + 1]);
			float rt = ts[i + 1];

			// (t - lt) / (mt - lt) * leftL = (1 - (t - mt) / (rt - mt)) * rightL
			// (t - lt) lr = (rt - t) rr
			float lr = length(left - mid) / (mt - lt);
			float rr = length(right - mid) / (rt - mt);

			float episilon = (rt - lt) * 1e-2f;
			float nt = lr + rr != 0 ? glm::clamp((lt * lr + rt * rr) / (lr + rr), lt + episilon, rt - episilon) : 0.5f * (lt + rt);
			ts[i] = mt = glm::mix(mt, nt, learningRate);
			mid = f(mt);

			left = mid;
			mid = right;
			lt = mt;
			mt = rt;
		}
	}

	return ts;
}

Eigen::VectorXf Kitten::lbfgsMin(int numVars, std::function<float(Eigen::VectorXf)> f,
	std::function<Eigen::VectorXf(Eigen::VectorXf)> g,
	Eigen::VectorXf& guess, const float tol, const int m) {
	using namespace Eigen;

	MatrixXf s(numVars, m);
	MatrixXf y(numVars, m);
	VectorXf rho(m);
	VectorXf a(m);

	int buffSize = 0;
	int buffInd = 0;

	VectorXf lastX = guess;
	VectorXf lastGv = g(guess);
	VectorXf x = guess - 0.001f * lastGv;
	float fv = f(x);
	VectorXf gv;
	VectorXf z;
	int itr = 0;
	while (true) {
		itr++;

		const int ind = buffInd % m;

		// Find descent direction with lbfgs

		gv = g(x);
		float k = (gv - lastGv).dot(x - lastX);
		if (k <= 0) {
			x -= 0.001f * gv;
			if (gv.norm() < tol) break;
			continue;
		}

		s.col(ind) = x - lastX;
		y.col(ind) = gv - lastGv;
		rho[ind] = 1 / k;

		lastX = x;
		lastGv = gv;
		z = gv;
		buffSize = std::min(m, buffSize + 1);

		for (int i = 0; i < buffSize; i++) {
			const int oi = (buffInd + m - i) % m;
			a[oi] = rho[oi] * z.dot(s.col(oi));
			z -= a[oi] * y.col(oi);
		}

		float Y = y.col(ind).dot(y.col(ind));
		if (Y != 0) {
			Y = s.col(ind).dot(y.col(ind)) / Y;
			z *= Y;
		}
		for (int i = 0; i < buffSize; i++) {
			const int oi = (buffInd + m + i - buffSize + 1) % m;
			z += s.col(oi) * (a[oi] - rho[oi] * y.col(oi).dot(z));
		}

		z *= -1;
		// Perform backtracking line search
		float t = 0.5f * gv.dot(z);
		float a = 1.f;
		float v = f(x);
		while (f(x + a * z) > v + a * t)
			a *= 0.5f;

		// Update guess
		x += a * z;
		float newFv = f(x);
		float diff = newFv - fv;
		fv = newFv;

		buffInd++;

		if (gv.norm() < tol && abs(diff) < tol) break;
	}
	// printf("ITR: %d\n", itr);
	return x;
}

Eigen::VectorXd Kitten::cg(
	Eigen::SparseMatrix<double, Eigen::RowMajor>& A,
	Eigen::VectorXd& b,
	const double tol,
	const int itrLim) {

	// Initialize preconditioner
	VectorXd invDiag(A.rows());
#pragma omp parallel for schedule(static, 512)
	for (int i = 0; i < invDiag.size(); i++) {
		double v = A.coeff(i, i);
		invDiag[i] = abs(v) < 1e-10 ? 1 : 1 / v;
	}

	VectorXd x = VectorXd::Zero(invDiag.size());
	VectorXd r = b;
	VectorXd d = invDiag.array() * r.array();

	VectorXd q(x.size());
	VectorXd s(x.size());

	const int UPDATE_ITR = std::max(100, (int)sqrt(A.cols()));
	double rDotD[2] = { r.dot(d), 0 };
	const double relTol = tol * tol * rDotD[0];
	size_t itr = 1;
	for (; (itrLim < 0 || itr <= itrLim) && rDotD[0] > relTol; itr++) {
		q = A * d;
		double alpha = rDotD[0] / d.dot(q);
		x += alpha * d;

		if (itr % UPDATE_ITR == 0)
			r = b - A * x;
		else
			r -= alpha * q;

		s = invDiag.array() * r.array();
		rDotD[1] = rDotD[0];
		rDotD[0] = r.dot(s);
		double beta = rDotD[0] / rDotD[1];
		d = s + beta * d;
	}
	// printf("%zd\n", itr);

	return x;
}

Eigen::VectorXd Kitten::bccg(
	Eigen::SparseMatrix<double, Eigen::RowMajor>& A,
	Eigen::VectorXd& b,
	Eigen::VectorXd& lower,
	const double tol,
	const int itrLim) {

	// Initialize preconditioner
	VectorXd invDiag(A.rows());
	//VectorXd x(invDiag.size());
	VectorXd x;

	{
		ConjugateGradient<SparseMatrix<double, RowMajor>, Lower | Upper, DiagonalPreconditioner<double>> cg;
		cg.setTolerance(512 * tol);
		cg.compute(A);
		x = cg.solve(b);
	}

#pragma omp parallel for schedule(static, 512)
	for (int i = 0; i < x.size(); i++) {
		double v = A.coeff(i, i);
		invDiag[i] = abs(v) < 1e-10 ? 1 : 1 / v;
		x[i] = std::max(x[i], lower[i]);
	}

	VectorXd r_tilde = b - A * x;

	// Initialize bounded set
	bool* boundSet = new bool[x.size()];
#pragma omp parallel for schedule(static, 1024)
	for (int i = 0; i < x.size(); i++)
		boundSet[i] = x[i] <= lower[i] && r_tilde[i] < 0;

	// Initialize the rest of CG
	VectorXd d = invDiag.array() * r_tilde.array();
	VectorXd r(x.size());
#pragma omp parallel for schedule(static, 1024)
	for (int i = 0; i < r.size(); i++)
		r[i] = boundSet[i] ? 0 : r_tilde[i];

	VectorXd q(x.size());
	VectorXd s(x.size());

	double rDotD[2] = { r.dot(d), 0 };
	const double relTol = tol * tol * rDotD[0];
	bool projected = false;
	bool haveUnreleased = false;
	bool boundsChanged = false;
	int itrSinceRes = 0;
	const int UPDATE_ITR = std::max(100, (int)sqrt(A.cols()));

	size_t itr = 1;
	for (; (itrLim < 0 || itr <= itrLim) && (rDotD[0] > relTol || boundsChanged || haveUnreleased); itr++, itrSinceRes++) {
		q = A * d;
		double alpha = rDotD[0] / d.dot(q);
		x += alpha * d;

		if (projected || itrSinceRes >= UPDATE_ITR) {
			r_tilde = b - A * x;
			itrSinceRes = 0;
		}
		else
			r_tilde -= alpha * q;

		// Update bounded set and projected
		boundsChanged = projected = haveUnreleased = false;
#pragma omp parallel for schedule(static, 1024)
		for (int i = 0; i < x.size(); i++) {
			bool bounded = x[i] <= lower[i] && r_tilde[i] < 0;

			if (itr % 64 == 1) {
				// Use the full bound update method
				if (boundSet[i] != bounded) {
					boundSet[i] = bounded;
					boundsChanged = true;
				}
			}
			else if (bounded && !boundSet[i]) {
				// We only want to bound things and not release too often to prevent slow convergence due to oscillations
				boundSet[i] = true;
				boundsChanged = true;
			}
			else if (boundSet[i] != bounded)
				haveUnreleased = true; // We dont want to exit before all the unreleased stuff is released

			if (x[i] < lower[i]) {
				x[i] = lower[i];
				projected = true;
			}

			r[i] = boundSet[i] ? 0 : r_tilde[i];
		}
		rDotD[1] = rDotD[0];

		if (boundsChanged) {
			d = invDiag.array() * r.array();
			rDotD[0] = r.dot(d);
		}
		else {
			s = invDiag.array() * r.array();
			rDotD[0] = r.dot(s);
			d = s + (rDotD[0] / rDotD[1]) * d;
		}
	}
	// printf("%zd\n", itr);
	delete[] boundSet;
	return x;
}

Eigen::VectorXd Kitten::rbccg(Eigen::SparseMatrix<double,
	Eigen::RowMajor>& A,
	Eigen::VectorXd& b,
	Eigen::VectorXd& lower,
	Eigen::VectorXd& shift,
	const double regAlpha,
	const double tol,
	const int itrLim) {
	using namespace Eigen;
	// Initialize preconditioner
	VectorXd invDiag(A.rows());
	VectorXd x(invDiag.size());
	VectorXd sb = regAlpha * shift + b;

#pragma omp parallel for schedule(static, 512)
	for (int i = 0; i < x.size(); i++) {
		double v = A.coeff(i, i) + regAlpha;
		invDiag[i] = abs(v) < 1e-10 ? 1 : 1 / v;
		x[i] = std::max(0., lower[i]);
	}

	VectorXd r_tilde = sb - A * x - regAlpha * x;

	// Initialize bounded set
	bool* boundSet = new bool[x.size()];
#pragma omp parallel for schedule(static, 1024)
	for (int i = 0; i < x.size(); i++)
		boundSet[i] = x[i] <= lower[i] && r_tilde[i] < 0;

	// Initialize the rest of CG
	VectorXd d = invDiag.array() * r_tilde.array();
	VectorXd r(x.size());
#pragma omp parallel for schedule(static, 1024)
	for (int i = 0; i < r.size(); i++)
		r[i] = boundSet[i] ? 0 : r_tilde[i];

	VectorXd q(x.size());
	VectorXd s(x.size());

	double rDotD[2] = { r.dot(d), 0 };
	const double relTol = tol * tol * rDotD[0];
	bool projected = false;
	bool haveUnreleased = false;
	bool boundsChanged = false;
	int itrSinceRes = 0;
	const int UPDATE_ITR = std::max(100, (int)sqrt(A.cols()));

	size_t itr = 1;
	for (; (itrLim < 0 || itr <= itrLim) && (rDotD[0] > relTol || boundsChanged || haveUnreleased); itr++, itrSinceRes++) {
		q = A * d + regAlpha * d;
		double alpha = rDotD[0] / d.dot(q);
		x += alpha * d;

		if (projected || itrSinceRes >= UPDATE_ITR) {
			r_tilde = sb - A * x - regAlpha * x;
			itrSinceRes = 0;
		}
		else
			r_tilde -= alpha * q;

		// Update bounded set and projected
		boundsChanged = projected = haveUnreleased = false;
#pragma omp parallel for schedule(static, 1024)
		for (int i = 0; i < x.size(); i++) {
			bool bounded = x[i] <= lower[i] && r_tilde[i] < 0;

			if (itr % 64 == 1) {
				// Use the full bound update method
				if (boundSet[i] != bounded) {
					boundSet[i] = bounded;
					boundsChanged = true;
				}
			}
			else if (bounded && !boundSet[i]) {
				// We only want to bound things and not release too often to prevent slow convergence due to oscillations
				boundSet[i] = true;
				boundsChanged = true;
			}
			else if (boundSet[i] != bounded)
				haveUnreleased = true; // We dont want to exit before all the unreleased stuff is released

			if (x[i] < lower[i]) {
				x[i] = lower[i];
				projected = true;
			}

			r[i] = boundSet[i] ? 0 : r_tilde[i];
		}
		rDotD[1] = rDotD[0];

		if (boundsChanged) {
			d = invDiag.array() * r.array();
			rDotD[0] = r.dot(d);
		}
		else {
			s = invDiag.array() * r.array();
			rDotD[0] = r.dot(s);
			d = s + (rDotD[0] / rDotD[1]) * d;
		}
	}
	// printf("%zd\n", itr);
	delete[] boundSet;
	return x;
}

VectorXd Kitten::ebccg(
	SparseMatrix<double>& A,
	VectorXd& b,
	VectorXd& lower,
	VectorXd& upper,
	const double tol,
	const int itrLim,
	const int k) {

	VectorXd x(b.size());
	bool* boundSet = new bool[b.size()];

#pragma omp parallel for schedule(static, 4096)
	for (int i = 0; i < b.size(); i++) {
		x[i] = std::max(lower[i], std::min(upper[i], 0.));
		boundSet[i] = false;
	}

	VectorXd g = A * x - b;
	VectorXd p(b.size());
	VectorXd q(b.size());
	VectorXd r[2]{ VectorXd(b.size()), VectorXd(b.size()) };

	int cur = 0;
	double tolSquared = pow2(tol);
	double lastRDotR;

	for (size_t itr = 0; itr < itrLim; itr++) {
		bool updateBounds = itr % k == 0;

#pragma omp parallel for schedule(static, 4096)
		for (int i = 0; i < b.size(); i++)
			r[cur][i] = boundSet[i] ? 0 : -g[i];

		double rDotR = r[cur].dot(r[cur]);

		// Polak Ribiere CG
		// Unlike plain BCCG, we don't care if the bounded set changed last loop
		double beta = itr ? (rDotR - r[cur].dot(r[1 - cur])) / lastRDotR : 0;
		if (beta > 0) p = r[cur] + beta * p;
		else p = r[cur];

		// Compute optimal step size
		q = A * p;
		double alpha = r[cur].dot(p) / p.dot(q);

		if (updateBounds) {
			double update = 0;
			double norm = 0;
			bool boundSetChanged = false;

#pragma omp parallel
			{
				double l_update = 0;
				double l_norm = 0;
				bool l_boundSetChanged = false;

#pragma omp for schedule(static, 2048)
				for (int i = 0; i < b.size(); i++) {
					// Update and accumulate delta
					double nx = x[i] + alpha * p[i];
					l_update += pow2(x[i] - nx);
					l_norm += pow2(x[i]);

					// Check if ~x != x
					if (nx < lower[i]) {
						nx = lower[i];
						l_boundSetChanged = true;
					}
					if (nx > upper[i]) {
						nx = upper[i];
						l_boundSetChanged = true;
					}
					x[i] = nx;

					// Check if B^k != B^{k-1}
					bool bound = (nx == lower[i] && g[i] > 0) || (nx == upper[i] && g[i] < 0);

					// Update B^k
					if (bound != boundSet[i]) {
						boundSet[i] = bound;
						l_boundSetChanged = true;
					}
				}
#pragma omp critical
				{
					update += l_update;
					norm += l_norm;
					boundSetChanged |= l_boundSetChanged;
				}
			}

			if (boundSetChanged)
				g = A * x - b;
			else {
				if (update < norm * tolSquared) break;
				g = g + alpha * q;
			}
		}
		else {
			x += alpha * p;
			g = g + alpha * q;
		}

		lastRDotR = rDotR;
		cur = 1 - cur;
	}

	delete[] boundSet;
	return x;
}