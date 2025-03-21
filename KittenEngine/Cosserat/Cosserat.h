#pragma once
#include <glm/glm.hpp>
#include <KittenEngine/includes/KittenEngine.h>
#include <KittenEngine/includes/modules/ComputeBuffer.h>
#include <KittenEngine/includes/modules/Rotor.h>
#include <functional>
#include <vector>

#define COSSERAT_THREAD_BATCH (90)
#define COSSERAT_THREAD_BATCH_STATIC (128)

namespace Cosserat {
	using namespace glm;
	using namespace Kitten;

	struct Seg {
		Rotor q;
		ivec2 i;
		float kss;		// Stretch stiffness. Negative if pinned.
		float l;
	};

	struct Vert {
		vec3 pos;
		float invMass;
	};

	struct RestAngle {
		Rotor qRest;	// qRest = q0^-1 * q1
		ivec2 i;		// Segment ids
		float kbt;		// Bending stiffness
	};

	struct MetaData {
		/// Simulation parameters
		float h = 1 / 120.f;	// Current time step size
		vec3 gravity;			// Gravity
		float damping = 1e-7;	// Damping
		float drag = 0.1;		// Drag

		float mu = 0.2;			// Friction
		float radius = 0.01f;	// Rod radius

		int useProportionalDamping = 1;	// If enabled, segments of all stiffnesses will have the same effective damping
		int numItr = 4;			// Number of iterations per step
		float lastH;			// Last time step size
		float time = 0;			// Current time

		float segStrainThreshold = INFINITY;	// Strain threshold for segment destruction
		float angleStrainThreshold = INFINITY;	// Strain threshold for angle destruction
	};

	enum class IntegrationMethod {
		StableCosserat,			// Stable Cosserat
		StableCosseratExact,	// Stable Cosserat exact
		VBD,			// Vertex Block Descent
		LCVBD,			// VBD with Linearized Constraint
		PSDPLCVBD,		// VBD with Linearized Constraint and Positive Semi-Definite Projection
		XPBD,			// Extended Position Based Dynamics
		NUM_METHODS		// Number of methods
	};

	class Sim {
	public:
		/// Data

		vector<Seg> segs;
		vector<Vert> verts;
		vector<RestAngle> restAngles;

		vector<vec3> lastVels;
		vector<vec3> vels;

		vector<vec3> xpbdSegWs;			// XPBD angular velocity
		vector<float> lambdaTradGammas;	// Fixed point lamabda stored as gamma for our method.

		// Index tables

		vector<int> vertSegIndices;		// Prefix sum of the number of segments connected to each vertex
		vector<int> vertSegIds;			// Segment ids connected to each vertex

		vector<int> segAngleIndices;	// Prefix sum of the number of rest angles connected to each segment
		vector<int> segAngleIds;		// Rest angle ids connected to each segment

		/// Simulation settings

		float maxH = 1 / 120.f;			// Maximum time step size allowed

		MetaData meta;
		int stepCounter = 0;			// A counter of how many substeps have been taken

		bool multithreaded = false;

		// The integration method to use. May be changed at any time
		IntegrationMethod method = IntegrationMethod::StableCosserat;

		/// Render settings

		float meshPhongThreshold = 0.5f;// Phong sin angle for mesh export
		float meshHeightLen = 0.025f;	// Mesh export seg sample length
		int meshRadialSegs = 12;		// Mesh export seg radial segments
		float meshJunctionPenalty = 0;	// Phong penalty for junctions (vertices with more than 2 segments)
		int meshExcludeInlinePhong = 0;	// Exclude inline vertices from phong breaking

		float renderDensity = 1.f;		// Fake density used for estimating rod volume
		float renderRadius = 0.01f;		// Radius of the rods
		float renderByDensity = 0.f;	// A blend of whether to render the rods scaled by vertex mass

		// Set this to a function that returns the external force on a vertex
		// Arguemnts are the simulation, the vertex and the simulation time
		std::function<vec3(Sim*, Vert&, float)> externalForce = nullptr;

	private:
		ComputeBuffer* vertBuffer = nullptr;
		ComputeBuffer* segBuffer = nullptr;
		Mesh* renderMesh = nullptr;

		vector<vec3> dxs;				// Displacements for VBD
		vector<vec3> xpbdSegLambdas;	// XPBD lagrange multipliers for segment constraints
		vector<vec4> xpbdAngleLambdas;	// XPBD lagrange multipliers for rest angle constraints
		vector<Rotor> xpbdLastQs;		// XPBD last orientations
		vector<int> vertexLastSegment;	// Last segment connected to a vertex. Only used for initialization

		void startIterate();
		void iterate();
		void endIterate();

		void startIterateXPBD();
		void iterateVertVBD();
		void iterateVertXPBD();
		void endIterateXPBD();

		void iterateSegLambda();
		void iterateSegLambdaExact();
		void iterateSegLCVBD();
		void iterateSegPSDPLCVBD();
		void iterateSegVBD();
		void iterateSegXPBD();

		void uploadGPU();

	public:
		Sim();
		~Sim();

		// CALL THIS after adding all vertices, segments and rest angles.
		// Initializes the simulation and generates all the index tables
		void init();

		// Add a vertex to the simulation
		void addVertex(vec3 pos);
		// Add a segment linking two vertices
		void addSegment(ivec2 i, float kss = 1e2f, float density = 1.f);
		// Add a rest angle between two segments for bending/twisting
		void addRestAngle(ivec2 i, float kbt = 1.f);

		// Pin a vertex in place
		void pinVertex(int i);
		// Pin segment in place
		void pinSegment(int i);

		// Deform or destory any segment/angles that doesnt pass a predicate.
		// The function should return true if the segment/angle should be kept.
		// Changes to rest shape can be done directly in the function through reference.
		// The predicate provides the simulation, the segment/angle and the strain.
		// Pass a nullptr to skip either segments or angles checks.
		void deformDestroy(
			std::function<bool(Sim*, Seg&, vec3)> segPredicate,
			std::function<bool(Sim*, RestAngle&, vec4)> anglePredicate);

		// Removes any segments or angles that exceed a strain threshold
		void destroyByStrainMagnitude(float segStrainThreshold = INFINITY, float angleStrainThreshold = INFINITY);

		// Sets a wind field defined on position and time
		void setWindFieldExternalForce(std::function<vec3(vec3, float)> windField, float dragCoeff = 1e-1f);

		// Uses a single time step to advance the simulation
		void step(float dt);
		// Advances the simulation until the given time using multiple time steps
		float advance(float dt);

		void render();

		// Debugging and convergence checkers
		void performSelfCheck();
		float getPotentialEnergy();	// Gets the total internal potential energy of the system
		float getIncrementalPotential();
		void checkSegConvergence();	// Checks convergence
		void checkSegConvergenceOverTime();	// Checks error over time
		void checkExactLambdaConvergence();	// Checks convergence of exact lambda

		void exportMesh(const char* path);

		// Save States

		void loadFromBinary(const char* path, bool loadOnlyGeometry = false);
		void saveToBinary(const char* path, bool saveOnlyGeomtry = false);
	};
}