// Reference of the host-side playback of the inflow table, standalone (no solver code).
// build: clang++ -O3 -std=c++17 inflow_reference.cpp -o inflow_reference
// usage: inflow_reference TABLE.bin STEP UREF OUT.f32   -> writes u'(y, z) of that solver step, 3*nz*ny float32, component-major,
//        cell a = y + ny*z, in units of U_H (lattice velocity = UREF * value), and prints the time of one update.
// The struct, read_inflow() and inflow_plane() are meant to be copied into setup.cpp unchanged.
// Format and method: make_inflow_table.py (write_table() and the docstring). The file is little-endian; so is every machine we use.
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

struct InflowTable {
	int32_t version, nb, ny, nz, n, taper; // nb = cells across b, n = number of modes
	int64_t seed;
	double u_c; // convection speed in U_H (information only)
	std::vector<double> omega, psi; // [n] rad per b/U_H, rad
	std::vector<float> U, wy, wz; // [nz] mean profile in U_H, [ny] and [nz] taper weights
	std::vector<float> B_re, B_im; // [3][nz][n]
	std::vector<float> Ey_re, Ey_im; // [n][ny]
	double sum_w; // sum of wy*wz over the plane
	std::vector<float> cs, sn; // work arrays [n]
};

static bool read_inflow(const char* path, InflowTable& T) {
	FILE* f = fopen(path, "rb");
	if(!f) return false;
	char magic[8];
	bool ok = fread(magic, 1, 8, f)==8 && memcmp(magic, "SCINFLW1", 8)==0;
	ok = ok && fread(&T.version, 4, 6, f)==6 && fread(&T.seed, 8, 1, f)==1 && fread(&T.u_c, 8, 1, f)==1;
	ok = ok && T.version==1 && T.n>0 && T.ny>0 && T.nz>0;
	if(ok) {
		const size_t n = (size_t)T.n, ny = (size_t)T.ny, nz = (size_t)T.nz;
		T.omega.resize(n); T.psi.resize(n); T.U.resize(nz); T.wy.resize(ny); T.wz.resize(nz);
		T.B_re.resize(3*nz*n); T.B_im.resize(3*nz*n); T.Ey_re.resize(n*ny); T.Ey_im.resize(n*ny);
		ok = fread(T.omega.data(), 8, n, f)==n && fread(T.psi.data(), 8, n, f)==n
		  && fread(T.U.data(), 4, nz, f)==nz && fread(T.wy.data(), 4, ny, f)==ny && fread(T.wz.data(), 4, nz, f)==nz
		  && fread(T.B_re.data(), 4, 3*nz*n, f)==3*nz*n && fread(T.B_im.data(), 4, 3*nz*n, f)==3*nz*n
		  && fread(T.Ey_re.data(), 4, n*ny, f)==n*ny && fread(T.Ey_im.data(), 4, n*ny, f)==n*ny;
		ok = ok && fgetc(f)==EOF; // nothing may follow
		double sy = 0.0, sz = 0.0;
		for(size_t y=0; y<ny; y++) sy += (double)T.wy[y];
		for(size_t z=0; z<nz; z++) sz += (double)T.wz[z];
		T.sum_w = sy*sz;
		T.cs.resize(n); T.sn.resize(n);
	}
	fclose(f);
	return ok;
}

// u'(y, z) at solver step `step` (may be fractional) into out[3*nz*ny], component-major, a = y + ny*z, in U_H.
// u_ref is the lattice velocity of U_H (UREF). The solver's dt is dx*u_ref/u_si, so the time in b/U_H is step*u_ref/nb.
static void inflow_plane(InflowTable& T, const double step, const double u_ref, float* out) {
	const int n = T.n, ny = T.ny, nz = T.nz;
	const double time = step*u_ref/(double)T.nb;
	for(int m=0; m<n; m++) { // phases in double: at 120000 steps they reach 1e4 rad
		const double ph = T.psi[m]-T.omega[m]*time;
		T.cs[m] = (float)cos(ph); T.sn[m] = (float)sin(ph);
	}
	for(int i=0; i<3; i++) for(int z=0; z<nz; z++) { // rows are independent: this loop can run in parallel
		const float* br = &T.B_re[((size_t)i*nz+z)*n];
		const float* bi = &T.B_im[((size_t)i*nz+z)*n];
		float* row = out+((size_t)i*nz+z)*ny;
		for(int y=0; y<ny; y++) row[y] = 0.0f;
		for(int m=0; m<n; m++) {
			const float cr = br[m]*T.cs[m]-bi[m]*T.sn[m]; // B*exp(i*phase), real and imaginary part
			const float ci = br[m]*T.sn[m]+bi[m]*T.cs[m];
			const float* er = &T.Ey_re[(size_t)m*ny];
			const float* ei = &T.Ey_im[(size_t)m*ny];
			for(int y=0; y<ny; y++) row[y] += cr*er[y]-ci*ei[y];
		}
		for(int y=0; y<ny; y++) row[y] *= T.wy[y];
	}
	double sum = 0.0; // remove the mean of u over the plane, weighted with the taper: no net flux fluctuation
	for(size_t a=0; a<(size_t)ny*nz; a++) sum += (double)out[a];
	const float c = (float)(sum/T.sum_w);
	for(int z=0; z<nz; z++) for(int y=0; y<ny; y++) out[(size_t)z*ny+y] -= c*T.wz[z]*T.wy[y];
}

int main(int argc, char** argv) {
	if(argc<5) { fprintf(stderr, "usage: inflow_reference TABLE.bin STEP UREF OUT.f32\n"); return 2; }
	InflowTable T;
	if(!read_inflow(argv[1], T)) { fprintf(stderr, "cannot read %s\n", argv[1]); return 1; }
	const double step = atof(argv[2]), u_ref = atof(argv[3]);
	std::vector<float> out((size_t)3*T.nz*T.ny);
	inflow_plane(T, step, u_ref, out.data());
	FILE* f = fopen(argv[4], "wb");
	if(!f || fwrite(out.data(), 4, out.size(), f)!=out.size()) { fprintf(stderr, "cannot write %s\n", argv[4]); return 1; }
	fclose(f);
	std::vector<float> tmp(out.size());
	const int reps = 20;
	const auto t0 = std::chrono::steady_clock::now();
	for(int r=0; r<reps; r++) inflow_plane(T, step+10.0*r, u_ref, tmp.data());
	const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now()-t0).count()/reps;
	double sum = 0.0;
	for(size_t a=0; a<(size_t)T.ny*T.nz; a++) sum += (double)out[a];
	printf("b%d: %d modes, %d x %d cells, one update %.2f ms on one thread, plane mean of u' %.2e U_H\n", T.nb, T.n, T.ny, T.nz, ms, sum/((double)T.ny*(T.nz-1)));
	return 0;
}
