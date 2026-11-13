// Altered from FluidX3D commit b85d64e (v3.8) for a wind study of Warsaw Spire and Łucka street, 2026-11; see README.md of the published repository.
#include "setup.hpp"
#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <future>

// Wind over a city: the solid is a 2.5D height map (cells above flat ground) plus an optional STL, the wind comes from the west (+x).
// Environment: HM HX HY NZ DX_M STL UREF ZREF_M ALPHA NU T_SI SAMPLE_SI FRAMES OUT PX PY ZPED ZHI, cameras QA QB QD QF QX / WA WB WD WF / SD SF,
// box VX0 VX1 VY0 VY1 VZ1, animation ANIM_T ANIM_SI; FX_U_MAX and FX_Q are read in lbm.cpp. Every variable is parsed strictly and checked in main_setup().
// Added for the wind-tunnel check (AIJ case A, docs/VALIDATION.md), all off unless set: XST (x-normal planes for second moments), PROFILE (tabulated inflow), VOL_TAGGED (box mean
// kept per frame), FX_CS (Smagorinsky constant, read in lbm.cpp), Z0_M (roughness of the wall-stress ground, only in the GROUND_WALL_MODEL build),
// and three hooks for the tests only: PERIODIC_XY, WIND_DEG, ZPROF_EVERY.
// Added for the turbulent inflow of that check, all off unless set: INFLOW (mode table of the turbulent inflow, played back at the plane x = 0; it also carries
// the mean profile of every fixed face and of the initial field), INFLOW_K (steps between two inlet planes, default 10), VOL_M2 (standard deviations in the 3D box),
// and two hooks for the playback test only: INFLOW_STEP0, INFLOW_DUMP.
// Added before the scored runs of that check, all off unless set and none of them touching the flow: SPIN_SI (seconds before the averaging starts), ZPEDS (near-ground
// planes on more layers), RENDER (0: no rendered frames), ANIM_SECTION (u and w on the plane y = PY per film frame), and the test hook FLAGS_DUMP.
static void fail(const string& msg) { print_error("wind setup: "+msg); } // print_error exits
static bool parse_float(const char* s, float& v) { char* e=nullptr; errno=0; v=(float)strtod(s, &e); return e!=s&&*e=='\0'&&errno==0; }
static bool parse_uint(const char* s, uint& v) { char* e=nullptr; errno=0; const unsigned long x=strtoul(s, &e, 10); v=(uint)x; return s[0]>='0'&&s[0]<='9'&&*e=='\0'&&errno==0&&x<=0xFFFFFFFFul; }
static uint env_uint(const char* k, const uint d) { const char* v = getenv(k); uint r=d; if(v&&!parse_uint(v, r)) fail(string(k)+"=\""+v+"\" is not a non-negative integer"); return r; }
static float env_float(const char* k, const float d) { const char* v = getenv(k); float r=d; if(v&&!parse_float(v, r)) fail(string(k)+"=\""+v+"\" is not a number"); return r; }
static string env_str(const char* k, const string& d) { const char* v = getenv(k); return v ? string(v) : d; }
static void write_f32(const string& path, const vector<float>& a) { std::ofstream f(path, std::ios::binary); f.write((const char*)a.data(), (std::streamsize)(a.size()*sizeof(float))); }
static bool parse_uint_list(const char* s, vector<uint>& v) { // "a,b,c": every item a non-negative integer, no empty item
	v.clear(); string item;
	for(const char* p=s; ; p++) {
		if(*p!=','&&*p!='\0') { item += *p; continue; }
		uint x = 0u; if(!parse_uint(item.c_str(), x)) return false;
		v.push_back(x); item.clear();
		if(*p=='\0') return true;
	}
}
static void read_profile(const string& path, vector<float>& zs, vector<float>& us) { // rows "z_m U": height above the ground in metres and speed as a fraction of UREF; a line starting with '#' is a comment
	std::ifstream f(path); if(!f) fail("cannot open PROFILE=\""+path+"\"");
	string line; uint row = 0u;
	while(std::getline(f, line)) {
		row++;
		if(line.empty()||line[0]=='#') continue;
		const char* a = line.c_str(); char* e = nullptr; errno = 0;
		const double z = strtod(a, &e); const char* b = e; const double u = strtod(b, &e);
		const bool two = b!=a&&e!=b; // both numbers were read
		while(*e==' '||*e=='\t'||*e=='\r') e++;
		if(!two||*e!='\0'||errno!=0||!(z>=0.0&&z<1.0E9)||!(u>=0.0&&u<1.0E9)||(!zs.empty()&&!((float)z>zs.back()))) fail("PROFILE row "+to_string(row)+" \""+line+"\" must be two numbers \"z_m U\", z_m >= 0 and rising, U >= 0");
		zs.push_back((float)z); us.push_back((float)u);
	}
	if(zs.size()<2u) fail("PROFILE needs at least two rows");
}

// Turbulent inflow: playback of a table of random Fourier modes. Method and units: validation/scripts/make_inflow_table.py; file format: write_table() there.
// The struct, read_inflow() and inflow_plane() follow the reference validation/scripts/inflow_reference.cpp (checked against the Python generator); the only change is
// that the rows of a plane are computed by several threads. Everything in the table is in units of the building width b and of the speed at building height U_H.
struct InflowTable {
	int32_t version=0, nb=0, ny=0, nz=0, n=0, taper=0; // nb = cells across b, n = number of modes
	int64_t seed = 0;
	double u_c = 0.0; // convection speed in U_H (information only)
	vector<double> omega, psi; // [n] rad per b/U_H, rad
	vector<float> U, wy, wz; // [nz] mean profile in U_H, [ny] and [nz] taper weights
	vector<float> B_re, B_im; // [3][nz][n]
	vector<float> Ey_re, Ey_im; // [n][ny]
	double sum_w = 0.0; // sum of wy*wz over the plane
	vector<float> cs, sn; // work arrays [n]
};
static bool read_inflow(const char* path, InflowTable& T) { // false when the file is missing, is not a version-1 table, is shorter or longer than its header says
	FILE* f = fopen(path, "rb");
	if(!f) return false;
	char magic[8];
	bool ok = fread(magic, 1, 8, f)==8&&memcmp(magic, "SCINFLW1", 8)==0;
	ok = ok&&fread(&T.version, 4, 6, f)==6&&fread(&T.seed, 8, 1, f)==1&&fread(&T.u_c, 8, 1, f)==1; // six int32 in a row: version, nb, ny, nz, n, taper
	ok = ok&&T.version==1&&T.nb>0&&T.n>0&&T.ny>0&&T.nz>0&&T.n<=65536&&T.ny<=65536&&T.nz<=65536; // the upper bounds keep a damaged header from asking for absurd memory
	if(ok) {
		const size_t n = (size_t)T.n, ny = (size_t)T.ny, nz = (size_t)T.nz;
		T.omega.resize(n); T.psi.resize(n); T.U.resize(nz); T.wy.resize(ny); T.wz.resize(nz);
		T.B_re.resize(3*nz*n); T.B_im.resize(3*nz*n); T.Ey_re.resize(n*ny); T.Ey_im.resize(n*ny);
		ok = fread(T.omega.data(), 8, n, f)==n&&fread(T.psi.data(), 8, n, f)==n
		  &&fread(T.U.data(), 4, nz, f)==nz&&fread(T.wy.data(), 4, ny, f)==ny&&fread(T.wz.data(), 4, nz, f)==nz
		  &&fread(T.B_re.data(), 4, 3*nz*n, f)==3*nz*n&&fread(T.B_im.data(), 4, 3*nz*n, f)==3*nz*n
		  &&fread(T.Ey_re.data(), 4, n*ny, f)==n*ny&&fread(T.Ey_im.data(), 4, n*ny, f)==n*ny;
		ok = ok&&fgetc(f)==EOF; // nothing may follow
		double sy = 0.0, sz = 0.0;
		for(size_t y=0; y<ny; y++) sy += (double)T.wy[y];
		for(size_t z=0; z<nz; z++) sz += (double)T.wz[z];
		T.sum_w = sy*sz;
		T.cs.resize(n); T.sn.resize(n);
	}
	fclose(f);
	return ok;
}
// u'(y, z) at solver step `step` into out[3*nz*ny], component-major, cell a = y + ny*z, in U_H. u_ref is the lattice velocity of U_H (UREF).
// The solver's dt is dx*u_ref/u_si, so the time in b/U_H is step*u_ref/nb.
static void inflow_plane(InflowTable& T, const double step, const double u_ref, float* out) {
	const int n = T.n, ny = T.ny, nz = T.nz;
	const double time = step*u_ref/(double)T.nb;
	for(int m=0; m<n; m++) { // phases in double: at 120000 steps they reach 1e4 rad
		const double ph = T.psi[m]-T.omega[m]*time;
		T.cs[m] = (float)cos(ph); T.sn[m] = (float)sin(ph);
	}
	const uint threads = max(min((uint)thread::hardware_concurrency(), 8u), 1u); // a rented container reports the cores of its host, not its own share: more threads than this only cost their start-up
	parallel_for((ulong)(3*nz), threads, [&](ulong r) { // row r = component*nz + z; rows are independent
		const float* br = &T.B_re[(size_t)r*n];
		const float* bi = &T.B_im[(size_t)r*n];
		float* row = out+(size_t)r*ny;
		for(int y=0; y<ny; y++) row[y] = 0.0f;
		for(int m=0; m<n; m++) {
			const float cr = br[m]*T.cs[m]-bi[m]*T.sn[m]; // B*exp(i*phase), real and imaginary part
			const float ci = br[m]*T.sn[m]+bi[m]*T.cs[m];
			const float* er = &T.Ey_re[(size_t)m*ny];
			const float* ei = &T.Ey_im[(size_t)m*ny];
			for(int y=0; y<ny; y++) row[y] += cr*er[y]-ci*ei[y];
		}
		for(int y=0; y<ny; y++) row[y] *= T.wy[y];
	});
	double sum = 0.0; // remove the mean of u over the plane, weighted with the taper: no net flux fluctuation
	for(size_t a=0; a<(size_t)ny*nz; a++) sum += (double)out[a];
	const float c = (float)(sum/T.sum_w);
	for(int z=0; z<nz; z++) for(int y=0; y<ny; y++) out[(size_t)z*ny+y] -= c*T.wz[z]*T.wy[y];
}

void main_setup() {
	const string hm = env_str("HM", "");
	const uint HX = env_uint("HX", 1000u), HY = env_uint("HY", 600u), Nz = env_uint("NZ", 200u);
	const float dx = env_float("DX_M", 3.0f); // metres per cell
	const float u_ref = env_float("UREF", 0.05f); // lattice velocity at z_ref
	const float t_si = env_float("T_SI", 600.0f), sample_si = env_float("SAMPLE_SI", 1.0f);
	if(HX==0u||HY==0u||Nz==0u) fail("HX, HY and NZ must be > 0");
	if(!(dx>0.0f)||!(u_ref>0.0f)||!(t_si>0.0f)) fail("DX_M, UREF and T_SI must be > 0");
	if(getenv("FX_U_MAX")&&!(env_float("FX_U_MAX", 1.0f)>0.0f)) fail("FX_U_MAX must be > 0"); // lbm.cpp reads it when the kernels are built
	(void)env_float("FX_Q", 0.0f); // only checks that FX_Q is a number; lbm.cpp reads it
	const float z_ref = env_float("ZREF_M", 100.0f)/dx, alpha = env_float("ALPHA", 0.24f); // power law, PN-EN 1991-1-4 NA terrain IV exponent
	const float nu = env_float("NU", 2.0E-8f); // air at full scale is ~1E-8 in lattice units; SUBGRID (Smagorinsky) keeps it stable
	const float u_si = 10.0f; // assumed wind speed at z_ref [m/s], only for converting time and speed
	const float dt = dx*u_ref/u_si; // seconds per time step
	const ulong steps = (ulong)(t_si/dt);
	const ulong sample = max((ulong)(sample_si/dt), 1ull); // averaging interval in steps; the real period is sample*dt
	const uint frames = env_uint("FRAMES", 4u); // render sets during the run, the last one at the end
	const uint px = env_uint("PX", 633u), py = env_uint("PY", 300u); // landmark cell (the tallest building): the cameras aim at it and the section passes through y = py
	const uint zped = env_uint("ZPED", 1u); // near-ground slice, cells above the ground layer
	const uint zhi = env_uint("ZHI", (uint)(100.0f/dx+0.5f)); // slice at about 100 m
	string out = env_str("OUT", get_exe_path()+"export/"); if(out.empty()||out.back()!='/') out += "/";
	const string stl = env_str("STL", ""); // optional solid as a binary STL in cell-index coordinates: cell (x, y, z) is centred at (x, y, z) and the origin is the domain corner, so a mesh drawn in domain metres is divided by DX_M and shifted by -0.5 cell
	const bool vol_set = getenv("VX0")!=nullptr||getenv("VX1")!=nullptr||getenv("VY0")!=nullptr||getenv("VY1")!=nullptr||getenv("VZ1")!=nullptr; // optional 3D box (cells VX0..VX1-1, VY0..VY1-1, 0..VZ1-1): mean and last velocity for streamlines, mean density for pressure
	const uint vx0=env_uint("VX0", 0u), vx1=env_uint("VX1", 0u), vy0=env_uint("VY0", 0u), vy1=env_uint("VY1", 0u), vz1=env_uint("VZ1", 0u);
	vector<uint> xst; // x-normal planes (cell indices) on which sums of u, v, w, uu, vv, ww, uw are accumulated over the averaging window
	if(getenv("XST")&&!parse_uint_list(getenv("XST"), xst)) fail(string("XST=\"")+getenv("XST")+"\" is not a comma-separated list of non-negative integers");
	for(const uint x : xst) if(x>=HX) fail("every XST plane must be < HX");
	vector<float> prof_z, prof_u; // tabulated inflow profile; empty = power law
	if(getenv("PROFILE")) read_profile(getenv("PROFILE"), prof_z, prof_u);
	const uint vol_tagged = env_uint("VOL_TAGGED", 0u); // 1: keep the box mean of every frame (vol_mean-<frame>.f32), so the settling check needs no copy made while the solver runs
	if(vol_tagged>1u) fail("VOL_TAGGED must be 0 or 1");
	if(getenv("FX_CS")&&!(env_float("FX_CS", 0.0f)>0.0f&&env_float("FX_CS", 0.0f)<=0.5f)) fail("FX_CS must be > 0 and <= 0.5"); // Smagorinsky constant; lbm.cpp reads it when the kernels are built
	const uint periodic_xy = env_uint("PERIODIC_XY", 0u); // test hook: 1 leaves the four side faces periodic instead of fixed-profile faces
	const float wind_deg = env_float("WIND_DEG", 0.0f); // test hook: turns the inflow about z by this angle (0 = from the west, along +x)
	const uint zprof_every = env_uint("ZPROF_EVERY", 0u); // test hook: every so many steps record the sums of density and momentum of every layer z
	if(periodic_xy>1u) fail("PERIODIC_XY must be 0 or 1");
	if(!(wind_deg>-360.0f&&wind_deg<360.0f)) fail("WIND_DEG must be between -360 and 360");
	const uint vol_m2 = env_uint("VOL_M2", 0u); // 1: also keep the standard deviations of u, v, w and the covariance of u and w in the 3D box (vol_std.f32)
	if(vol_m2>1u) fail("VOL_M2 must be 0 or 1");
	if(vol_m2==1u&&!vol_set) fail("VOL_M2=1 needs the 3D box (VX0 VX1 VY0 VY1 VZ1)");
	const float spin_si = env_float("SPIN_SI", -1.0f); // seconds of wind before the averaging starts; unset: the second half of the run is averaged
	if(getenv("SPIN_SI")&&!(spin_si>=0.0f&&spin_si<t_si)) fail("SPIN_SI must be >= 0 and below T_SI");
	vector<uint> zpeds; // near-ground layers that get the planes of ZPED under their own names (hz<k>_mean.f32, hz<k>_mean_ux.f32, hz<k>_mean_uy.f32, hz<k>_std.f32, hz<k>_mean-<frame>.f32)
	if(getenv("ZPEDS")&&!parse_uint_list(getenv("ZPEDS"), zpeds)) fail(string("ZPEDS=\"")+getenv("ZPEDS")+"\" is not a comma-separated list of non-negative integers");
	for(ulong i=0ull; i<zpeds.size(); i++) { if(zpeds[i]==0u||zpeds[i]>=Nz) fail("every ZPEDS layer must be >= 1 and < NZ"); for(ulong j=0ull; j<i; j++) if(zpeds[j]==zpeds[i]) fail("ZPEDS names a layer twice"); }
	const uint render_on = env_uint("RENDER", 1u); // 0: no rendered frames (q3d, q3d_wide, side, the film image); every plane is still written
	if(render_on>1u) fail("RENDER must be 0 or 1");
	const uint anim_section = env_uint("ANIM_SECTION", 0u); // 1: every film frame also writes u and w on the plane y = PY (anim/sec_u-<k>.f32, anim/sec_w-<k>.f32, nx*nz float32 each, x fastest)
	if(anim_section>1u) fail("ANIM_SECTION must be 0 or 1");
	if(anim_section==1u&&!(env_float("ANIM_T", 0.0f)>0.0f)) fail("ANIM_SECTION=1 needs a film: ANIM_T > 0");
	const uint flags_dump = env_uint("FLAGS_DUMP", 0u); // test hook: 1 writes the flag field as the device holds it after the initialisation (planes/flags.u8, one byte per cell)
	if(flags_dump>1u) fail("FLAGS_DUMP must be 0 or 1");
	const bool inflow_on = getenv("INFLOW")!=nullptr; // turbulent inflow: the mode table is played back at the plane x = 0
	InflowTable inflow;
	const uint inflow_k = env_uint("INFLOW_K", 10u); // steps between two inlet planes; in between the inlet velocity is interpolated linearly in time
	const ulong inflow_step0 = (ulong)env_uint("INFLOW_STEP0", 0u); // test hook: playback step of solver step 0, so that a late plane can be checked in a short run
	vector<uint> inflow_dump; // test hook: playback steps after which the inlet velocity is read back from the device and written to planes/inlet-<step>.f32
	if(getenv("INFLOW_DUMP")&&!parse_uint_list(getenv("INFLOW_DUMP"), inflow_dump)) fail(string("INFLOW_DUMP=\"")+getenv("INFLOW_DUMP")+"\" is not a comma-separated list of non-negative integers");
	if(!inflow_on&&(getenv("INFLOW_K")||getenv("INFLOW_STEP0")||getenv("INFLOW_DUMP"))) fail("INFLOW_K, INFLOW_STEP0 and INFLOW_DUMP need INFLOW");
	if(inflow_on) {
		if(inflow_k==0u) fail("INFLOW_K must be > 0");
		if(getenv("PROFILE")) fail("INFLOW carries its own mean profile: unset PROFILE");
		if(periodic_xy!=0u||wind_deg!=0.0f) fail("INFLOW needs the inlet at x = 0: PERIODIC_XY and WIND_DEG must be 0");
		if(!read_inflow(getenv("INFLOW"), inflow)) fail(string("INFLOW=\"")+getenv("INFLOW")+"\" cannot be read as a mode table (magic SCINFLW1, version 1, exact size)");
		if((uint)inflow.ny!=HY||(uint)inflow.nz!=Nz) fail("the INFLOW table is made for a plane of "+to_string((uint)inflow.ny)+" x "+to_string((uint)inflow.nz)+" cells, the grid has HY x NZ = "+to_string(HY)+" x "+to_string(Nz));
		for(const uint s : inflow_dump) if((ulong)s<inflow_step0) fail("every INFLOW_DUMP step must be >= INFLOW_STEP0");
	}
#ifdef GROUND_WALL_MODEL
	const float z0 = env_float("Z0_M", 0.0f); // roughness length of the ground in metres
	if(!(z0>0.0f&&z0<0.5f*dx)) fail("this build has the wall-stress ground: Z0_M must be > 0 and below the first cell centre (DX_M/2)");
	fx_ground_cd = sq(0.4f/log(0.5f*dx/z0)); // log law at the first cell centre z1 = dx/2: u_tau = kappa*U(z1)/ln(z1/z0), kappa = 0.4, tau = rho*u_tau^2 = rho*Cd*U(z1)^2
	if(!(fx_ground_cd*u_ref<0.1f)) fail("Z0_M gives a wall stress too large for one time step"); // the cell must keep most of its momentum
#else // GROUND_WALL_MODEL
	if(getenv("Z0_M")) fail("Z0_M is set but this build has no wall-stress ground (build with -DGROUND_WALL_MODEL)");
#endif // GROUND_WALL_MODEL
	if(frames==0u) fail("FRAMES must be > 0");
	if(zped>=Nz||zhi>=Nz) fail("ZPED and ZHI must be < NZ");
	if(px>=HX||py>=HY) fail("PX must be < HX and PY must be < HY");
	if(vol_set&&!(vx0<vx1&&vx1<=HX&&vy0<vy1&&vy1<=HY&&vz1>0u&&vz1<=Nz)) fail("the box needs VX0 < VX1 <= HX, VY0 < VY1 <= HY and 0 < VZ1 <= NZ");

	vector<ushort> H((ulong)HX*HY);
	{ std::ifstream f(hm, std::ios::binary); f.read((char*)H.data(), (std::streamsize)(H.size()*sizeof(ushort))); if(!f) fail("cannot read "+to_string((ulong)(H.size()*sizeof(ushort)))+" bytes of height map HM=\""+hm+"\""); }

	LBM lbm(HX, HY, Nz, nu);
	const uint Nx = lbm.get_Nx(), Ny = lbm.get_Ny();
	const ulong N = lbm.get_N(), NxNy = (ulong)Nx*Ny;
	const auto profile = [&](const float z) -> float { // wall sits at z=0.5
		if(inflow_on) return u_ref*inflow.U[min((uint)z, Nz-1u)]; // the table holds the mean speed of every layer
		if(prof_z.empty()) return u_ref*pow(fmax(z-0.5f, 0.5f)/z_ref, alpha);
		const float zm = fmax(z-0.5f, 0.5f)*dx; // height of the cell centre above the ground in metres
		if(zm<=prof_z.front()) return u_ref*prof_u.front(); // outside the table the end values hold
		if(zm>=prof_z.back()) return u_ref*prof_u.back();
		const ulong k = (ulong)(std::upper_bound(prof_z.begin(), prof_z.end(), zm)-prof_z.begin()); // prof_z[k-1] <= zm < prof_z[k]
		return u_ref*(prof_u[k-1ull]+(prof_u[k]-prof_u[k-1ull])*(zm-prof_z[k-1ull])/(prof_z[k]-prof_z[k-1ull])); // linear between the rows
	};
	const float wind_cos = cos(wind_deg*pif/180.0f), wind_sin = sin(wind_deg*pif/180.0f);
	if(stl!="") { // voxelize the solid first (reads flags back to the host), then the height map and the inflow around it
		Mesh* mesh = read_stl(stl, 1.0f);
		lbm.voxelize_mesh_on_device(mesh, TYPE_S);
		ulong solid = 0ull; for(ulong n=0ull; n<N; n++) solid += (lbm.flags[n]&TYPE_S) ? 1ull : 0ull;
		print_info("STL "+stl+": "+to_string(mesh->triangle_number)+" triangles, "+to_string(solid)+" solid cells");
	}
	parallel_for(N, [&](ulong n) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		if(lbm.flags[n]&TYPE_S) return; // the STL solid
		if(z<(uint)H[(ulong)y*HX+x]) { lbm.flags[n] = TYPE_S; return; } // ground layer and buildings
		lbm.u.x[n] = profile((float)z);
		if(wind_deg!=0.0f) { lbm.u.y[n] = lbm.u.x[n]*wind_sin; lbm.u.x[n] *= wind_cos; }
		if((periodic_xy==0u&&(x==0u||x==Nx-1u||y==0u||y==Ny-1u))||z==Nz-1u) lbm.flags[n] = TYPE_E; // inflow, outflow, sides, top
#ifdef GROUND_WALL_MODEL
		else if(z==1u&&H[(ulong)y*HX+x]==1u) lbm.flags[n] = TYPE_X; // first air cell above open ground: the kernel applies the wall stress there; roofs keep bounce-back (at a roof edge a link from below comes from a fluid cell)
#endif // GROUND_WALL_MODEL
	});
#ifdef GROUND_WALL_MODEL
	print_info("wall-stress ground: z0 = "+to_string(z0, 6u)+" m, first cell centre at "+to_string(0.5f*dx, 3u)+" m, Cd = "+to_string(fx_ground_cd, 6u));
	{ // count the marked cells and check what apply_ground_stress() relies on: the cell below a marked cell and its four side neighbours in the ground layer are solid, so the four upward diagonal populations of the marked cell are its own, bounced back
		ulong marked = 0ull, bad = 0ull;
		for(ulong i=0ull; i<NxNy; i++) if(lbm.flags[NxNy+i]&TYPE_X) { // marks exist in the layer z = 1 only
			marked++;
			const ulong x = i%(ulong)Nx, y = i/(ulong)Nx;
			const ulong below[5] = { i, (x+1ull)%Nx+y*Nx, (x+Nx-1ull)%Nx+y*Nx, x+((y+1ull)%Ny)*Nx, x+((y+Ny-1ull)%Ny)*Nx }; // neighbours wrap as in the kernel
			for(const ulong b : below) if((lbm.flags[b]&(TYPE_S|TYPE_E))!=TYPE_S) { bad++; break; }
		}
		if(bad>0ull) fail(to_string(bad)+" cells marked for the wall stress have a cell of the ground layer below or beside them that is not solid (a height map value of 0?)");
		print_info("wall-stress ground: "+to_string(marked)+" cells of the first air layer marked, the ground layer below each of them is solid");
	}
#endif // GROUND_WALL_MODEL
	if(!prof_z.empty()) print_info("inflow: tabulated profile, "+to_string((ulong)prof_z.size())+" rows from PROFILE");
	// turbulent inflow: plane j is the fluctuation at playback step j*inflow_k. The device keeps two planes (slot j%2 holds plane j); before a time step the kernel
	// set_inlet blends them. The host computes the plane after next in a second thread while the device runs, and uploads it when the playback reaches an update step.
	const ulong Ax = (ulong)Ny*Nz; // cells of the inlet plane
	vector<float> inflow_stage; // one plane in U_H, the layout of inflow_plane()
	ulong inflow_j_dev = 0ull; // highest plane index on the device
	double inflow_wait_s = 0.0; ulong inflow_uploads = 0ull; // time the solver stood still because the next plane was not ready, and the number of planes uploaded during the run
	std::future<void> inflow_next; // computes plane inflow_j_dev+1 into inflow_stage
	const auto inflow_launch = [&](const ulong j) { // start computing plane j into inflow_stage; the caller joins it (inflow_next.get()) before it touches inflow_stage again
		const double step = (double)(j*(ulong)inflow_k);
		inflow_next = std::async(std::launch::async, [&inflow, &inflow_stage, step, u_ref]() { inflow_plane(inflow, step, (double)u_ref, inflow_stage.data()); });
	};
	const auto inflow_store = [&](const ulong j) { // inflow_stage (plane j) -> lattice units with the mean profile -> host copy of slot j%2
		float* slot = lbm.lbm_domain[0]->inlet.data()+(j%2ull)*3ull*Ax;
		for(ulong a=0ull; a<Ax; a++) slot[a] = u_ref*(inflow.U[a/(ulong)Ny]+inflow_stage[a]);
		for(ulong a=Ax; a<3ull*Ax; a++) slot[a] = u_ref*inflow_stage[a];
	};
	if(inflow_on) {
		if(lbm.get_D()!=1u) fail("INFLOW works on a single device only");
		lbm.lbm_domain[0]->enable_inlet(inflow_k, inflow_step0);
		inflow_stage.resize(3ull*Ax);
		const ulong j0 = inflow_step0/(ulong)inflow_k;
		for(ulong j=j0; j<=j0+1ull; j++) { inflow_plane(inflow, (double)(j*(ulong)inflow_k), (double)u_ref, inflow_stage.data()); inflow_store(j); }
		lbm.lbm_domain[0]->inlet.write_to_device(); // both planes
		inflow_j_dev = j0+1ull;
		inflow_launch(inflow_j_dev+1ull);
		print_info("inflow: turbulent, "+to_string((uint)inflow.n)+" modes from INFLOW (table for "+to_string((uint)inflow.nb)+" cells across b), a new plane every "+to_string(inflow_k)+" steps"+(inflow_step0>0ull ? ", playback starts at step "+to_string(inflow_step0) : ""));
	}
	const auto inflow_advance = [&]() { // make sure the plane after the current interval is on the device; call before running steps
		const ulong j = (lbm.get_t()+inflow_step0)/(ulong)inflow_k; // the coming steps blend planes j and j+1
		while(inflow_j_dev<j+1ull) {
			{ Clock wait; inflow_next.get(); inflow_wait_s += wait.stop(); inflow_uploads++; } // inflow_stage now holds plane inflow_j_dev+1; the device is idle between two run() calls
			inflow_j_dev++;
			inflow_store(inflow_j_dev);
			lbm.lbm_domain[0]->inlet.write_to_device((inflow_j_dev%2ull)*3ull*Ax, 3ull*Ax); // only the slot that changed
			inflow_launch(inflow_j_dev+1ull);
		}
	};
	const auto inflow_dump_plane = [&](const ulong sp) { // test hook: what the device holds at x = 0 after the time step of playback step sp (lattice units, mean profile included; 3*nz*ny float32, component-major, cell a = y + ny*z)
		Memory<float>& u = lbm.lbm_domain[0]->u;
		u.read_from_device(); // a no-op on zero-copy devices
		const float* d = u.data();
		vector<float> a(3ull*Ax);
		for(ulong k=0ull; k<3ull; k++) for(ulong i=0ull; i<Ax; i++) a[k*Ax+i] = d[k*N+i*(ulong)Nx];
		create_folder(out+"planes/");
		write_f32(out+"planes/inlet-"+to_string(sp)+".f32", a);
	};
	if(!xst.empty()) { // the planes are gathered on the device at every sample, so only they cross to the host
		if(lbm.get_D()!=1u) fail("XST works on a single device only");
		lbm.lbm_domain[0]->enable_xst(xst);
	}
	print_info("wind setup: "+to_string(Nx)+"x"+to_string(Ny)+"x"+to_string(Nz)+" cells, dx = "+to_string(dx, 1u)+" m, dt = "+to_string(dt, 4u)+" s, steps = "+to_string(steps));

	// cameras: azimuth a from +x towards +y, elevation b; the camera looks back at the target
	const float3 c = float3(0.5f*(float)Nx, 0.5f*(float)Ny, 0.5f*(float)Nz);
	const auto look_at = [&](const float3& target, const float a, const float b, const float dist, const float fov) {
		const float ar = a*pif/180.0f, br = b*pif/180.0f;
		const float3 rz = float3(cos(ar)*cos(br), sin(ar)*cos(br), sin(br));
		lbm.graphics.set_camera_free(target-c+dist*rz, a, b, fov);
	};
	const float3 landmark = float3((float)px, (float)py, 40.0f/dx); // camera target: the landmark cell, 40 m above ground
	const auto render = [&](const string& tag) {
		lbm.graphics.slice_mode = 0;
		lbm.graphics.visualization_modes = VIS_FLAG_SURFACE|VIS_Q_CRITERION;
		look_at(landmark+float3(env_float("QX", 60.0f)/dx, 0.0f, 0.0f), env_float("QA", -115.0f), env_float("QB", 28.0f), env_float("QD", 1500.0f)/dx, env_float("QF", 40.0f));
		lbm.graphics.write_frame(out+"q3d/", "q3d-"+tag);
		look_at(landmark+float3(-150.0f/dx, 0.0f, 0.0f), env_float("WA", -150.0f), env_float("WB", 32.0f), env_float("WD", 3400.0f)/dx, env_float("WF", 40.0f));
		lbm.graphics.write_frame(out+"q3d_wide/", "q3d_wide-"+tag);
		lbm.graphics.visualization_modes = VIS_FIELD; // section only, solids left blank like a drawn section
		lbm.graphics.slice_mode = 2; lbm.graphics.slice_y = (int)py;
		look_at(landmark+float3(150.0f/dx, 0.0f, 60.0f/dx), -90.0f, 0.0f, env_float("SD", 1700.0f)/dx, env_float("SF", 40.0f));
		lbm.graphics.write_frame(out+"side/", "side-"+tag);
	};

	// planes read from the device: speed |u| in lattice units
	vector<float> hz_inst(NxNy), hz_ux(NxNy, 0.0f), hz_uy(NxNy, 0.0f), hi_inst(NxNy);
	vector<double> hz_mean(NxNy, 0.0), hi_mean(NxNy, 0.0), hz_sq(NxNy, 0.0), hi_sq(NxNy, 0.0); // sums and sums of squares of the speed in double: the standard deviation is a small difference of two large sums
	vector<float> vy_inst((ulong)Nx*Nz), vy_mean((ulong)Nx*Nz, 0.0f);
	vector<double> zp_mean(NxNy*zpeds.size(), 0.0), zp_sq(NxNy*zpeds.size(), 0.0); vector<float> zp_ux(NxNy*zpeds.size(), 0.0f), zp_uy(NxNy*zpeds.size(), 0.0f); // ZPEDS: per listed layer the four sums kept for ZPED, with the same types and the same arithmetic
	uint samples = 0u;
	const ulong NyNz = (ulong)Ny*Nz;
	vector<double> xsum(7ull*NyNz*xst.size(), 0.0); // per plane: sums of u, v, w, uu, vv, ww, uw (lattice units), each Ny*Nz with y fastest; double because a variance is a small difference of two large sums
	vector<double> zprof; uint zprof_records = 0u; // test hook: per record and layer z the sums of rho, rho*ux, rho*uy, rho*uz over the cells that are neither solid nor fixed
	const bool vol = vol_set;
	const ulong vnx=vol ? (ulong)(vx1-vx0) : 0ull, vny=vol ? (ulong)(vy1-vy0) : 0ull, vnz=(ulong)vz1, vn=vnx*vny*vnz;
	vector<float> vmean(3ull*vn, 0.0f), vrho(vn, 0.0f); // mean velocity and mean density (pressure) in the box
	vector<double> vsum2(vol_m2==1u ? 7ull*vn : 0ull, 0.0); // VOL_M2: sums of u, v, w, uu, vv, ww, uw per box cell (lattice units); double because a variance is a small difference of two large sums
	uint vsamples = 0u;
	const auto read_box = [&](const bool with_rho) { // bring the box to the host: the slab z < VZ1 is contiguous in n = x + Nx*(y + Ny*z); a no-op on zero-copy devices, required on discrete GPUs
		Memory<float>& u = lbm.lbm_domain[0]->u;
		for(uint k=0u; k<3u; k++) u.read_from_device((ulong)k*N, vnz*NxNy);
		if(with_rho) lbm.lbm_domain[0]->rho.read_from_device(0ull, vnz*NxNy);
	};
	const auto read_planes = [&](const bool accumulate) {
		lbm.update_fields(); // without INTERACTIVE_GRAPHICS the velocity field is not refreshed every step
		Memory<float>& u = lbm.lbm_domain[0]->u;
		for(uint k=0u; k<3u; k++) for(const uint zz : {zped, zhi}) u.read_from_device((ulong)k*N+(ulong)zz*NxNy, NxNy);
		for(uint k=0u; k<3u; k++) for(uint z=0u; z<Nz; z++) u.read_from_device((ulong)k*N+(ulong)z*NxNy+(ulong)py*Nx, Nx);
		if(accumulate) for(uint k=0u; k<3u; k++) for(const uint zz : zpeds) u.read_from_device((ulong)k*N+(ulong)zz*NxNy, NxNy);
		const float* d = u.data();
		if(accumulate&&!xst.empty()) { // the planes are not contiguous in the velocity field: a kernel copies them into a compact buffer and only that crosses to the host (until 2026-11-10 the whole field was read here; the numbers are the same)
			lbm.lbm_domain[0]->gather_planes();
			const float* g = lbm.lbm_domain[0]->xst_buf.data(); const ulong T3 = (ulong)xst.size()*NyNz; // component k of cell i of plane s: g[k*T3+s*NyNz+i]
			parallel_for(T3, max(min((uint)thread::hardware_concurrency(), 8u), 1u), [&](ulong j) { const ulong s = j/NyNz, i = j%NyNz; double* q = xsum.data()+7ull*NyNz*s;
				const double ux=g[j], uy=g[T3+j], uz=g[2ull*T3+j];
				q[i] += ux; q[NyNz+i] += uy; q[2ull*NyNz+i] += uz; q[3ull*NyNz+i] += ux*ux; q[4ull*NyNz+i] += uy*uy; q[5ull*NyNz+i] += uz*uz; q[6ull*NyNz+i] += ux*uz; });
		}
		for(ulong i=0ull; i<NxNy; i++) {
			const ulong a = (ulong)zped*NxNy+i, b = (ulong)zhi*NxNy+i;
			const float ux=d[a], uy=d[N+a], uz=d[2ull*N+a];
			hz_inst[i] = sqrt(ux*ux+uy*uy+uz*uz);
			hi_inst[i] = sqrt(sq(d[b])+sq(d[N+b])+sq(d[2ull*N+b]));
			if(accumulate) { const double qz=hz_inst[i], qi=hi_inst[i]; hz_mean[i] += qz; hz_sq[i] += qz*qz; hi_mean[i] += qi; hi_sq[i] += qi*qi; hz_ux[i] += ux; hz_uy[i] += uy; }
		}
		if(accumulate) for(ulong l=0ull; l<zpeds.size(); l++) { // the sums of ZPED once more for every listed layer
			double* m = zp_mean.data()+l*NxNy; double* m2 = zp_sq.data()+l*NxNy; float* sx = zp_ux.data()+l*NxNy; float* sy = zp_uy.data()+l*NxNy;
			for(ulong i=0ull; i<NxNy; i++) {
				const ulong a = (ulong)zpeds[l]*NxNy+i;
				const float ux=d[a], uy=d[N+a], uz=d[2ull*N+a];
				const float speed = sqrt(ux*ux+uy*uy+uz*uz);
				const double qz = speed; m[i] += qz; m2[i] += qz*qz; sx[i] += ux; sy[i] += uy;
			}
		}
		for(uint z=0u; z<Nz; z++) for(uint x=0u; x<Nx; x++) {
			const ulong a = (ulong)z*NxNy+(ulong)py*Nx+x, j = (ulong)z*Nx+x;
			vy_inst[j] = sqrt(sq(d[a])+sq(d[N+a])+sq(d[2ull*N+a]));
			if(accumulate) vy_mean[j] += vy_inst[j];
		}
		if(accumulate&&vol) {
			read_box(true);
			const float* r = lbm.lbm_domain[0]->rho.data();
			parallel_for(vn, [&](ulong i) { const ulong x=vx0+i%vnx, y=vy0+(i/vnx)%vny, z=i/(vnx*vny), n=x+y*(ulong)Nx+z*NxNy;
				vmean[i] += d[n]; vmean[vn+i] += d[N+n]; vmean[2ull*vn+i] += d[2ull*N+n]; vrho[i] += r[n];
				if(vol_m2==1u) { const double ux=d[n], uy=d[N+n], uz=d[2ull*N+n]; double* q = vsum2.data()+i;
					q[0] += ux; q[vn] += uy; q[2ull*vn] += uz; q[3ull*vn] += ux*ux; q[4ull*vn] += uy*uy; q[5ull*vn] += uz*uz; q[6ull*vn] += ux*uz; } });
			vsamples++;
		}
		if(accumulate) samples++;
	};
	const ulong t_avg = getenv("SPIN_SI") ? (ulong)(spin_si/dt+0.5f) : steps/2ull; // samples taken at t > t_avg are averaged: after SPIN_SI seconds, or the second half of the run
	const float anim_t = env_float("ANIM_T", 0.0f), anim_si = env_float("ANIM_SI", 1.0f);
	const ulong t_anim = anim_t>0.0f ? steps-min(steps, (ulong)(anim_t/dt)) : steps+1ull, anim_every = max((ulong)(anim_si/dt), 1ull);
	const auto write_planes = [&](const string& tag) {
		const string p = out+"planes/";
		create_folder(p);
		write_f32(p+"hz_inst-"+tag+".f32", hz_inst); write_f32(p+"hi_inst-"+tag+".f32", hi_inst); write_f32(p+"vy_inst-"+tag+".f32", vy_inst);
		{ std::ofstream ms(p+"samples-"+tag+".txt"); ms << samples << "\n"; } // 0 for stages that end before the averaging window starts
		if(samples>0u) {
			const double s = 1.0/(double)samples; const float sf = 1.0f/(float)samples;
			vector<float> a(NxNy), d(NxNy), sz(NxNy), si(NxNy), b(hz_ux), cc(hz_uy), e(vy_mean);
			for(ulong i=0ull; i<NxNy; i++) {
				const double mz=hz_mean[i]*s, mi=hi_mean[i]*s;
				a[i] = (float)mz; d[i] = (float)mi;
				sz[i] = (float)sqrt(fmax(hz_sq[i]*s-mz*mz, 0.0)); si[i] = (float)sqrt(fmax(hi_sq[i]*s-mi*mi, 0.0)); // standard deviation of the instantaneous speed over the averaging window
			}
			for(float& v : b) v *= sf; for(float& v : cc) v *= sf; for(float& v : e) v *= sf;
			write_f32(p+"hz_mean.f32", a); write_f32(p+"hz_mean_ux.f32", b); write_f32(p+"hz_mean_uy.f32", cc); write_f32(p+"hi_mean.f32", d); write_f32(p+"vy_mean.f32", e);
			write_f32(p+"hz_mean-"+tag+".f32", a); write_f32(p+"hi_mean-"+tag+".f32", d); // tagged copies: the mean so far, to check that it has settled
			write_f32(p+"hz_std.f32", sz); write_f32(p+"hi_std.f32", si);
			for(ulong l=0ull; l<zpeds.size(); l++) { // ZPEDS: the planes of ZPED for layer k as hz<k>_...
				const string hk = p+"hz"+to_string(zpeds[l]);
				vector<float> za(NxNy), zs(NxNy), zb(zp_ux.begin()+(long)(l*NxNy), zp_ux.begin()+(long)((l+1ull)*NxNy)), zc(zp_uy.begin()+(long)(l*NxNy), zp_uy.begin()+(long)((l+1ull)*NxNy));
				for(ulong i=0ull; i<NxNy; i++) {
					const double mz = zp_mean[l*NxNy+i]*s;
					za[i] = (float)mz; zs[i] = (float)sqrt(fmax(zp_sq[l*NxNy+i]*s-mz*mz, 0.0));
				}
				for(float& v : zb) v *= sf; for(float& v : zc) v *= sf;
				write_f32(hk+"_mean.f32", za); write_f32(hk+"_mean_ux.f32", zb); write_f32(hk+"_mean_uy.f32", zc); write_f32(hk+"_std.f32", zs); write_f32(hk+"_mean-"+tag+".f32", za);
			}
			if(!xst.empty()) { // the sums so far, one file per frame: differences of two frames give the statistics of the time between them
				{ std::ofstream xf(p+"xst-"+tag+".f64", std::ios::binary); xf.write((const char*)xsum.data(), (std::streamsize)(xsum.size()*sizeof(double))); }
				std::ofstream xm(p+"xst_meta-"+tag+".txt");
				xm << "planes"; for(const uint x : xst) xm << " " << x;
				xm << "\nny " << Ny << "\nnz " << Nz << "\nsamples " << samples << "\nlayout per plane 7 arrays of float64 sums (u v w uu vv ww uw, lattice units), each ny*nz with y fastest then z\n";
			}
		}
		if(vol&&vsamples>0u) {
			read_box(false); // the last velocity in the box
			const float sv = 1.0f/(float)vsamples; vector<float> vm(vmean), vi(3ull*vn);
			for(float& v : vm) v *= sv;
			const float* d = lbm.lbm_domain[0]->u.data();
			parallel_for(vn, [&](ulong i) { const ulong x=vx0+i%vnx, y=vy0+(i/vnx)%vny, z=i/(vnx*vny), n=x+y*(ulong)Nx+z*NxNy; vi[i]=d[n]; vi[vn+i]=d[N+n]; vi[2ull*vn+i]=d[2ull*N+n]; });
			vector<float> vr(vrho); for(float& v : vr) v *= sv;
			write_f32(p+"vol_mean.f32", vm); write_f32(p+"vol_inst-"+tag+".f32", vi); write_f32(p+"vol_rho_mean.f32", vr);
			if(vol_tagged==1u) write_f32(p+"vol_mean-"+tag+".f32", vm); // the mean so far, kept per frame
			if(vol_m2==1u) { // standard deviations of u, v, w and the covariance of u and w over the samples so far: four arrays of the box, lattice units
				vector<float> vs(4ull*vn); const double s2 = 1.0/(double)vsamples;
				parallel_for(vn, [&](ulong i) { const double* q = vsum2.data()+i; const double mu=q[0]*s2, mv=q[vn]*s2, mw=q[2ull*vn]*s2;
					vs[i] = (float)sqrt(fmax(q[3ull*vn]*s2-mu*mu, 0.0)); vs[vn+i] = (float)sqrt(fmax(q[4ull*vn]*s2-mv*mv, 0.0)); vs[2ull*vn+i] = (float)sqrt(fmax(q[5ull*vn]*s2-mw*mw, 0.0));
					vs[3ull*vn+i] = (float)(q[6ull*vn]*s2-mu*mw); });
				write_f32(p+"vol_std.f32", vs);
			}
			for(const string& name : vol_tagged==1u ? vector<string>{"vol_meta.txt", "vol_meta-"+tag+".txt"} : vector<string>{"vol_meta.txt"}) {
				std::ofstream vmeta(p+name);
				vmeta << "x0 " << vx0 << "\ny0 " << vy0 << "\nnx " << vnx << "\nny " << vny << "\nnz " << vnz << "\nsamples " << vsamples << "\nlayout component-major, x fastest, then y, then z\n";
			}
		}
		if(zprof_records>0u) {
			{ std::ofstream zf(p+"zprof.f64", std::ios::binary); zf.write((const char*)zprof.data(), (std::streamsize)(zprof.size()*sizeof(double))); }
			std::ofstream zm(p+"zprof_meta.txt");
			zm << "nz " << Nz << "\nrecords " << zprof_records << "\nevery_steps " << zprof_every << "\nlayout per record and layer z four float64: sums of rho, rho*ux, rho*uy, rho*uz over cells that are neither solid nor fixed (lattice units); record k is taken after k*every_steps steps\n";
		}
		std::ofstream m(p+"meta.txt");
		m << "nx " << Nx << "\nny " << Ny << "\nnz " << Nz << "\ndx_m " << dx << "\nu_ref_lattice " << u_ref << "\nu_ref_si " << u_si << "\nz_ref_m " << z_ref*dx
		  << "\nzped_cell " << zped << "\nzhi_cell " << zhi << "\npy " << py << "\npx " << px << "\nsamples " << samples << "\nt_si " << (float)lbm.get_t()*dt
		  << "\nsample_steps " << sample << "\nsample_period_s " << (float)sample*dt << "\navg_start_s " << (float)t_avg*dt << "\navg_end_s " << (float)steps*dt
		  << "\nanim_start_s " << (anim_t>0.0f ? (float)t_anim*dt : 0.0f) << "\nanim_period_s " << (anim_t>0.0f ? (float)anim_every*dt : 0.0f) << "\n";
		if(!zpeds.empty()) { m << "zpeds_cells "; for(ulong l=0ull; l<zpeds.size(); l++) m << (l>0ull ? "," : "") << zpeds[l]; m << "\n"; } // only with ZPEDS: without it the file is what it was
	};

	uint anim_k = 0u;
	const auto anim_frame = [&]() {
		read_planes(false);
		create_folder(out+"anim/");
		char name[32]; snprintf(name, sizeof(name), "hi-%04u.f32", anim_k);
		write_f32(out+"anim/"+string(name), hi_inst);
		if(anim_section==1u) { // u and w on the plane y = PY, from the rows read_planes() has just brought to the host
			const float* d = lbm.lbm_domain[0]->u.data(); vector<float> su((ulong)Nx*Nz), sw((ulong)Nx*Nz);
			for(uint z=0u; z<Nz; z++) for(uint x=0u; x<Nx; x++) { const ulong a = (ulong)z*NxNy+(ulong)py*Nx+x, j = (ulong)z*Nx+x; su[j] = d[a]; sw[j] = d[2ull*N+a]; }
			snprintf(name, sizeof(name), "sec_u-%04u.f32", anim_k); write_f32(out+"anim/"+string(name), su);
			snprintf(name, sizeof(name), "sec_w-%04u.f32", anim_k); write_f32(out+"anim/"+string(name), sw);
		}
		if(render_on==1u) {
			lbm.graphics.slice_mode = 0; lbm.graphics.visualization_modes = VIS_FLAG_SURFACE|VIS_Q_CRITERION;
			look_at(landmark+float3(env_float("QX", 60.0f)/dx, 0.0f, 0.0f), env_float("QA", -115.0f), env_float("QB", 28.0f), env_float("QD", 1500.0f)/dx, env_float("QF", 40.0f));
			snprintf(name, sizeof(name), "q-%04u", anim_k);
			lbm.graphics.write_frame(out+"anim/", string(name));
		}
		anim_k++;
	};
	const auto record_zprof = [&]() { // the fields as the next time step will see them (after streaming and the boundary conditions, before the collision)
		lbm.update_fields();
		lbm.lbm_domain[0]->u.read_from_device(); lbm.lbm_domain[0]->rho.read_from_device(); // no-ops on zero-copy devices
		const float* d = lbm.lbm_domain[0]->u.data(); const float* r = lbm.lbm_domain[0]->rho.data();
		const ulong base = zprof.size(); zprof.resize(base+4ull*Nz, 0.0);
		for(ulong n=0ull; n<N; n++) {
			if(lbm.flags[n]&(TYPE_S|TYPE_E)) continue;
			const ulong o = base+4ull*(n/NxNy); const double rr = r[n];
			zprof[o] += rr; zprof[o+1ull] += rr*d[n]; zprof[o+2ull] += rr*d[N+n]; zprof[o+3ull] += rr*d[2ull*N+n];
		}
		zprof_records++;
	};
	lbm.run(0u);
	if(flags_dump==1u) { // test hook: the flags as the device holds them once the simulation is initialised
		lbm.flags.read_from_device();
		vector<uchar> fl(N); for(ulong n=0ull; n<N; n++) fl[n] = lbm.flags[n];
		create_folder(out+"planes/");
		std::ofstream ff(out+"planes/flags.u8", std::ios::binary); ff.write((const char*)fl.data(), (std::streamsize)N);
	}
	if(zprof_every>0u) record_zprof();
	Clock clock;
	for(uint k=1u; k<=frames; k++) {
		const ulong target = steps*k/frames;
		while(lbm.get_t()<target) {
			ulong chunk = min(sample, target-lbm.get_t());
			if(zprof_every>0u) chunk = min(chunk, (ulong)zprof_every-lbm.get_t()%(ulong)zprof_every); // stop at every record step
			if(!inflow_on) lbm.run(chunk);
			else for(const ulong stop=lbm.get_t()+chunk; lbm.get_t()<stop;) { // the same steps in parts that end at every inlet update, so the samples below are taken at the same steps as without INFLOW
				inflow_advance(); // the plane this interval blends towards must be on the device
				const ulong sp = lbm.get_t()+inflow_step0; // playback step of the next time step
				ulong part = min(stop-lbm.get_t(), (ulong)inflow_k-sp%(ulong)inflow_k); // stop at every update step
				for(const uint s : inflow_dump) if((ulong)s>=sp) part = min(part, (ulong)s+1ull-sp); // and right after every step that is dumped
				lbm.run(part);
				for(const uint s : inflow_dump) if((ulong)s+1ull==lbm.get_t()+inflow_step0) inflow_dump_plane((ulong)s);
			}
			if(zprof_every>0u&&lbm.get_t()%(ulong)zprof_every==0ull) record_zprof();
			if(lbm.get_t()>t_avg) read_planes(true);
			if(lbm.get_t()>=t_anim&&(lbm.get_t()-t_anim)/anim_every>=(ulong)anim_k) anim_frame();
		}
		print_info("t = "+to_string(lbm.get_t())+" steps ("+to_string((float)lbm.get_t()*dt, 0u)+" s of wind), wall time "+to_string(clock.stop(), 0u)+" s, samples "+to_string(samples));
		read_planes(false);
		write_planes(to_string(k));
		if(render_on==1u) render(to_string(k));
	}
	if(inflow_on) print_info("inflow: "+to_string(inflow_uploads)+" planes uploaded during the run; the solver waited "+to_string((float)inflow_wait_s, 2u)+" s for them in total ("+to_string(1000.0f*(float)inflow_wait_s/(float)max(inflow_uploads, 1ull), 2u)+" ms per plane)");
}
