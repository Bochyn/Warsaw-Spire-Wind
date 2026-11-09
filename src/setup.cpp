// Altered from FluidX3D commit b85d64e (v3.8) for a wind study of Warsaw Spire and Łucka street, 2026-11; see README.md of the published repository.
#include "setup.hpp"
#include <cerrno>
#include <cstdlib>
#include <fstream>

// Wind over a city: the solid is a 2.5D height map (cells above flat ground) plus an optional STL, the wind comes from the west (+x).
// Environment: HM HX HY NZ DX_M STL UREF ZREF_M ALPHA NU T_SI SAMPLE_SI FRAMES OUT PX PY ZPED ZHI, cameras QA QB QD QF QX / WA WB WD WF / SD SF,
// box VX0 VX1 VY0 VY1 VZ1, animation ANIM_T ANIM_SI; FX_U_MAX and FX_Q are read in lbm.cpp. Every variable is parsed strictly and checked in main_setup().
static void fail(const string& msg) { print_error("wind setup: "+msg); } // print_error exits
static bool parse_float(const char* s, float& v) { char* e=nullptr; errno=0; v=(float)strtod(s, &e); return e!=s&&*e=='\0'&&errno==0; }
static bool parse_uint(const char* s, uint& v) { char* e=nullptr; errno=0; const unsigned long x=strtoul(s, &e, 10); v=(uint)x; return s[0]>='0'&&s[0]<='9'&&*e=='\0'&&errno==0&&x<=0xFFFFFFFFul; }
static uint env_uint(const char* k, const uint d) { const char* v = getenv(k); uint r=d; if(v&&!parse_uint(v, r)) fail(string(k)+"=\""+v+"\" is not a non-negative integer"); return r; }
static float env_float(const char* k, const float d) { const char* v = getenv(k); float r=d; if(v&&!parse_float(v, r)) fail(string(k)+"=\""+v+"\" is not a number"); return r; }
static string env_str(const char* k, const string& d) { const char* v = getenv(k); return v ? string(v) : d; }
static void write_f32(const string& path, const vector<float>& a) { std::ofstream f(path, std::ios::binary); f.write((const char*)a.data(), (std::streamsize)(a.size()*sizeof(float))); }

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
	if(frames==0u) fail("FRAMES must be > 0");
	if(zped>=Nz||zhi>=Nz) fail("ZPED and ZHI must be < NZ");
	if(px>=HX||py>=HY) fail("PX must be < HX and PY must be < HY");
	if(vol_set&&!(vx0<vx1&&vx1<=HX&&vy0<vy1&&vy1<=HY&&vz1>0u&&vz1<=Nz)) fail("the box needs VX0 < VX1 <= HX, VY0 < VY1 <= HY and 0 < VZ1 <= NZ");

	vector<ushort> H((ulong)HX*HY);
	{ std::ifstream f(hm, std::ios::binary); f.read((char*)H.data(), (std::streamsize)(H.size()*sizeof(ushort))); if(!f) fail("cannot read "+to_string((ulong)(H.size()*sizeof(ushort)))+" bytes of height map HM=\""+hm+"\""); }

	LBM lbm(HX, HY, Nz, nu);
	const uint Nx = lbm.get_Nx(), Ny = lbm.get_Ny();
	const ulong N = lbm.get_N(), NxNy = (ulong)Nx*Ny;
	const auto profile = [&](const float z) { return u_ref*pow(fmax(z-0.5f, 0.5f)/z_ref, alpha); }; // wall sits at z=0.5
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
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==Nz-1u) lbm.flags[n] = TYPE_E; // inflow, outflow, sides, top
	});
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
	uint samples = 0u;
	const bool vol = vol_set;
	const ulong vnx=vol ? (ulong)(vx1-vx0) : 0ull, vny=vol ? (ulong)(vy1-vy0) : 0ull, vnz=(ulong)vz1, vn=vnx*vny*vnz;
	vector<float> vmean(3ull*vn, 0.0f), vrho(vn, 0.0f); // mean velocity and mean density (pressure) in the box
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
		const float* d = u.data();
		for(ulong i=0ull; i<NxNy; i++) {
			const ulong a = (ulong)zped*NxNy+i, b = (ulong)zhi*NxNy+i;
			const float ux=d[a], uy=d[N+a], uz=d[2ull*N+a];
			hz_inst[i] = sqrt(ux*ux+uy*uy+uz*uz);
			hi_inst[i] = sqrt(sq(d[b])+sq(d[N+b])+sq(d[2ull*N+b]));
			if(accumulate) { const double qz=hz_inst[i], qi=hi_inst[i]; hz_mean[i] += qz; hz_sq[i] += qz*qz; hi_mean[i] += qi; hi_sq[i] += qi*qi; hz_ux[i] += ux; hz_uy[i] += uy; }
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
				vmean[i] += d[n]; vmean[vn+i] += d[N+n]; vmean[2ull*vn+i] += d[2ull*N+n]; vrho[i] += r[n]; });
			vsamples++;
		}
		if(accumulate) samples++;
	};
	const ulong t_avg = steps/2ull; // samples taken at t > t_avg are averaged: the second half of the run
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
		}
		if(vol&&vsamples>0u) {
			read_box(false); // the last velocity in the box
			const float sv = 1.0f/(float)vsamples; vector<float> vm(vmean), vi(3ull*vn);
			for(float& v : vm) v *= sv;
			const float* d = lbm.lbm_domain[0]->u.data();
			parallel_for(vn, [&](ulong i) { const ulong x=vx0+i%vnx, y=vy0+(i/vnx)%vny, z=i/(vnx*vny), n=x+y*(ulong)Nx+z*NxNy; vi[i]=d[n]; vi[vn+i]=d[N+n]; vi[2ull*vn+i]=d[2ull*N+n]; });
			vector<float> vr(vrho); for(float& v : vr) v *= sv;
			write_f32(p+"vol_mean.f32", vm); write_f32(p+"vol_inst-"+tag+".f32", vi); write_f32(p+"vol_rho_mean.f32", vr);
			std::ofstream vmeta(p+"vol_meta.txt");
			vmeta << "x0 " << vx0 << "\ny0 " << vy0 << "\nnx " << vnx << "\nny " << vny << "\nnz " << vnz << "\nsamples " << vsamples << "\nlayout component-major, x fastest, then y, then z\n";
		}
		std::ofstream m(p+"meta.txt");
		m << "nx " << Nx << "\nny " << Ny << "\nnz " << Nz << "\ndx_m " << dx << "\nu_ref_lattice " << u_ref << "\nu_ref_si " << u_si << "\nz_ref_m " << z_ref*dx
		  << "\nzped_cell " << zped << "\nzhi_cell " << zhi << "\npy " << py << "\npx " << px << "\nsamples " << samples << "\nt_si " << (float)lbm.get_t()*dt
		  << "\nsample_steps " << sample << "\nsample_period_s " << (float)sample*dt << "\navg_start_s " << (float)t_avg*dt << "\navg_end_s " << (float)steps*dt
		  << "\nanim_start_s " << (anim_t>0.0f ? (float)t_anim*dt : 0.0f) << "\nanim_period_s " << (anim_t>0.0f ? (float)anim_every*dt : 0.0f) << "\n";
	};

	uint anim_k = 0u;
	const auto anim_frame = [&]() {
		read_planes(false);
		create_folder(out+"anim/");
		char name[32]; snprintf(name, sizeof(name), "hi-%04u.f32", anim_k);
		write_f32(out+"anim/"+string(name), hi_inst);
		lbm.graphics.slice_mode = 0; lbm.graphics.visualization_modes = VIS_FLAG_SURFACE|VIS_Q_CRITERION;
		look_at(landmark+float3(env_float("QX", 60.0f)/dx, 0.0f, 0.0f), env_float("QA", -115.0f), env_float("QB", 28.0f), env_float("QD", 1500.0f)/dx, env_float("QF", 40.0f));
		snprintf(name, sizeof(name), "q-%04u", anim_k);
		lbm.graphics.write_frame(out+"anim/", string(name));
		anim_k++;
	};
	lbm.run(0u);
	Clock clock;
	for(uint k=1u; k<=frames; k++) {
		const ulong target = steps*k/frames;
		while(lbm.get_t()<target) {
			lbm.run(min(sample, target-lbm.get_t()));
			if(lbm.get_t()>t_avg) read_planes(true);
			if(lbm.get_t()>=t_anim&&(lbm.get_t()-t_anim)/anim_every>=(ulong)anim_k) anim_frame();
		}
		print_info("t = "+to_string(lbm.get_t())+" steps ("+to_string((float)lbm.get_t()*dt, 0u)+" s of wind), wall time "+to_string(clock.stop(), 0u)+" s, samples "+to_string(samples));
		read_planes(false);
		write_planes(to_string(k));
		render(to_string(k));
	}
}
