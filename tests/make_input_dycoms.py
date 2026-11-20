# Builds data/input_dycoms.nc for the CI tests, in the format and on the grid of data/input_marine_StCu.nc.
# Run from tests/:  python make_input_dycoms.py
# DYCOMS-II RF01 nocturnal stratocumulus (Stevens et al. 2005, MWR 133, 1443): z_i 840 m;
# below: theta_l 289.0 K, q_t 9.0 g/kg; above: theta_l 297.5 + (z-z_i)^(1/3), q_t 1.5 g/kg;
# p_s 1017.8 hPa, SST 292.5 K, SHF 15 W m-2, LHF 115 W m-2, u 7, v -5.5 m/s.
# Longwave heating from the RF01 parameterized net flux on the initial cloud (F0 70, F1 22 W m-2,
# kappa 85 m2 kg-1, D 3.75e-6 s-1 in the free-troposphere term), held fixed in time as RTHRATEN.
# Above z_i q_t is capped at 50 % RH so the 28 km column does not saturate aloft.
import netCDF4 as nc, numpy as np, shutil, os
D = os.path.dirname(os.path.abspath(__file__))
src, dst = f"{D}/data/input_marine_StCu.nc", f"{D}/data/input_dycoms.nc"
shutil.copyfile(src, dst)
f = nc.Dataset(dst, 'r+')
f.TITLE = ' DYCOMS-II RF01 single column (Stevens et al. 2005, MWR 133, 1443) on the input_marine_StCu.nc grid'
f.history = ('initial profile, surface fluxes and longwave heating (RTHRATEN, held fixed) replaced by the DYCOMS-II RF01 case '
             'of Stevens et al. (2005) using make_input_dycoms.py; other fields and attributes are from input_marine_StCu.nc\n' + f.history)
Rd, cp, Rv, g, Lv, p0 = 287.0, 1004.0, 461.5, 9.81, 2.5e6, 1.0e5
dz = np.array(f['dz'][1].ravel(), float); nz = dz.size
zw = np.concatenate(([0.], np.cumsum(dz))); zm = zw[:-1] + 0.5*dz
zi = 840.
thl = np.where(zm < zi, 289.0, 297.5 + np.maximum(zm-zi, 0.)**(1./3.))
qt = np.where(zm < zi, 9.0e-3, 1.5e-3)
u = 7.0 + 0*zm; v = -5.5 + 0*zm
def qsat(T, p):
    es = 611.2*np.exp(17.67*(T-273.15)/(T-29.65))
    return 0.622*es/(p-es)
ps = 101780.0
p = np.zeros(nz); th = thl.copy(); qv = qt.copy(); qc = np.zeros(nz)
pbot = ps
for k in range(nz):
    for it in range(30):
        pk = pbot*np.exp(-g*(0.5*dz[k])/(Rd*th[k]*(pbot/p0)**(Rd/cp)*(1+0.61*qv[k]-qc[k])))
        ex = (pk/p0)**(Rd/cp)
        T = thl[k]*ex
        qs = qsat(T, pk)
        if zm[k] > zi: qt[k] = min(1.5e-3, 0.5*qs)   # free troposphere: RF01 1.5 g/kg, capped at 50 % RH aloft
        if qt[k] > qs:
            for n in range(20):
                qs = qsat(T, pk); dqs = qs*17.67*243.5/((T-29.65)**2)
                fT = T - thl[k]*ex - Lv/cp*(qt[k]-qs); T -= fT/(1+Lv/cp*dqs)
            qc[k] = qt[k]-qsat(T, pk); qv[k] = qt[k]-qc[k]
        else:
            qc[k] = 0.; qv[k] = qt[k]
        th[k] = T/ex
    p[k] = pk
    pbot = pk*np.exp(-g*(0.5*dz[k])/(Rd*T*(1+0.61*qv[k]-qc[k])))
ex = (p/p0)**(Rd/cp); T = th*ex
rho = p/(Rd*T*(1+0.61*qv-qc))
# RF01 net longwave flux on the w faces: F = F0 exp(-Q(z,top)) + F1 exp(-Q(0,z)) + free-troposphere term
lw = rho*qc*dz                                   # liquid water path per layer, kg m-2
Qab = np.concatenate((np.cumsum(lw[::-1])[::-1], [0.]))   # LWP above each face
Qbe = np.concatenate(([0.], np.cumsum(lw)))                # LWP below each face
kap, F0, F1, Dd = 85., 70., 22., 3.75e-6
rhoi = np.interp(zi, zm, rho)
dzi = np.maximum(zw-zi, 0.)
F = F0*np.exp(-kap*Qab) + F1*np.exp(-kap*Qbe) + rhoi*cp*Dd*(0.25*dzi**(4./3.) + zi*dzi**(1./3.))*(zw > zi)
dTdt = -(F[1:]-F[:-1])/(rho*cp*dz)               # K/s
rthr = dTdt/ex                                    # theta tendency, K/s
nt = f.dimensions['Time'].size
def put3(vn, a):
    for t in range(nt): f[vn][t,:,0,0] = a
for vn,a in (("dz",dz),("TH",th),("QVAPOR",qv),("QCLOUD",qc),("QICE",0*qc),("TK",T),("PTOTAL",p),("EXNER",ex),
            ("RHO",rho),("W_mass",0*qc),("EXCH_H_mass",0*qc),("EXCH_M_mass",0*qc),("RTHRATEN",rthr)):
    put3(vn, a)
if "CLDFRA_BL" in f.variables: put3("CLDFRA_BL", (qc > 0).astype(float))
for t in range(nt):
    f['U_mass'][t,:] = u; f['V_mass'][t,:] = v
rho0 = rho[0]
for vn,a in (("PSFC",ps),("TSK",292.5),("QSFC",qsat(292.5,ps)),("UST",0.25),("HFX",15.),
            ("QFX",115./Lv),("WSPD",np.hypot(7.,5.5)),("ZNT",2.0e-4),("XLAND",2.0),("UOCE",0.),("VOCE",0.),
            ("FLHC",rho0*cp*1.1e-3*np.hypot(7.,5.5))):
    f[vn][:] = a
f.close()
kc = np.where(qc > 0)[0]
print(f"dycoms: nz {nz}; cloud k {kc.min()+1}..{kc.max()+1} (z {zw[kc.min()]:.0f}-{zw[kc.max()+1]:.0f} m), qc max {qc.max()*1e3:.3f} g/kg, LWP {lw.sum()*1e3:.1f} g m-2")
for k in range(max(kc.max()-3,0), min(kc.max()+3,nz)): print(f"  k{k+1} zm {zm[k]:.0f} dz {dz[k]:.0f} qc {qc[k]*1e3:.3f} rthraten {rthr[k]:.2e} K/s ({rthr[k]*86400:.1f} K/day)")
