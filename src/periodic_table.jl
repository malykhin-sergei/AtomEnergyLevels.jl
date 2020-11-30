# Sober Physicists Don't Find Giraffes Hiding In Kitchens
atomic_shell = Dict(0 => 'S', 1 => 'P', 2 => 'D', 3 => 'F',
                    4 => 'G', 5 => 'H', 6 => 'I', 7 => 'K')

# Electronic configurations of the elements, taken from
# https://math.nist.gov/DFTdata/atomdata/configuration.html
atomic_electron_configuration = (
H  = conf_enc("1s1"),
He = conf_enc("[He]"),
Li = conf_enc("[He] 2s1"),
Be = conf_enc("[He] 2s2"),
B  = conf_enc("[He] 2s2 2p1"),
C  = conf_enc("[He] 2s2 2p2"),
N  = conf_enc("[He] 2s2 2p3"),
O  = conf_enc("[He] 2s2 2p4"),
F  = conf_enc("[He] 2s2 2p5"),
Ne = conf_enc("[Ne]"),
Na = conf_enc("[Ne] 3s1"),
Mg = conf_enc("[Ne] 3s2"),
Al = conf_enc("[Ne] 3s2 3p1"),
Si = conf_enc("[Ne] 3s2 3p2"),
P  = conf_enc("[Ne] 3s2 3p3"),
S  = conf_enc("[Ne] 3s2 3p4"),
Cl = conf_enc("[Ne] 3s2 3p5"),
Ar = conf_enc("[Ar]"),
K  = conf_enc("[Ar] 4s1"),
Ca = conf_enc("[Ar] 4s2"),
Sc = conf_enc("[Ar] 3d1 4s2"),
Ti = conf_enc("[Ar] 3d2 4s2"),
V  = conf_enc("[Ar] 3d3 4s2"),
Cr = conf_enc("[Ar] 3d5 4s1"),
Mn = conf_enc("[Ar] 3d5 4s2"),
Fe = conf_enc("[Ar] 3d6 4s2"),
Co = conf_enc("[Ar] 3d7 4s2"),
Ni = conf_enc("[Ar] 3d8 4s2"),
Cu = conf_enc("[Ar] 3d10 4s1"),
Zn = conf_enc("[Ar] 3d10 4s2"),
Ga = conf_enc("[Ar] 3d10 4s2 4p1"),
Ge = conf_enc("[Ar] 3d10 4s2 4p2"),
As = conf_enc("[Ar] 3d10 4s2 4p3"),
Se = conf_enc("[Ar] 3d10 4s2 4p4"),
Br = conf_enc("[Ar] 3d10 4s2 4p5"),
Kr = conf_enc("[Kr]"),
Rb = conf_enc("[Kr] 5s1"),
Sr = conf_enc("[Kr] 5s2"),
Y  = conf_enc("[Kr] 4d1 5s2"),
Zr = conf_enc("[Kr] 4d2 5s2"),
Nb = conf_enc("[Kr] 4d4 5s1"),
Mo = conf_enc("[Kr] 4d5 5s1"),
Tc = conf_enc("[Kr] 4d5 5s2"),
Ru = conf_enc("[Kr] 4d7 5s1"),
Rh = conf_enc("[Kr] 4d8 5s1"),
Pd = conf_enc("[Kr] 4d10"),
Ag = conf_enc("[Kr] 4d10 5s1"),
Cd = conf_enc("[Kr] 4d10 5s2"),
In = conf_enc("[Kr] 4d10 5s2 5p1"),
Sn = conf_enc("[Kr] 4d10 5s2 5p2"),
Sb = conf_enc("[Kr] 4d10 5s2 5p3"),
Te = conf_enc("[Kr] 4d10 5s2 5p4"),
I  = conf_enc("[Kr] 4d10 5s2 5p5"),
Xe = conf_enc("[Xe]"),
Cs = conf_enc("[Xe] 6s1"),
Ba = conf_enc("[Xe] 6s2"),
La = conf_enc("[Xe] 5d1 6s2"),
Ce = conf_enc("[Xe] 4f1 5d1 6s2"),
Pr = conf_enc("[Xe] 4f3 6s2"),
Nd = conf_enc("[Xe] 4f4 6s2"),
Pm = conf_enc("[Xe] 4f5 6s2"),
Sm = conf_enc("[Xe] 4f6 6s2"),
Eu = conf_enc("[Xe] 4f7 6s2"),
Gd = conf_enc("[Xe] 4f7 5d1 6s2"),
Tb = conf_enc("[Xe] 4f9 6s2"),
Dy = conf_enc("[Xe] 4f10 6s2"),
Ho = conf_enc("[Xe] 4f11 6s2"),
Er = conf_enc("[Xe] 4f12 6s2"),
Tm = conf_enc("[Xe] 4f13 6s2"),
Yb = conf_enc("[Xe] 4f14 6s2"),
Lu = conf_enc("[Xe] 4f14 5d1 6s2"),
Hf = conf_enc("[Xe] 4f14 5d2 6s2"),
Ta = conf_enc("[Xe] 4f14 5d3 6s2"),
W  = conf_enc("[Xe] 4f14 5d4 6s2"),
Re = conf_enc("[Xe] 4f14 5d5 6s2"),
Os = conf_enc("[Xe] 4f14 5d6 6s2"),
Ir = conf_enc("[Xe] 4f14 5d7 6s2"),
Pt = conf_enc("[Xe] 4f14 5d9 6s1"),
Au = conf_enc("[Xe] 4f14 5d10 6s1"),
Hg = conf_enc("[Xe] 4f14 5d10 6s2"),
Tl = conf_enc("[Xe] 4f14 5d10 6s2 6p1"),
Pb = conf_enc("[Xe] 4f14 5d10 6s2 6p2"),
Bi = conf_enc("[Xe] 4f14 5d10 6s2 6p3"),
Po = conf_enc("[Xe] 4f14 5d10 6s2 6p4"),
At = conf_enc("[Xe] 4f14 5d10 6s2 6p5"),
Rn = conf_enc("[Rn]"),
Fr = conf_enc("[Rn] 7s1"),
Ra = conf_enc("[Rn] 7s2"),
Ac = conf_enc("[Rn] 6d1 7s2"),
Th = conf_enc("[Rn] 6d2 7s2"),
Pa = conf_enc("[Rn] 5f2 6d1 7s2"),
 U = conf_enc("[Rn] 5f3 6d1 7s2"));
