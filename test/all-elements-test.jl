# This code calculates the total energies for all atoms Z = 1..92
# in the LDA(SVWN) DFT approximation and compares them with
# the NIST data available at https://www.nist.gov/document/dftdata-targz

ENIST = [   -0.445671;     -2.834836;     -7.335195;    -14.447209;    -24.344198;
           -37.425749;    -54.025016;    -74.473077;    -99.099648;   -128.233481;
          -161.440060;   -199.139406;   -241.315573;   -288.198397;   -339.946219;
          -396.716081;   -458.664179;   -525.946195;   -598.200590;   -675.742283;
          -758.679275;   -847.277216;   -941.678904;  -1042.030238;  -1148.449372;
         -1261.093056;  -1380.091264;  -1505.580197;  -1637.785861;  -1776.573850;
         -1921.846456;  -2073.807332;  -2232.534978;  -2398.111440;  -2570.620700;
         -2750.147940;  -2936.337293;  -3129.453161;  -3329.520604;  -3536.737751;
         -3751.196175;  -3973.013235;  -4202.188857;  -4438.981228;  -4683.301031;
         -4935.368406;  -5195.031215;  -5462.390982;  -5737.309064;  -6019.953353;
         -6310.376268;  -6608.631413;  -6914.773092;  -7228.856107;  -7550.557710;
         -7880.111578;  -8217.575230;  -8563.360285;  -8917.664369;  -9280.311037;
         -9651.484134; -10031.259090; -10419.710775; -10816.653877; -11222.941975;
        -11637.869664; -12061.770549; -12494.718304; -12936.786494; -13388.048594;
        -13848.230375; -14317.493720; -14795.884202; -15283.448822; -15780.236024;
        -16286.295408; -16801.677471; -17326.576377; -17860.790943; -18404.274220;
        -18956.957627; -19518.993145; -20090.414449; -20671.256539; -21261.555215;
        -21861.346869; -22470.319655; -23088.688083; -23716.496952; -24353.832231;
        -25001.291382; -25658.417889];

using AtomEnergyLevels
using Test, Printf

function test_all()
  N = length(ENIST)
  E = zeros(N) 
  for at_number=1:N
    element = keys(atomic_electron_configuration)[at_number]
    @info "Calculating $element:"
    E[at_number], ρ, ψ = lda(Z = at_number, β = 0.3)
    p = sortperm(ψ[3,:])
    ψ = ψ[:, p]
    @info @sprintf("Energy levels for the %s atom are:", element)
    for (i, ϵᵢ) in enumerate(ψ[3,:])
      l, nᵣ, nᵢ = ψ[1, i], ψ[2, i], ψ[4, i]
      n = nᵣ + l + 1
      @info @sprintf("\t%i%s\t(%4.1f)\t%14.6f", n, atomic_shell[l+1], nᵢ, ϵᵢ)
    end
    @info @sprintf("ΔE = %12.6f", ENIST[at_number]-E[at_number])
  end
  @info "== RESULTS SUMMARY =="
  for Z=1:N
    element = keys(atomic_electron_configuration)[Z]
    @info @sprintf("%2s: ΔE = %+e", element, ENIST[Z]-E[Z])
  end  
end

#=
[ Info: == RESULTS SUMMARY ==
[ Info:  H: ΔE = -4.825849e-07
[ Info: He: ΔE = -3.773020e-07
[ Info: Li: ΔE = +1.885772e-07
[ Info: Be: ΔE = +4.739539e-07
[ Info:  B: ΔE = +1.026720e-07
[ Info:  C: ΔE = -4.641473e-07
[ Info:  N: ΔE = +1.392992e-07
[ Info:  O: ΔE = -1.956511e-07
[ Info:  F: ΔE = +1.081396e-07
[ Info: Ne: ΔE = +2.685419e-07
[ Info: Na: ΔE = +3.087938e-07
[ Info: Mg: ΔE = +3.143902e-07
[ Info: Al: ΔE = +4.059358e-07
[ Info: Si: ΔE = -3.976903e-07
[ Info:  P: ΔE = -4.026344e-08
[ Info:  S: ΔE = -4.999658e-07
[ Info: Cl: ΔE = +4.915568e-07
[ Info: Ar: ΔE = -8.318011e-08
[ Info:  K: ΔE = -2.904915e-07
[ Info: Ca: ΔE = -3.886527e-07
[ Info: Sc: ΔE = +3.623470e-07
[ Info: Ti: ΔE = -1.119326e-07
[ Info:  V: ΔE = +3.112017e-07
[ Info: Cr: ΔE = -4.147864e-08
[ Info: Mn: ΔE = +4.454487e-07
[ Info: Fe: ΔE = -1.634728e-07
[ Info: Co: ΔE = +1.986248e-07
[ Info: Ni: ΔE = -3.038651e-07
[ Info: Cu: ΔE = -1.411345e-07
[ Info: Zn: ΔE = -3.286452e-07
[ Info: Ga: ΔE = -1.312649e-07
[ Info: Ge: ΔE = -6.353412e-08
[ Info: As: ΔE = +4.522440e-07
[ Info: Se: ΔE = +1.660305e-07
[ Info: Br: ΔE = -1.467338e-07
[ Info: Kr: ΔE = +4.020594e-07
[ Info: Rb: ΔE = +3.240161e-07
[ Info: Sr: ΔE = +3.584414e-07
[ Info:  Y: ΔE = -2.323027e-07
[ Info: Zr: ΔE = -1.045250e-07
[ Info: Nb: ΔE = -2.323941e-08
[ Info: Mo: ΔE = +2.575384e-07
[ Info: Tc: ΔE = +2.610768e-07
[ Info: Ru: ΔE = +1.888957e-07
[ Info: Rh: ΔE = -5.547354e-07
[ Info: Pd: ΔE = -3.359319e-07
[ Info: Ag: ΔE = -3.742452e-07
[ Info: Cd: ΔE = -4.515005e-08
[ Info: In: ΔE = -5.496395e-07
[ Info: Sn: ΔE = -1.398103e-07
[ Info: Sb: ΔE = -1.468561e-07
[ Info: Te: ΔE = +2.745874e-07
[ Info:  I: ΔE = -2.020943e-07
[ Info: Xe: ΔE = -5.579741e-07
[ Info: Cs: ΔE = -3.086016e-07
[ Info: Ba: ΔE = -7.467042e-08
[ Info: La: ΔE = -2.223151e-07
[ Info: Ce: ΔE = -1.469998e-07
[ Info: Pr: ΔE = -2.569850e-07
[ Info: Nd: ΔE = -3.546229e-07
[ Info: Pm: ΔE = +2.059005e-07
[ Info: Sm: ΔE = -4.220838e-07
[ Info: Eu: ΔE = -4.312333e-07
[ Info: Gd: ΔE = +1.065309e-07
[ Info: Tb: ΔE = -2.803372e-07
[ Info: Dy: ΔE = +3.004407e-07
[ Info: Ho: ΔE = +1.368790e-07
[ Info: Er: ΔE = -1.903045e-07
[ Info: Tm: ΔE = -3.280365e-07
[ Info: Yb: ΔE = +1.912740e-07
[ Info: Lu: ΔE = -5.482034e-07
[ Info: Hf: ΔE = +2.027809e-08
[ Info: Ta: ΔE = -3.143068e-07
[ Info:  W: ΔE = -1.371791e-07
[ Info: Re: ΔE = +2.109100e-07
[ Info: Os: ΔE = +2.784636e-07
[ Info: Ir: ΔE = -7.087801e-07
[ Info: Pt: ΔE = -3.284476e-07
[ Info: Au: ΔE = -5.660950e-07
[ Info: Hg: ΔE = -2.111010e-07
[ Info: Tl: ΔE = +1.368026e-07
[ Info: Pb: ΔE = -3.709865e-07
[ Info: Bi: ΔE = -4.356298e-07
[ Info: Po: ΔE = +1.854714e-07
[ Info: At: ΔE = -1.899498e-07
[ Info: Rn: ΔE = -2.945744e-07
[ Info: Fr: ΔE = -6.667688e-07
[ Info: Ra: ΔE = -2.192173e-07
[ Info: Ac: ΔE = -7.238996e-07
[ Info: Th: ΔE = -1.004701e-07
[ Info: Pa: ΔE = -9.101495e-08
[ Info:  U: ΔE = -3.833484e-07
=#
    