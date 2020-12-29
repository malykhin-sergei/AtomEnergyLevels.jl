# This code calculates the total energies for all atoms Z = 1..92
# in the LDA(SVWN) DFT approximation and compares them with
# the NIST data available at https://www.nist.gov/document/dftdata-targz

const ENIST = [-0.445671;     -2.834836;     -7.335195;    -14.447209;    -24.344198;
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

function test_all(N = length(ENIST))
  E = zeros(N)
  for at_number=1:N
    element = keys(atomic_electron_configuration)[at_number]
    @info "Calculating $element:" 
    results = lda(at_number, -35:0.1:20, δn = 1e-8)
    E[at_number] = results.energy.total
    @info @sprintf("Energy levels for the %s atom are:", element)
    for (quantum_numbers, ψ) in sort(collect(results.orbitals), by = x -> last(x).ϵᵢ)
      nᵣ, l = quantum_numbers 
      n = nᵣ + l + 1
      @info @sprintf("\t%i%s\t(%4.1f)\t%14.6f", n, atomic_shell[l], ψ.nᵢ, ψ.ϵᵢ)
    end
    @info @sprintf("ΔE = %12.6f", ENIST[at_number]-E[at_number])
  end
  @info "== RESULTS SUMMARY =="
  for Z=1:N
    element = keys(atomic_electron_configuration)[Z]
    @info @sprintf("%2s: ΔE = %+e", element, ENIST[Z]-E[Z])
  end  
end

test_all()
