using Tsys
using Test

# These J1939-6342 (aka 1934-638) flux models are from different epochs, so
# returned flux is expected to be similar, but not identical.
model_1e6 = [-30.7667, 26.4908, -7.0977, 0.605334] # From MeerKAT, ν1=1e6
model_1e9 = [1.186E+00, 1.191E-01, -1.191E+00]     # From ATCA, ν1=1e9

@testset "Tsys.jl" begin
    @test model_flux(model_1e6) == 10^model_1e6[1]
    @test model_flux(model_1e9) == 10^model_1e9[1]
end
