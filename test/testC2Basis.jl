@testset "Testing the C2 Basis" begin

import ResponseManneville.C2BasisDefinition

N = 16
B = C2Basis(N)
@test B.p == [i/N for i in 0:N]

for i in 1:N+1
    ϕ, ϕ′ = ResponseManneville.C2BasisDefinition.basis_element_with_der(B, i)
    v = zeros(N+1)
    v[i] = 1
    @test maximum(abs.(mid.(ϕ.(Interval.(B.p))) - v)) < 2^-15
    @test maximum(abs.(mid.(ϕ′.(Interval.(B.p))) - zeros(N+1))) < 2^-15 
end 

for i in 1:N+1
    ν, ν′ = ResponseManneville.C2BasisDefinition.basis_element_with_der(B, i+N+1)
    v = zeros(N+1)
    v[i] = 1
    @test maximum(abs.(mid.(ν′.(Interval.(B.p))) - v)) < 2^-15
    @test maximum(abs.(mid.(ν.(Interval.(B.p))) - zeros(N+1))) < 2^-15 
end 


end;