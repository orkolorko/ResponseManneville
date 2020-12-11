@testset "Testing the Dynamic" begin

# We start by testing the coordinate change
import ResponseManneville.InducedLSVMapDefinition: CoordinateChange, InvCoordinateChange

@test CoordinateChange(0.5) == 0
@test CoordinateChange(1) == 1
@test InvCoordinateChange(0) == 0.5
@test CoordinateChange(1) == 1

# We start now the map
import ResponseManneville.InducedLSVMapDefinition: ApproxInducedLSV, derleft, derderleft 
D = ApproxInducedLSV(0.5, 10)
@test D.α == 0.5
@test derleft(D, 0) == 1
@test derderleft(D, 0) == Inf

using TaylorSeries, IntervalArithmetic

T(x) = ResponseManneville.DynamicDefinition.plottable(D, x)
x = Taylor1([0.6, 1], 2)
y = T(x)

# we test the Taylor Series, the last branch is linear with derivative 2
@test y[1] == 2
@test y[2] == 0

# we test now the preimwithder_derder function
test_x = Base.rand(Float64, 10)

for t in test_x
    for k in 1:ResponseManneville.DynamicDefinition.nbranches(D)
        # preimwithder_derder computes the preimage through the k-th branch
        # with tight enclosures for the derivative and the second derivative
        val = ResponseManneville.InducedLSVMapDefinition.preimwithder_derder(D, k, t, 0.000001)
        val_z = T(Taylor1([mid(val[1]), 1], 2)) #this is computed through TaylorSeries
        @test abs(val_z[0]-t)<2^-15 && abs(val_z[1]-mid(Interval(val[2])))<2^-15 && abs(2*val_z[2]-mid(Interval(val[3])))<2^-15   
    end
end

end;