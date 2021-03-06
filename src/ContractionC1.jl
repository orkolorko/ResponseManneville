using .C2BasisDefinition: C1, C1Norm



function opnormbound(N::Type{<:C1}, M, B)
    n = size(M, 2)

    est = 0

    for i in 1:n
        v = zeros(n)
        v[i] = 1
        z = C1Norm(B, v)
        est+= C1Norm(B, M*v)/z
    end
    return est.hi
end

function norms_of_powers_basis(B::Basis, m::Integer, Q::DiscretizedOperator, f::AbstractArray;
    normv0::Real=-1., #used as "missing" value
    normQ::Real=-1.,
    normE::Real=-1.,
    normEF::Real=-1.,
    normIEF::Real=-1.,
    normN::Real=-1.)

    @assert eltype(f) <: Interval
    T = typeof(zero(eltype(Q.L)).hi) # gets "Float64" from Q.L
    n = size(Q.L, 1)
    M = mid.(Q.L)
    R = radius.(Q.L)
    #δ = opnormbound(N, R)
    #γz = gamma(T, max_nonzeros_per_row(Q.L))
    #γn = gamma(T, n+3) # not n+2 like in the paper, because we wish to allow for f to be the result of rounding
    #ϵ = zero(T)
    midf = mid.(f)

    norm(v) = C2BasisDefinition.C1Norm(B, v)

    norms = zeros(m)

    S = zeros((n, m))


    for v in AverageZero(B)
        norm_0 = norm(v)
        v/= norm_0.lo
        for i in 1:m
            w = M*v 
            v = w - Q.e * (midf*w)[1]
            S[:, i]+= abs.(v)        
        end
    end

    return [norm(S[:, i]) for i in 1:m]
    
end








