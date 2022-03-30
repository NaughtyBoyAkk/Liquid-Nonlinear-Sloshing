###########################################################################
# Gauss-Kronrod integration-weight computation for arbitrary floating-point
# types and precision, implemented based on the description in:
#
#    Dirk P. Laurie, "Calculation of Gauss-Kronrod quadrature rules,"
#    Mathematics of Computation, vol. 66, no. 219, pp. 1133-1145 (1997).
#
# for the Kronrod rule, and for the Gauss rule from the description in
#
#    Lloyd N. Trefethen and David Bau, Numerical Linear Algebra (SIAM, 1997).
#
# Arbitrary-precision eigenvalue (eignewt & eigpoly) and eigenvector
# (eigvec1) routines are written by SGJ, independent of the above sources.
#
# Since we only implement Gauss-Kronrod rules for the unit weight function,
# the diagonals of the Jacobi matrices are zero and certain things simplify
# compared to the general case of an arbitrary weight function.

# Given a symmetric tridiagonal matrix H with H[i,i] = 0 and
# H[i-1,i] = H[i,i-1] = b[i-1], compute p(z) = det(z I - H) and its
# derivative p'(z), returning (p,p').
function eigpoly(b,z,m=length(b)+1)
    d1 = z
    d1deriv = d2 = one(z)
    d2deriv = zero(z)
    for i = 2:m
        b2 = b[i-1]^2
        d = z * d1 - b2 * d2
        dderiv = d1 + z * d1deriv - b2 * d2deriv
        d2 = d1
        d1 = d
        d2deriv = d1deriv
        d1deriv = dderiv
    end
    return (d1, d1deriv)
end

# compute the n smallest eigenvalues of the symmetric tridiagonal matrix H
# (defined from b as in eigpoly) using a Newton iteration
# on det(H - lambda I).  Unlike eig, handles BigFloat.
function eignewt(b,m,n)
    # get initial guess from eig on Float64 matrix
    H = SymTridiagonal(zeros(m), Float64[ b[i] for i in 1:m-1 ])
    lambda0 = sort(eigvals(H))

    lambda = Array{eltype(b)}(undef,n)
    for i = 1:n
        lambda[i] = lambda0[i]
        for k = 1:1000
            (p,pderiv) = eigpoly(b,lambda[i],m)
            lambda[i] = (lamold = lambda[i]) - p / pderiv
            if abs(lambda[i] - lamold) < 10 * eps(lambda[i]) * abs(lambda[i])
                break
            end
        end
        # do one final Newton iteration for luck and profit:
        (p,pderiv) = eigpoly(b,lambda[i],m)
        lambda[i] = lambda[i] - p / pderiv
    end
    return lambda
end

# given an eigenvalue z and the matrix H(b) from above, return
# the corresponding eigenvector, normalized to 1.
using LinearAlgebra
function eigvec1(b,z::Number,m=length(b)+1)
    # "cheat" and use the fact that our eigenvector v must have a
    # nonzero first entries (since it is a quadrature weight), so we
    # can set v[1] = 1 to solve for the rest of the components:.
    v = Array{eltype(b)}(undef,m)
    v[1] = 1
    if m > 1
        s = v[1]
        v[2] = z * v[1] / b[1]
        s += v[2]^2
        for i = 3:m
            v[i] = - (b[i-2]*v[i-2] - z*v[i-1]) / b[i-1]
            s += v[i]^2
        end
        rmul!(v, 1 / sqrt(s))
    end
    return v
end

"""
gauss([T,] N)

Return a pair `(x, w)` of `N` quadrature points `x[i]` and weights `w[i]` to
integrate functions on the interval `(-1, 1)`,  i.e. `sum(w .* f.(x))`
approximates the integral.  Uses the method described in Trefethen &
Bau, Numerical Linear Algebra, to find the `N`-point Gaussian quadrature
in O(`N`²) operations.

`T` is an optional parameter specifying the floating-point type, defaulting
to `Float64`. Arbitrary precision (`BigFloat`) is also supported.
"""
# function gauss{T<:AbstractFloat}(::Type{T}, N::Integer)
function gauss(::Type{T}, N::Integer) where {T<:AbstractFloat}
    if N < 1
        throw(ArgumentError("Gauss rules require positive order"))
    end
    o = one(T)
    b = T[ n / sqrt(4n^2 - o) for n = 1:N-1 ]
    x = eignewt(b,N,N)
    w = T[ 2*eigvec1(b,x[i])[1]^2 for i = 1:N ]
    return (x, w)
end

gauss(N::Integer) = gauss(Float64, N)
