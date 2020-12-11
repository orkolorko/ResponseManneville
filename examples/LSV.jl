using ResponseManneville
using ValidatedNumerics

D = ApproxInducedLSV(0.5, 10)
B = C2Basis(1024)
Q = DiscretizedOperator(B, D)