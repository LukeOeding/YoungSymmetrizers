-- Young symmetrizers applied to cp factorizations

-- input: (1) the CP factorization of a tensor as factor matrices.
-- eg. A, B, C so that S = \sum_i a_i\otimes b_i \otimes c_i, with x_i the i-th column of X.
-- (2) a filling of a multi-tableau
-- eg T = {{{0,1,2},{3}}, {{0,1,3},{2}}, {{0,2,3},{1}}}
-- output: The product of determinants that (should) symmetrize to the value of the Young symmetrizer.

restart
tab = {{{0,1,2},{3}}, {{0,1,3},{2}}, {{0,2,3},{1}}}
A = matrix{{1,0,0,0},{0,1,0,0}, {0,0,1,0}, {0,0,0,1}} 
B = matrix{{1,0,0,0},{0,1,0,0}, {0,0,1,0}, {0,0,0,1}} 
C = matrix{{1,0,0,0},{0,1,0,0}, {0,0,1,0}, {0,0,0,1}} 
factorMats ={A,B,C}

t = sum(rank source A, r-> apply(rank(target A), i->apply(rank(target B), j-> apply(rank(target C), k-> A_(i,r)*B_(j,r)*C_(k,r)))))
matrix\t
A = matrix{{1,0,0,0,1},{0,1,0,0,1}, {0,0,1,0,1}, {0,0,0,1,1}} 
B = matrix{{1,0,0,0,1},{0,1,0,0,1}, {0,0,1,0,1}, {0,0,0,1,1}} 
C = matrix{{1,0,0,0,1},{0,1,0,0,1}, {0,0,1,0,1}, {0,0,0,1,1}} 
factorMats ={A,B,C}
tab

A = matrix{{1,0,0,1},{0,1,0,1}, {0,0,1,1}} 
B = matrix{{1,0,0,1},{0,1,0,1}, {0,0,1,1}} 
C = matrix{{1,0,0,1},{0,1,0,1}, {0,0,1,1}} 
factorMats ={A,B,C}
tab
tab = {{{0,1,2},{3}}, {{0,1,3},{2}}, {{0,2,3},{1}}}
t = sum(rank source A, r-> apply(rank(target A), i->apply(rank(target B), j-> apply(rank(target C), k-> A_(i,r)*B_(j,r)*C_(k,r)))))
matrix\t

-- for higher degree d and an expression of rank R, choose d rank-1 terms with replacement. This is a sequence in {0..R-1}^d.
-- Then use this sequence as a rule for how to choose the columns of the matrices. 
-- Note that if you're choosing for the first tableau, you don't choose repeats for columns, otherwise you're obviously zero.
-- The choices also have to be consistent for the other factors.

mtabs = {{{0,1,2},{3,4,5}}, {{0,2,4},{1,3,5}}, {{0,1,3},{2,4,5}}}
choices = {0,1,2,0,1,2} -- pick degree-many entries from the vectors in the A-factor, with replacement. 
apply(mtabs, tab -> apply(tab, row ->  apply(row, t-> choices_t)))
-- note that there are repeats in the third factor

apply(#factorMats, i->  apply(tab#i, t0 -> (factorMats#i)_t0^{0..#t0-1}))

YSCP = (Tabs, FacMats) -> (
fA = 
?
)