--symbolic (block) Laplace expansion of determinants

--------------------
-- expand the determinant as a sum of determinants, but the determinants are symoblic until requested
--------------------
--------------------
-- example
--------------------
restart
--------------------
-- make a ring with the variables that factor the way we're expecting:
--------------------
R = QQ; for s in {4,3,2} do ( R= R[flatten apply(subsets(4,s), j-> apply(subsets(4,s), i-> a_(i,j)))]); 
-- a_{ss,tt} = the determinant with rows ss and columns tt
R = R[a_(0,0)..a_(3,3)]
for i to 3 do for j to 3 do a_({i},{j}) = a_(i,j) -- get rid of the braces at the bottom level
-- look at these expressions - they factor the variables with the longer indices
(a_(0,0)+a_(1,1))*a_({1,2},{1,3})
(a_(1,1))*(a_({1,2},{1,3}) + a_({1,3},{1,3}))
(a_(0,0)+a_(1,1))*(a_({1,2},{1,3}) + a_({1,3},{1,3}))
(a_(0,0)+a_(1,1))*(a_({1,2},{1,3}) + a_({1,3},{1,3}))*(a_({0,1,2},{0,1,3}) + a_({0,1,3},{0,1,3}))

--- look at what happens if we revers the order we append the variables:
restart
R0 = QQ; R1 = R0[a_(0,0)..a_(3,3)]
for s to 2 do print peek value ("R"|(toString s))
R2
for s in {2,3,4} do ( toString(concatenate("R",toString(s)))= value("R"|toString(s-1))[flatten apply(subsets(4,s), j-> apply(subsets(4,s), i-> a_(i,j)))]); 
for i to 3 do for j to 3 do a_({i},{j}) = a_(i,j) -- get rid of the braces at the bottom level
-- look at these expressions - they factor the variables with the longer indices
(a_(0,0)+a_(1,1))*a_({1,2},{1,3})
(a_(1,1))*(a_({1,2},{1,3}) + a_({1,3},{1,3}))
(a_(0,0)+a_(1,1))*(a_({1,2},{1,3}) + a_({1,3},{1,3}))
(a_(0,0)+a_(1,1))*(a_({1,2},{1,3}) + a_({1,3},{1,3}))*(a_({0,1,2},{0,1,3}) + a_({0,1,3},{0,1,3}))
-- the factoring order reverses, but also the intermediate rings don't have names, so the display is a mess.


A = matrix apply(4, i-> apply(4, j-> a_(i,j)))
eval2 = flatten apply(subsets(4,2), i-> apply(subsets(4,2), j-> a_(i,j) => det(A^i_j)))
eval3 = flatten apply(subsets(4,3), i-> apply(subsets(4,3), j-> a_(i,j) => det(A^i_j)))

setComplement:= (I,n) ->  toList(set(0..(n-1)) - set I);
-- expand on rows i
LaplaceRowExpand:= method();
LaplaceRowExpand(List, Matrix) := (L, Atmp)->(
n := rank source Atmp;
k := length L;
sum(subsets(n,k), K -> (-1)^((sum L) - (sum K))*a_(L,K)*a_(setComplement(L,n), setComplement(K,n)) )
) 

 f = sub(LaplaceRowExpand({1,3},A), R)
 det A - sub(f, eval2)
-- expand on col i
LaplaceColExpand:= method();
LaplaceColExpand(List, Matrix) := (L, Atmp)->(
n := rank source Atmp;
k := length L;
sum(subsets(n,k), K -> (-1)^((sum L) - (sum K))*a_(K,L)*a_(setComplement(K,n), setComplement(L,n)) )
) 

f = sub(LaplaceColExpand({1,3},A), R)
det A - sub(f, eval2)

f = sub(LaplaceRowExpand({0,1,2},A), R)
det A - sub(f, eval3)
last baseName a_({1,2},{1,2})

LaplaceRowExpand({0,1,2},A)*LaplaceRowExpand({0,1,3},A)
