
-- expand on rows i
-- Laplace Row Expansion
--------------------

LaplaceRowExpand:= method();
LaplaceRowExpand(List, IndexedVariableTable, List, List) := (L, a, P, Q)->(
--- expand on rows L (usual indexing), the matrix has symbols a, row indices P (any indices), and column indices Q (any indices)
-- L must be a subset of P
k := length L;
n := length P;
print matrix(apply(P, p-> apply(Q, q-> a_(p,q))));
sum(subsets(Q,k), K -> (-1)^((sum L) - (sum K))*a_(L,K)*a_(setComplement(L,P), setComplement(K,Q)) )
) 

------------------
-- Laplace diff for matrices:
------------------
LaplaceDiff:=method();
LaplaceDiff(ZZ,ZZ,List,List) :=(i,j,I,J) -> (-1)^(locateInList(i,I)- locateInList(j,J))*a_(i,j)*a_(setComplement({i},I), setComplement({j},J))
LaplaceDiff(ZZ,ZZ,RingElement) :=(i,j,v) -> (IJ := last baseName v; LaplaceDiff(i,j,IJ#0,IJ#1))

-- differentiate with respect to the a_(i,j)-entry of a matrix
-- whose rows are I and columns J
end 
--------------------
-- example
--------------------

--------------------
-- make a ring with the variables that factor the way we're expecting:
--------------------
-- this makes the rings, but forgets the names...
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
-- you can name rings with subscripts in their name, but this might interfere with using subscript to get generators
R0 = QQ; 
R1 = R0[a_(0,0)..a_(3,3)]
R2 = R1[flatten apply(subsets(4,2), j-> apply(subsets(4,2), i-> a_(i,j)))]
R3 = R2[flatten apply(subsets(4,3), j-> apply(subsets(4,3), i-> a_(i,j)))]
--------------------
R4 = R3[flatten apply(subsets(4,4), j-> apply(subsets(4,4), i-> a_(i,j)))]
for i to 3 do for j to 3 do a_({i},{j}) = a_(i,j) -- get rid of the braces at the bottom level
-- look at these expressions - they factor the variables with the longer indices
(a_(0,0)+a_(1,1))*a_({1,2},{1,3})
(a_(1,1))*(a_({1,2},{1,3}) + a_({1,3},{1,3}))
(a_(0,0)+a_(1,1))*(a_({1,2},{1,3}) + a_({1,3},{1,3}))
(a_(0,0)+a_(1,1))*(a_({1,2},{1,3}) + a_({1,3},{1,3}))*(a_({0,1,2},{0,1,3}) + a_({0,1,3},{0,1,3}))
-- now see how the block Laplace expansion works
A = matrix apply(4, i-> apply(4, j-> a_(i,j)))
-- evaluation maps that compute the actual determinants in place of their placeholders
for i to 3 do for j to 3 do a_({i},{j}) = a_(i,j) -- get rid of the braces at the bottom level, we don't need a map for this one
eval2 = flatten apply(subsets(4,2), i-> apply(subsets(4,2), j-> a_(i,j) => det(A^i_j))) 
eval3 = flatten apply(subsets(4,3), i-> apply(subsets(4,3), j-> a_(i,j) => det(A^i_j)))
eval4 = flatten apply(subsets(4,4), i-> apply(subsets(4,4), j-> a_(i,j) => det(A^i_j)))

setComplement:=method();
setComplement (List, ZZ) := (I,n) ->  toList(set(0..(n-1)) - set I);
setComplement (List, List):= (I,P) -> if isSubset(set I, set P) then toList(set P - set I) else print("error, not a subset");
setComplement({1,2},{1,2,3})
setComplement({1,2},3)
-- expand on rows i
LaplaceRowExpand:= method();
LaplaceRowExpand(List, IndexedVariableTable, List, List) := (L, a, P, Q)->(
--- expand on rows L (usual indexing), the matrix has symbols a, row indices P (any indices), and column indices Q (any indices)
-- L must be a subset of P
k := length L;
n := length P;
print matrix(apply(P, p-> apply(Q, q-> a_(p,q))));
sum(subsets(Q,k), K -> (-1)^((sum L) - (sum K))*a_(L,K)*a_(setComplement(L,P), setComplement(K,Q)) )
) 
f = LaplaceRowExpand({1,3},a,{0,1,2,3},{0,1,2,3})
det A - sub(f, eval2)
-- can also take expansion of a matrix with different row/column indices
-- 
f = LaplaceRowExpand({1,3},a,{1,2,3},{0,1,3})

------------------
-- Laplace diff for matrices:
------------------
-- assume I is a list with no repeats, find the position of i in I:
locateInList:=method();
locateInList (ZZ,List):= (i,I) -> for s to length I-1 do if I_s ==i then return s else continue;
LaplaceDiff:=method();
LaplaceDiff(ZZ,ZZ,List,List, IndexedVariableTable) :=(i,j,I,J,aa) -> (-1)^(locateInList(i,I)- locateInList(j,J))*aa_(i,j)*aa_(setComplement({i},I), setComplement({j},J))
LaplaceDiff(1,2,{0,1,2},{0,2,3})
LaplaceDiff(ZZ,ZZ,RingElement) :=(i,j,v) -> (IJ := last baseName v; LaplaceDiff(i,j,IJ#0,IJ#1))
LaplaceDiff(1,2,a_({0,1,2},{0,2,3}))
LaplaceDiff(1,2,a_({0,1,2,3},{0,1,2,3}))

-- differentiate with respect to the a_(i,j)-entry of a matrix
-- whose rows are I and columns J


-- old method -- note it doesn't actually depend on the actual matrix...
-- expand on col i
LaplaceColExpand:= method();
LaplaceColExpand(List, Matrix) := (L, Atmp)->(
n := rank source Atmp;
k := length L;
sum(subsets(n,k), K -> (-1)^((sum L) - (sum K))*a_(K,L)*a_(setComplement(K,n), setComplement(L,n)) )
) 

f = LaplaceColExpand({1,3},A)
det A - sub(f, eval2)

f = LaplaceRowExpand({0,1,2},A)
det A - sub(f, eval3)
