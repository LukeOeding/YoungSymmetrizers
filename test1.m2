---Young Symmetrizers  -- with tensor product rings, and one partial evaluation trick
-- the polynomial is the 3x3x3 degree 9 Strassen Invariant
--restart
(dg, aa, bb, cc) = (9,2,2,2);-- degree of the invariant, the projective dimensions of the spaces
--Ra = QQ[a_(0,0)..a_(aa,dg-1)] -- ring on A times CC^d 
--Rb = QQ[b_(0,0)..b_(bb,dg-1)] -- ring on B times CC^d 
--Rc = QQ[c_(0,0)..c_(cc,dg-1)] -- ring on C times CC^d 
--Rx = QQ[x_(0,0,0)..x_(aa,bb,cc)] -- ring on A \otimes B\otimes C
--R = Ra**Rb**Rc**Rx  -- big ring -- Later we compare to subsequent ring extensions instead of tensor product
--R0 = QQ[a_(0,0)..a_(aa,dg-1), b_(0,0)..b_(bb,dg-1),c_(0,0)..c_(cc,dg-1)] 
--R = R0[x_(0,0,0)..x_(aa,bb,cc)]
R0 = QQ[x_(0,0,0)..x_(aa,bb,cc)]
R = R0[a_(0,0)..a_(aa,dg-1), b_(0,0)..b_(bb,dg-1),c_(0,0)..c_(cc,dg-1)] 

-- the following matrices should have their column indices be determined by the 
--  columns of some triple of tableaux
A1= matrix apply(3,i -> apply({0,1,2}, j-> a_(i,j) ));
A2= matrix apply(3,i -> apply({3,4,5}, j-> a_(i,j) ));
A3= matrix apply(3,i -> apply({6,7,8}, j-> a_(i,j) ));

B1= matrix apply(3,i -> apply({0,3,6}, j-> b_(i,j) ));
B2= matrix apply(3,i -> apply({1,4,7}, j-> b_(i,j) ));
B3= matrix apply(3,i -> apply({2,5,8}, j-> b_(i,j) ));

C1= matrix apply(3,i -> apply({0,3,5}, j-> c_(i,j) ));
C2= matrix apply(3,i -> apply({1,4,6}, j-> c_(i,j) ));
C3= matrix apply(3,i -> apply({2,7,8}, j-> c_(i,j) ));

-- do one part at a time because we know that the determinant of the other A'matrices 
--  are just coming along for the ride
tmp = ()->(
time T1 =det(A1)*det(B1)*det(B2)*det(B3)*det(C1)*det(C2)*det(C3); -- 1.56s
F1=T1; time for d from 0 to 2 do F1= sum(cc+1,k-> sum(bb+1, j-> sum(aa+1,i->x_(i,j,k)*diff(a_(i,d)*b_(j,d)*c_(k,d),F1)))); --20.9s
time F12 = det(A2)*F1;--2.1s
time for d from 3 to 5 do F12= sum(cc+1,k-> sum(bb+1, j-> sum(aa+1,i->x_(i,j,k)*diff(a_(i,d)*b_(j,d)*c_(k,d),F12)))); --76s
time F123 = det(A3)*F12; --1.8s
time for d from 6 to 8 do F123= sum(cc+1,k-> sum(bb+1, j-> sum(aa+1,i->x_(i,j,k)*diff(a_(i,d)*b_(j,d)*c_(k,d),F123)))); --42s
F123
)
elapsedTime strassen9 = tmp(); 
-- 157s (using tensor product of rings)
-- 140s (using all at same level)
-- 150s (abc's first, then x's R[abc][x])
-- 148s (x's first, then abc's R[x][abc]