---Young Symmetrizers  -- with tensor product rings, and one partial evaluation trick
restart
-- the folowing code produces a degree 9 polynomial in the 27 variables x_(i,j,k),
-- using auxillary variables a_(i,d), b_(j,d), c_(k,d)
-- the polynomial is the 3x3x3 degree 9 Strassen Invariant
(dg, aa, bb, cc) = (9,2,2,2);-- degree of the invariant, the projective dimensions of the spaces
R = QQ;  for i from 0 to dg-1 do(R= R[a_(0,dg-i-1)..a_(aa,dg-i-1), b_(0,dg-i-1)..b_(bb,dg-i-1), c_(0,dg-i-1)..c_(cc,dg-i-1)];) -- ring on C times CC^d 
R = R[x_(0,0,0)..x_(aa,bb,cc)] -- ring on A \otimes B\otimes C
--R = Ra_(aa)**Rb_(bb)**Rc_(cc)**Rx  -- big ring -- Later we compare to subsequent ring extensions instead of tensor product
A1= matrix apply(3,i -> apply({0,1,2}, j-> a_(i,j) ));A2= matrix apply(3,i -> apply({3,4,5}, j-> a_(i,j) ));A3= matrix apply(3,i -> apply({6,7,8}, j-> a_(i,j) ));
B1= matrix apply(3,i -> apply({0,3,6}, j-> b_(i,j) ));B2= matrix apply(3,i -> apply({1,4,7}, j-> b_(i,j) ));B3= matrix apply(3,i -> apply({2,5,8}, j-> b_(i,j) ));
C1= matrix apply(3,i -> apply({0,3,5}, j-> c_(i,j) ));C2= matrix apply(3,i -> apply({1,4,6}, j-> c_(i,j) ));C3= matrix apply(3,i -> apply({2,7,8}, j-> c_(i,j) ));
myList = flatten flatten apply(aa+1, i-> apply(bb+1, j-> apply(cc+1, k-> {i,j,k}  )));
tmp = ()->(
time T1 =det(A1)*det(B1)*det(B2)*det(B3)*det(C1)*det(C2)*det(C3); -- 1.38s
F=T1; time for d from 0 to 2 do F = sum(myList, L-> x_(L#0, L#1, L#2)*diff(a_(L#0,d)*b_(L#1,d)*c_(L#2,d),F));--27s
time F = det(A2)*F;--6.2.3s
time for d from 3 to 5 do F = sum(myList, L-> x_(L#0, L#1, L#2)*diff(a_(L#0,d)*b_(L#1,d)*c_(L#2,d),F)); -- 96.6s
time F = det(A3)*F; --2.1s
time for d from 6 to 8 do F = sum(myList, L-> x_(L#0, L#1, L#2)*diff(a_(L#0,d)*b_(L#1,d)*c_(L#2,d),F)); -- 130s
F
)
elapsedTime strassen9 = tmp();  --195s
#terms F
det(A1)*det(B1)