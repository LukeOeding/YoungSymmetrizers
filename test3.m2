---Young Symmetrizers  -- with tensor product rings, and one partial evaluation trick
restart
(dg, aa, bb, cc) = (9,2,2,2);-- degree of the invariant, the projective dimensions of the spaces
R0 = QQ[x_(0,0,0)..x_(aa,bb,cc)];
Ra = R0[a_(0,0)..a_(aa,dg-1)] -- ring on A times CC^d 
Rb = Ra[b_(0,0)..b_(bb,dg-1)] -- ring on B times CC^d 
Rc = Rb[c_(0,0)..c_(cc,dg-1)] -- ring on C times CC^d 
A1= matrix apply(3,i -> apply({0,1,2}, j-> a_(i,j) ));
A2= matrix apply(3,i -> apply({3,4,5}, j-> a_(i,j) ));
A3= matrix apply(3,i -> apply({6,7,8}, j-> a_(i,j) ));
B1= matrix apply(3,i -> apply({0,3,6}, j-> b_(i,j) ));
B2= matrix apply(3,i -> apply({1,4,7}, j-> b_(i,j) ));
B3= matrix apply(3,i -> apply({2,5,8}, j-> b_(i,j) ));
C1= matrix apply(3,i -> apply({0,3,5}, j-> c_(i,j) ));
C2= matrix apply(3,i -> apply({1,4,6}, j-> c_(i,j) ));
C3= matrix apply(3,i -> apply({2,7,8}, j-> c_(i,j) ));
tmp = ()->(
time T1 =det(A1)*det(B1)*det(B2)*det(B3)*det(C1)*det(C2)*det(C3); -- 1.5s
F1=T1; time for d from 0 to 2 do F1= sum(cc+1,k-> sum(bb+1, j-> sum(aa+1,i->x_(i,j,k)*coefficient(a_(i,d)*b_(j,d)*c_(k,d),F1)))); --20.4s
time F12 = det(A2)*F1;--2.3s
time for d from 3 to 5 do F12= sum(cc+1,k-> sum(bb+1, j-> sum(aa+1,i->x_(i,j,k)*coefficient(a_(i,d)*b_(j,d)*c_(k,d),F12)))); --78s
time F123 = det(A3)*F12; --2.4s
time for d from 6 to 8 do F123= sum(cc+1,k-> sum(bb+1, j-> sum(aa+1,i->x_(i,j,k)*coefficient(a_(i,d)*b_(j,d)*c_(k,d),F123)))); --44s
F123
)
elapsedTime strassen9 = tmp(); -- 160s
-- what happens when you re-use variables? -- doesn't seem to matter...
tmp = ()->(
time F =det(A1)*det(B1)*det(B2)*det(B3)*det(C1)*det(C2)*det(C3); -- 1.5s
time for d from 0 to 2 do F= sum(cc+1,k-> sum(bb+1, j-> sum(aa+1,i->x_(i,j,k)*coefficient(a_(i,d)*b_(j,d)*c_(k,d),F)))); --20.4s
time F = det(A2)*F;--2.3s
time for d from 3 to 5 do F= sum(cc+1,k-> sum(bb+1, j-> sum(aa+1,i->x_(i,j,k)*coefficient(a_(i,d)*b_(j,d)*c_(k,d),F)))); --78s
time F = det(A3)*F; --2.4s
time for d from 6 to 8 do F= sum(cc+1,k-> sum(bb+1, j-> sum(aa+1,i->x_(i,j,k)*coefficient(a_(i,d)*b_(j,d)*c_(k,d),F)))); --44s
F
)
elapsedTime strassen9 = tmp(); -- 161s
-- what happens if you add in the B'dets one at a time?
tmp = ()->(
time F =det(A1)*det(B1)*det(C1)*det(C2)*det(C3); -- 0.04s
time for d from 0 to 0 do F= sum(cc+1,k-> sum(bb+1, j-> sum(aa+1,i->x_(i,j,k)*coefficient(a_(i,d)*b_(j,d)*c_(k,d),F)))); --0.17s
time F = det(B2)*F; -- 0.05
time for d from 1 to 1 do F= sum(cc+1,k-> sum(bb+1, j-> sum(aa+1,i->x_(i,j,k)*coefficient(a_(i,d)*b_(j,d)*c_(k,d),F)))); --1s
time F = det(B3)*F; --0.34
time for d from 2 to 2 do F= sum(cc+1,k-> sum(bb+1, j-> sum(aa+1,i->x_(i,j,k)*coefficient(a_(i,d)*b_(j,d)*c_(k,d),F)))); -- 6.19
time F = det(A2)*F;--2.3s
time for d from 3 to 5 do F= sum(cc+1,k-> sum(bb+1, j-> sum(aa+1,i->x_(i,j,k)*coefficient(a_(i,d)*b_(j,d)*c_(k,d),F)))); --78s
time F = det(A3)*F; --2.5s
time for d from 6 to 8 do F= sum(cc+1,k-> sum(bb+1, j-> sum(aa+1,i->x_(i,j,k)*coefficient(a_(i,d)*b_(j,d)*c_(k,d),F)))); --45.7s
F
)
elapsedTime strassen9 = tmp(); -- 149s

-- what happens if you add in the B'dets and the C'dets one at a time?
tmp = ()->(
time F =det(A1)*det(B1)*det(C1); -- 0.001s
time for d from 0 to 0 do F= sum(cc+1,k-> sum(bb+1, j-> sum(aa+1,i->x_(i,j,k)*coefficient(a_(i,d)*b_(j,d)*c_(k,d),F)))); --0.008s
time F = det(B2)*det(C2)*F; -- 0.017
time for d from 1 to 1 do F= sum(cc+1,k-> sum(bb+1, j-> sum(aa+1,i->x_(i,j,k)*coefficient(a_(i,d)*b_(j,d)*c_(k,d),F)))); --0.15s
time F = det(B3)*det(C3)*F; --0.44
time for d from 2 to 2 do F= sum(cc+1,k-> sum(bb+1, j-> sum(aa+1,i->x_(i,j,k)*coefficient(a_(i,d)*b_(j,d)*c_(k,d),F)))); -- 6.41
time F = det(A2)*F;--2.3s
time for d from 3 to 5 do F= sum(cc+1,k-> sum(bb+1, j-> sum(aa+1,i->x_(i,j,k)*coefficient(a_(i,d)*b_(j,d)*c_(k,d),F)))); --79s
time F = det(A3)*F; --2.75s
time for d from 6 to 8 do F= sum(cc+1,k-> sum(bb+1, j-> sum(aa+1,i->x_(i,j,k)*coefficient(a_(i,d)*b_(j,d)*c_(k,d),F)))); --44s
F
)
elapsedTime strassen9 = tmp(); -- 147.6s