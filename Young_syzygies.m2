---Young Symmetrizers  -- trying to compute syzygies, try the Koszul complex
restart
(dg, aa, bb) = (9,2,2);-- degree of the invariant, the projective dimensions of the spaces
R0 = QQ[x_(0,0)..x_(aa,bb)] -- ring on A \otimes B
Ra = R0[a_(0,0)..a_(aa,dg-1)] -- ring on A times CC^d 
Rb = Ra[b_(0,0)..b_(bb,dg-1)] -- ring on B times CC^d 

-- suppose we already know that the bi-partition {{1},{1}} is the first module M, so we compute the second one by computing all the possible syzygy modules in the tensor product M ** R

-- two modules show up:
A1= matrix apply(2,i -> apply({0,1}, j-> a_(i,j) ));
B1= matrix apply(2,i -> apply({0,1}, j-> b_(i,j) ));
T1 = det(A1)*det(B1)
-- full evaluation:
val = map(R0, R0, random(R0^1,R0^(#gens R0)))
F1=T1; for d from 0 to 0 do F0= sum(bb+1, j-> sum(aa+1,i-> (val x_(i,j))*diff(a_(i,d)*b_(j,d),F1)));
for d from 1 to 1 do F0= sum(bb+1, j-> sum(aa+1,i->  x_(i,j)*diff(a_(i,d)*b_(j,d),F0)));
factor F0

F1=T1; for d from 1 to 1 do F1= sum(bb+1, j-> sum(aa+1,i-> (val x_(i,j))*diff(a_(i,d)*b_(j,d),F1)));
for d from 0 to 0 do F1= sum(bb+1, j-> sum(aa+1,i->  x_(i,j)*diff(a_(i,d)*b_(j,d),F1)));
factor F1
F1 - F0


A1= matrix apply(1,i -> apply({0}, j-> a_(i,j) ));
A2= matrix apply(1,i -> apply({1}, j-> a_(i,j) ));
B1= matrix apply(1,i -> apply({0}, j-> b_(i,j) ));
B2= matrix apply(1,i -> apply({1}, j-> b_(i,j) ));



-- do one part at a time because we know that the determinant of the other A'matrices 
--  are just coming along for the ride
tmp = ()->(
time T1 =det(A1)*det(B1)*det(B2)*det(B3);
F1=T1; time for d from 0 to 2 do F1= sum(cc+1,k-> sum(bb+1, j-> sum(aa+1,i->x_(i,j,k)*diff(a_(i,d)*b_(j,d)*c_(k,d),F1)))); --20.9s
time F12 = det(A2)*F1;--2.3s
time for d from 3 to 5 do F12= sum(cc+1,k-> sum(bb+1, j-> sum(aa+1,i->x_(i,j,k)*diff(a_(i,d)*b_(j,d)*c_(k,d),F12)))); --78s
time F123 = det(A3)*F12; --2.3s
time for d from 6 to 8 do F123= sum(cc+1,k-> sum(bb+1, j-> sum(aa+1,i->x_(i,j,k)*diff(a_(i,d)*b_(j,d)*c_(k,d),F123)))); --45s
F123
)
elapsedTime strassen9 = tmp();--164s
 strassen9