--test7 the idea here is to make an auxillary ring with partial progress - 2 factors at a time:
restart
(dg, aa, bb, cc) = (9,2,2,2);-- degree of the invariant, the projective dimensions of the spaces
Rx = QQ[x_(0,0,0)..x_(aa,bb,cc)] -- ring on A \otimes B\otimes C
Rz = Rx[z_(0,0,0)..z_(aa,bb,dg-1)] -- intermediate ring
Ra = Rz[a_(0,0)..a_(aa,dg-1)] -- ring on A times CC^d 
Rb = Ra[b_(0,0)..b_(bb,dg-1)] -- ring on B times CC^d 
Rc = Rb[c_(0,0)..c_(cc,dg-1)] -- ring on C times CC^d 
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

tmp = ()->(
time T1 =det(A1)*det(B1)*det(B2)*det(B3);
F1=T1; time for d from 0 to 2 do F1=  sum(bb+1, j-> sum(aa+1,i->z_(i,j,d)*diff(a_(i,d)*b_(j,d),F1))); 
time F12 = det(A2)*F1;
time for d from 3 to 5 do F12=  sum(bb+1, j-> sum(aa+1,i->z_(i,j,d)*diff(a_(i,d)*b_(j,d),F12))); 
time F123 = det(A3)*F12; 
time for d from 6 to 8 do F123= sum(bb+1, j-> sum(aa+1,i->z_(i,j,d)*diff(a_(i,d)*b_(j,d),F123))); 
time F =det(C1)*F123;
time for d in {0,3,5} do F=  sum(aa+1, i-> sum(bb+1, j-> sum(cc+1,k->x_(i,j,k)*diff(c_(k,d)*z_(i,j,d),F)))); 
time F =det(C2)*F;
time for d in {1,4,6} do F=  sum(aa+1, i-> sum(bb+1, j-> sum(cc+1,k->x_(i,j,k)*diff(c_(k,d)*z_(i,j,d),F)))); 
time F =det(C3)*F;
time for d in {2,7,8} do F=  sum(aa+1, i-> sum(bb+1, j-> sum(cc+1,k->x_(i,j,k)*diff(c_(k,d)*z_(i,j,d),F)))); 
F)
elapsedTime strassen9 = tmp(); -- 147s

#terms strassen9