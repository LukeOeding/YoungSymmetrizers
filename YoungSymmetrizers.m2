---Young Symmetrizers
restart
-- the folowing code produces a degree 4 polynomial in the 8 variables x_(i,j,k),
-- using auxillary variables a_(i,d), b_(j,d), c_(k,d)

-- the polynomial is the 2x2x2 hyperdeterminant


Ra = QQ[a_(0,0)..a_(1,3)] -- ring on  A times CC^d 
Rb = QQ[b_(0,0)..b_(1,3)] -- ring on  B times CC^d 
Rc = QQ[c_(0,0)..c_(1,3)] -- ring on  C times CC^d 
Rx = QQ[x_(0,0,0)..x_(1,1,1)] -- ring on A \otimes B\otimes C
R = Ra**Rb**Rc**Rx  -- big ring


-- the following matrices should have their column indices be determined by the 
--  columns of some triple of tableaux, and ideally there should be a function that makes 
--  them dynamically
A1= matrix apply(2,i -> apply({0,1}, j-> a_(i,j) ))
A2= matrix apply(2,i -> apply({2,3}, j-> a_(i,j) ))
B1= matrix apply(2,i -> apply({0,1}, j-> b_(i,j) ))
B2= matrix apply(2,i -> apply({2,3}, j-> b_(i,j) ))
C1= matrix apply(2,i -> apply({0,2}, j-> c_(i,j) ))
C2= matrix apply(2,i -> apply({1,3}, j-> c_(i,j) ))

-- this should be held as a product and not expanded for larger problems
T =  det(A1)*det(A2)*det(B1)*det(B2)*det(C1)*det(C2)


-- this loop builds the poylnomial in Rx, one factor at a time.  This should also be held as a product.
F=T;
for d from 0 to 3 do F= sum(2,k-> sum(2, j-> sum(2,i->x_(i,j,k)*
	    diff(a_(i,d)*b_(j,d)*c_(k,d),F))));
F

-- there should also be an evaluation option, instead of multiplying by each x_I
 -- evauate at each stage.

---Young Symmetrizers  -- try II
restart
-- the folowing code produces a degree 4 polynomial in the variables x_(i,j,k),
-- using auxillary variables a_(i,d), b_(j,d), c_(k,d)

-- the polynomial is the 3x3x3 degree 4 Strassen polynomial(s)
xx = 4
aa = 2
bb = 2 
cc = 2

Ra = QQ[a_(0,0)..a_(aa,xx)] -- ring on  A times CC^d 
Rb = QQ[b_(0,0)..b_(bb,xx)] -- ring on  B times CC^d 
Rc = QQ[c_(0,0)..c_(cc,xx)] -- ring on  C times CC^d 
Rx = QQ[x_(0,0,0)..x_(aa,bb,cc)] -- ring on A \otimes B\otimes C
R = Ra**Rb**Rc**Rx  -- big ring


-- the following matrices should have their column indices be determined by the 
--  columns of some triple of tableaux, and ideally there should be a function that makes 
--  them dynamically
A1= matrix apply(3,i -> apply({0,1,2}, j-> a_(i,j) ))
A2= matrix apply(1,i -> apply({3}, j-> a_(i,j) ))

B1= matrix apply(3,i -> apply({0,1,3}, j-> b_(i,j) ))
B2= matrix apply(1,i -> apply({2}, j-> b_(i,j) ))

C1= matrix apply(3,i -> apply({0,2,3}, j-> c_(i,j) ))
C2= matrix apply(1,i -> apply({1}, j-> c_(i,j) ))

-- this should be held as a product and not expanded for larger problems
T = det(A1)*det(A2)*det(B1)*det(B2)*det(C1)*det(C2);

-- this loop builds the poylnomial in Rx, one factor at a time.  This should also be held as a product.
F=T;
for d from 0 to xx-1 do F= sum(cc+1,k-> sum(bb+1, j-> sum(aa+1,i->x_(i,j,k)*contract(a_(i,d)*b_(j,d)*c_(k,d),F))));

uco = F-> gcd unique apply( unique flatten entries (coefficients F)_1, p-> abs sub(p,QQ))
F = F/((uco F)*(leadCoefficient F))
-- there should also be an evaluation option, instead of multiplying by each x_I
 -- evauate at each stage.

--- next step --- construct Lie-Algebra operators

lie = (F,m,p,q)->(
     if m ==1 then tmp = sum(bb+1,j-> sum(cc+1,k->  x_(q,j,k) *diff(x_(p,j,k) ,F) ));
     if m ==2 then tmp = sum(aa+1,i-> sum(cc+1,k->  x_(i,q,k) *diff(x_(i,p,k) ,F) ));
     if m ==3 then tmp = sum(aa+1,i-> sum(bb+1,j->  x_(i,j,q) *diff(x_(i,j,p) ,F) ));
if tmp != 0 then tmp = tmp/uco(tmp);
tmp
);
F
lie(F,1,1,0)
-- apply the Lie algebra many times in every factor, eventually get a basis of the module
L = {F};
for m from 1 to 3 do for f in L do  for j from 0 to bb do ( for k from 0 to cc do( tmp = lie(f,m,j,k); if tmp!=0 and not member(tmp,L) and not member(-tmp,L) then L = L|{tmp};))
length L


---Young Symmetrizers  -- try III
restart
-- the folowing code produces a degree 4 polynomial in the 8 variables x_(i,j,k),
-- using auxillary variables a_(i,d), b_(j,d), c_(k,d)

-- the polynomial is the 3x3x3 degree 9 Strassen Invariant
xx = 9
aa = 2
bb = 2 
cc = 2

Ra = QQ[a_(0,0)..a_(aa,xx)] -- ring on  A times CC^d 
Rb = QQ[b_(0,0)..b_(bb,xx)] -- ring on  B times CC^d 
Rc = QQ[c_(0,0)..c_(cc,xx)] -- ring on  C times CC^d 
Rx = QQ[x_(0,0,0)..x_(aa,bb,cc)] -- ring on A \otimes B\otimes C
R = Ra**Rb**Rc**Rx  -- big ring


-- the following matrices should have their column indices be determined by the 
--  columns of some triple of tableaux, and ideally there should be a function that makes 
--  them dynamically
A1= matrix apply(3,i -> apply({0,1,2}, j-> a_(i,j) ))
A2= matrix apply(3,i -> apply({3,4,5}, j-> a_(i,j) ))
A3= matrix apply(3,i -> apply({6,7,8}, j-> a_(i,j) ))

B1= matrix apply(3,i -> apply({0,3,6}, j-> b_(i,j) ))
B2= matrix apply(3,i -> apply({1,4,7}, j-> b_(i,j) ))
B3= matrix apply(3,i -> apply({2,5,8}, j-> b_(i,j) ))

C1= matrix apply(3,i -> apply({0,4,8}, j-> c_(i,j) ))
C2= matrix apply(3,i -> apply({1,5,6}, j-> c_(i,j) ))
C3= matrix apply(3,i -> apply({2,3,7}, j-> c_(i,j) ))

-- this should be held as a product and not expanded for larger problems
time T = det(A1)*det(A2)*det(A3)*det(B1)*det(B2)*det(B3)*det(C1)*det(C2)*det(C3);
--T = hold det(A1)*det(A2)*det(A3)*det(B1)*det(B2)*det(B3)*det(C1)*det(C2)*det(C3);

-- this loop builds the poylnomial in Rx, one factor at a time.  This should also be held as a product.
F=T
for d from 0 to xx do F= sum(cc+1,k-> sum(bb+1, j-> sum(cc+1,i->x_(i,j,k)*diff(a_(i,d)*b_(j,d)*c_(k,d),F))));
F

-- there should also be an evaluation option, instead of multiplying by each x_I
 -- evauate at each stage.


restart
-- the folowing code produces a degree 6 polynomial in the 32 variables x_(i,j,k,l,m),
-- using auxillary variables a_(i,p), b_(j,p), c_(k,p), d_(l,p), e_(m,p)

xx = 6
n = 1


Ra = QQ[a_(0,0)..a_(n,xx)] -- ring on  A times CC^d 
Rb = QQ[b_(0,0)..b_(n,xx)] -- ring on  B times CC^d 
Rc = QQ[c_(0,0)..c_(n,xx)] -- ring on  C times CC^d 
Rd = QQ[d_(0,0)..d_(n,xx)] -- ring on  C times CC^d 
Re = QQ[e_(0,0)..e_(n,xx)] -- ring on  C times CC^d 


Rx = QQ[x_(0,0,0,0,0)..x_(n,n,n,n,n)] -- ring on A \otimes B\otimes C
R = Ra**Rb**Rc**Rd**Re**Rx  -- big ring

--[1, 2, 3, 4, 5, 6], [1, 2, 3, 5, 4, 6], [1, 3, 2, 4, 5, 6], [1, 3, 2, 5, 4, 6], [1, 4, 2, 5, 3, 6]

-- the following matrices should have their column indices be determined by the 
--  columns of some triple of tableaux, and ideally there should be a function that makes 
--  them dynamically
A1= matrix apply(3,i -> apply({0,1,2}, j-> a_(i,j) ))
A2= matrix apply(3,i -> apply({3,4,5}, j-> a_(i,j) ))
A3= matrix apply(3,i -> apply({6,7,8}, j-> a_(i,j) ))

B1= matrix apply(3,i -> apply({0,3,6}, j-> b_(i,j) ))
B2= matrix apply(3,i -> apply({1,4,7}, j-> b_(i,j) ))
B3= matrix apply(3,i -> apply({2,5,8}, j-> b_(i,j) ))

C1= matrix apply(3,i -> apply({0,4,8}, j-> c_(i,j) ))
C2= matrix apply(3,i -> apply({1,5,6}, j-> c_(i,j) ))
C3= matrix apply(3,i -> apply({2,3,7}, j-> c_(i,j) ))

-- this should be held as a product and not expanded for larger problems
time T = det(A1)*det(A2)*det(A3)*det(B1)*det(B2)*det(B3)*det(C1)*det(C2)*det(C3);
--T = hold det(A1)*det(A2)*det(A3)*det(B1)*det(B2)*det(B3)*det(C1)*det(C2)*det(C3);

-- this loop builds the poylnomial in Rx, one factor at a time.  This should also be held as a product.
F=T
for d from 0 to xx do F= sum(cc+1,k-> sum(bb+1, j-> sum(cc+1,i->x_(i,j,k)*diff(a_(i,d)*b_(j,d)*c_(k,d),F))));
F

-- there should also be an evaluation option, instead of multiplying by each x_I
 -- evauate at each stage.


--------------------------------
--------------------------------
restart
-- the polynomial is a 3x3x3 degree 6 invariant
xx = 5
aa = 2
bb = 2 
cc = 2

Ra = QQ[a_(0,0)..a_(aa,xx)] -- ring on  A times CC^d 
Rb = QQ[b_(0,0)..b_(bb,xx)] -- ring on  B times CC^d 
Rc = QQ[c_(0,0)..c_(cc,xx)] -- ring on  C times CC^d 
Rx = QQ[x_(0,0,0)..x_(aa,bb,cc)] -- ring on A \otimes B\otimes C
R = Ra**Rb**Rc**Rx  -- big ring


-- the following matrices should have their column indices be determined by the 
--  columns of some triple of tableaux, and ideally there should be a function that makes 
--  them dynamically
A1= matrix apply(3,i -> apply({0,1,2}, j-> a_(i,j) ))
A2= matrix apply(3,i -> apply({3,4,5}, j-> a_(i,j) ))


B1= matrix apply(3,i -> apply({0,2,4}, j-> b_(i,j) ))
B2= matrix apply(3,i -> apply({1,3,5}, j-> b_(i,j) ))

C1= matrix apply(3,i -> apply({0,1,3}, j-> c_(i,j) ))
C2= matrix apply(3,i -> apply({2,4,5}, j-> c_(i,j) ))

-- this should be held as a product and not expanded for larger problems
time T = det(A1)*det(A2)*det(B1)*det(B2)*det(C1)*det(C2);
F=T;
size T
time for d from 0 to xx do( print(size F); F= sum(cc+1,k-> sum(bb+1, j-> sum(cc+1,i->x_(i,j,k)*diff(a_(i,d)*b_(j,d)*c_(k,d),F))));)
size F
toString F

