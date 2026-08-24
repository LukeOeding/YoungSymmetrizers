---Young Symmetrizers for 6 factors
--- tableaux found previously... 
if class(checks)===Symbol then checks=true;
syt2 = {{{1,2}}};
syt4 = {{{1,2},{3,4}},{{1,3},{2,4}}};
syt6 = {{{1,2},{3,4},{5,6}}, {{1,2},{3,5},{4,6}},{{1,3},{2,4},{5,6}}, {{1,3},{2,5},{4,6}}, {{1,4},{2,5},{3,6}}};
fill8={{1,2,3,4,5,6,7,8},{1,2,3,4,5,7,6,8},{1,2,3,5,4,6,7,8},{1,2,3,5,4,7,6,8},{1,2,3,6,4,7,5,8},{1,3,2,4,5,6,7,8},{1,3,2,4,5,7,6,8},{1,3,2,5,4,6,7,8},{1,3,2,5,4,7,6,8},{1,3,2,6,4,7,5,8},{1,4,2,5,3,6,7,8},{1,4,2,5,3,7,6,8},{1,4,2,6,3,7,5,8},{1,5,2,6,3,7,4,8}};
m8to2row = L -> {{L#0,L#1}, {L#2,L#3}, {L#4,L#5}, {L#6,L#7}};
syt8 = m8to2row\fill8;
goodTabList8 = {{0,0,2,9,12},{0,0,2,12,9},{0,0,9,2,12},{0,0,9,12,2},{0,0,12,2,9},{0,2,0,9,12},{0,2,9,0,12},{0,2,9,12,0},{0,9,0,2,12},{0,9,0,12,2},{0,9,2,0,12},{0,9,2,12,0},{0,9,12,0,2},{0,9,12,2,0},{0,12,0,2,9},{2,0,0,9,12},{2,0,9,0,12},{2,0,9,12,0},{2,9,0,0,12},{2,9,0,12,0},{2,9,12,0,0},{9,0,0,2,12},{9,0,0,12,2},{9,0,2,0,12},{9,0,2,12,0},{9,0,12,0,2},{9,12,0,0,2},{12,0,0,2,9},{0,2,3,12,13},{0,2,3,13,12},{0,2,12,3,13},{0,2,12,13,3},{0,2,13,3,12},{0,3,2,12,13},{0,0,12,12,12},{0,3,6,9,13}};
mgTabList8 = {{0,0,2,9,12},{0,0,9,2,12},{0,0,9,12,2},{0,2,0,9,12},{0,2,9,0,12},{0,2,9,12,0},{0,9,0,2,12},{0,9,0,12,2},{0,9,2,0,12},{0,9,2,12,0},{0,9,12,0,2},{0,9,12,2,0},{2,0,0,9,12},{2,0,9,0,12},{2,0,9,12,0},{2,9,0,0,12},{2,9,0,12,0},{2,9,12,0,0},{9,0,0,2,12},{9,0,2,0,12},{0,3,6,9,13}};
fill10= {{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, {1, 2, 3, 4, 5, 6, 7, 9, 8, 10}, {1, 2, 3, 4, 5, 7, 6, 8, 9, 10}, {1, 2, 3, 4, 5, 7, 6, 9, 8, 10}, {1, 2, 3, 4, 5, 8, 6, 9, 7, 10}, {1, 2, 3, 5, 4, 6, 7, 8, 9, 10}, {1, 2, 3, 5, 4, 6, 7, 9, 8, 10}, {1, 2, 3, 5, 4, 7, 6, 8, 9, 10}, {1, 2, 3, 5, 4, 7, 6, 9, 8, 10}, {1, 2, 3, 5, 4, 8, 6, 9, 7, 10}, {1, 2, 3, 6, 4, 7, 5, 8, 9, 10}, {1, 2, 3, 6, 4, 7, 5, 9, 8, 10}, {1, 2, 3, 6, 4, 8, 5, 9, 7, 10}, {1, 2, 3, 7, 4, 8, 5, 9, 6, 10}, {1, 3, 2, 4, 5, 6, 7, 8, 9, 10}, {1, 3, 2, 4, 5, 6, 7, 9, 8, 10}, {1, 3, 2, 4, 5, 7, 6, 8, 9, 10}, {1, 3, 2, 4, 5, 7, 6, 9, 8, 10}, {1, 3, 2, 4, 5, 8, 6, 9, 7, 10}, {1, 3, 2, 5, 4, 6, 7, 8, 9, 10}, {1, 3, 2, 5, 4, 6, 7, 9, 8, 10}, {1, 3, 2, 5, 4, 7, 6, 8, 9, 10}, {1, 3, 2, 5, 4, 7, 6, 9, 8, 10}, {1, 3, 2, 5, 4, 8, 6, 9, 7, 10}, {1, 3, 2, 6, 4, 7, 5, 8, 9, 10}, {1, 3, 2, 6, 4, 7, 5, 9, 8, 10}, {1, 3, 2, 6, 4, 8, 5, 9, 7, 10}, {1, 3, 2, 7, 4, 8, 5, 9, 6, 10}, {1, 4, 2, 5, 3, 6, 7, 8, 9, 10}, {1, 4, 2, 5, 3, 6, 7, 9, 8, 10}, {1, 4, 2, 5, 3, 7, 6, 8, 9, 10}, {1, 4, 2, 5, 3, 7, 6, 9, 8, 10}, {1, 4, 2, 5, 3, 8, 6, 9, 7, 10}, {1, 4, 2, 6, 3, 7, 5, 8, 9, 10}, {1, 4, 2, 6, 3, 7, 5, 9, 8, 10}, {1, 4, 2, 6, 3, 8, 5, 9, 7, 10}, {1, 4, 2, 7, 3, 8, 5, 9, 6, 10}, {1, 5, 2, 6, 3, 7, 4, 8, 9, 10}, {1, 5, 2, 6, 3, 7, 4, 9, 8, 10}, {1, 5, 2, 6, 3, 8, 4, 9, 7, 10}, {1, 5, 2, 7, 3, 8, 4, 9, 6, 10}, {1, 6, 2, 7, 3, 8, 4, 9, 5, 10}};
m10to2row = L -> {{L#0,L#1}, {L#2,L#3}, {L#4,L#5}, {L#6,L#7}, {L#8,L#9}};
syt10 = m10to2row\fill10;
m12to2row = L -> {{L#0,L#1}, {L#2,L#3}, {L#4,L#5}, {L#6,L#7}, {L#8,L#9}, {L#10,L#11}};
psyt12 = {{1,2,3,4,5,6,7,8,9,10,11,12},{1,2,3,4,5,6,7,8,9,11,10,12},{1,2,3,4,5,6,7,9,8,10,11,12},{1,2,3,4,5,6,7,9,8,11,10,12},{1,2,3,4,5,6,7,10,8,11,9,12},{1,2,3,4,5,7,6,8,9,10,11,12},{1,2,3,4,5,7,6,8,9,11,10,12},{1,2,3,4,5,7,6,9,8,10,11,12},{1,2,3,4,5,7,6,9,8,11,10,12},{1,2,3,4,5,7,6,10,8,11,9,12},{1,2,3,4,5,8,6,9,7,10,11,12},{1,2,3,4,5,8,6,9,7,11,10,12},{1,2,3,4,5,8,6,10,7,11,9,12},{1,2,3,4,5,9,6,10,7,11,8,12},{1,2,3,5,4,6,7,8,9,10,11,12},{1,2,3,5,4,6,7,8,9,11,10,12},{1,2,3,5,4,6,7,9,8,10,11,12},{1,2,3,5,4,6,7,9,8,11,10,12},{1,2,3,5,4,6,7,10,8,11,9,12},{1,2,3,5,4,7,6,8,9,10,11,12},{1,2,3,5,4,7,6,8,9,11,10,12},{1,2,3,5,4,7,6,9,8,10,11,12},{1,2,3,5,4,7,6,9,8,11,10,12},{1,2,3,5,4,7,6,10,8,11,9,12},{1,2,3,5,4,8,6,9,7,10,11,12},{1,2,3,5,4,8,6,9,7,11,10,12},{1,2,3,5,4,8,6,10,7,11,9,12},{1,2,3,5,4,9,6,10,7,11,8,12},{1,2,3,6,4,7,5,8,9,10,11,12},{1,2,3,6,4,7,5,8,9,11,10,12},{1,2,3,6,4,7,5,9,8,10,11,12},{1,2,3,6,4,7,5,9,8,11,10,12},{1,2,3,6,4,7,5,10,8,11,9,12},{1,2,3,6,4,8,5,9,7,10,11,12},{1,2,3,6,4,8,5,9,7,11,10,12},{1,2,3,6,4,8,5,10,7,11,9,12},{1,2,3,6,4,9,5,10,7,11,8,12},{1,2,3,7,4,8,5,9,6,10,11,12},{1,2,3,7,4,8,5,9,6,11,10,12},{1,2,3,7,4,8,5,10,6,11,9,12},{1,2,3,7,4,9,5,10,6,11,8,12},{1,2,3,8,4,9,5,10,6,11,7,12},{1,3,2,4,5,6,7,8,9,10,11,12},{1,3,2,4,5,6,7,8,9,11,10,12},{1,3,2,4,5,6,7,9,8,10,11,12},{1,3,2,4,5,6,7,9,8,11,10,12},{1,3,2,4,5,6,7,10,8,11,9,12},{1,3,2,4,5,7,6,8,9,10,11,12},{1,3,2,4,5,7,6,8,9,11,10,12},{1,3,2,4,5,7,6,9,8,10,11,12},{1,3,2,4,5,7,6,9,8,11,10,12},{1,3,2,4,5,7,6,10,8,11,9,12},{1,3,2,4,5,8,6,9,7,10,11,12},{1,3,2,4,5,8,6,9,7,11,10,12},{1,3,2,4,5,8,6,10,7,11,9,12},{1,3,2,4,5,9,6,10,7,11,8,12},{1,3,2,5,4,6,7,8,9,10,11,12},{1,3,2,5,4,6,7,8,9,11,10,12},{1,3,2,5,4,6,7,9,8,10,11,12},{1,3,2,5,4,6,7,9,8,11,10,12},{1,3,2,5,4,6,7,10,8,11,9,12},{1,3,2,5,4,7,6,8,9,10,11,12},{1,3,2,5,4,7,6,8,9,11,10,12},{1,3,2,5,4,7,6,9,8,10,11,12},{1,3,2,5,4,7,6,9,8,11,10,12},{1,3,2,5,4,7,6,10,8,11,9,12},{1,3,2,5,4,8,6,9,7,10,11,12},{1,3,2,5,4,8,6,9,7,11,10,12},{1,3,2,5,4,8,6,10,7,11,9,12},{1,3,2,5,4,9,6,10,7,11,8,12},{1,3,2,6,4,7,5,8,9,10,11,12},{1,3,2,6,4,7,5,8,9,11,10,12},{1,3,2,6,4,7,5,9,8,10,11,12},{1,3,2,6,4,7,5,9,8,11,10,12},{1,3,2,6,4,7,5,10,8,11,9,12},{1,3,2,6,4,8,5,9,7,10,11,12},{1,3,2,6,4,8,5,9,7,11,10,12},{1,3,2,6,4,8,5,10,7,11,9,12},{1,3,2,6,4,9,5,10,7,11,8,12},{1,3,2,7,4,8,5,9,6,10,11,12},{1,3,2,7,4,8,5,9,6,11,10,12},{1,3,2,7,4,8,5,10,6,11,9,12},{1,3,2,7,4,9,5,10,6,11,8,12},{1,3,2,8,4,9,5,10,6,11,7,12},{1,4,2,5,3,6,7,8,9,10,11,12},{1,4,2,5,3,6,7,8,9,11,10,12},{1,4,2,5,3,6,7,9,8,10,11,12},{1,4,2,5,3,6,7,9,8,11,10,12},{1,4,2,5,3,6,7,10,8,11,9,12},{1,4,2,5,3,7,6,8,9,10,11,12},{1,4,2,5,3,7,6,8,9,11,10,12},{1,4,2,5,3,7,6,9,8,10,11,12},{1,4,2,5,3,7,6,9,8,11,10,12},{1,4,2,5,3,7,6,10,8,11,9,12},{1,4,2,5,3,8,6,9,7,10,11,12},{1,4,2,5,3,8,6,9,7,11,10,12},{1,4,2,5,3,8,6,10,7,11,9,12},{1,4,2,5,3,9,6,10,7,11,8,12},{1,4,2,6,3,7,5,8,9,10,11,12},{1,4,2,6,3,7,5,8,9,11,10,12},{1,4,2,6,3,7,5,9,8,10,11,12},{1,4,2,6,3,7,5,9,8,11,10,12},{1,4,2,6,3,7,5,10,8,11,9,12},{1,4,2,6,3,8,5,9,7,10,11,12},{1,4,2,6,3,8,5,9,7,11,10,12},{1,4,2,6,3,8,5,10,7,11,9,12},{1,4,2,6,3,9,5,10,7,11,8,12},{1,4,2,7,3,8,5,9,6,10,11,12},{1,4,2,7,3,8,5,9,6,11,10,12},{1,4,2,7,3,8,5,10,6,11,9,12},{1,4,2,7,3,9,5,10,6,11,8,12},{1,4,2,8,3,9,5,10,6,11,7,12},{1,5,2,6,3,7,4,8,9,10,11,12},{1,5,2,6,3,7,4,8,9,11,10,12},{1,5,2,6,3,7,4,9,8,10,11,12},{1,5,2,6,3,7,4,9,8,11,10,12},{1,5,2,6,3,7,4,10,8,11,9,12},{1,5,2,6,3,8,4,9,7,10,11,12},{1,5,2,6,3,8,4,9,7,11,10,12},{1,5,2,6,3,8,4,10,7,11,9,12},{1,5,2,6,3,9,4,10,7,11,8,12},{1,5,2,7,3,8,4,9,6,10,11,12},{1,5,2,7,3,8,4,9,6,11,10,12},{1,5,2,7,3,8,4,10,6,11,9,12},{1,5,2,7,3,9,4,10,6,11,8,12},{1,5,2,8,3,9,4,10,6,11,7,12},{1,6,2,7,3,8,4,9,5,10,11,12},{1,6,2,7,3,8,4,9,5,11,10,12},{1,6,2,7,3,8,4,10,5,11,9,12},{1,6,2,7,3,9,4,10,5,11,8,12},{1,6,2,8,3,9,4,10,5,11,7,12},{1,7,2,8,3,9,4,10,5,11,6,12}};
syt12 = m12to2row\ psyt12;
if class(KK)===Symbol then KK=QQ;
--KK = QQ;
X = KK[x_(0,0,0,0,0,0)..x_(1,1,1,1,1,1)];
primaryInvariants = {};
rndPt = () -> map(KK, X, random(KK^1, KK^(numgens X))) 
rp = rndPt();
perm = P-> map(X,X, apply(toList( (0,0,0,0,0,0)..(1,1,1,1,1,1)), ll -> x_(toSequence apply(P, p-> ll#p))));
-- cycle = map(X,X, flatten flatten flatten flatten flatten apply(2, i1-> apply(2, i2->apply(2, i3->apply(2, i4->apply(2, i5-> apply(2, i6-> x_(i2,i3,i4,i5,i6,i1) ) ) ) )) ));
-- cycle2 = map(X,X, flatten flatten flatten flatten flatten apply(2, i1-> apply(2, i2->apply(2, i3->apply(2, i4->apply(2, i5-> apply(2, i6-> x_(i3,i4,i5,i6,i1,i2) ) ) ) )) ));
-- cycle3 = map(X,X, flatten flatten flatten flatten flatten apply(2, i1-> apply(2, i2->apply(2, i3->apply(2, i4->apply(2, i5-> apply(2, i6-> x_(i4,i5,i6,i1,i2,i3) ) ) ) )) ));
-- cycle4 = map(X,X, flatten flatten flatten flatten flatten apply(2, i1-> apply(2, i2->apply(2, i3->apply(2, i4->apply(2, i5-> apply(2, i6-> x_(i5,i6,i1,i2,i3,i4) ) ) ) )) ));
-- cycle5 = map(X,X, flatten flatten flatten flatten flatten apply(2, i1-> apply(2, i2->apply(2, i3->apply(2, i4->apply(2, i5-> apply(2, i6-> x_(i6,i1,i2,i3,i4,i5) ) ) ) )) ));
-- swap = map(X,X, flatten flatten flatten flatten flatten apply(2, i1-> apply(2, i2->apply(2, i3->apply(2, i4->apply(2, i5-> apply(2, i6-> x_(i2,i1,i3,i4,i5,i6) ) ) ) )) ));

list2mats = (L,a) -> matrix apply(2, i-> apply( L, j-> a_(i,j-1)));
tab2poly = (L)->(
 ta:=L#0; tb:=L#1;tc:=L#2;td:=L#3;te:=L#4;tf:=L#5;
 dg := length flatten ta;
 R := X[a_(0,0)..a_(1,dg-1),b_(0,0)..b_(1,dg-1),c_(0,0)..c_(1,dg-1),d_(0,0)..d_(1,dg-1),e_(0,0)..e_(1,dg-1),f_(0,0)..f_(1,dg-1), ];
 A := apply(ta, L -> list2mats(L, a));
 B := apply(tb, L -> list2mats(L, b));
 C := apply(tc, L -> list2mats(L, c));
 D := apply(td, L -> list2mats(L, d));
 E := apply(te, L -> list2mats(L, e));
 F := apply(tf, L -> list2mats(L, f));
 usedA := {};usedB := {};usedC := {};usedD := {};usedE := {};usedF := {};
 Fp :=1_R; T :=1_R;
 for z to dg-1 do time ( T=product(
 (for i to sub(dg/2-1,ZZ) list if isMember(z + 1,ta#i) and not isMember(i,usedA) then (usedA = usedA|{i}; det A_i) else continue)|
 (for i to sub(dg/2-1,ZZ) list if isMember(z + 1,tb#i) and not isMember(i,usedB) then (usedB = usedB|{i}; det B_i) else continue)|
 (for i to sub(dg/2-1,ZZ) list if isMember(z + 1,tc#i) and not isMember(i,usedC) then (usedC = usedC|{i}; det C_i) else continue)|
 (for i to sub(dg/2-1,ZZ) list if isMember(z + 1,td#i) and not isMember(i,usedD) then (usedD = usedD|{i}; det D_i) else continue)|
 (for i to sub(dg/2-1,ZZ) list if isMember(z + 1,te#i) and not isMember(i,usedE) then (usedE = usedE|{i}; det E_i) else continue)|
 (for i to sub(dg/2-1,ZZ) list if isMember(z + 1,tf#i) and not isMember(i,usedF) then (usedF = usedF|{i}; det F_i) else continue));
 Fp = Fp*T; -- you could insert a time before this...
 Fp= sum(2,n->sum(2,m-> sum(2, l-> sum(2,k-> sum(2, j-> sum(2,i->x_(i,j,k,l,m,n)*diff(a_(i,z)*b_(j,z)*c_(k,z)*d_(l,z)*e_(m,z)*f_(n,z),Fp)))))));
 ); -- after many checks (parallel sum, making auxillary arrays, using diff over a tensor, etc), this seems to be the most efficient way to execute the sum.
 sub(Fp,X)
);

tab2polyEval = (L,rules)->(
 ta:=L#0; tb:=L#1;tc:=L#2;td:=L#3;te:=L#4;tf:=L#5;
 dg := length flatten ta;
 R := X[a_(0,0)..a_(1,dg-1),b_(0,0)..b_(1,dg-1),c_(0,0)..c_(1,dg-1),d_(0,0)..d_(1,dg-1),e_(0,0)..e_(1,dg-1),f_(0,0)..f_(1,dg-1), ];
 A := apply(ta, L -> list2mats(L, a));
 B := apply(tb, L -> list2mats(L, b));
 C := apply(tc, L -> list2mats(L, c));
 D := apply(td, L -> list2mats(L, d));
 E := apply(te, L -> list2mats(L, e));
 F := apply(tf, L -> list2mats(L, f));
 usedA := {};usedB := {};usedC := {};usedD := {};usedE := {};usedF := {};
 Fp :=1_R; T :=1_R; -- took out "time" for the z-for loop
 for z to dg-1 do  ( T=product(
 (for i to sub(dg/2-1,ZZ) list if isMember(z + 1,ta#i) and not isMember(i,usedA) then (usedA = usedA|{i}; det A_i) else continue)|
 (for i to sub(dg/2-1,ZZ) list if isMember(z + 1,tb#i) and not isMember(i,usedB) then (usedB = usedB|{i}; det B_i) else continue)|
 (for i to sub(dg/2-1,ZZ) list if isMember(z + 1,tc#i) and not isMember(i,usedC) then (usedC = usedC|{i}; det C_i) else continue)|
 (for i to sub(dg/2-1,ZZ) list if isMember(z + 1,td#i) and not isMember(i,usedD) then (usedD = usedD|{i}; det D_i) else continue)|
 (for i to sub(dg/2-1,ZZ) list if isMember(z + 1,te#i) and not isMember(i,usedE) then (usedE = usedE|{i}; det E_i) else continue)|
 (for i to sub(dg/2-1,ZZ) list if isMember(z + 1,tf#i) and not isMember(i,usedF) then (usedF = usedF|{i}; det F_i) else continue));
 Fp = Fp*T; -- you could insert a time before this...
 Fp= sum(2,n->sum(2,m-> sum(2, l-> sum(2,k-> sum(2, j-> sum(2,i->sub(x_(i,j,k,l,m,n),rules)*diff(a_(i,z)*b_(j,z)*c_(k,z)*d_(l,z)*e_(m,z)*f_(n,z),Fp)))))));
 ); -- after many checks (parallel sum, making auxillary arrays, using diff over a tensor, etc), this seems to be the most efficient way to execute the sum.
 sub(Fp,X)
);

-- tab2polyEval = (L,rules)->(
--  ta:=L#0; tb:=L#1;tc:=L#2;td:=L#3;te:=L#4;tf:=L#5;
--  dg := length flatten ta;
--  R := X[a_(0,0)..a_(1,dg-1),b_(0,0)..b_(1,dg-1),c_(0,0)..c_(1,dg-1),d_(0,0)..d_(1,dg-1),e_(0,0)..e_(1,dg-1),f_(0,0)..f_(1,dg-1), ];
--  A := apply(ta, L -> list2mats(L, a));
--  B := apply(tb, L -> list2mats(L, b));
--  C := apply(tc, L -> list2mats(L, c));
--  D := apply(td, L -> list2mats(L, d));
--  E := apply(te, L -> list2mats(L, e));
--  F := apply(tf, L -> list2mats(L, f));
--  usedA := {};usedB := {};usedC := {};usedD := {};usedE := {};usedF := {};
--  Fp :=1_R; T :=1_R;
--  for z to dg-1 do ( T=product(
--  (for i to sub(dg/2-1,ZZ) list if isMember(z + 1,ta#i) and not isMember(i,usedA) then (usedA = usedA|{i}; det A_i) else continue)|
--  (for i to sub(dg/2-1,ZZ) list if isMember(z + 1,tb#i) and not isMember(i,usedB) then (usedB = usedB|{i}; det B_i) else continue)|
--  (for i to sub(dg/2-1,ZZ) list if isMember(z + 1,tc#i) and not isMember(i,usedC) then (usedC = usedC|{i}; det C_i) else continue)|
--  (for i to sub(dg/2-1,ZZ) list if isMember(z + 1,td#i) and not isMember(i,usedD) then (usedD = usedD|{i}; det D_i) else continue)|
--  (for i to sub(dg/2-1,ZZ) list if isMember(z + 1,te#i) and not isMember(i,usedE) then (usedE = usedE|{i}; det E_i) else continue)|
--  (for i to sub(dg/2-1,ZZ) list if isMember(z + 1,tf#i) and not isMember(i,usedF) then (usedF = usedF|{i}; det F_i) else continue));
--  Fp = Fp*T; -- you could insert a time before this...
--  Fp= sum(2,n->sum(2,m-> sum(2, l-> sum(2,k-> sum(2, j-> sum(2,i->sub(x_(i,j,k,l,m,n),rules)*diff(a_(i,z)*b_(j,z)*c_(k,z)*d_(l,z)*e_(m,z)*f_(n,z),Fp)))))));
--  ); -- after many checks (parallel sum, making auxillary arrays, using diff over a tensor, etc), this seems to be the most efficient way to execute the sum.
--  sub(Fp,X)
-- );

--- points 
gx = gens X
pseudoCartan = apply(2^5,i-> gx_(2^6-1-i) => gx_i)
randomLine = apply(2^6,i-> gx_(i) => random(KK)*gx_0 + random(KK)*gx_15 + random(KK)*gx_31)
randomPoint = apply(2^6,i-> gx_(i) => random(KK))
--- tableaux
end;



inv2 = tab2poly({syt2#0,syt2#0,syt2#0,syt2#0,syt2#0,syt2#0});
F2 = {inv2};
I2 = ideal inv2
--inv40 = tab2poly({syt4#0,syt4#0,syt4#0,syt4#0,syt4#0,syt4#0});
--factor inv40 -- this is the square of the degree 2 invariant
--inv4 = tab2poly({syt4#0,syt4#0,syt4#0,syt4#0,syt4#0,syt4#1});
--factor inv4 -- this is the square of the degree 2 invariant
-- so... this one is new, not in the ideal of the degree 2 invariant
inv42 = tab2poly({syt4#0,syt4#0,syt4#0,syt4#0,syt4#1,syt4#1}); 
f43 =  tab2poly({syt4#0,syt4#0,syt4#0,syt4#1,syt4#1,syt4#1}); 

pointLibrary = for i to 16 list apply(2^6,i-> gx_(i) => random(KK));
time newCol = matrix apply(pointLibrary, p ->  {tab2polyEval({syt4#0,syt4#0,syt4#0,syt4#0,syt4#1, syt4#1}, p)});  
M4 = newCol; rk=1;
for P in drop(uniquePermutations {0,0,0,0,1,1},1) do(
    time newCol = matrix apply(pointLibrary, p ->  {tab2polyEval(apply(6, i-> syt4#(P#i)), p)});  
    M4 = M4|newCol;
    print(toString(P)|" rank M4 = "|toString(rank M4)));
toString uniquePermutations {0,0,0,0,1,1}
M4
ker M4 -- notice that a signed sum of the first 10 multi-fillings give the remining 5.

myPerms42 ={{0,1,2,3,4,5}, {0,1,2,4,3,5}, {0,1,4,3,2,5}, {0,4,2,3,1,5}, {0,1,2,4,5,3}, {0,1,4,2,5,3}, {0,4,2,1,5,3}, {0,1,4,5,2,3}, {0,4,1,5,2,3}, {0,4,5,1,2,3}, {4,1,2,3,0,5}, {4,1,2,0,5,3}, {4,1,0,5,2,3}, {4,0,5,1,2,3}, {4,5,0,1,2,3}};
time deg4List = for P in myPerms42 list ( (perm P) inv42);
--apply(deg4List, ff -> #factor ff)
betti ideal deg4List
betti mingens ideal deg4List
betti mingens ideal deg4List_{0..9} -- notice that the first 10 are maximally independent
--factor sum deg4List
F4 ={};rk =0; pmList = for P in myPerms42 list ( tmpPerm = perm P; tmp = tmpPerm inv42; if rank source mingens ideal (F4|{tmp})>rk then (F4 = F4|{tmp}; rk=rk + 1;P ) else continue);
LL = {0,0,0,0,1,1};
netList  for P in pmList list apply(P, i->LL_i)
--time F42 = for P in myPerms42 list (tmpPerm = perm P; tmpPerm inv42);
I4 = ideal F4;
betti I4
betti mingens(I4+I2)
--rank source mingens I4
if checks ==true then time betti mingens (I2 + I4)
-- of course you can make other invariants, but they won't be independent.
--inv43 = tab2poly({syt4#0,syt4#0,syt4#0,syt4#1,syt4#1,syt4#1});
--myPerms33 = {{0,1,2,3,4,5},{0,1,3,2,4,5},{0,1,3,4,2,5},{0,1,3,4,4,2},{0,3,1,2,4,5},{0,3,1,4,2,5},{0,3,1,4,5,2},{0,3,4,1,2,5},{0,3,4,1,5,2},{0,3,4,5,1,2},{3,0,1,2,4,5},{3,0,1,4,2,5},{3,0,1,4,5,2},{3,0,4,1,2,5},{3,0,4,1,5,2},{3,0,4,5,1,2},{3,4,0,1,2,5},{3,4,0,1,4,2},{3,4,0,5,1,2},{3,4,5,0,1,2}}
--time F43 = for P in myPerms33 list (tmpPerm = perm P; tmpPerm inv43);
--I43 = ideal F43;
--betti mingens I43
--betti mingens (I2 + I42 + ideal(inv43))

-- these are the first 5 primary invariants - but let's do a randomized jacobian to make sure we think they're algebraically independent.
J4 = jacobian I4;
phi = rndPt();
rank phi J4 
if checks ==true then betti mingens (I2 + I4)
J24 = jacobian (I2 + I4);
phi = rndPt();
rank phi J24 

primaryInvariants = F2|F4;

-- note that we have 10 + 1 invariants in degree 4, which matches the hilbert function
--elapsedTime inv6 = tab2poly({syt6#0, syt6#1, syt6#2, syt6#3, syt6#4,syt6#4}); -- 26s 
--#unique permutations{0,1,2,3,4,4}
--#terms inv6 -- 68192
--elapsedTime inv60 = tab2poly({syt6#0, syt6#0, syt6#0, syt6#0, syt6#0,syt6#0}); -- 3.6s 
--#terms inv60 -- 5984, and it's the cube of the degree 2 invariant... 
-- note 000001 and similar produce the same polynomial.. 
--elapsedTime inv62 = tab2poly({syt6#0, syt6#0, syt6#0, syt6#0, syt6#1,syt6#1}); -- 3.6s 
--#terms inv62 -- 24416
--factor inv62 -- this is a quadratic times a quartic
--betti mingens (I2 + ideal inv62 )
netList syt6
elapsedTime inv62 = tab2poly({syt6#0, syt6#0, syt6#0, syt6#0, syt6#3,syt6#3}); -- 14s 
-- this is not the most efficient way to do this:
--time F62 = apply(unique permutations {0,0,0,0,3,3}, p-> tab2poly({syt6#(p#0),syt6#(p#1),syt6#(p#2),syt6#(p#3),syt6#(p#4),syt6#(p#5)}));
--myPerms ={{0,1,2,3,4,5}, {0,1,2,4,3,5}, {0,1,4,3,2,5}, {0,4,2,3,1,5}, {4,1,2,3,0,5}, {0,1,2,4,5,3}, {0,1,4,2,5,3}, {0,4,2,1,5,3}, {4,1,2,0,5,3}, {0,1,4,5,2,3}, {0,4,1,5,2,3}, {4,1,0,5,2,3}, {0,4,5,1,2,3}, {4,0,5,1,2,3}, {4,5,0,1,2,3}}
time F62 = for P in myPerms42 list (tmpPerm = perm P; tmpPerm inv62);
I62 = ideal F62;
if checks ==true then time betti mingens I62 -- 15, so they're all independent

if checks ==true then time betti mingens (I2 + I4 + I62) -- and they're independent from the lower degree terms
--elapsedTime inv62b = tab2poly({syt6#0, syt6#0, syt6#0, syt6#0, syt6#4,syt6#4}); -- this one is not independent from the other one
--time betti mingens (I62 + ideal inv62b)
--elapsedTime inv62b = tab2poly({syt6#2, syt6#2, syt6#2, syt6#2, syt6#3,syt6#3}); -- check for independence-- this one factors
--factor inv62b
--elapsedTime inv62b = tab2poly({syt6#2, syt6#2, syt6#2, syt6#2, syt6#4,syt6#4}); -- check for independence
--time betti mingens (I62 + ideal inv62b)
--elapsedTime inv62b = tab2poly({syt6#3, syt6#3, syt6#3, syt6#3, syt6#4,syt6#4}); -- check for independence
--time betti mingens (I62 + ideal inv62b)
--#factor inv62b -- it factors... 
elapsedTime inv63= tab2poly({syt6#0, syt6#0, syt6#0, syt6#4, syt6#4,syt6#4}); -- check for independence
if checks == true then time betti mingens (I62 + ideal inv63)
if checks == true then time betti mingens (I2 + I4 + I62 + ideal inv63) -- and it is independent from the lower degree terms
up = unique permutations {0,0,0,4,4,4}
myPerms33 = {{0,1,2,3,4,5},{0,1,3,2,4,5},{0,1,3,4,2,5},{0,1,3,4,4,2},{0,3,1,2,4,5},{0,3,1,4,2,5},{0,3,1,4,5,2},{0,3,4,1,2,5},{0,3,4,1,5,2},{0,3,4,5,1,2},{3,0,1,2,4,5},{3,0,1,4,2,5},{3,0,1,4,5,2},{3,0,4,1,2,5},{3,0,4,1,5,2},{3,0,4,5,1,2},{3,4,0,1,2,5},{3,4,0,1,4,2},{3,4,0,5,1,2},{3,4,5,0,1,2}}
time F63 = for P in myPerms33 list (tmpPerm = perm P; tmpPerm inv63);
F63 ={};rk =0;
pmList = for P in myPerms33 list ( tmpPerm = perm P; tmp = tmpPerm inv63; if rank source mingens ideal (F63|{tmp})>rk then (F63 = F63|{tmp}; rk=rk + 1;P ) else continue)
I63 = ideal F63;
LL = {0,0,0,1,1,1}
netList for P in pmList list apply(P, i->LL_i)

if checks ==true then time betti mingens I63 -- 10, so there's 10 independent among the 15
if checks ==true then time betti mingens (I2 + I4 + I62 + I63) -- and they're independent from the lower degree terms
-- here is one more invariant that is of a different type than the others
elapsedTime inv622 = tab2poly({syt6#0, syt6#0, syt6#2, syt6#2, syt6#4,syt6#4}); --
if checks ==true then time betti mingens (I62 + I63 + ideal inv622)
if checks ==true then time betti mingens (I2 + I4 + I62 + I63 + ideal inv622) 
-- we can guess the permutations of the fillings that will be independent.
myPerms222 = {{0,1,2,3,4,5}, {0,2,1,3,4,5}, {0,1,2,4,3,5}, {0,2,1,4,3,5}, {0,4,2,3,1,5}}
time F622 = for P in myPerms222 list (tmpPerm = perm P; tmpPerm inv622);
I622 = ideal F622;
if checks ==true then time betti mingens I622 
if checks ==true then time betti mingens (I62 + I63 + I622)
#F62
#F63
#F622
F6 = F62|F63|F622;
I6 = I62 + I63 + I622;
if checks ==true then time betti mingens (I2 + I4 + I6) 
betti(I2 + I4 + I6)
if checks ==true then J246 = jacobian (I2 + I4 + I6);
if checks ==true then phi J246;
if checks ==true then rank phi J246
-- you can also check the jacobian one step at a time if you think that some of the polynomials will not be algebraically independent
--bx = value toString transpose basis(1,X);
--phi = rndPt();
--MJ = phi diff(bx, F2#0); rk = 1;pmi = {F2#0}; for i to #F4-1 do (newCol = phi diff(bx, F4#i); tmpM = MJ|newCol; if rank tmpM>rk then (print(i); rk = rk + 1; pmi = pmi|{F4#i}; MJ =MJ|newCol; ) else continue)
--rank MJ 
--for i to #F6-1 do (newCol = phi diff(bx, F6#i); tmpM = MJ|newCol; if rank tmpM>rk then (print(i); rk = rk + 1; pmi = pmi|{F6#i}; MJ =MJ|newCol; ) else continue)
--rank MJ 
primaryInvariants = F2|F4|F6;
#primaryInvariants
end;
restart
allowableThreads = maxAllowableThreads
checks =false;
KK = ZZ/nextPrime(100);
--allowableThreads = 1--  maxAllowableThreads-1
input"/Users/lao0004/Documents/GitHub/YoungSymmetrizers/6factor.m2";
-- standard settings -- inv62 took 15s, betti mingens I2 + I4 + I62 took -- used 73.2505s (cpu); 18.8735s (thread); 0s (gc)
-- 14 threads: -- used 74.5685s (cpu); 16.8846s (thread); 0s (gc)
-- 1 thread: -- used 75.4438s (cpu); 17.2099s (thread); 0s (gc)
-- in d egree 8 it starts to take too long to evaluate and find the polynomials.
-- We can also get independence statements from point evaluations
netList syt8
-- to find independence we don't need to be able to take a jacobian, just evaluate at a lot of points.
-- to compute jacobians, we could hold all variables constant except one and evaluate them, then take the derivative with respect to the variable of interest...
-- try to work with 70 points for now
betti (I2+ I4 + I6)
pointLibrary = for i to 20 list apply(2^6,i-> gx_(i) => random(KK));
time M6 = matrix apply(pointLibrary, p -> apply(F2|F4|F6, f-> sub(f, p)));
-- check non-zero
time tab2polyEval({syt8#0,syt8#0,syt8#0,syt8#0,syt8#1,syt8#1}, randomPoint) -- about 6 seconds per point working over F_101
binomial(6,3)
pointLibrary = for i to 22 list apply(2^6,i-> gx_(i) => random(KK));
time newCol = matrix apply(pointLibrary, p ->  {tab2polyEval({syt8#0,syt8#0,syt8#0,syt8#0,syt8#1, syt8#1}, p)}); -- expect around 
M82 = newCol;
GoodFill82 ={{0,0,0,0,1,1}}; rk =1; 
for P in drop(unique permutations {0,0,0,0,1,1},1) do (print("testing "|toString(P)|" with current rank "|toString(rk));
    time newCol = matrix apply(pointLibrary, p ->  {tab2polyEval({syt8#(P#0),syt8#(P#1),syt8#(P#2),syt8#(P#3),syt8#(P#4), syt8#(P#5)}, p)}); 
    if rank (M82|newCol)>rk then (GoodFill82 = GoodFill82|{P}; M82 = M82|newCol; rk=rk + 1; ) else continue);
-- this finished on the server, and was much faster. 
--Fill8={{0,0,0,0,1,1},{0,0,0,1,0,1},{0,0,0,1,1,0},{0,0,1,0,0,1},{0,0,1,0,1,0},{0,0,1,1,0,0},{0,1,0,0,0,1},{0,1,0,0,1,0},{0,1,0,1,0,0},{0,1,1,0,0,0}}

-- checking other fillings of format {4,2}
rk = rank(M82)
for j from 2 to nops(syt8)-1 do (print("testing "|toString(syt8#i)|" with current rank "|toString(rk));
    time newCol = matrix apply(pointLibrary, p ->  {tab2polyEval({syt8#0,syt8#0,syt8#0,syt8#0,syt8#j, syt8#j}, p)}, p)}); 
    if rank (M82|newCol)>rk then (GoodFill82 = GoodFill82|{P}; M82 = M82|newCol; rk=rk + 1; ) else continue);



time newCol = matrix apply(pointLibrary, p ->  {tab2polyEval({syt8#0,syt8#0,syt8#0,syt8#1,syt8#1, syt8#1}, p)});
M83 = newCol;
GoodFill83 ={}; rk =0; 
for P in unique permutations {0,0,0,1,1,1} do (print("testing "|toString(P)|" with current rank "|toString(rk));
    time newCol = matrix apply(pointLibrary, p ->  {tab2polyEval(apply(6, i-> syt8#(P#i)), p)}); 
    if rank (M83|newCol)>rk then (GoodFill83 = GoodFill83|{P}; M83 = M83|newCol; rk=rk + 1; ) else continue);
--- above does some numerical experiments to try to determine which fillings are good ones to work with
rank(M82|M83)
apply(11, i-> rank(M82|M83_{i}))
M82|M83
rank((M82*matrix apply(10, i-> {1}))|(M83*matrix apply(11, i-> {1})))
time inv82 = tab2poly({syt8#0,syt8#0,syt8#0,syt8#0,syt8#1, syt8#1}); -- this takes about 242s on MBP.
#terms inv82
I = I2+ I4 + I6;
--time betti mingens (I+ideal inv8) -- we don't expect this to finish
--myPerms42 ={{0,1,2,3,4,5}, {0,1,2,4,3,5}, {0,1,4,3,2,5}, {0,4,2,3,1,5}, {4,1,2,3,0,5}, {0,1,2,4,5,3}, {0,1,4,2,5,3}, {0,4,2,1,5,3}, {4,1,2,0,5,3}, {0,1,4,5,2,3}, {0,4,1,5,2,3}, {4,1,0,5,2,3}, {0,4,5,1,2,3}, {4,0,5,1,2,3}, {4,5,0,1,2,3}};
time F8 = for P in myPerms42 list ((perm P) inv82);
--F8 ={};rk =0; pmList = for P in myPerms42 list ( tmpPerm = perm P; tmp = tmpPerm inv82; if rank source mingens ideal (F8|{tmp})>rk then (F8 = F8|{tmp}; rk=rk + 1;P ) else continue);
binomial(6,2)

restart
KK = ZZ
input"/Users/lao0004/Documents/GitHub/YoungSymmetrizers/6factor.m2";

time tab2polyEval({syt8#0,syt8#0,syt8#0,syt8#0,syt8#1,syt8#1}, randomPoint) 
-- on MBP over QQ: used 0.530551s (cpu); 0.305336s (thread); 0s (gc)
-- on MBP over ZZ/101 used 0.458976s (cpu); 0.269948s (thread); 0s (gc)

pointLibrary = for i to 69 list apply(2^6,i-> gx_(i) => random(KK));
time newCol = matrix apply(pointLibrary, p ->  {tab2polyEval({syt8#0,syt8#0,syt8#0,syt8#0,syt8#1, syt8#1}, p)}); 
-- on MBP over QQ: used 34.1592s (cpu); 20.3076s (thread); 0s (gc)
-- on MBP over ZZ/101: -- used 32.0425s (cpu); 18.3652s (thread); 0s (gc)

time inv82 = tab2poly({syt8#0,syt8#0,syt8#0,syt8#0,syt8#1, syt8#1}); 