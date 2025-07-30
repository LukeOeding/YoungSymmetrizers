---Young Symmetrizers for 5 factors
--restart -- input"5factor.m2"
--KK = ZZ/nextPrime(100); --KK = QQ;
KK = QQ;
X = KK[x_(0,0,0,0,0)..x_(1,1,1,1,1)];
primaryInvariants = {};
rndPt = () -> map(KK, X, random(KK^1, KK^(numgens X))) 
rp = rndPt();
cycle = map(X,X, flatten flatten flatten flatten  apply(2, i1-> apply(2, i2->apply(2, i3->apply(2, i4->apply(2, i5-> x_(i2,i3,i4,i5,i1) ) ) ) )) );
cycle2 = map(X,X, flatten flatten flatten flatten  apply(2, i1-> apply(2, i2->apply(2, i3->apply(2, i4->apply(2, i5-> x_(i3,i4,i5,i1,i2) ) ) ) )) );
cycle3 = map(X,X, flatten flatten flatten flatten  apply(2, i1-> apply(2, i2->apply(2, i3->apply(2, i4->apply(2, i5-> x_(i4,i5,i1,i2,i3) ) ) ) )) );
cycle4 = map(X,X, flatten flatten flatten flatten  apply(2, i1-> apply(2, i2->apply(2, i3->apply(2, i4->apply(2, i5-> x_(i5,i1,i2,i3,i4) ) ) ) )) );
swap = map(X,X, flatten flatten flatten flatten  apply(2, i1-> apply(2, i2->apply(2, i3->apply(2, i4->apply(2, i5-> x_(i2,i1,i3,i4,i5) ) ) ) )) );
list2mats = (L,a) -> matrix apply(2, i-> apply( L, j-> a_(i,j-1)));
tab2poly = (L)->(
    ta:=L#0; tb:=L#1;tc:=L#2;td:=L#3;te:=L#4;
    dg := length flatten ta;
    R := X[a_(0,0)..a_(1,dg-1),b_(0,0)..b_(1,dg-1),c_(0,0)..c_(1,dg-1),d_(0,0)..d_(1,dg-1),e_(0,0)..e_(1,dg-1)];
    A := apply(ta, L -> list2mats(L, a));
    B := apply(tb, L -> list2mats(L, b));
    C := apply(tc, L -> list2mats(L, c));
    D := apply(td, L -> list2mats(L, d));
    E := apply(te, L -> list2mats(L, e));
    usedA := {};usedB := {};usedC := {};usedD := {};usedE := {};
    F :=1_R; T :=1_R;
    for z to dg-1 do time ( T=product(
        (for i to sub(dg/2-1,ZZ) list if isMember(z+1,ta#i) and not isMember(i,usedA) then (usedA = usedA|{i}; det A_i) else continue)|
        (for i to sub(dg/2-1,ZZ) list if isMember(z+1,tb#i) and not isMember(i,usedB) then (usedB = usedB|{i}; det B_i) else continue)|
        (for i to sub(dg/2-1,ZZ) list if isMember(z+1,tc#i) and not isMember(i,usedC) then (usedC = usedC|{i}; det C_i) else continue)|
        (for i to sub(dg/2-1,ZZ) list if isMember(z+1,td#i) and not isMember(i,usedD) then (usedD = usedD|{i}; det D_i) else continue)|
        (for i to sub(dg/2-1,ZZ) list if isMember(z+1,te#i) and not isMember(i,usedE) then (usedE = usedE|{i}; det E_i) else continue));
        F = F*T;
        time F= sum(2,m-> sum(2, l-> sum(2,k-> sum(2, j-> sum(2,i->x_(i,j,k,l,m)*diff(a_(i,z)*b_(j,z)*c_(k,z)*d_(l,z)*e_(m,z),F))))));
--<<<<<<< HEAD
--=======
    ); -- after many checks (parallel sum, making auxillary arrays, using diff over a tensor, etc), this seems to be the most efficient way to execute the sum.
    sub(F,X)
);

tab2polyEval = (L, rules)->(
    ta:=L#0; tb:=L#1;tc:=L#2;td:=L#3;te:=L#4;
    dg := length flatten ta;
    R := X[a_(0,0)..a_(1,dg-1),b_(0,0)..b_(1,dg-1),c_(0,0)..c_(1,dg-1),d_(0,0)..d_(1,dg-1),e_(0,0)..e_(1,dg-1)];
    A := apply(ta, L -> list2mats(L, a));
    B := apply(tb, L -> list2mats(L, b));
    C := apply(tc, L -> list2mats(L, c));
    D := apply(td, L -> list2mats(L, d));
    E := apply(te, L -> list2mats(L, e));
    usedA := {};usedB := {};usedC := {};usedD := {};usedE := {};
    F :=1_R; T :=1_R;
    for z to dg-1 do ( T=product(
        (for i to sub(dg/2-1,ZZ) list if isMember(z+1,ta#i) and not isMember(i,usedA) then (usedA = usedA|{i}; det A_i) else continue)|
        (for i to sub(dg/2-1,ZZ) list if isMember(z+1,tb#i) and not isMember(i,usedB) then (usedB = usedB|{i}; det B_i) else continue)|
        (for i to sub(dg/2-1,ZZ) list if isMember(z+1,tc#i) and not isMember(i,usedC) then (usedC = usedC|{i}; det C_i) else continue)|
        (for i to sub(dg/2-1,ZZ) list if isMember(z+1,td#i) and not isMember(i,usedD) then (usedD = usedD|{i}; det D_i) else continue)|
        (for i to sub(dg/2-1,ZZ) list if isMember(z+1,te#i) and not isMember(i,usedE) then (usedE = usedE|{i}; det E_i) else continue));
        F = F*T;
        F= sum(2,m-> sum(2, l-> sum(2,k-> sum(2, j-> sum(2,i->sub(x_(i,j,k,l,m),rules)*diff(a_(i,z)*b_(j,z)*c_(k,z)*d_(l,z)*e_(m,z),F))))));
    ); -- after many checks (parallel sum, making auxillary arrays, using diff over a tensor, etc), this seems to be the most efficient way to execute the sum.
    sub(F,X)
);

-- this doesn't seem to find much of a speed up, likely because the first tableau is already in the standard order.
tab2polyList = (ta,tb,tc,td,te,L)->(
    dg := length flatten ta;
    R := X[a_(0,0)..a_(1,dg-1),b_(0,0)..b_(1,dg-1),c_(0,0)..c_(1,dg-1),d_(0,0)..d_(1,dg-1),e_(0,0)..e_(1,dg-1)];
    A := apply(ta, L -> list2mats(L, a));
    B := apply(tb, L -> list2mats(L, b));
    C := apply(tc, L -> list2mats(L, c));
    D := apply(td, L -> list2mats(L, d));
    E := apply(te, L -> list2mats(L, e));
    usedA := {};usedB := {};usedC := {};usedD := {};usedE := {};
    F :=1_R; T :=1_R;
    for z in L do ( T=product(
        (for i to sub(dg/2-1,ZZ) list if isMember(z+1,ta#i) and not isMember(i,usedA) then (usedA = usedA|{i}; det A_i) else continue)|
        (for i to sub(dg/2-1,ZZ) list if isMember(z+1,tb#i) and not isMember(i,usedB) then (usedB = usedB|{i}; det B_i) else continue)|
        (for i to sub(dg/2-1,ZZ) list if isMember(z+1,tc#i) and not isMember(i,usedC) then (usedC = usedC|{i}; det C_i) else continue)|
        (for i to sub(dg/2-1,ZZ) list if isMember(z+1,td#i) and not isMember(i,usedD) then (usedD = usedD|{i}; det D_i) else continue)|
        (for i to sub(dg/2-1,ZZ) list if isMember(z+1,te#i) and not isMember(i,usedE) then (usedE = usedE|{i}; det E_i) else continue));
        F = F*T;
        F= sum(2,m-> sum(2, l-> sum(2,k-> sum(2, j-> sum(2,i->x_(i,j,k,l,m)*diff(a_(i,z)*b_(j,z)*c_(k,z)*d_(l,z)*e_(m,z),F))))));
    ); -- should check to see if you make auxillary variables for the derivative wrt some of the abc's and split up the computation if that speeds things up more.
    sub(F,X)
);
-- reuse the variables, just set the pair to the first variable.

--- points 
gx = gens X
pseudoCartan = apply(16,i-> gx_(31-i) => gx_i)
randomLine = apply(32,i-> gx_(i) => random(KK)*gx_0 + random(KK)*gx_15+ random(KK)*gx_31)
randomPoint = apply(32,i-> gx_(i) => random(KK))
---  tableaux
syt4 = {{{1,2},{3,4}},{{1,3},{2,4}}}
f4 = tab2poly(syt4#0,syt4#0,syt4#0,syt4#0,syt4#1);
I4 = ideal(f4);
betti mingens I4 -- these are the first 5 primary invariants - but let's do a randomized jacobian to make sure we think they're algebraically independent.
time I4 = I4 + cycle I4 + cycle2 I4 + cycle3 I4 + cycle4 I4;
J4 = jacobian I4;
phi = rndPt();
rank phi J4 -- random rank 5 gives us some confidence that it is rank 5. -- I think: lower bounds are certain, upper bounds are probabilistic
--time ind4 = det J4^{0..4}; -- if this is not the zero polynomial, then we'd have confirmation that they're actually independent. However it is likely to take foreover. On the other hand, the zero polynomnial doesn't have non-zero values. So we could just evaluate it at random points and if it's not zero, then the rank is 5. That's what the randomized test does. 
--hilbertFunction(4,I4)
primaryInvariants = flatten entries gens I4;
elapsedTime inv6 = tab2poly(syt6#0, syt6#1, syt6#2, syt6#3, syt6#4); -- 1.09s (tensor product ring) vs .934s (flat ring)
I6 = ideal inv6;
I = I4+I6;
J46 = jacobian(I); -- to show that these polynomials are algebraically independent, you could have a general point in their intersection, then compute their Jacobian. 
rp = rndPt();
rank rp J46
primaryInvariants = primaryInvariants|{inv6};


elapsedTime inv8 = tab2poly(syt8#0,syt8#0,syt8#2,syt8#9,syt8#12); -- 20s
elapsedTime inv8tmp = tab2polySeq(syt8#0,syt8#0,syt8#2,syt8#9,syt8#12); -- 22s -- similar
inv8-inv8tmp
elapsedTime inv8b = tab2poly(syt8#0,syt8#2,syt8#3,syt8#12,syt8#13); -- 30s
elapsedTime inv8btmp = tab2polySeq(syt8#0,syt8#2,syt8#3,syt8#12,syt8#13); -- 28.9s
elapsedTime inv8c = tab2poly(syt8#0,syt8#0,syt8#12,syt8#12,syt8#12); -- 24s -- not independent of deg 4
elapsedTime inv8ctmp = tab2polySeq(syt8#0,syt8#0,syt8#12,syt8#12,syt8#12); -- 22.5s -- not independent of deg 4
elapsedTime inv8d = tab2poly(syt8#0,syt8#3,syt8#6,syt8#9,syt8#13); -- 39s
I8 = ideal{ inv8, inv8c, inv8b,inv8d};  -- leave out inv8c from the non-minimal generators. 
primaryInvariants = primaryInvariants|{inv8, inv8c, inv8d};
rank rp jacobian ideal primaryInvariants
PI8 = {inv8b , cycle inv8b , cycle2 inv8b , cycle3 inv8b , cycle4 inv8b};
rp = rndPt();
rank rp jacobian( ideal(PI8) + ideal primaryInvariants)
PI = ideal(PI8) + ideal primaryInvariants;
time betti mingens PI
primaryInvariants = primaryInvariants |PI8;

I = I4+I6+I8;
betti mingens (I)
betti I
J468 = jacobian(I); -- to show that these polynomials are algebraically independent, you could have a general point in their intersection, then compute their Jacobian. 
rp = rndPt();
rank rp J468
--
--J = jacobian(I4 + I6 + ideal{inv8, cycle inv8, cycle2 inv8, cycle3 inv8, cycle4 inv8}); 
--rp = rndPt();
--rank rp J

--time I8 = ideal(inv8, cycle inv8, cycle2 inv8, cycle3 inv8, cycle4 inv8); 
--time I8s = swap I8; betti I8
--time I8 = (I8+I8s); betti I8
time I8s = swap I8;
print betti (I8+I8s)
time I8 = I8 + I8s + cycle I8 + cycle2 I8 + cycle3 I8 + cycle4 I8  + cycle I8s + cycle2 I8s + cycle3 I8s + cycle4 I8s; 
time I8 = ideal mingens I8; betti I8 -- 25 polynomials when you only take inv8,inv8b,inv8d. 30 if you include inv8c
time I8s = swap I8; betti I8s
time I8 = ideal mingens (I8+I8s); betti I8 -- 36 poly's when you start from inv8, inv8b, inv8c, inv8d

I = I4+I6+I8;
betti I
tex betti mingens (I)
betti I
J468 = jacobian(I); -- to show that these polynomials are algebraically independent, you could have a general point in their intersection, then compute their Jacobian. 
rp = rndPt();
rank rp J468
--apply(9, i-> hilbertFunction(i, I))
netList\ transpose\ syt10
elapsedTime inv10b = tab2poly(syt10#0,syt10#17,syt10#24,syt10#37,syt10#40); -- approx 448s
elapsedTime inv10b = tab2polySeq(syt10#0,syt10#17,syt10#24,syt10#37,syt10#40); -- approx 448s
time rank rp jacobian( ideal(inv10b) + ideal primaryInvariants)
PI = ideal(inv10b) + ideal primaryInvariants;
betti  PI

I10 = ideal inv10b;
time I = ideal mingens (I+I10);
betti I
elapsedTime inv10 = tab2poly(syt10#0,syt10#0,syt10#4,syt10#34,syt10#40); -- approx 724s
I10 = I10 +  ideal inv10b;
I = I+I10;
time I = ideal mingens I;
betti I

-- try to get the rest of the generators:
I10 = I10 + I6*I4;
I10s = swap I10;
I10 = I10 + I10s;
time I10 = I10 + cycle I10 + cycle2 I10 + cycle3 I10 + cycle4 I10;
time I10 = ideal mingens I10; 
time betti I10
I10s = swap I10;
I10 = I10 + I10s;
--time hilbertFunction(10,I10);
time I10 = ideal mingens I10; 
time betti I10  

I = I+I10;
betti I
I = ideal mingens (I)
tex oo
--apply(11, i-> hilbertFunction(i, I))

m12to2row \last fill12


g12 = {{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}, {1, 3, 2, 4, 5, 6, 7, 8, 9, 11, 10, 12}, {1, 3, 2, 4, 5, 9, 6, 10, 7, 11, 8, 12}, {1, 5, 2, 7, 3, 8, 4, 9, 6, 10, 11, 12}, {1, 5, 2, 7, 3, 9, 4, 10, 6, 11, 8, 12}}
for gg in g12 list for ss to #syt12 -1 list if gg==(syt12#ss) then ss else continue 
g = "g"
elapsedTime inv12 = tab2poly(syt12#0, syt12#43, syt12#54, syt12#121, syt12#124); -- approx ?
-- used 0.003743s (cpu); 0.003737s (thread); 0s (gc)
 -- used 0.019548s (cpu); 0.0195507s (thread); 0s (gc)
 -- used 0.095837s (cpu); 0.0958415s (thread); 0s (gc)
 -- used 0.348665s (cpu); 0.348668s (thread); 0s (gc)
 -- used 6.92518s (cpu); 3.12194s (thread); 0s (gc)
 -- used 40.3507s (cpu); 23.272s (thread); 0s (gc)
 -- used 504.258s (cpu); 195.734s (thread); 0s (gc)
 -- used 546.026s (cpu); 383.372s (thread); 0s (gc)
 -- used 5026.21s (cpu); 2568.05s (thread); 0s (gc)
 -- used 9281.81s (cpu); 5324.85s (thread); 0s (gc)

time rank rp jacobian( ideal(inv12) + PI)
PI = PI + ideal(inv12) 
betti  PI

time rank rp jacobian( ideal(inv12) + ideal primaryInvariants)
PI = ideal(inv10b) + ideal primaryInvariants;
betti  PI
=======
---  tableaux found previously... 
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
goodTabList10 = {{0,0,4,34,40},{0,0,34,4,40},{0,0,34,40,4},{0,4,0,34,40},{0,4,34,0,40},{0,4,34,40,0},{0,34,0,4,40},{0,34,0,40,4},{0,34,40,0,4},{0,17,24,37,40}};
m12to2row = L -> {{L#0,L#1}, {L#2,L#3}, {L#4,L#5}, {L#6,L#7}, {L#8,L#9}, {L#10,L#11}};
psyt12 = {{1,2,3,4,5,6,7,8,9,10,11,12},{1,2,3,4,5,6,7,8,9,11,10,12},{1,2,3,4,5,6,7,9,8,10,11,12},{1,2,3,4,5,6,7,9,8,11,10,12},{1,2,3,4,5,6,7,10,8,11,9,12},{1,2,3,4,5,7,6,8,9,10,11,12},{1,2,3,4,5,7,6,8,9,11,10,12},{1,2,3,4,5,7,6,9,8,10,11,12},{1,2,3,4,5,7,6,9,8,11,10,12},{1,2,3,4,5,7,6,10,8,11,9,12},{1,2,3,4,5,8,6,9,7,10,11,12},{1,2,3,4,5,8,6,9,7,11,10,12},{1,2,3,4,5,8,6,10,7,11,9,12},{1,2,3,4,5,9,6,10,7,11,8,12},{1,2,3,5,4,6,7,8,9,10,11,12},{1,2,3,5,4,6,7,8,9,11,10,12},{1,2,3,5,4,6,7,9,8,10,11,12},{1,2,3,5,4,6,7,9,8,11,10,12},{1,2,3,5,4,6,7,10,8,11,9,12},{1,2,3,5,4,7,6,8,9,10,11,12},{1,2,3,5,4,7,6,8,9,11,10,12},{1,2,3,5,4,7,6,9,8,10,11,12},{1,2,3,5,4,7,6,9,8,11,10,12},{1,2,3,5,4,7,6,10,8,11,9,12},{1,2,3,5,4,8,6,9,7,10,11,12},{1,2,3,5,4,8,6,9,7,11,10,12},{1,2,3,5,4,8,6,10,7,11,9,12},{1,2,3,5,4,9,6,10,7,11,8,12},{1,2,3,6,4,7,5,8,9,10,11,12},{1,2,3,6,4,7,5,8,9,11,10,12},{1,2,3,6,4,7,5,9,8,10,11,12},{1,2,3,6,4,7,5,9,8,11,10,12},{1,2,3,6,4,7,5,10,8,11,9,12},{1,2,3,6,4,8,5,9,7,10,11,12},{1,2,3,6,4,8,5,9,7,11,10,12},{1,2,3,6,4,8,5,10,7,11,9,12},{1,2,3,6,4,9,5,10,7,11,8,12},{1,2,3,7,4,8,5,9,6,10,11,12},{1,2,3,7,4,8,5,9,6,11,10,12},{1,2,3,7,4,8,5,10,6,11,9,12},{1,2,3,7,4,9,5,10,6,11,8,12},{1,2,3,8,4,9,5,10,6,11,7,12},{1,3,2,4,5,6,7,8,9,10,11,12},{1,3,2,4,5,6,7,8,9,11,10,12},{1,3,2,4,5,6,7,9,8,10,11,12},{1,3,2,4,5,6,7,9,8,11,10,12},{1,3,2,4,5,6,7,10,8,11,9,12},{1,3,2,4,5,7,6,8,9,10,11,12},{1,3,2,4,5,7,6,8,9,11,10,12},{1,3,2,4,5,7,6,9,8,10,11,12},{1,3,2,4,5,7,6,9,8,11,10,12},{1,3,2,4,5,7,6,10,8,11,9,12},{1,3,2,4,5,8,6,9,7,10,11,12},{1,3,2,4,5,8,6,9,7,11,10,12},{1,3,2,4,5,8,6,10,7,11,9,12},{1,3,2,4,5,9,6,10,7,11,8,12},{1,3,2,5,4,6,7,8,9,10,11,12},{1,3,2,5,4,6,7,8,9,11,10,12},{1,3,2,5,4,6,7,9,8,10,11,12},{1,3,2,5,4,6,7,9,8,11,10,12},{1,3,2,5,4,6,7,10,8,11,9,12},{1,3,2,5,4,7,6,8,9,10,11,12},{1,3,2,5,4,7,6,8,9,11,10,12},{1,3,2,5,4,7,6,9,8,10,11,12},{1,3,2,5,4,7,6,9,8,11,10,12},{1,3,2,5,4,7,6,10,8,11,9,12},{1,3,2,5,4,8,6,9,7,10,11,12},{1,3,2,5,4,8,6,9,7,11,10,12},{1,3,2,5,4,8,6,10,7,11,9,12},{1,3,2,5,4,9,6,10,7,11,8,12},{1,3,2,6,4,7,5,8,9,10,11,12},{1,3,2,6,4,7,5,8,9,11,10,12},{1,3,2,6,4,7,5,9,8,10,11,12},{1,3,2,6,4,7,5,9,8,11,10,12},{1,3,2,6,4,7,5,10,8,11,9,12},{1,3,2,6,4,8,5,9,7,10,11,12},{1,3,2,6,4,8,5,9,7,11,10,12},{1,3,2,6,4,8,5,10,7,11,9,12},{1,3,2,6,4,9,5,10,7,11,8,12},{1,3,2,7,4,8,5,9,6,10,11,12},{1,3,2,7,4,8,5,9,6,11,10,12},{1,3,2,7,4,8,5,10,6,11,9,12},{1,3,2,7,4,9,5,10,6,11,8,12},{1,3,2,8,4,9,5,10,6,11,7,12},{1,4,2,5,3,6,7,8,9,10,11,12},{1,4,2,5,3,6,7,8,9,11,10,12},{1,4,2,5,3,6,7,9,8,10,11,12},{1,4,2,5,3,6,7,9,8,11,10,12},{1,4,2,5,3,6,7,10,8,11,9,12},{1,4,2,5,3,7,6,8,9,10,11,12},{1,4,2,5,3,7,6,8,9,11,10,12},{1,4,2,5,3,7,6,9,8,10,11,12},{1,4,2,5,3,7,6,9,8,11,10,12},{1,4,2,5,3,7,6,10,8,11,9,12},{1,4,2,5,3,8,6,9,7,10,11,12},{1,4,2,5,3,8,6,9,7,11,10,12},{1,4,2,5,3,8,6,10,7,11,9,12},{1,4,2,5,3,9,6,10,7,11,8,12},{1,4,2,6,3,7,5,8,9,10,11,12},{1,4,2,6,3,7,5,8,9,11,10,12},{1,4,2,6,3,7,5,9,8,10,11,12},{1,4,2,6,3,7,5,9,8,11,10,12},{1,4,2,6,3,7,5,10,8,11,9,12},{1,4,2,6,3,8,5,9,7,10,11,12},{1,4,2,6,3,8,5,9,7,11,10,12},{1,4,2,6,3,8,5,10,7,11,9,12},{1,4,2,6,3,9,5,10,7,11,8,12},{1,4,2,7,3,8,5,9,6,10,11,12},{1,4,2,7,3,8,5,9,6,11,10,12},{1,4,2,7,3,8,5,10,6,11,9,12},{1,4,2,7,3,9,5,10,6,11,8,12},{1,4,2,8,3,9,5,10,6,11,7,12},{1,5,2,6,3,7,4,8,9,10,11,12},{1,5,2,6,3,7,4,8,9,11,10,12},{1,5,2,6,3,7,4,9,8,10,11,12},{1,5,2,6,3,7,4,9,8,11,10,12},{1,5,2,6,3,7,4,10,8,11,9,12},{1,5,2,6,3,8,4,9,7,10,11,12},{1,5,2,6,3,8,4,9,7,11,10,12},{1,5,2,6,3,8,4,10,7,11,9,12},{1,5,2,6,3,9,4,10,7,11,8,12},{1,5,2,7,3,8,4,9,6,10,11,12},{1,5,2,7,3,8,4,9,6,11,10,12},{1,5,2,7,3,8,4,10,6,11,9,12},{1,5,2,7,3,9,4,10,6,11,8,12},{1,5,2,8,3,9,4,10,6,11,7,12},{1,6,2,7,3,8,4,9,5,10,11,12},{1,6,2,7,3,8,4,9,5,11,10,12},{1,6,2,7,3,8,4,10,5,11,9,12},{1,6,2,7,3,9,4,10,5,11,8,12},{1,6,2,8,3,9,4,10,5,11,7,12},{1,7,2,8,3,9,4,10,5,11,6,12}};
syt12 = m12to2row\ psyt12;
goodTabList12={{0,6,35,113,117},{0,6,35,117,113},{0,6,113,35,117},{0,6,113,117,35},{0,6,117,35,113},{0,6,117,113,35},{0,35,6,113,117},{0,35,6,117,113},{0,35,113,6,117},{0,35,113,117,6},{0,35,117,6,113},{0,35,117,113,6},{0,113,6,35,117},{0,113,6,117,35},{0,113,35,6,117},{0,113,35,117,6},{0,113,117,6,35},{0,113,117,35,6},{0,117,6,35,113},{0,117,6,113,35},{0,117,35,6,113},{0,117,35,113,6},{0,117,113,6,35},{0,117,113,35,6},{6,0,35,113,117},{6,0,35,117,113},{6,0,113,35,117},{6,0,113,117,35},{6,0,117,35,113},{6,0,117,113,35},{6,35,0,113,117},{6,35,0,117,113},{6,35,113,0,117},{6,35,113,117,0},{6,35,117,0,113},{6,35,117,113,0},{6,113,0,35,117},{6,113,0,117,35},{6,113,35,0,117},{6,113,35,117,0},{6,113,117,0,35},{6,113,117,35,0},{6,117,0,35,113},{6,117,0,113,35},{6,117,35,0,113},{6,117,35,113,0},{6,117,113,0,35},{6,117,113,35,0},{35,0,6,113,117},{35,0,6,117,113},{35,0,113,6,117},{35,0,113,117,6},{35,0,117,6,113},{35,0,117,113,6},{35,6,0,113,117},{35,6,0,117,113},{35,6,113,0,117},{35,6,113,117,0},{35,6,117,0,113},{35,6,117,113,0},{35,113,0,6,117},{35,113,0,117,6},{35,113,6,0,117},{35,113,6,117,0},{35,113,117,0,6},{35,113,117,6,0},{35,117,0,6,113},{35,117,0,113,6},{35,117,6,0,113},{35,117,6,113,0},{35,117,113,0,6},{113,0,6,35,117},{113,0,6,117,35},{113,0,35,6,117},{113,0,35,117,6},{113,0,117,6,35},{113,0,117,35,6},{113,6,0,35,117},{113,6,0,117,35},{113,6,35,0,117},{113,6,35,117,0},{113,6,117,0,35},{113,6,117,35,0},{113,35,0,6,117},{113,35,0,117,6},{113,35,6,0,117},{113,35,6,117,0},{113,35,117,0,6},{113,35,117,6,0},{113,117,0,6,35},{113,117,0,35,6},{113,117,6,0,35},{113,117,6,35,0},{113,117,35,0,6},{117,0,6,35,113},{117,0,6,113,35},{117,0,35,6,113},{117,0,35,113,6},{117,0,113,6,35},{117,0,113,35,6},{117,6,0,35,113},{117,6,0,113,35},{117,6,35,0,113},{117,6,35,113,0},{117,6,113,0,35},{117,35,0,6,113},{117,35,0,113,6},{117,35,113,0,6},{117,113,0,6,35},{117,113,0,35,6},{117,113,6,0,35},{0,5,32,36,91},{0,5,32,91,36},{0,5,36,32,91},{0,5,36,91,32},{0,5,91,32,36},{0,5,91,36,32},{0,32,5,36,91},{0,32,5,91,36},{0,32,36,5,91},{0,32,36,91,5},{0,32,91,5,36},{0,32,91,36,5},{0,36,5,32,91},{0,36,5,91,32},{0,36,32,5,91},{0,36,32,91,5},{0,36,91,5,32},{0,36,91,32,5},{0,91,5,32,36},{0,91,5,36,32},{0,91,32,5,36},{5,0,32,36,91},{5,0,32,91,36},{5,0,36,32,91},{5,0,36,91,32},{5,0,91,32,36},{5,32,0,36,91},{5,32,0,91,36},{5,32,36,0,91},{5,32,36,91,0},{5,32,91,0,36},{5,36,0,32,91},{5,36,0,91,32},{5,36,91,0,32},{5,91,0,32,36},{32,0,5,36,91},{32,0,5,91,36},{32,0,36,5,91},{32,0,36,91,5},{32,0,91,5,36},{32,5,0,36,91},{32,5,36,0,91},{32,36,0,5,91},{32,36,0,91,5},{32,36,91,0,5},{32,91,0,5,36},{36,0,5,32,91},{36,0,5,91,32},{36,0,32,5,91},{36,0,32,91,5},{36,0,91,5,32},{36,91,0,5,32},{91,0,5,32,36},{0,6,64,65,129},{0,6,64,129,65},{0,6,65,64,129},{0,6,65,129,64},{0,6,129,64,65},{0,64,6,65,129},{0,64,6,129,65},{0,64,65,6,129},{0,64,65,129,6},{0,64,129,6,65},{0,64,129,65,6},{0,65,6,64,129},{0,65,6,129,64},{0,65,129,6,64},{0,129,6,64,65},{0,129,64,6,65},{0,129,64,65,6},{6,0,64,65,129},{6,0,64,129,65},{6,0,65,64,129},{6,0,65,129,64},{6,0,129,64,65},{6,64,0,65,129},{6,64,0,129,65},{6,64,65,0,129},{6,64,65,129,0},{6,64,129,0,65},{6,65,0,64,129},{6,65,0,129,64},{6,65,129,0,64},{6,129,0,64,65},{64,0,6,65,129},{64,0,6,129,65},{64,0,65,6,129},{64,0,65,129,6},{64,0,129,6,65},{64,6,0,65,129},{64,6,65,0,129},{64,65,0,6,129},{64,129,0,6,65},{65,0,6,64,129},{129,0,6,64,65},{0,14,87,102,110},{0,14,87,110,102},{0,14,102,87,110},{0,14,102,110,87},{0,14,110,87,102},{0,87,14,102,110},{0,87,102,14,110},{0,87,102,110,14},{0,102,14,87,110},{0,110,14,87,102},{14,0,87,102,110},{87,0,14,102,110},{102,0,14,87,110},{0,41,78,115,116},{0,43,55,121,124},{0,17,34,72,119},{0,2,36,89,110},{0,3,42,48,116},{0,3,13,15,110},{0,1,12,44,46},{0,5,39,40,42},{0,27,58,80,93}};
minList12 = sort unique (sort\goodTabList12);
fill12 = apply(goodTabList12, L-> apply(L, l ->syt12#l));
end;
-->>>>>>> tmp
