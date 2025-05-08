---Young Symmetrizers  -- with tensor product rings, and one partial evaluation trick
-- the polynomial is the 3x3x3 degree 9 Strassen Invariant
--restart -- input"5factor.m2"
KK = ZZ/101;
X = KK[x_(0,0,0,0,0)..x_(1,1,1,1,1)]
primaryInvariants = {};
rndPt = () -> map(KK, X, random(KK^1, KK^(numgens X))) 
rp = rndPt()
cycle = map(X,X, flatten flatten flatten flatten  apply(2, i1-> apply(2, i2->apply(2, i3->apply(2, i4->apply(2, i5-> x_(i2,i3,i4,i5,i1) ) ) ) )) );
cycle2 = map(X,X, flatten flatten flatten flatten  apply(2, i1-> apply(2, i2->apply(2, i3->apply(2, i4->apply(2, i5-> x_(i3,i4,i5,i1,i2) ) ) ) )) );
cycle3 = map(X,X, flatten flatten flatten flatten  apply(2, i1-> apply(2, i2->apply(2, i3->apply(2, i4->apply(2, i5-> x_(i4,i5,i1,i2,i3) ) ) ) )) );
cycle4 = map(X,X, flatten flatten flatten flatten  apply(2, i1-> apply(2, i2->apply(2, i3->apply(2, i4->apply(2, i5-> x_(i5,i1,i2,i3,i4) ) ) ) )) );
swap = map(X,X, flatten flatten flatten flatten  apply(2, i1-> apply(2, i2->apply(2, i3->apply(2, i4->apply(2, i5-> x_(i2,i1,i3,i4,i5) ) ) ) )) );
list2mats = (L,a) -> matrix apply(2, i-> apply( L, j-> a_(i,j-1)))
tab2poly = (ta,tb,tc,td,te)->(
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
    ); -- should check to see if you make auxillary variables for the derivative wrt some of the abc's and split up the computation if that speeds things up more.
    sub(F,X)
)
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
