---Young Symmetrizers for 4 factors
-- !work in progress!
--restart 
--KK = ZZ/nextPrime(100); --KK = QQ;
--restart
KK = QQ;n=4;
dims = {2,2,2,2}
X = KK[x_(0,0,0,0)..x_(dims_0-1,dims_1-1,dims_2-1,dims_3-1)];
perm = P-> map(X,X, apply(toList( (0,0,0,0)..(1,1,1,1)), ll -> x_(toSequence apply(P, p-> ll#p))));
sgn = p-> det((id_(KK^(max(p)+1)))_p^(sort p))
rndPt = () -> map(KK, X, random(KK^1, KK^(numgens X))) 
rp = rndPt();
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
goodTabList10 = {{0,0,4,34,40},{0,0,34,4,40},{0,0,34,40,4},{0,4,0,34,40},{0,4,34,0,40},{0,4,34,40,0},{0,34,0,4,40},{0,34,0,40,4},{0,34,40,0,4},{0,17,24,37,40}};
m12to2row = L -> {{L#0,L#1}, {L#2,L#3}, {L#4,L#5}, {L#6,L#7}, {L#8,L#9}, {L#10,L#11}};
psyt12 = {{1,2,3,4,5,6,7,8,9,10,11,12},{1,2,3,4,5,6,7,8,9,11,10,12},{1,2,3,4,5,6,7,9,8,10,11,12},{1,2,3,4,5,6,7,9,8,11,10,12},{1,2,3,4,5,6,7,10,8,11,9,12},{1,2,3,4,5,7,6,8,9,10,11,12},{1,2,3,4,5,7,6,8,9,11,10,12},{1,2,3,4,5,7,6,9,8,10,11,12},{1,2,3,4,5,7,6,9,8,11,10,12},{1,2,3,4,5,7,6,10,8,11,9,12},{1,2,3,4,5,8,6,9,7,10,11,12},{1,2,3,4,5,8,6,9,7,11,10,12},{1,2,3,4,5,8,6,10,7,11,9,12},{1,2,3,4,5,9,6,10,7,11,8,12},{1,2,3,5,4,6,7,8,9,10,11,12},{1,2,3,5,4,6,7,8,9,11,10,12},{1,2,3,5,4,6,7,9,8,10,11,12},{1,2,3,5,4,6,7,9,8,11,10,12},{1,2,3,5,4,6,7,10,8,11,9,12},{1,2,3,5,4,7,6,8,9,10,11,12},{1,2,3,5,4,7,6,8,9,11,10,12},{1,2,3,5,4,7,6,9,8,10,11,12},{1,2,3,5,4,7,6,9,8,11,10,12},{1,2,3,5,4,7,6,10,8,11,9,12},{1,2,3,5,4,8,6,9,7,10,11,12},{1,2,3,5,4,8,6,9,7,11,10,12},{1,2,3,5,4,8,6,10,7,11,9,12},{1,2,3,5,4,9,6,10,7,11,8,12},{1,2,3,6,4,7,5,8,9,10,11,12},{1,2,3,6,4,7,5,8,9,11,10,12},{1,2,3,6,4,7,5,9,8,10,11,12},{1,2,3,6,4,7,5,9,8,11,10,12},{1,2,3,6,4,7,5,10,8,11,9,12},{1,2,3,6,4,8,5,9,7,10,11,12},{1,2,3,6,4,8,5,9,7,11,10,12},{1,2,3,6,4,8,5,10,7,11,9,12},{1,2,3,6,4,9,5,10,7,11,8,12},{1,2,3,7,4,8,5,9,6,10,11,12},{1,2,3,7,4,8,5,9,6,11,10,12},{1,2,3,7,4,8,5,10,6,11,9,12},{1,2,3,7,4,9,5,10,6,11,8,12},{1,2,3,8,4,9,5,10,6,11,7,12},{1,3,2,4,5,6,7,8,9,10,11,12},{1,3,2,4,5,6,7,8,9,11,10,12},{1,3,2,4,5,6,7,9,8,10,11,12},{1,3,2,4,5,6,7,9,8,11,10,12},{1,3,2,4,5,6,7,10,8,11,9,12},{1,3,2,4,5,7,6,8,9,10,11,12},{1,3,2,4,5,7,6,8,9,11,10,12},{1,3,2,4,5,7,6,9,8,10,11,12},{1,3,2,4,5,7,6,9,8,11,10,12},{1,3,2,4,5,7,6,10,8,11,9,12},{1,3,2,4,5,8,6,9,7,10,11,12},{1,3,2,4,5,8,6,9,7,11,10,12},{1,3,2,4,5,8,6,10,7,11,9,12},{1,3,2,4,5,9,6,10,7,11,8,12},{1,3,2,5,4,6,7,8,9,10,11,12},{1,3,2,5,4,6,7,8,9,11,10,12},{1,3,2,5,4,6,7,9,8,10,11,12},{1,3,2,5,4,6,7,9,8,11,10,12},{1,3,2,5,4,6,7,10,8,11,9,12},{1,3,2,5,4,7,6,8,9,10,11,12},{1,3,2,5,4,7,6,8,9,11,10,12},{1,3,2,5,4,7,6,9,8,10,11,12},{1,3,2,5,4,7,6,9,8,11,10,12},{1,3,2,5,4,7,6,10,8,11,9,12},{1,3,2,5,4,8,6,9,7,10,11,12},{1,3,2,5,4,8,6,9,7,11,10,12},{1,3,2,5,4,8,6,10,7,11,9,12},{1,3,2,5,4,9,6,10,7,11,8,12},{1,3,2,6,4,7,5,8,9,10,11,12},{1,3,2,6,4,7,5,8,9,11,10,12},{1,3,2,6,4,7,5,9,8,10,11,12},{1,3,2,6,4,7,5,9,8,11,10,12},{1,3,2,6,4,7,5,10,8,11,9,12},{1,3,2,6,4,8,5,9,7,10,11,12},{1,3,2,6,4,8,5,9,7,11,10,12},{1,3,2,6,4,8,5,10,7,11,9,12},{1,3,2,6,4,9,5,10,7,11,8,12},{1,3,2,7,4,8,5,9,6,10,11,12},{1,3,2,7,4,8,5,9,6,11,10,12},{1,3,2,7,4,8,5,10,6,11,9,12},{1,3,2,7,4,9,5,10,6,11,8,12},{1,3,2,8,4,9,5,10,6,11,7,12},{1,4,2,5,3,6,7,8,9,10,11,12},{1,4,2,5,3,6,7,8,9,11,10,12},{1,4,2,5,3,6,7,9,8,10,11,12},{1,4,2,5,3,6,7,9,8,11,10,12},{1,4,2,5,3,6,7,10,8,11,9,12},{1,4,2,5,3,7,6,8,9,10,11,12},{1,4,2,5,3,7,6,8,9,11,10,12},{1,4,2,5,3,7,6,9,8,10,11,12},{1,4,2,5,3,7,6,9,8,11,10,12},{1,4,2,5,3,7,6,10,8,11,9,12},{1,4,2,5,3,8,6,9,7,10,11,12},{1,4,2,5,3,8,6,9,7,11,10,12},{1,4,2,5,3,8,6,10,7,11,9,12},{1,4,2,5,3,9,6,10,7,11,8,12},{1,4,2,6,3,7,5,8,9,10,11,12},{1,4,2,6,3,7,5,8,9,11,10,12},{1,4,2,6,3,7,5,9,8,10,11,12},{1,4,2,6,3,7,5,9,8,11,10,12},{1,4,2,6,3,7,5,10,8,11,9,12},{1,4,2,6,3,8,5,9,7,10,11,12},{1,4,2,6,3,8,5,9,7,11,10,12},{1,4,2,6,3,8,5,10,7,11,9,12},{1,4,2,6,3,9,5,10,7,11,8,12},{1,4,2,7,3,8,5,9,6,10,11,12},{1,4,2,7,3,8,5,9,6,11,10,12},{1,4,2,7,3,8,5,10,6,11,9,12},{1,4,2,7,3,9,5,10,6,11,8,12},{1,4,2,8,3,9,5,10,6,11,7,12},{1,5,2,6,3,7,4,8,9,10,11,12},{1,5,2,6,3,7,4,8,9,11,10,12},{1,5,2,6,3,7,4,9,8,10,11,12},{1,5,2,6,3,7,4,9,8,11,10,12},{1,5,2,6,3,7,4,10,8,11,9,12},{1,5,2,6,3,8,4,9,7,10,11,12},{1,5,2,6,3,8,4,9,7,11,10,12},{1,5,2,6,3,8,4,10,7,11,9,12},{1,5,2,6,3,9,4,10,7,11,8,12},{1,5,2,7,3,8,4,9,6,10,11,12},{1,5,2,7,3,8,4,9,6,11,10,12},{1,5,2,7,3,8,4,10,6,11,9,12},{1,5,2,7,3,9,4,10,6,11,8,12},{1,5,2,8,3,9,4,10,6,11,7,12},{1,6,2,7,3,8,4,9,5,10,11,12},{1,6,2,7,3,8,4,9,5,11,10,12},{1,6,2,7,3,8,4,10,5,11,9,12},{1,6,2,7,3,9,4,10,5,11,8,12},{1,6,2,8,3,9,4,10,5,11,7,12},{1,7,2,8,3,9,4,10,5,11,6,12}};
syt12 = m12to2row\ psyt12;


list2mats = (L,a) -> matrix apply(#L, i-> apply( L, j-> a_(i,j-1)));
--list2mats = (L,a) -> matrix apply(2, i-> apply( L, j-> a_(i,j-1)));

tab2poly = (L)->(
    ta:=L#0; tb:=L#1;tc:=L#2;td:=L#3;
    dg := length flatten ta;    
    R := X[a_(0,0)..a_(dims_0-1,dg-1),b_(0,0)..b_(dims_1-1,dg-1),c_(0,0)..c_(dims_2-1,dg-1),d_(0,0)..d_(dims_3-1,dg-1)];
    A := apply(ta, L -> list2mats(L, a));
    B := apply(tb, L -> list2mats(L, b));
    C := apply(tc, L -> list2mats(L, c));
    D := apply(td, L -> list2mats(L, d));
    usedA := {};usedB := {};usedC := {};usedD := {};
    F :=1_R; T :=1_R;
    -- check that these loops end at the right place... dg/2-1?
    for z to dg-1 do ( T=product(
        (for i to sub(dg/2-1,ZZ) list if isMember(z+1,ta#i) and not isMember(i,usedA) then (usedA = usedA|{i}; det A_i) else continue)|
        (for i to sub(dg/2-1,ZZ) list if isMember(z+1,tb#i) and not isMember(i,usedB) then (usedB = usedB|{i}; det B_i) else continue)|
        (for i to sub(dg/2-1,ZZ) list if isMember(z+1,tc#i) and not isMember(i,usedC) then (usedC = usedC|{i}; det C_i) else continue)|
        (for i to sub(dg/2-1,ZZ) list if isMember(z+1,td#i) and not isMember(i,usedD) then (usedD = usedD|{i}; det D_i) else continue));
        F = F*T;
         F= sum(dims_0, l-> sum(dims_1,k-> sum(dims_2, j-> sum(dims_3,i->x_(i,j,k,l)*diff(a_(i,z)*b_(j,z)*c_(k,z)*d_(l,z),F)))));
    ); -- after many checks (parallel sum, making auxillary arrays, using diff over a tensor, etc), this seems to be the most efficient way to execute the sum.
    sub(F,X)
);

tab2polyEval = (L, rules)->(
    ta:=L#0; tb:=L#1;tc:=L#2;td:=L#3;
    dg := length flatten ta;
    R := X[a_(0,0)..a_(dims_0-1,dg-1),b_(0,0)..b_(dims_1-1,dg-1),c_(0,0)..c_(dims_2-1,dg-1),d_(0,0)..d_(dims_3-1,dg-1)];
    A := apply(ta, L -> list2mats(L, a));
    B := apply(tb, L -> list2mats(L, b));
    C := apply(tc, L -> list2mats(L, c));
    D := apply(td, L -> list2mats(L, d));
    usedA := {};usedB := {};usedC := {};usedD := {};
    F :=1_R; T :=1_R;
    for z to dg-1 do ( T=product(
        (for i to sub(dg/2-1,ZZ) list if isMember(z+1,ta#i) and not isMember(i,usedA) then (usedA = usedA|{i}; det A_i) else continue)|
        (for i to sub(dg/2-1,ZZ) list if isMember(z+1,tb#i) and not isMember(i,usedB) then (usedB = usedB|{i}; det B_i) else continue)|
        (for i to sub(dg/2-1,ZZ) list if isMember(z+1,tc#i) and not isMember(i,usedC) then (usedC = usedC|{i}; det C_i) else continue)|
        (for i to sub(dg/2-1,ZZ) list if isMember(z+1,td#i) and not isMember(i,usedD) then (usedD = usedD|{i}; det D_i) else continue));
        F = F*T;
        F= sum(2,m-> sum(2, l-> sum(2,k-> sum(2, j-> sum(2,i->sub(x_(i,j,k,l,m),rules)*diff(a_(i,z)*b_(j,z)*c_(k,z)*d_(l,z)*e_(m,z),F))))));
    ); -- after many checks (parallel sum, making auxillary arrays, using diff over a tensor, etc), this seems to be the most efficient way to execute the sum.
    sub(F,X)
);

-- this doesn't seem to find much of a speed up, likely because the first tableau is already in the standard order.
-- tab2polyList = (ta,tb,tc,td,te,L)->(
--     dg := length flatten ta;
--     R := X[a_(0,0)..a_(1,dg-1),b_(0,0)..b_(1,dg-1),c_(0,0)..c_(1,dg-1),d_(0,0)..d_(1,dg-1),e_(0,0)..e_(1,dg-1)];
--     A := apply(ta, L -> list2mats(L, a));
--     B := apply(tb, L -> list2mats(L, b));
--     C := apply(tc, L -> list2mats(L, c));
--     D := apply(td, L -> list2mats(L, d));
--     E := apply(te, L -> list2mats(L, e));
--     usedA := {};usedB := {};usedC := {};usedD := {};usedE := {};
--     F :=1_R; T :=1_R;
--     for z in L do ( T=product(
--         (for i to sub(dg/2-1,ZZ) list if isMember(z+1,ta#i) and not isMember(i,usedA) then (usedA = usedA|{i}; det A_i) else continue)|
--         (for i to sub(dg/2-1,ZZ) list if isMember(z+1,tb#i) and not isMember(i,usedB) then (usedB = usedB|{i}; det B_i) else continue)|
--         (for i to sub(dg/2-1,ZZ) list if isMember(z+1,tc#i) and not isMember(i,usedC) then (usedC = usedC|{i}; det C_i) else continue)|
--         (for i to sub(dg/2-1,ZZ) list if isMember(z+1,td#i) and not isMember(i,usedD) then (usedD = usedD|{i}; det D_i) else continue)|
--         (for i to sub(dg/2-1,ZZ) list if isMember(z+1,te#i) and not isMember(i,usedE) then (usedE = usedE|{i}; det E_i) else continue));
--         F = F*T;
--         F= sum(2,m-> sum(2, l-> sum(2,k-> sum(2, j-> sum(2,i->x_(i,j,k,l,m)*diff(a_(i,z)*b_(j,z)*c_(k,z)*d_(l,z)*e_(m,z),F))))));
--     ); -- should check to see if you make auxillary variables for the derivative wrt some of the abc's and split up the computation if that speeds things up more.
--     sub(F,X)
-- );
-- reuse the variables, just set the pair to the first variable.

--- points 
gx = gens X
---  tableaux
f2 = tab2poly(apply(4, i-> syt2#0));
f2
P = {0,0,0,1}
tgt = uniquePermutations P;
ct = 0; myPerms = flatten for i to #tgt-1 list if ct == i then for xx in permutations(n) list if ct == i and P_xx == tgt#i then (ct = ct+1; xx ) else continue
#myPerms
f4a = tab2poly({syt4#0,syt4#0,syt4#0,syt4#1});
factor f4a
betti mingens ideal {f2,f4a}
time F4big = unique apply(myPerms, p-> (perm p) f4a);
#F4big
betti mingens ideal F4big
f4b = tab2poly({syt4#0,syt4#0,syt4#1,syt4#1});
f4b
P = {0,0,1,1}
tgt = uniquePermutations P;
ct = 0; myPerms = flatten for i to #tgt-1 list if ct == i then for xx in permutations(n) list if ct == i and P_xx == tgt#i then (ct = ct+1; xx ) else continue
#myPerms

time F4big = unique apply(myPerms, p-> (perm p) f4b);
#F4big
betti mingens ideal F4big
netList (factor\F4big)

betti mingens ideal({f2} |F4big)

#syt6
f6 = tab2poly({syt6#0,syt6#1,syt6#2,syt6#3});
f6 = sum(permutations 4, p-> (perm p) f6);
F6 = unique flatten apply(permutations 5, p -> tab2poly(apply(n, i-> syt6#(p#i))));
--F6 = {sum(permutations n, p -> (perm p) f6)}
--F6 = unique apply(permutations n, p -> (perm p) f6);
#F6
betti mingens ideal F6_{0..3}
betti mingens ideal F6
F6 = F6_{0..3};
input"/Users/lao0004/Documents/GitHub/YoungSymmetrizers/YoungSpecht.m2"

myTab = {{0,1,2,3}}; V6 = clambda_myTab \ F6;
betti mingens ideal V6
myTab = {{0,1,2},{3}}; V631 = clambda_myTab \ F6;
betti mingens ideal V631
myTab = {{0,1},{2,3}}; V622 = clambda_myTab \ F6;
betti mingens ideal V622
myTab = {{0,2},{1,3}}; V622 = V622|(clambda_myTab \ F6);
betti mingens ideal V622

myTab = {{0,1},{2},{3}}; V6211 = clambda_myTab \ F6;
betti mingens ideal V6211
myTab = {{0,2},{1},{3}}; V6211 = clambda_myTab \ F6;
betti mingens ideal V6211
myTab = {{0,3},{2},{1}}; V6211 = clambda_myTab \ F6;
betti mingens ideal V6211


myTab = {{0},{1},{2},{3}}; V61111 = clambda_myTab \ F6;
betti mingens ideal V61111

betti mingens ideal F6
sum(permutations n, p -> sgn(p)*(perm p) f6)

--betti mingens ideal{f2, f4b, f4c,f6}


loadPackage"SpechtModule"
n=4;
charTable = characterTable  n
CT = matrix charTable#values
-- partition to conjugacy class rep:
-- takes a partition as a list and produces 
 -- permutation that is composed of disjoint cycles,
 --  in order, with the cycle in each case being the
 --  next cycle of the integers in that part of the list
cyc = L-> L_{1..#L-1,0}
partToClassRep = L ->(T := toList(0..(sum L)); S:={0}|apply(#L,i-> sum(L_{0..i}));
  flatten apply(#L, l ->cyc  T_{S#l..S#(l+1)-1} )
)


reps = apply(partitions n, p -> partToClassRep toList  p)

G = permutations(n);
con = h->  sort unique apply(G, g-> ((inversePermutation(g))_h)_g)
conjClasses = (unique (con\reps));
orders = apply(conjClasses, c-> #c) -- sizes of conjugacy classes.
ip = (av,bv) -> 1/n!*trace( av*diagonalMatrix(orders)*transpose bv)

computeChar = polList ->(
    pts:= apply(#polList, i-> rndPt());
    makeMg := myPerm -> sub(matrix apply(pts_{0..#polList-1},pt-> apply(polList, f-> (pt*myPerm) f)), baseRing X) ; 
    print("computing the representation matrices");
    Mgs := for rep in reps list time makeMg(perm rep);
    M := last Mgs; 
    apply(Mgs, MMg->trace(M^(-1)*MMg ))
)
decomposeCharacter = mychar ->    (
    dc := apply(#partitions n , i-> ip(matrix({mychar}), CT^{i}));
    print(netList {toList\(partitions n), dc});
    dc)
decomposeCharacter computeChar({f2})
decomposeCharacter computeChar(F4big)
decomposeCharacter computeChar({f2^2})
decomposeCharacter computeChar(F6_{0..3})

I = ideal ({f2}|F4big_{0,1,2}|F6);
betti mingens I
time r = res I;
betti r
end;
-->>>>>>> tmp
restart
input"/Users/lao0004/Documents/GitHub/YoungSymmetrizers/4factor.m2"

loadPackage"SchurRings"
n=4
S = schurRing(s,n,GroupActing => "Sn")
rep = 2*s_(4) + s_{2,2}
exteriorPower(2, s_{2,2})
symmetricPower(2, s_{2,2})
netList schurResolution(rep, {s_n}, DegreeLimit =>4)
betti r
symmetricPower(2,s_n)
rep = 2*s_(4) + s_{2,2}
netList schurResolution(rep, {s_n}, DegreeLimit =>4)
-- it's not clear how to include the degree for each representation
-- but you can put the degrees back in by computing the exterior powers by hand
--
r26 = res ideal {f2,f6};
betti r26

netList schurResolution(2*s_(4), {s_n}, DegreeLimit =>4)

netList schurResolution(2*s_(4), {s_n}, DegreeLimit =>4)

splice{1..5}
-- the invariant ring for 4 qubits looks like this in even degrees:
M = {s_4, s_4 +  s_{2,2}, 2* s_4  + s_{2,2}}
-- minimal generators
rep = {s_4, s_{2,2}, s_4 }
netList schurResolution(sum rep, M, DegreeLimit =>4)

rep = s_n + s_{n-1,1}
M = {s_n}|splice{n:rep}
