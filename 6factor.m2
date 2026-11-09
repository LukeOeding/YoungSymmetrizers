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
if class(KK)===Symbol then KK=ZZ/101;
--KK = QQ;
X = KK[x_(0,0,0,0,0,0)..x_(1,1,1,1,1,1)];
primaryInvariants = {};
rndPt = () -> map(KK, X, random(KK^1, KK^(numgens X))) 
rp = rndPt();
perm = P-> map(X,X, apply(toList( (0,0,0,0,0,0)..(1,1,1,1,1,1)), ll -> x_(toSequence apply(P, p-> ll#p))));
list2mats = (L,a) -> matrix apply(2, i-> apply( L, j-> a_(i,j-1)));
tab2poly = (L)->(
    ta:=L#0; tb:=L#1;tc:=L#2;td:=L#3;te:=L#4;tf:=L#5;
    dg := length flatten ta;
    R := X[a_(0,0)..a_(1,dg-1),b_(0,0)..b_(1,dg-1),c_(0,0)..c_(1,dg-1),d_(0,0)..d_(1,dg-1),e_(0,0)..e_(1,dg-1),f_(0,0)..f_(1,dg-1), ];
    A := apply(ta, L -> list2mats(L, a)); B := apply(tb, L -> list2mats(L, b));
    C := apply(tc, L -> list2mats(L, c)); D := apply(td, L -> list2mats(L, d));
    E := apply(te, L -> list2mats(L, e)); F := apply(tf, L -> list2mats(L, f));
    usedA := {};usedB := {};usedC := {};usedD := {}; usedE := {}; usedF := {};
    Fp :=1_R; T :=1_R;
    for z to dg-1 do (
        T=product(
            (for i to sub(dg/2-1,ZZ) list if isMember(z + 1,ta#i) and not isMember(i,usedA) then (usedA = usedA|{i}; det A_i) else continue)|
            (for i to sub(dg/2-1,ZZ) list if isMember(z + 1,tb#i) and not isMember(i,usedB) then (usedB = usedB|{i}; det B_i) else continue)|
            (for i to sub(dg/2-1,ZZ) list if isMember(z + 1,tc#i) and not isMember(i,usedC) then (usedC = usedC|{i}; det C_i) else continue)|
            (for i to sub(dg/2-1,ZZ) list if isMember(z + 1,td#i) and not isMember(i,usedD) then (usedD = usedD|{i}; det D_i) else continue)|
            (for i to sub(dg/2-1,ZZ) list if isMember(z + 1,te#i) and not isMember(i,usedE) then (usedE = usedE|{i}; det E_i) else continue)|
            (for i to sub(dg/2-1,ZZ) list if isMember(z + 1,tf#i) and not isMember(i,usedF) then (usedF = usedF|{i}; det F_i) else continue)
        );
    Fp = Fp*T; -- you could insert a time before this...
 Fp= sum(2,n->sum(2,m-> sum(2, l-> sum(2,k-> sum(2, j-> sum(2,i->x_(i,j,k,l,m,n)*diff(a_(i,z)*b_(j,z)*c_(k,z)*d_(l,z)*e_(m,z)*f_(n,z),Fp)))))));
 ); -- after many checks (parallel sum, making auxillary arrays, using diff over a tensor, etc), this seems to be the most efficient way to execute the sum.
 sub(Fp,X)
);

-- this doesn't work faster than just differentiating...
-- tab2polyML = (L)->(
--     ta:=L#0; tb:=L#1;tc:=L#2;td:=L#3;te:=L#4;tf:=L#5;
--     dg := length flatten ta;
--     R := X[a_(0,0)..a_(1,dg-1),b_(0,0)..b_(1,dg-1),c_(0,0)..c_(1,dg-1),d_(0,0)..d_(1,dg-1),e_(0,0)..e_(1,dg-1),f_(0,0)..f_(1,dg-1), ];
--     A := apply(ta, L -> list2mats(L, a)); B := apply(tb, L -> list2mats(L, b));
--     C := apply(tc, L -> list2mats(L, c)); D := apply(td, L -> list2mats(L, d));
--     E := apply(te, L -> list2mats(L, e)); F := apply(tf, L -> list2mats(L, f));
--     usedA := {};usedB := {};usedC := {};usedD := {}; usedE := {}; usedF := {};
--     Fp :=1_R; T :=1_R;
--     for z to dg-1 do (
--         T=product(
--             (for i to sub(dg/2-1,ZZ) list if isMember(z + 1,ta#i) and not isMember(i,usedA) then (usedA = usedA|{i}; det A_i) else continue)|
--             (for i to sub(dg/2-1,ZZ) list if isMember(z + 1,tb#i) and not isMember(i,usedB) then (usedB = usedB|{i}; det B_i) else continue)|
--             (for i to sub(dg/2-1,ZZ) list if isMember(z + 1,tc#i) and not isMember(i,usedC) then (usedC = usedC|{i}; det C_i) else continue)|
--             (for i to sub(dg/2-1,ZZ) list if isMember(z + 1,td#i) and not isMember(i,usedD) then (usedD = usedD|{i}; det D_i) else continue)|
--             (for i to sub(dg/2-1,ZZ) list if isMember(z + 1,te#i) and not isMember(i,usedE) then (usedE = usedE|{i}; det E_i) else continue)|
--             (for i to sub(dg/2-1,ZZ) list if isMember(z + 1,tf#i) and not isMember(i,usedF) then (usedF = usedF|{i}; det F_i) else continue)
--         );
--     Fp = Fp*T; -- you could insert a time before this...
--  Fp= sum(2,n->sum(2,m-> sum(2, l-> sum(2,k-> sum(2, j-> sum(2,i->x_(i,j,k,l,m,n)*
--  (sub(Fp,{a_(i,z) =>0, b_(j,z) =>0, c_(k,z) =>0, d_(l,z) =>0, e_(m,z) =>0, f_(n,z)  =>0}) - 
--   sub(Fp,{a_(i,z) =>1, b_(j,z) =>1, c_(k,z) =>1, d_(l,z) =>1, e_(m,z) =>1, f_(n,z)  =>1}))
--  ))))));
--  ); -- after many checks (parallel sum, making auxillary arrays, using diff over a tensor, etc), this seems to be the most efficient way to execute the sum.
--  sub(Fp,X)
-- );

tab2polyEvalRules = (L,rules)->(
    ta:=L#0; tb:=L#1;tc:=L#2;td:=L#3;te:=L#4;tf:=L#5;
    dg := length flatten ta;
    R := X[a_(0,0)..a_(1,dg-1),b_(0,0)..b_(1,dg-1),c_(0,0)..c_(1,dg-1),d_(0,0)..d_(1,dg-1),e_(0,0)..e_(1,dg-1),f_(0,0)..f_(1,dg-1), ];
    A := apply(ta, L -> list2mats(L, a)); B := apply(tb, L -> list2mats(L, b));
    C := apply(tc, L -> list2mats(L, c)); D := apply(td, L -> list2mats(L, d));
    E := apply(te, L -> list2mats(L, e)); F := apply(tf, L -> list2mats(L, f));
    usedA := {}; usedB := {}; usedC := {}; usedD := {}; usedE := {}; usedF := {};
    Fp :=1_R; T :=1_R;
    for z to dg-1 do ( 
        T=product(
            (for i to sub(dg/2-1,ZZ) list if isMember(z + 1,ta#i) and not isMember(i,usedA) then (usedA = usedA|{i}; det A_i) else continue)|
            (for i to sub(dg/2-1,ZZ) list if isMember(z + 1,tb#i) and not isMember(i,usedB) then (usedB = usedB|{i}; det B_i) else continue)|
            (for i to sub(dg/2-1,ZZ) list if isMember(z + 1,tc#i) and not isMember(i,usedC) then (usedC = usedC|{i}; det C_i) else continue)|
            (for i to sub(dg/2-1,ZZ) list if isMember(z + 1,td#i) and not isMember(i,usedD) then (usedD = usedD|{i}; det D_i) else continue)|
            (for i to sub(dg/2-1,ZZ) list if isMember(z + 1,te#i) and not isMember(i,usedE) then (usedE = usedE|{i}; det E_i) else continue)|
            (for i to sub(dg/2-1,ZZ) list if isMember(z + 1,tf#i) and not isMember(i,usedF) then (usedF = usedF|{i}; det F_i) else continue)
        );
        Fp = Fp*T; -- you could insert a time before this...
        Fp= sum(2,n->sum(2,m-> sum(2, l-> sum(2,k-> sum(2, j-> sum(2,i->sub(x_(i,j,k,l,m,n),rules)*diff(a_(i,z)*b_(j,z)*c_(k,z)*d_(l,z)*e_(m,z)*f_(n,z),Fp)))))));
    ); -- after many checks (parallel sum, making auxillary arrays, using diff over a tensor, etc), this seems to be the most efficient way to execute the sum.
    sub(Fp,X)
);

tab2polyEval = (L,vec)->(
    a :=local a; b:=local b; c:= local c; d := local d; e := local e; f:= local f; x:= local x;
    X := KK[x_(0,0,0,0,0,0)..x_(1,1,1,1,1,1)];
    gx:= gens X;
    rules := apply(2^6,i-> gx_i => vec_i);
    ta:=L#0; tb:=L#1; tc:=L#2; td:=L#3; te:=L#4; tf:=L#5;
    dg := length flatten ta;
    R := X[a_(0,0)..a_(1,dg-1),b_(0,0)..b_(1,dg-1),c_(0,0)..c_(1,dg-1),d_(0,0)..d_(1,dg-1),e_(0,0)..e_(1,dg-1),f_(0,0)..f_(1,dg-1), ];
    A := apply(ta, L -> list2mats(L, a)); B := apply(tb, L -> list2mats(L, b));
    C := apply(tc, L -> list2mats(L, c)); D := apply(td, L -> list2mats(L, d));
    E := apply(te, L -> list2mats(L, e)); F := apply(tf, L -> list2mats(L, f));
    usedA := {}; usedB := {}; usedC := {}; usedD := {}; usedE := {}; usedF := {};
    Fp :=1_R; T :=1_R; -- took out "time" for the z-for loop
    for z to dg-1 do (
        T=product(
            (for i to sub(dg/2-1,ZZ) list if isMember(z + 1,ta#i) and not isMember(i,usedA) then (usedA = usedA|{i}; det A_i) else continue)|
            (for i to sub(dg/2-1,ZZ) list if isMember(z + 1,tb#i) and not isMember(i,usedB) then (usedB = usedB|{i}; det B_i) else continue)|
            (for i to sub(dg/2-1,ZZ) list if isMember(z + 1,tc#i) and not isMember(i,usedC) then (usedC = usedC|{i}; det C_i) else continue)|
            (for i to sub(dg/2-1,ZZ) list if isMember(z + 1,td#i) and not isMember(i,usedD) then (usedD = usedD|{i}; det D_i) else continue)|
            (for i to sub(dg/2-1,ZZ) list if isMember(z + 1,te#i) and not isMember(i,usedE) then (usedE = usedE|{i}; det E_i) else continue)|
            (for i to sub(dg/2-1,ZZ) list if isMember(z + 1,tf#i) and not isMember(i,usedF) then (usedF = usedF|{i}; det F_i) else continue));
        Fp = Fp*T; -- you could insert a time before this...
        Fp= sum(2,n->sum(2,m-> sum(2, l-> sum(2,k-> sum(2, j-> sum(2,i->sub(x_(i,j,k,l,m,n),rules)*diff(a_(i,z)*b_(j,z)*c_(k,z)*d_(l,z)*e_(m,z)*f_(n,z),Fp)))))));
    ); 
    sub(Fp,KK)
);


--- points 
gx = gens X
--pseudoCartan = apply(2^5,i-> gx_(2^6-1-i) => gx_i)
--randomLine = apply(2^6,i-> gx_(i) => random(KK)*gx_0 + random(KK)*gx_15 + random(KK)*gx_31)
randomPoint = apply(2^6,i-> gx_(i) => random(KK))
randomVec = apply(2^6,i-> random(KK))
--- tableaux
--- timing on MBP with maxAllowableThreads
--allowableThreads = maxAllowableThreads -- 2nd time we don't use this, compare times
--allowableThreads = 1--  maxAllowableThreads-1
--allowableThreads = maxAllowableThreads
--elapsedTime parallelApply(1..5, j-> (R:=QQ[Variables =>10]; betti mingens ideal apply(20, i-> random(5, R))))
time inv2 = tab2poly({syt2#0,syt2#0,syt2#0,syt2#0,syt2#0,syt2#0});
-- used 0.017699s (cpu); 0.0177007s (thread); 0s (gc)
-- used 0.020683s (cpu); 0.0206809s (thread); 0s (gc)
F2 = {inv2};
I2 = ideal inv2;
--inv40 = tab2poly({syt4#0,syt4#0,syt4#0,syt4#0,syt4#0,syt4#0});
--factor inv40 -- this is the square of the degree 2 invariant
--inv4 = tab2poly({syt4#0,syt4#0,syt4#0,syt4#0,syt4#0,syt4#1});
--factor inv4 -- this is the square of the degree 2 invariant
-- so... this one is new, not in the ideal of the degree 2 invariant
time inv42 = tab2poly({syt4#0,syt4#0,syt4#0,syt4#0,syt4#1,syt4#1});
-- used 0.354864s (cpu); 0.28645s (thread); 0s (gc)
-- used 0.381847s (cpu); 0.310352s (thread); 0s (gc)
time f43 =  tab2poly({syt4#0,syt4#0,syt4#0,syt4#1,syt4#1,syt4#1});
-- used 0.321175s (cpu); 0.29434s (thread); 0s (gc)
-- used 0.332638s (cpu); 0.301291s (thread); 0s (gc)

-- vecLibrary = for i to 16 list apply(2^6,i-> random(KK));
-- --time newCol = matrix apply(pointLibrary, p ->  {tab2polyEval({syt4#0,syt4#0,syt4#0,syt4#0,syt4#1, syt4#1}, p)});  
-- elapsedTime newCol = matrix apply(vecLibrary, p ->  {tab2polyEval({syt4#0,syt4#0,syt4#0,syt4#0,syt4#1, syt4#1}, p)});  
-- --elapsedTime newCol = matrix parallelApply(vecLibrary, p ->  {tab2polyEvalVec({syt4#0,syt4#0,syt4#0,syt4#0,syt4#1, syt4#1}, p)});  
-- -- used 4.03423s (cpu); 2.494s (thread); 0s (gc)
-- -- used 4.19182s (cpu); 2.59258s (thread); 0s (gc)

-- M4 = newCol; rk=1;
-- elapsedTime for P in drop(uniquePermutations {0,0,0,0,1,1},1) do(
--     time newCol = matrix apply(vecLibrary, p ->  {tab2polyEvalVec(apply(6, i-> syt4#(P#i)), p)});  
--     M4 = M4|newCol;
--     print(toString(P)|" rank M4 = "|toString(rank M4)));
--  -- used 58.0856s (cpu); 35.621s (thread); 0s (gc) over QQ, 36s over ZZ/101
--  -- used 59.9313s (cpu); 37.342s (thread); 0s (gc)
-- elapsedTime M4 = matrix apply(vecLibrary, p -> apply(uniquePermutations {0,0,0,0,1,1}, P->  tab2polyEvalVec(apply(6, k-> syt4#(P#k)), p)));  
-- -- 40.0761s elapsed on MBP over QQ, 39s over ZZ/101
-- rank M4
-- ker M4
-- --elapsedTime M4tmp = matrix apply(vecLibrary, p -> parallelApply(uniquePermutations {0,0,0,0,1,1}, P->  tab2polyEvalVec(apply(6, k-> syt4#(P#k)), p)));  
-- -- 100.532s elapsed on MBP over QQ
-- --elapsedTime M4tmp = transpose matrix parallelApply(uniquePermutations {0,0,0,0,1,1}, P-> apply(vecLibrary, p -> tab2polyEvalVec(apply(6, k-> syt4#(P#k)), p)));  
--  -- 97.1781s elapsed on MBP over QQ
-- -- this just stops and doesn't actually compute anything... maybe it can't allocate enough resources to do the computation?
-- --elapsedTime M4tmp = matrix parallelApply(vecLibrary, p -> parallelApply(uniquePermutations {0,0,0,0,1,1}, P->  tab2polyEval(apply(6, k-> syt4#(P#k)), p)));  
myPerms42 ={{0,1,2,3,4,5}, {0,1,2,4,3,5}, {0,1,4,3,2,5}, {0,4,2,3,1,5}, {0,1,2,4,5,3}, {0,1,4,2,5,3}, {0,4,2,1,5,3}, {0,1,4,5,2,3}, {0,4,1,5,2,3}, {0,4,5,1,2,3}, {4,1,2,3,0,5}, {4,1,2,0,5,3}, {4,1,0,5,2,3}, {4,0,5,1,2,3}, {4,5,0,1,2,3}};
time deg4List = for P in myPerms42 list ( (perm P) inv42);
-- used 0.069648s (cpu); 0.0696382s (thread); 0s (gc)
-- used 0.054087s (cpu); 0.0540874s (thread); 0s (gc)
--apply(deg4List, ff -> #factor ff)
betti ideal deg4List
--time betti mingens ideal deg4List
-- used 0.057871s (cpu); 0.0578629s (thread); 0s (gc)
 -- used 0.049638s (cpu); 0.0496364s (thread); 0s (gc)
-- here it's much faster to just make the polynomial and then permute the variables, then compute mingens, than to com
F4 = deg4List_{0..9};
I4 = ideal F4;

betti mingens I4 -- notice that the first 10 are maximally independent
-- this is a sequential way to do this:
--factor sum deg4List
--F4 ={};rk =0; time pmList = for P in myPerms42 list ( tmpPerm = perm P; tmp = tmpPerm inv42; if rank source mingens ideal (F4|{tmp})>rk then (F4 = F4|{tmp}; rk=rk + 1;P ) else continue);
 -- used 0.375883s (cpu); 0.314275s (thread); 0s (gc)
 -- used 0.435511s (cpu); 0.329417s (thread); 0s (gc)
--LL = {0,0,0,0,1,1};
--netList  for P in pmList list apply(P, i->LL_i)
--time F42 = for P in myPerms42 list (tmpPerm = perm P; tmpPerm inv42);
--betti mingens(I4 + ideal f43)
--betti mingens(I4 + ideal f43+I2) -- this shows that f43 is in the ideal generated by the degree 2 and the other degree 4 parts
--time betti mingens(I4+I2)
-- used 0.050143s (cpu); 0.0501329s (thread); 0s (gc)
-- used 0.057863s (cpu); 0.0578575s (thread); 0s (gc)
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
J24 = jacobian (I2 + I4);
phi = rndPt();
rank phi J24 

primaryInvariants = F2|F4;
#primaryInvariants
-- note that we have 10 + 1 invariants in degree 4, which matches the hilbert function
--elapsedTime inv6 = tab2poly({syt6#0, syt6#1, syt6#2, syt6#3, syt6#4,syt6#4}); -- 26s 
--#unique permutations{0,1,2,3,4,4}
--#terms inv6 -- 68192
--elapsedTime inv60 = tab2poly({syt6#0, syt6#0, syt6#0, syt6#0, syt6#0,syt6#0}); -- 3.6s 
-- time inv61 = tab2poly({syt6#0, syt6#0, syt6#0, syt6#0, syt6#0,syt6#4}); 
-- factor inv61
--#terms inv61 -- 5984, and it's the cube of the degree 2 invariant... 
-- note 000001 and similar produce the same polynomial.. 
-- time tmp = tab2poly({syt6#0, syt6#0, syt6#0, syt6#0, syt6#1,syt6#1}); -- 3.6s 
-- #terms tmp -- 24416
-- apply(factor tmp, xx-> degree value xx) -- this is a quadratic times a quartic
-- time tmp = tab2poly({syt6#0, syt6#0, syt6#0, syt6#0, syt6#1,syt6#2}); -- this is irreducible degree 6
-- #terms tmp
-- apply(factor tmp, xx-> degree value xx)

--betti mingens (I2 + ideal inv62 )
--netList syt6
time inv62 = tab2poly(apply({0,0,0,0,3,3}, p-> syt6#p)); 
-- used 17.5454s (cpu); 14.3367s (thread); 0s (gc)
-- used 17.6273s (cpu); 14.8133s (thread); 0s (gc)

-- this is not the most efficient way to do this:
--time F62 = apply(unique permutations {0,0,0,0,3,3}, p-> tab2poly({syt6#(p#0),syt6#(p#1),syt6#(p#2),syt6#(p#3),syt6#(p#4),syt6#(p#5)}));
--myPerms ={{0,1,2,3,4,5}, {0,1,2,4,3,5}, {0,1,4,3,2,5}, {0,4,2,3,1,5}, {4,1,2,3,0,5}, {0,1,2,4,5,3}, {0,1,4,2,5,3}, {0,4,2,1,5,3}, {4,1,2,0,5,3}, {0,1,4,5,2,3}, {0,4,1,5,2,3}, {4,1,0,5,2,3}, {0,4,5,1,2,3}, {4,0,5,1,2,3}, {4,5,0,1,2,3}}
time F62 = for P in myPerms42 list (tmpPerm = perm P; tmpPerm inv62);
-- used 2.06377s (cpu); 1.5739s (thread); 0s (gc)
-- used 2.19184s (cpu); 1.62901s (thread); 0s (gc)
I62 = ideal F62;
if checks ==true then time betti mingens I62 -- 15, so they're all independent
-- used 6.92331s (cpu); 4.7452s (thread); 0s (gc)
-- used 6.69836s (cpu); 4.58144s (thread); 0s (gc)

if checks ==true then time betti mingens (I2 + I4 + I62) -- and they're independent from the lower degree terms
-- used 20.0459s (cpu); 12.6877s (thread); 0s (gc)
 -- used 27.8768s (cpu); 15.1133s (thread); 0s (gc)
--elapsedTime inv62b = tab2poly({syt6#0, syt6#0, syt6#0, syt6#0, syt6#4,syt6#4}); -- this one is not independent from the other one
--time betti mingens (I62 + ideal inv62b)
--elapsedTime inv62b = tab2poly({syt6#2, syt6#2, syt6#2, syt6#2, syt6#3,syt6#3}); -- check for independence-- this one factors
--factor inv62b
--elapsedTime inv62b = tab2poly({syt6#2, syt6#2, syt6#2, syt6#2, syt6#4,syt6#4}); -- check for independence
--time betti mingens (I62 + ideal inv62b)
--elapsedTime inv62b = tab2poly({syt6#3, syt6#3, syt6#3, syt6#3, syt6#4,syt6#4}); -- check for independence
--time betti mingens (I62 + ideal inv62b)
--#factor inv62b -- it factors... 
time inv63 = tab2poly(apply({0,0,0,4,4,4}, p-> syt6#p)); 
--time inv63= tab2poly({syt6#0, syt6#0, syt6#0, syt6#4, syt6#4,syt6#4}); -- check for independence
-- used 22.1626s (cpu); 16.3867s (thread); 0s (gc)
-- used 26.8649s (cpu); 17.0764s (thread); 0s (gc)

if checks == true then time betti mingens (I62 + ideal inv63)
-- used 7.17916s (cpu); 5.02327s (thread); 0s (gc)
-- used 9.97143s (cpu); 5.15708s (thread); 0s (gc)

if checks == true then time betti mingens (I2 + I4 + I62 + ideal inv63) -- and it is independent from the lower degree terms
-- used 22.0298s (cpu); 18.4967s (thread); 0s (gc)
-- used 28.0628s (cpu); 19.2103s (thread); 0s (gc)

up = unique permutations {0,0,0,4,4,4}
myPerms33 = {{0,1,2,3,4,5},{0,1,3,2,4,5},{0,1,3,4,2,5},{0,1,3,4,5,2},{0,3,1,2,4,5},{0,3,1,4,2,5},{0,3,1,4,5,2},{0,3,4,1,2,5},{0,3,4,1,5,2},{0,3,4,5,1,2},{3,0,1,2,4,5},{3,0,1,4,2,5},{3,0,1,4,5,2},{3,0,4,1,2,5},{3,0,4,1,5,2},{3,0,4,5,1,2},{3,4,0,1,2,5},{3,4,0,1,5,2},{3,4,0,5,1,2},{3,4,5,0,1,2}}
--for i to 19 list for j from i+1 to 19 list if 1== #unique {(perm myPerms33#i) inv63, (perm myPerms33#j) inv63 } then {i,j} else continue
--{{{0,19}},{{1,18}},{{2,17}},{{3,16}},{{4,15}},{{5,14}},{{6,13}},{{7,12}},{{8,11}},{{9,10}},{},{},{},{},{},{},{},{},{},{}}
--for i to 19 list for j from i+1 to 19 list if 1== rank source mingens ideal {(perm myPerms33#i) inv63, (perm myPerms33#j) inv63 } then {i,j} else continue
--
time F63 = for P in myPerms33_{0..9} list ((perm P) inv63);
-- used 5.70057s (cpu); 3.26402s (thread); 0s (gc)
-- used 2.0796s (cpu); 1.30475s (thread); 0s (gc)
#(unique F63)
I63 = ideal F63; -- 
-- F63 ={};rk =0;
-- pmList = for P in myPerms33 list ( tmpPerm = perm P; tmp = tmpPerm inv63; if rank source mingens ideal (F63|{tmp})>rk then (F63 = F63|{tmp}; rk=rk + 1;P ) else continue)

-- LL = {0,0,0,1,1,1}
-- netList for P in pmList list apply(P, i->LL_i)

if checks ==true then time betti mingens I63 -- 10, so there's 10 independent among the 15
if checks ==true then time betti mingens (I2 + I4 + I62 + I63) -- and they're independent from the lower degree terms
-- here is one more invariant that is of a different type than the others
time inv622 = tab2poly(apply({0,0,2,2,4,4}, p-> syt6#p)); 
 -- used 26.3111s (cpu); 15.8824s (thread); 0s (gc)
--elapsedTime inv622 = tab2poly({syt6#0, syt6#0, syt6#2, syt6#2, syt6#4,syt6#4}); --
if checks ==true then time betti mingens (I62 + I63 + ideal inv622)
if checks ==true then time betti mingens (I2 + I4 + I62 + I63 + ideal inv622) 
-- we can guess the permutations of the fillings that will be independent.
myPerms222 = {{0,1,2,3,4,5}, {0,2,1,3,4,5}, {0,1,2,4,3,5}, {0,2,1,4,3,5}, {0,4,2,3,1,5}}
LL = {0,0,2,2,4,4};
netList for P in myPerms222 list apply(P, i->LL_i)

time F622 = for P in myPerms222 list (perm P) inv622;
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
--restart
--allowableThreads = maxAllowableThreads
checks =false;
KK = ZZ/nextPrime(100);
--allowableThreads = 1--  maxAllowableThreads-1
input"/Users/lao0004/Documents/GitHub/YoungSymmetrizers/6factor.m2";
-- standard settings -- inv62 took 15s, betti mingens I2 + I4 + I62 took -- used 73.2505s (cpu); 18.8735s (thread); 0s (gc)
-- 14 threads: -- used 74.5685s (cpu); 16.8846s (thread); 0s (gc)
-- 1 thread: -- used 75.4438s (cpu); 17.2099s (thread); 0s (gc)
-- in d egree 8 it starts to take too long to evaluate and find the polynomials.
-- We can also get independence statements from point evaluations
-- to find independence we don't need to be able to take a jacobian, just evaluate at a lot of points.
-- to compute jacobians, we could hold all variables constant except one and evaluate them, then take the derivative with respect to the variable of interest...
-- try to work with 70 points for now
-- betti (I2+ I4 + I6)
-- pointLibrary = for i to 20 list apply(2^6,i-> gx_i => random(KK));
-- time M6 = matrix apply(pointLibrary, p -> apply(F2|F4|F6, f-> sub(f, p)));
-- rank M6
-- check non-zero
time tab2polyEval(apply({0,0,0,0,1,1}, p-> syt8#p), randomVec) -- about 0.56 seconds per point working over F_101

vecLibrary = for i to 16 list apply(2^6,i-> random(KK));
myUP42= uniquePermutations {0,0,0,0,1,1};
for2List = flatten apply(15, i-> apply(15, j-> (i,j) ))
Mm = new MutableMatrix from matrix(apply(15, i-> apply(15, j-> 0_KK)));
elapsedTime tmp = apply(apply(for2List, p -> schedule(p-> Mm_(p#0,p#1) = tab2polyEval(apply(6, k-> syt8#((myUP42#(p#0))#k)), vecLibrary#(p#1)), p)), taskResult);
 -- 71.6543s elapsed
allowableThreads = maxAllowableThreads
elapsedTime tmp = apply(apply(for2List, p -> schedule(p-> Mm_(p#0,p#1) = tab2polyEval(apply(6, k-> syt8#((myUP42#(p#0))#k)), vecLibrary#(p#1)), p)), taskResult);
-- changing to maxAllowableThreads actually increased the evaluation time!
 -- 162.955s elapsed
elapsedTime newCol = matrix apply(vecLibrary, p ->  {tab2polyEval(apply({0,0,0,0,1,1}, i-> syt8#i), p)}); 
 -- 5.64512s elapsed

binomial(6,2)
M82 = newCol;
GoodFill82 ={{0,0,0,0,1,1}}; rk =1; 
elapsedTime for P in drop(unique permutations {0,0,0,0,1,1},1) do (print("testing "|toString(P)|" with current rank "|toString(rk));
    time newCol = matrix apply(vecLibrary, p ->  {tab2polyEval(apply(6, i-> syt8#(P#i)), p)}); 
    if rank (M82|newCol)>rk then (GoodFill82 = GoodFill82|{P}; M82 = M82|newCol; rk=rk + 1; ) else continue);




time inv82 = tab2poly(apply({0,0,0,0,1,1}, p-> syt8#p)); -- this takes about 242s on MBP -- 
-- used 671.965s (cpu); 196.656s (thread); 0s (gc)
-- used 443.682s (cpu); 172.586s (thread); 0s (gc)
 -- when using julia, the tab2poly function here takes about 9.6s, which is about 20x faster than this!
#terms inv82
time F82 = for P in myPerms42 list ((perm P) inv82);
-- used 170.021s (cpu); 25.7427s (thread); 0s (gc)
-- used 27.6193s (cpu); 15.3379s (thread); 0s (gc)
I82 = ideal F82;
--#F82
--time betti mingens I82
-- used 238.253s (cpu); 43.2273s (thread); 0s (gc) -- only 10 independent, let's pick the first 10 like before
I82 = ideal F82_{0..9};
time betti mingens I82 -- there's 10 independent!
-- used 258.072s (cpu); 40.9351s (thread); 0s (gc)

-- you could try to numerically determine which are independent if the mingens command doesn't work quickly:
-- 
--time vl = matrix for j to 17 list {tab2polyEval(apply({0,0,0,0,1,1}, i-> syt8#i), vecLibrary_j)}
-- 
time vl = matrix for v in vecLibrary list {tab2polyEval(apply({0,0,0,0,1,1}, i-> syt8#i), v)}
tab2polyEval(apply({0,0,0,0,1,1}, i-> syt8#i), vecLibrary_0)
-- used 18.1877s (cpu); 6.67367s (thread); 0s (gc)

 -- used 18.1406s (cpu); 6.744s (thread); 0s (gc)
 -- 5.77566s elapsed
-- elapsedTime newCol = matrix parallelApply(vecLibrary, p ->  {tab2polyEval(apply({0,0,0,0,1,1}, i-> syt8#i), p)}); 
-- 10.3762s elapsed
-- expect around 18 times .6 s if it is done in serial
-- this used all cores (even not having entered the command --allowableThreads = maxAllowableThreads)
allowableThreads
M82 = newCol;
GoodFill82 ={{0,0,0,0,1,1}}; rk =1; 
elapsedTime for P in drop(unique permutations {0,0,0,0,1,1},1) do (print("testing "|toString(P)|" with current rank "|toString(rk));
    time newCol = matrix apply(vecLibrary, p ->  {tab2polyEval(apply(6, i-> syt8#(P#i)), p)}); 
    if rank (M82|newCol)>rk then (GoodFill82 = GoodFill82|{P}; M82 = M82|newCol; rk=rk + 1; ) else continue);
-- this finished on the server, and was much faster. 
--Fill8={{0,0,0,0,1,1},{0,0,0,1,0,1},{0,0,0,1,1,0},{0,0,1,0,0,1},{0,0,1,0,1,0},{0,0,1,1,0,0},{0,1,0,0,0,1},{0,1,0,0,1,0},{0,1,0,1,0,0},{0,1,1,0,0,0}}
-- this confirms that the first 10 are independent,
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