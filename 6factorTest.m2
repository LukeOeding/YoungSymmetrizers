KK=ZZ/101;
syt8 = {{{1, 2}, {3, 4}, {5, 6}, {7, 8}}, {{1, 2}, {3, 4}, {5, 7}, {6, 8}}, {{1, 2}, {3, 5}, {4, 6}, {7, 8}}, {{1, 2}, {3, 5}, {4, 7}, {6, 8}}, {{1, 2}, {3,
        6}, {4, 7}, {5, 8}}, {{1, 3}, {2, 4}, {5, 6}, {7, 8}}, {{1, 3}, {2, 4}, {5, 7}, {6, 8}}, {{1, 3}, {2, 5}, {4, 6}, {7, 8}}, {{1, 3}, {2, 5}, {4, 7},
        {6, 8}}, {{1, 3}, {2, 6}, {4, 7}, {5, 8}}, {{1, 4}, {2, 5}, {3, 6}, {7, 8}}, {{1, 4}, {2, 5}, {3, 7}, {6, 8}}, {{1, 4}, {2, 6}, {3, 7}, {5, 8}}, {{1,
        5}, {2, 6}, {3, 7}, {4, 8}}};
X = KK[x_(0,0,0,0,0,0)..x_(1,1,1,1,1,1)];
gx = gens X;
rndPt = () -> map(KK, X, random(KK^1, KK^(numgens X)));
randomPoint = apply(2^6,i-> gx_(i) => random(KK));
randomVec = apply(2^6,i->  random(KK));
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
 for z to dg-1 do ( T=product(
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
tab2polyEval = (L,vec)->(
    a :=local a; b:=local b; c:= local c; d := local d; e := local e; f:= local f; x:= local x;
X := KK[x_(0,0,0,0,0,0)..x_(1,1,1,1,1,1)];
gx:= gens X;
rules := apply(2^6,i-> gx_(i) => vec_i);
 ta:=L#0; tb:=L#1; tc:=L#2; td:=L#3; te:=L#4; tf:=L#5;
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
 ); 
 sub(Fp,KK)
);
--end;
--restart;
--elapsedTime inv82 = tab2poly({syt8#0,syt8#0,syt8#0,syt8#0,syt8#1, syt8#1}); -- this takes about 242s on MBP, -- 351.488s elapsed on server.
pointLibrary = for i to 69 list apply(2^6,i-> gx_(i) => random(KK));
vecLibrary = for i to 69 list apply(2^6,i-> random(KK));
end;
restart
input"/Users/lao0004/Documents/GitHub/YoungSymmetrizers/6factorTest.m2";
vecLibrary = for i to 14 list apply(2^6,i-> random(KK));
elapsedTime tab2polyEval({syt8#0,syt8#0,syt8#0,syt8#0,syt8#1,syt8#1}, randomVec) 
-- used 0.490086s (cpu); 0.297022s (thread); 0s (gc) on MBP over ZZ/101
-- used 0.806229s (cpu); 0.598669s (thread); 0s (gc) on server over ZZ/101
allowableThreads = maxAllowableThreads
elapsedTime newCol = matrix parallelApply(vecLibrary, p ->  {tab2polyEval({syt8#0,syt8#0,syt8#0,syt8#0,syt8#1, syt8#1}, p)});  
 -- for 15 evaluations on MBP this takes 9.7s
 -- Saturn: -- 5.41572s elapsed
elapsedTime newColS = matrix apply(vecLibrary, p ->  {tab2polyEval({syt8#0,syt8#0,syt8#0,syt8#0,syt8#1, syt8#1}, p)});  
-- for 15 evaluations on MBP this takes 3.8s -- so faster to do this NOT in parallel.
-- on Saturn: -- 8.5933s elapsed
-- used 33.0477s (cpu); 19.6729s (thread); 0s (gc) on MBP over ZZ/101
 -- used 64.4684s (cpu); 42.6878s (thread); 0s (gc) on server over ZZ/101
allowableThreads = maxAllowableThreads

vecLibrary = for i to 14 list apply(2^6,i-> random(KK));
elapsedTime M82 = matrix apply(vecLibrary, p -> apply(uniquePermutations {0,0,0,0,1,1}, P->  tab2polyEval(apply(6, k-> syt8#(P#k)), p)));  
-- MBP:
-- 16.4892s elapsed for 4 rows
-- 75.1828s elapsed for 16 rows
-- 89.0724s elapsed for 21 rows
-- Saturn Server:


elapsedTime M82 = matrix parallelApply(vecLibrary, p -> parallelApply(uniquePermutations {0,0,0,0,1,1}, P->  tab2polyEval(apply(6, k-> syt8#(P#k)), p)));  
-- on MBP:
-- 32.8969s elapsed for 4 rows
-- 59.6232s elapsed for 15 rows
-- on Saturn:
-- 129.162s elapsed
-- if the vecLibrary is too big, then the process never starts, and you have to crash the system to start over
elapsedTime M82 = matrix apply(uniquePermutations {0,0,0,0,1,1}, P->  parallelApply(vecLibrary, p ->tab2polyEval(apply(6, k-> syt8#(P#k)), p)));  
-- on MBP:
-- 140.437s elapsed
-- on Saturn:
-- 83.8917s elapsed
Mm = new MutableMatrix from matrix(apply(15, i-> apply(15, j-> 0_KK)));
myUP42= uniquePermutations {0,0,0,0,1,1};
for2List = flatten apply(15, i-> apply(15, j-> (i,j) ))
-- elapsedTime parallelApply(for2List, p-> tab2polyEval(apply(6, k-> syt8#((myUP42#(p#0))#k)), vecLibrary#(p#1)))
-- -- on MBP:

-- -- on Saturn:
-- elapsedTime parallelApply(for2List, p-> tab2polyEval(apply(6, k-> syt8#((myUP42#(p#0))#k)), vecLibrary#(p#1)) , Strategy=>"raw")

Mm = new MutableMatrix from matrix(apply(15, i-> apply(15, j-> 0_KK)));
elapsedTime apply(apply(for2List, p -> schedule(p-> Mm_(p#0,p#1) = tab2polyEval(apply(6, k-> syt8#((myUP42#(p#0))#k)), vecLibrary#(p#1)), p)), taskResult);
rank matrix Mm
-- -- on MBP:
-- 153.327s elapsed
-- -- on Saturn:
 -- 89.8823s elapsed

myUP42= uniquePermutations {0,0,0,0,1,1};
for2List = flatten apply(15, i-> apply(15, j-> (i,j) ))
Mm = new MutableMatrix from matrix(apply(15, i-> apply(15, j-> 0_KK)));
elapsedTime tmp = apply(apply(for2List, p -> schedule(p-> Mm_(p#0,p#1) = tab2polyEval(apply(6, k-> syt8#((myUP42#(p#0))#k)), vecLibrary#(p#1)), p)), taskResult);

for i to 14 do (
elapsedTime apply(apply(15, j -> schedule(p-> Mm_(i,j) = tab2polyEval(apply(6, k-> syt8#((myUP42#(j))#k)), vecLibrary#(i)), p)), taskResult);
)
elapsedTime apply(apply(for2List, p -> schedule(p-> Mm_(p#0,p#1) = tab2polyEval(apply(6, k-> syt8#((myUP42#(p#0))#k)), vecLibrary#(p#1)), p)), taskResult);
rank matrix Mm


Mm = new MutableMatrix from matrix(apply(15, i-> apply(15, j-> 0_KK)));
elapsedTime for p in for2List do schedule(p-> Mm_(p#0,p#1) = tab2polyEval(apply(6, k-> syt8#((myUP42#(p#0))#k)), vecLibrary#(p#1)), p);
rank matrix Mm


elapsedTime M82 = matrix apply(uniquePermutations {0,0,0,0,1,1}, P->  parallelApply(vecLibrary, p ->tab2polyEval(apply(6, k-> syt8#(P#k)), p),Strategy=>"raw"));  
-- on MBP:
 -- -- 143.098s elapsed (39.9004s when was that?) elapsed for 15 vectors in vecLibrary
-- on Saturn:
 -- 89.3629s elapsed
elapsedTime newCol = matrix parallelApply(vecLibrary, p ->  {tab2polyEval({syt8#0,syt8#0,syt8#0,syt8#0,syt8#1, syt8#1}, p)});  
-- on MBP:
-- 9.33746s elapsed
-- on Saturn:
 -- 6.28434s elapsed
elapsedTime newCol = matrix parallelApply(vecLibrary, p ->  {tab2polyEval({syt8#0,syt8#0,syt8#0,syt8#0,syt8#1, syt8#1}, p)});  
-- on MBP:
-- 9.33746s elapsed
-- on Saturn:
 -- 6.28434s elapsed

elapsedTime newCol = matrix apply(vecLibrary, p ->  {tab2polyEval({syt8#0,syt8#0,syt8#0,syt8#0,syt8#1, syt8#1}, p)});  
-- on MBP:
-- 4.26847s elapsed
-- on Saturn:
 -- 9.02302s elapsed
M8 = newCol;
GoodFill82 ={{0,0,0,0,1,1}}; rk =1; 
elapsedTime for P in drop(uniquePermutations {0,0,0,0,1,1},1) do (print("testing "|toString(P)|" with current rank "|toString(rk));
    elapsedTime newCol = matrix parallelApply(vecLibrary, p ->  {tab2polyEval(apply(6, k-> syt8#(P#k)), p)});  
    if rank (M8|newCol)>rk then (GoodFill82 = GoodFill82|{P}; M8 = M8|newCol; rk=rk + 1; ) else continue);
-- this takes approximately 15x the time it takes for each main computation (20s or 40s)
-- on MBP with "apply"
-- 58.1771s elapsed
-- on Saturn with "apply"
 -- 126.365s elapsed
 --  -- 78.3278s elapsed with "parallelApply"
time newCol = matrix apply(vecLibrary, p ->  {tab2polyEval({syt8#0,syt8#0,syt8#0,syt8#0,syt8#3, syt8#3}, p)});  
rank(M8|newCol)
M8 = M8|newCol;
GoodFill82 =GoodFill82|{{0,0,0,0,3,3}}; rk =11; 
for P in drop(uniquePermutations {0,0,0,0,3,3},1) do (print("testing "|toString(P)|" with current rank "|toString(rk));
    time newCol = matrix apply(pointLibrary, p ->  {tab2polyEval(apply(6, k-> syt8#(P#k)), p)}); 
    if rank (M8|newCol)>rk then (GoodFill82 = GoodFill82|{P}; M8 = M8|newCol; rk=rk + 1; ) else continue);

time newCol = matrix apply(pointLibrary, p ->  {tab2polyEval({syt8#0,syt8#0,syt8#0,syt8#0,syt8#6, syt8#6}, p)});  
rank(M8|newCol)
M8 = M8|newCol;
GoodFill82 =GoodFill82|{{0,0,0,0,6,6}}; rk =26; 
for P in drop(uniquePermutations {0,0,0,0,6,6},1) do (print("testing "|toString(P)|" with current rank "|toString(rk));
    time newCol = matrix apply(pointLibrary, p ->  {tab2polyEval(apply(6, k-> syt8#(P#k)), p)}); 
    if rank (M8|newCol)>rk then (GoodFill82 = GoodFill82|{P}; M8 = M8|newCol; rk=rk + 1; ) else continue);

time newCol = matrix apply(pointLibrary, p ->  {tab2polyEval({syt8#0,syt8#0,syt8#0,syt8#0,syt8#8, syt8#8}, p)});  
rank(M8|newCol)
M8 = M8|newCol;
GoodFill82 =GoodFill82|{{0,0,0,0,8,8}}; rk =rank(M8|newCol); 
for P in drop(uniquePermutations {0,0,0,0,8,8},1) do (print("testing "|toString(P)|" with current rank "|toString(rk));
    time newCol = matrix apply(pointLibrary, p ->  {tab2polyEval(apply(6, k-> syt8#(P#k)), p)}); 
    if rank (M8|newCol)>rk then (GoodFill82 = GoodFill82|{P}; M8 = M8|newCol; rk=rk + 1; ) else continue);

rank(M8)
time newCol = matrix apply(pointLibrary, p ->  {tab2polyEval({syt8#0,syt8#0,syt8#0,syt8#0,syt8#9, syt8#9}, p)});  
rank(M8|newCol)
M8 = M8|newCol;
GoodFill82 =GoodFill82|{{0,0,0,0,9,9}}; rk =rank(M8|newCol); 
for P in drop(uniquePermutations {0,0,0,0,9,9},1) do (print("testing "|toString(P)|" with current rank "|toString(rk));
    time newCol = matrix apply(pointLibrary, p ->  {tab2polyEval(apply(6, k-> syt8#(P#k)), p)}); 
    if rank (M8|newCol)>rk then (GoodFill82 = GoodFill82|{P}; M8 = M8|newCol; rk=rk + 1; ) else continue);

netList GoodFill82
for i from 10 to 13 list(
print("testing "|toString(i));
time newCol = matrix apply(pointLibrary, p ->  {tab2polyEval({syt8#0,syt8#0,syt8#0,syt8#0,syt8#i, syt8#i}, p)});  
rtmp = rank (M8|newCol);
print("found new rank: "|toString(rtmp));
rtmp
)

GoodFill822 ={}; rk =0; 
for P in uniquePermutations {0,0,0,0,2,2} do (print("testing "|toString(P)|" with current rank "|toString(rk));
    time newCol = matrix apply(pointLibrary, p ->  {tab2polyEval({syt8#(P#0),syt8#(P#1),syt8#(P#2),syt8#(P#3),syt8#(P#4), syt8#(P#5)}, p)}); 
    if rank (M82|newCol)>rk then (GoodFill822 = GoodFill822|{P}; M82 = M82|newCol; rk=rk + 1; ) else continue);


 time inv82 = tab2poly({syt8#0,syt8#0,syt8#0,syt8#0,syt8#1, syt8#1}); 
 -- used 364.065s (cpu); 154.037s (thread); 0s (gc) on MBP over ZZ/101
  -- used 318.983s (cpu); 286.287s (thread); 0s (gc) on server over ZZ/101
end;
  input"6factorTest.m2";
