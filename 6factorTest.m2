KK=ZZ/101;
syt8 = {{{1, 2}, {3, 4}, {5, 6}, {7, 8}}, {{1, 2}, {3, 4}, {5, 7}, {6, 8}}, {{1, 2}, {3, 5}, {4, 6}, {7, 8}}, {{1, 2}, {3, 5}, {4, 7}, {6, 8}}, {{1, 2}, {3,
        6}, {4, 7}, {5, 8}}, {{1, 3}, {2, 4}, {5, 6}, {7, 8}}, {{1, 3}, {2, 4}, {5, 7}, {6, 8}}, {{1, 3}, {2, 5}, {4, 6}, {7, 8}}, {{1, 3}, {2, 5}, {4, 7},
        {6, 8}}, {{1, 3}, {2, 6}, {4, 7}, {5, 8}}, {{1, 4}, {2, 5}, {3, 6}, {7, 8}}, {{1, 4}, {2, 5}, {3, 7}, {6, 8}}, {{1, 4}, {2, 6}, {3, 7}, {5, 8}}, {{1,
        5}, {2, 6}, {3, 7}, {4, 8}}};
X = KK[x_(0,0,0,0,0,0)..x_(1,1,1,1,1,1)];
gx = gens X;
rndPt = () -> map(KK, X, random(KK^1, KK^(numgens X)));
randomPoint = apply(2^6,i-> gx_(i) => random(KK));
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
 ); 
 sub(Fp,KK)
);
--end;
--restart;
--elapsedTime inv82 = tab2poly({syt8#0,syt8#0,syt8#0,syt8#0,syt8#1, syt8#1}); -- this takes about 242s on MBP, -- 351.488s elapsed on server.
pointLibrary = for i to 69 list apply(2^6,i-> gx_(i) => random(KK));
end;
input"6factorTest.m2";
time tab2polyEval({syt8#0,syt8#0,syt8#0,syt8#0,syt8#1,syt8#1}, randomPoint) 
-- used 0.490086s (cpu); 0.297022s (thread); 0s (gc) on MBP over ZZ/101
-- used 0.806229s (cpu); 0.598669s (thread); 0s (gc) on server over ZZ/101
time newCol = matrix apply(pointLibrary, p ->  {tab2polyEval({syt8#0,syt8#0,syt8#0,syt8#0,syt8#1, syt8#1}, p)});  
-- used 33.0477s (cpu); 19.6729s (thread); 0s (gc) on MBP over ZZ/101
 -- used 64.4684s (cpu); 42.6878s (thread); 0s (gc) on server over ZZ/101
M8 = newCol;
GoodFill82 ={{0,0,0,0,1,1}}; rk =1; 
for P in drop(unique permutations {0,0,0,0,1,1},1) do (print("testing "|toString(P)|" with current rank "|toString(rk));
    time newCol = matrix apply(pointLibrary, p ->  {tab2polyEval({syt8#(P#0),syt8#(P#1),syt8#(P#2),syt8#(P#3),syt8#(P#4), syt8#(P#5)}, p)}); 
    if rank (M8|newCol)>rk then (GoodFill82 = GoodFill82|{P}; M8 = M8|newCol; rk=rk + 1; ) else continue);
-- this takes approximately 15x the time it takes for each main computation (20s or 40s)
time newCol = matrix apply(pointLibrary, p ->  {tab2polyEval({syt8#0,syt8#0,syt8#0,syt8#0,syt8#3, syt8#3}, p)});  
rank(M8|newCol)
M8 = M8|newCol;
GoodFill82 =GoodFill82|{{0,0,0,0,3,3}}; rk =11; 
for P in drop(unique permutations {0,0,0,0,3,3},1) do (print("testing "|toString(P)|" with current rank "|toString(rk));
    time newCol = matrix apply(pointLibrary, p ->  {tab2polyEval(apply(6, k-> syt8#(P#k)), p)}); 
    if rank (M8|newCol)>rk then (GoodFill82 = GoodFill82|{P}; M8 = M8|newCol; rk=rk + 1; ) else continue);

time newCol = matrix apply(pointLibrary, p ->  {tab2polyEval({syt8#0,syt8#0,syt8#0,syt8#0,syt8#6, syt8#6}, p)});  
rank(M8|newCol)
M8 = M8|newCol;
GoodFill82 =GoodFill82|{{0,0,0,0,6,6}}; rk =26; 
for P in drop(unique permutations {0,0,0,0,6,6},1) do (print("testing "|toString(P)|" with current rank "|toString(rk));
    time newCol = matrix apply(pointLibrary, p ->  {tab2polyEval(apply(6, k-> syt8#(P#k)), p)}); 
    if rank (M8|newCol)>rk then (GoodFill82 = GoodFill82|{P}; M8 = M8|newCol; rk=rk + 1; ) else continue);

time newCol = matrix apply(pointLibrary, p ->  {tab2polyEval({syt8#0,syt8#0,syt8#0,syt8#0,syt8#8, syt8#8}, p)});  
rank(M8|newCol)
M8 = M8|newCol;
GoodFill82 =GoodFill82|{{0,0,0,0,8,8}}; rk =rank(M8|newCol); 
for P in drop(unique permutations {0,0,0,0,8,8},1) do (print("testing "|toString(P)|" with current rank "|toString(rk));
    time newCol = matrix apply(pointLibrary, p ->  {tab2polyEval(apply(6, k-> syt8#(P#k)), p)}); 
    if rank (M8|newCol)>rk then (GoodFill82 = GoodFill82|{P}; M8 = M8|newCol; rk=rk + 1; ) else continue);

rank(M8)
time newCol = matrix apply(pointLibrary, p ->  {tab2polyEval({syt8#0,syt8#0,syt8#0,syt8#0,syt8#9, syt8#9}, p)});  
rank(M8|newCol)
M8 = M8|newCol;
GoodFill82 =GoodFill82|{{0,0,0,0,9,9}}; rk =rank(M8|newCol); 
for P in drop(unique permutations {0,0,0,0,9,9},1) do (print("testing "|toString(P)|" with current rank "|toString(rk));
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
for P in unique permutations {0,0,0,0,2,2} do (print("testing "|toString(P)|" with current rank "|toString(rk));
    time newCol = matrix apply(pointLibrary, p ->  {tab2polyEval({syt8#(P#0),syt8#(P#1),syt8#(P#2),syt8#(P#3),syt8#(P#4), syt8#(P#5)}, p)}); 
    if rank (M82|newCol)>rk then (GoodFill822 = GoodFill822|{P}; M82 = M82|newCol; rk=rk + 1; ) else continue);


 time inv82 = tab2poly({syt8#0,syt8#0,syt8#0,syt8#0,syt8#1, syt8#1}); 
 -- used 364.065s (cpu); 154.037s (thread); 0s (gc) on MBP over ZZ/101
  -- used 318.983s (cpu); 286.287s (thread); 0s (gc) on server over ZZ/101
end;
  input"6factorTest.m2";
