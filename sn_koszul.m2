-- compute the Koszul complex with the S_n action
-- we can compute the modules using SchurRings
--restart
n=5
R = QQ[x_0..x_(n-1)]
r = res ideal basis_1 R
apply(n, i-> rank(r.dd_i))
-- Buchsbaum-Eisenbud criterion for exactness is that the ranks of consecutive maps must equal the rank of the module:
apply(n, i-> rank(r.dd_i) + rank(r.dd_(i+1))  )
apply(n, i-> rank(source r.dd_i))
-- in the equivariant case, the B-E criterion says that every G-module that occurs in an R module must occur either before or after also, and all the equivariant maps are full rank.

-- non-equivariant version of Koszul via point evaluations:
-- take functions and turn them into column vectors by evaluating at random points:
rp = ()-> map(R,R, random(R^1, R^(numgens R)))

u = sum gens R
v = flatten entries ((basis_1 R)* gens ker diff(basis_1 R , u))

points = apply(5, i-> rp())
fn2vec = f-> matrix apply(points, pt-> {pt(f)})
fn2vec u
vecCat = vecs -> (m:=vecs#0; for i from 1 to #vecs-1 do m = m|vecs#i; m)
vecCat (fn2vec\v)
vvecs = fn2vec\v
-- first map is U+ V to exteriorPower(2, U+V) = U**V + exteriorPower(2,V)
-- second map is exteriorPower(2, U+V) = U**V + exteriorPower(2,V) to exteriorPower(3,U+ V)
for i from 0 to #v-1 list for j from i+1 to #v-1 list 
rank vecCat vvecs
rank transpose matrix apply(subsets( vvecs,2) , ss-> flatten entries exteriorPower(2, ss#0|ss#1))
rank transpose matrix apply(subsets( vvecs,3) , ss-> flatten entries exteriorPower(3, ss#0|ss#1|ss#2))
rank transpose matrix apply(subsets( vvecs,4) , ss-> flatten entries exteriorPower(4, ss#0|ss#1|ss#2|ss#3))
apply(subsets( vvecs,2) , ss-> maxMinors(ss#0|ss#1))
subsets(splice {0..3},2)

-- this gives you the modules:
loadPackage"SchurRings"
S = schurRing(s,n,GroupActing => "Sn")
rep = s_(n) + s_{n-1,1}
sr = flatten schurResolution(rep, {s_n}, DegreeLimit =>4)
apply(sr, t->dim t#1)


gens R
r = res ideal apply(2, i-> random(2,R))

mingens ideal flatten entries transpose r.dd_1 == mingens ideal flatten entries r.dd_2

restart
-- the idea is to use the Koszul complex to check if you have a regular sequence.
Rx = QQ[x_0..x_2]
rx = res ideal basis_1 Rx
Ry = QQ[y_0..y_3]
f0 = random(2, Ry)
f1 = random(3, Ry)
f2 = f0*f1;
-- f2 = random(5, Ry)
time fr = res coker matrix{{f0,f1,f2}}; -- this is the actual resolution
betti fr
fr.dd_0; rank oo
fr.dd_1; rank oo
fr.dd_2; rank oo
-- these are the BE-conditions:
apply(3, i-> rank(fr.dd_i) + rank(fr.dd_(i+1))  )
apply(3, i-> rank(source fr.dd_i))

Hfr = HH fr;
Hfr_0 == coker matrix{ {f0,f1,f2}}
Hfr_1 == 0 -- this is zero
Hfr_2 == 0 -- this is zero
-- now produce the Koszul complex by hand:
R = Rx**Ry
L = x_0*sub(f0,R) + x_1*sub(f1,R) + x_2*sub(f2,R);
ry1 = sub(diff(rx.dd_1,L), Ry);
ry2 = sub(diff(rx.dd_2,L), Ry);
ry3 = sub(diff(rx.dd_3,L), Ry);
-- check differentials square to zero:
ry1*ry2
ry2*ry3
-- 
ry1; rank ry1
ry2; rank ry2
ry3; rank ry3
-- 
rank ry1
rank ry1 + rank ry2
rank ry2 + rank ry3
rank ry3
MC = {ry1,ry2, ry3};
C = toChainComplex MC;
-- the Buchsbaum-Eisenbud criterion seems to be satisfied:
apply(4, i-> rank(C.dd_i) + rank(C.dd_(i+1))  )
apply(4, i-> rank(source C.dd_i))

-- make the koszul complex with built-in function:
C = koszul matrix{{f0,f1,f2}};
betti C
rank source mingens ideal C.dd_1 -- so not minimal generators
apply(4, i-> rank(C.dd_i) + rank(C.dd_(i+1))  )
apply(4, i-> rank(source C.dd_i))

mC = minimize C;
betti mC
betti C
-- now check the homology:
H = prune HH C
H_0 == coker matrix{ {f0,f1,f2}}
H_1 == 0 -- this is false when we make f2 = f0*f1 and true when f2 is random.
H_2 == 0

loadPackage"Depth"
isRegularSequence {f0,f1,f2}

rp = map(QQ, Rx, random(QQ^1,QQ^3))
rxp = rp rx
rank rxp.dd_1 , rank rx.dd_1 
rank rxp.dd_2 , rank rx.dd_2
rank rxp.dd_3 , rank rx.dd_3 

rpy = map(QQ, Ry, random(QQ^1,QQ^4))
rC = rpy C
rank rC.dd_1 , rank C.dd_1 
rank rC.dd_2 , rank C.dd_2
rank rC.dd_3 , rank C.dd_3 

restart
R = QQ[x_0..x_2]
K = koszul(matrix {{x_0,x_1, x_0^3}})
rank K.dd_1 
ker K.dd_1
image K.dd_2 
rank K.dd_2 
rank K.dd_3 
K.dd_1
K.dd_2
apply(4, i-> rank(K.dd_i) + rank(K.dd_(i+1))  )
apply(4, i-> rank(source K.dd_i))


H = prune HH K
H_0
H_1


K = koszul(matrix {{x_0,x_1, x_2}})
H = prune HH K
H_0
H_1
gens H_1


restart
R = QQ[x,y1,y2, z, Degrees=>{2,4,4,6}]
Kc = koszul(matrix {{x,x^3}})
Hc = prune HH Kc
Hc_0
Hc_1
Kc.dd
ker Kc.dd_0 
 image Kc.dd_1
ker Kc.dd_1
 image Kc.dd_2


Kb = koszul(matrix {{x,y1,y2,x^3}})
Hb = prune HH Kb
Hb_0
Hb_1
ker Kb.dd_0
target Kb.dd_2
Hb_0
(ker Kb.dd_0 / image Kb.dd_1)
mingens (ker Kb.dd_1 / image Kb.dd_2)

ker Kb.dd_1
image Kb.dd_2

Ka = koszul(matrix {{x,y1,y2,z}})
Ha = prune HH Ka
Ha_0
Ha_1
gens H_0
Kb.dd_1


rank K.dd_1 
ker K.dd_1
image K.dd_2 
rank K.dd_2 
rank K.dd_3 
K.dd_1
K.dd_2
apply(4, i-> rank(K.dd_i) + rank(K.dd_(i+1))  )
apply(4, i-> rank(source K.dd_i))


restart
loadPackage"Depth"
-- check the idea if any sequence of homogeneous minimal generators of length = codim is a regular sequence in a graded ring.
R = QQ[x_(0,0)..x_(2,2)]
X = genericMatrix(R,R_0,3,3)
I = minors(2, X)
codim I
G = mingens I
betti G

-- L = new HashTable from for S in subsets(9,4) list(
-- S=> isRegularSequence(flatten entries G_S)
-- );
-- net L

-- T = for S in subsets(9,4) list(
-- H = prune HH koszul(G_S);
-- (H_1 == 0) == isRegularSequence(flatten entries G_S)
-- );
-- tally T

myMinors = new HashTable from flatten for row in subsets(3,2) list for col in subsets(3,2) list(
{row, col} => det(X^row_col)
);
net myMinors
test = new HashTable from apply(subsets(keys myMinors, 3), ll -> ll=> isRegularSequence(apply(ll, l -> myMinors#l)));
net test
netList subsets(keys myMinors, 3)

successes = sort keys selectPairs(test, (k, v) -> v == true)
failures =  sort keys selectPairs(test, (k, v) -> v == false)

for fa in failures list codim ideal apply(fa, f -> myMinors#f)
pds = for fa in failures list primaryDecomposition ideal apply(fa, f -> myMinors#f);
net tally for pd in pds list codim \pd
tally for fa in successes list codim ideal apply(fa, f -> myMinors#f)
test#{{{0,1},{0,1}}, {{0,1},{0,2}}, {{0,2},{0,2}}}
test#{{{0,1},{0,1}}, {{1,2},{1,2}}, {{0,2},{0,2}}}

time test = new HashTable from apply(subsets(keys myMinors, 4), ll -> ll=> isRegularSequence(apply(ll, l -> myMinors#l)));
 
successes = sort keys selectPairs(test, (k, v) -> v == true);
failures =  sort keys selectPairs(test, (k, v) -> v == false);

tally for fa in failures list codim ideal apply(fa, f -> myMinors#f)
pds = for fa in failures list primaryDecomposition ideal apply(fa, f -> myMinors#f);
net tally for pd in pds list codim \pd
for fa in successes list codim ideal apply(fa, f -> myMinors#f)

for fa in successes list #primaryDecomposition ideal apply(fa, f -> myMinors#f)
tally oo


-- try again for 2x4 matrix
R = QQ[x_(0,0)..x_(1,3)]
rp = map(QQ, R, random(QQ^1, QQ^8))
X = genericMatrix(R,R_0,2,4)
I = minors(2, X)
codim I
-- G = mingens I
-- betti G, betti I
-- r = res I
-- net r.dd
-- betti r
-- syz gens I, gens I
myMinors = new HashTable from flatten for row in subsets(2,2) list for col in subsets(4,2) list(
{row, col} => det(X^row_col)
);

net for ss in subsets(keys myMinors, 3) list(
K1 = koszul(1,matrix {apply(ss, s -> myMinors#s)});
K2 = koszul(2,matrix {apply(ss, s -> myMinors#s)});
{mingens(ker K1 / image K2), isRegularSequence apply(ss, s -> myMinors#s)}
)
test2 = new HashTable from for ll in subsets(keys myMinors, 3) list( 
H = prune HH koszul matrix{ apply(ll, l -> myMinors#l)};
 ll=> H_1);
net test2

test = new HashTable from apply(subsets(keys myMinors, 3), ll -> ll=> isRegularSequence(apply(ll, l -> myMinors#l)));

testJ = new HashTable from apply(subsets(keys myMinors, 3), ll -> ll=> rank rp jacobian ideal apply(ll, l -> myMinors#l))
tally apply(keys test, key -> {test#key, testJ#key})

testcd = new HashTable from apply(subsets(keys myMinors, 3), ll -> ll=> codim ideal apply(ll, l -> myMinors#l))
tally apply(keys test, key -> {test#key, testJ#key, testcd#key})


successes = sort keys selectPairs(test, (k, v) -> v == true)
failures =  sort keys selectPairs(test, (k, v) -> v == false)

for fa in failures list codim ideal apply(fa, f -> myMinors#f)
for fa in successes list codim ideal apply(fa, f -> myMinors#f)
pds = for fa in failures list primaryDecomposition ideal apply(fa, f -> myMinors#f);
for pd in pds list codim \pd
(failures)#1
(successes)#0
tally for fa in successes list codim \primaryDecomposition ideal apply(fa, f -> myMinors#f)


-- T = for S in subsets(9,4) list(
-- H = prune HH koszul(G_S);
-- (H_1 == 0) == isRegularSequence(flatten entries G_S)
-- );
-- tally T
