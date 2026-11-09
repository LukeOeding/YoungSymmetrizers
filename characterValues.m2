-- compute characters by using point evaluations
--restart
input"/Users/lao0004/Documents/GitHub/YoungSymmetrizers/5factor.m2"
loadPackage"SpechtModule"
charTable = characterTable 5
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
partToClassRep({3,1,1,1})
(partitions 5)#3
partToClassRep({3,1,1})

reps = apply(partitions 5, p -> partToClassRep toList  p)

G = permutations(5);
con = h->  sort unique apply(G, g-> ((inversePermutation(g))_h)_g)
conjClasses = (unique (con\reps));
orders = apply(conjClasses, c-> #c) -- sizes of conjugacy classes.
ip = (av,bv) -> 1/120*trace( av*diagonalMatrix(orders)*transpose bv)

computeChar = polList ->(
    pts:= apply(#polList, i-> rndPt());
    makeMg := myPerm -> sub(matrix apply(pts,pt-> apply(polList, f-> (pt*myPerm) f)), baseRing X) ; 
    print("computing the representation matrices");
    Mgs := for rep in reps list time makeMg(perm rep);
    M := last Mgs; 
    apply(Mgs, MMg->trace(M^(-1)*MMg ))
)
computeCharProd = polProdList ->(
    pts:= apply(#polProdList, i-> rndPt());
    makeMg := myPerm -> sub(matrix apply(pts,pt-> apply(polProdList, pf-> value(product(toList pf, f-> (pt*myPerm) f)))), baseRing X) ; 
    print("computing the representation matrices");
    Mgs := for rep in reps list time makeMg(perm rep);
    M := last Mgs; 
    apply(Mgs, MMg->trace(M^(-1)*MMg ))
)
decomposeCharacter = mychar ->    (
    dc := apply(7, i-> ip(matrix({mychar}), CT^{i}));
    print(netList {toList\(partitions 5), dc});
    dc)
end;    
-- character of the degree 4 invariants
betti ideal F4
char4 = computeChar(F4)
decomposeCharacter(char4)

-- character for the deg 6 invariant
F6 ={inv6};
char6 = computeChar(F6)
decomposeCharacter(char6)

-- two ways to compute an invariant minimial generating set
-- 1) put in new polynomials one at a time, but take an invariant complement to the orbit closure?
-- 2) project to the invariant subspaces inside an orbit, and check if each invariant subspace increases the number of minimal generators


--goodTabList8 = {{0,0,2,9,12},{0,0,2,12,9},{0,0,9,2,12},{0,0,9,12,2},{0,0,12,2,9},{0,2,0,9,12},{0,2,9,0,12},{0,2,9,12,0},{0,9,0,2,12},{0,9,0,12,2},{0,9,2,0,12},{0,9,2,12,0},{0,9,12,0,2},{0,9,12,2,0},{0,12,0,2,9},{2,0,0,9,12},{2,0,9,0,12},{2,0,9,12,0},{2,9,0,0,12},{2,9,0,12,0},{2,9,12,0,0},{9,0,0,2,12},{9,0,0,12,2},{9,0,2,0,12},{9,0,2,12,0},{9,0,12,0,2},{9,12,0,0,2},{12,0,0,2,9},{0,2,3,12,13},{0,2,3,13,12},{0,2,12,3,13},{0,2,12,13,3},{0,2,13,3,12},{0,3,2,12,13},{0,0,12,12,12},{0,3,6,9,13}};
--print tally (sort\goodTabList8)
-- #goodTabList8
-- pts = apply(#goodTabList8, i-> rndPt());
-- time inv8v = tab2polyMap(apply(5, i-> syt8#(goodTabList8#0#i)), pts#0);

--time M = sub(matrix apply(pts,pt-> apply(goodTabList8, Q -> tab2polyMap(apply(5, i-> syt8#(Q#i)), pt))), baseRing X) ; 
-- expect 125s 
-- In degree 8 we have invariants that are products of degree 4 invariants:
F44 = flatten apply(5, i-> apply(5-i, j-> F4#i*F4#(i+j)));
time char44 = computeChar(F44)
dchar44 = decomposeCharacter(char44)

-- then we add in invariants that are independent of those we already have.
P = {0,0,2,9,12};
tgt = mgTabList8_{0..19} -- just grab the known contributors to the mingens
--tgt =goodTabList8_{0..27}; -- use this if you want to just have an independent set of mingens
--tgt = uniquePermutations P;
ct = 0; myPerms = flatten for i to #tgt-1 list if ct == i then for xx in permutations({0,1,2,3,4}) list if ct == i and P_xx == tgt#i then (ct = ct+1; xx ) else continue

elapsedTime inv8 = tab2poly(apply(5, i-> syt8#(P#i))); -- time = 14s
time F8big = unique apply(myPerms, p-> (perm p) inv8); -- time if 1.2s if you only use the good permutations,  4s if you use the 60 permutations of 0,0,2,9,12
#F8big -- we will find that this spanns a 28 dimensional invariant subspace, and if we take a minimal spanning set, that will still span an invariant space

-- point evaluation method takes about 19.7s, all the other steps are essentially free.
pts = apply(#F8big, i-> rndPt());
elapsedTime Mb = sub(matrix apply(pts_{0..#F8big-1},pt-> apply(F8big, f-> pt f)), baseRing X) ; 
rankList = apply(#F8big, i-> rank(Mb_{0..i-1}^{0..i-1}))
jumps = {0}|for i to #rankList-1 list if rankList#(i) > rankList#(i-1) then i else continue;
rank(Mb_jumps^jumps)
F8a = F8big_jumps;

betti ideal F8a
time betti mingens ideal F8a -- verify that they're independent
-- compute character of the entire module (F8 is only invariant after modding out by F44):
time char448a = computeChar(F44|F8a)
decomposeCharacter(char448a)
last mgTabList8
-- what happens with the next filling?
unique(sort\ goodTabList8)

P = {0,3,6,9,13};
inv8d = tab2poly(apply(5, i-> syt8#(P#i))); -- time = 50s
F8d = {inv8d};
betti mingens ideal (F44|F8a|F8d)
time char448ad = computeChar(F44|F8a|F8d)
dcharr8 = decomposeCharacter(char448ad)
dcharr8 - dchar44
F8 = (F44|F8a|F8d);
betti ideal F8
I8 = ideal F8;
I = I4 + I6 + I8 + I10;
betti I
--- of course, if you want to discover the spanning vectors, you can go like this:
--netList goodTabList8
P = {0,2,3,12,13};
--tgt = delete(last mgTabList8, mgTabList8)

tgt = uniquePermutations P;
ct = 0; myPerms = flatten for i to #tgt-1 list if ct == i then for xx in permutations({0,1,2,3,4}) list if ct == i and P_xx == tgt#i then (ct = ct+1; xx ) else continue
#myPerms
elapsedTime inv8b = tab2poly(apply(5, i-> syt8#(P#i))); -- time = 30s
time F8bigb = unique apply(myPerms, p-> (perm p) inv8b); -- time 76s
#F8bigb -- we will find that this spanns a 29? dimensional invariant subspace, and if we take a minimal spanning set, that will still span an invariant space
betti mingens ideal F8bigb
-- point evaluation method takes most of the time, the rest of the steps are essentially free

-- if you know the dimension in advance, you can make a smaller matrix to compute the rank...
pts = apply(#F8bigb, i-> rndPt());
--Mb = sub(matrix apply(pts_{0..#F8bigb-1},pt-> apply(F8bigb, f-> pt f)), baseRing X) ; 
Mb = sub(matrix apply(pts_{0..29-1},pt-> apply(F8bigb, f-> pt f)), baseRing X) ; 
-- takes 85s?
Mb;
rankList = apply(120, i-> rank(Mb_{0..i}));
rankList
jumps = {0}|for i to #rankList-1 list if rankList#(i) > rankList#(i-1) then i else continue;
jumps
rank Mb
apply(29, i-> rank Mb_(jumps_{0..i}))
rank(Mb_jumps) -- this might not have full rank?
F8b = F8bigb_jumps;
time char8b = computeChar(F8b)
decomposeCharacter(char8a)
decomposeCharacter(char8b)
-- if F8|F8b spans an invariant space, and F8 spans an invariant space,
-- then a subset will span their span, and if the whole space is invariant we should get
-- a basis of the invariant space...

betti ideal F44
betti mingens ideal F44
F8big = F44|F8;
pts = apply(#F8big, i-> rndPt());
Mb = sub(matrix apply(pts_{0..#F8big-1},pt-> apply(F8big, f-> pt f)), baseRing X) ; 

rankList = apply(#F8big, i-> rank(Mb_{0..i}));
rankList
jumps = {0}|for i to #rankList-1 list if rankList#(i) > rankList#(i-1) then i else continue;
jumps
rank Mb
apply(34, i-> rank Mb_(jumps_{0..i}))
rank(Mb_jumps) -- this might not have full rank?
F448 = F8big_jumps;
time char448 = computeChar(F448)
decomposeCharacter(char448)


F8bigab = F448|F8b;
pts = apply(#F8bigab, i-> rndPt());
-- we know what the dimension is supposed to be... 34
Mb = sub(matrix apply(pts_{0..34-1},pt-> apply(F8bigab, f-> pt f)), baseRing X) ; 
-- takes 85s?
rank Mb -- this is 34, so it says that the orbit of this polynomial inv8b doesn't add to the list of minimal generators
rankList = apply(#F8bigab, i-> rank(Mb_{0..i}));
jumps = {0}|for i to #rankList-1 list if rankList#(i) > rankList#(i-1) then i else continue;
rank Mb
F8ab = F8bigab_jumps;
time char8ab = computeChar(F8ab)
decomposeCharacter(char8ab) -- this confirms that you don't get anything new here...


-- what happens with the next filling?
P = {0,0,12,12,12};
tgt = uniquePermutations P;
ct = 0; myPerms = flatten for i to #tgt-1 list if ct == i then for xx in permutations({0,1,2,3,4}) list if ct == i and P_xx == tgt#i then (ct = ct+1; xx ) else continue
#myPerms
elapsedTime inv8c = tab2poly(apply(5, i-> syt8#(P#i))); -- time = 50s
time F8bigc = unique apply(myPerms, p-> (perm p) inv8c); -- time .2s
#F8bigc

-- point evaluation method takes most of the time, the rest of the steps are essentially free

pts = apply(#F8bigc, i-> rndPt());
Mc = sub(matrix apply(pts_{0..#F8bigc-1},pt-> apply(F8bigc, f-> pt f)), baseRing X) ; 
rankList = apply(#F8bigc+1, i-> rank(Mc_{0..i-1}^{0..i-1}))
jumps = for i to #F8bigc-1 list if rankList#(i+1) > rankList#i then i else continue;
jumps
rank(Mc_jumps^jumps) -- this might not have full rank?
F8c = F8bigc_jumps;
betti ideal F8bigc
time betti mingens ideal F8bigc

time char8c = computeChar(F8c)
decomposeCharacter(char8c) -- this confirms that you don't get anything new here...


F8big = F448|F8bigc;
pts = apply(#F8big, i-> rndPt());
Mb = sub(matrix apply(pts_{0..#F8big-1},pt-> apply(F8big, f-> pt f)), baseRing X) ; 
rank Mb

rankList = apply(#F8big, i-> rank(Mb_{0..i}));
jumps = {0}|for i to #rankList-1 list if rankList#(i) > rankList#(i-1) then i else continue;
jumps
rank Mb
F8abc = F8big_jumps;
time char8abc = computeChar(F8abc)
decomposeCharacter(char8abc) -- this confirms that you don't get anything new here...



time betti mingens (ideal F8+ideal F8b + ideal F8c ) 
-- takes about 46s - get a 35 dimensional space from their span
time betti mingens (I4 + ideal F8+ideal F8b + ideal F8c)
 -- takes about 50s -- and it tells us that the polynomials from this new module are not mingens after including the deg4 terms
-- compute the character anyway:
makeMgc = myPerm -> sub(matrix apply(pts_{0..#F8c-1},pt-> apply(F8c, f-> (pt*myPerm) f)), baseRing X) ; 
Mgsc = for rep in reps list time makeMgc(perm rep);
-- maybe this takes about 20 minutes with YS, or with polynomials 2s for each step (over ZZ/101) or 8s over ZZ/nextPrime(10000) per iteration
-- the way that the partitions list out, the one corresponding to the identity is last:
M = last Mgsc; mychi = apply(Mgsc, MMg->trace(M^(-1)*MMg ))
apply(7, i-> ip(matrix({mychi}), CT^{i}))
-- character {0,0,1,1,2,4,10}
-- decomposition {1,1,1,0,0,0,0}
netList {toList\(partitions 5),  apply(7, i-> ip(matrix({mychi}), CT^{i}))}



P = {0,3,6,9,13};
tgt = uniquePermutations P;
ct = 0; myPerms = flatten for i to #tgt-1 list if ct == i then for xx in permutations({0,1,2,3,4}) list if ct == i and P_xx == tgt#i then (ct = ct+1; xx ) else continue
#myPerms
elapsedTime inv8d = tab2poly(apply(5, i-> syt8#(P#i))); -- time = 50s
time F8bigd = unique apply(myPerms, p-> (perm p) inv8d); -- time 
#F8bigd


-- just add in the last invariant to fill out the basis of the degree 8 invariants:
F8abcd = F8abc|{inv8d};
betti mingens ideal F8abcd
pts = apply(#F8abcd, i-> rndPt());
Mb = sub(matrix apply(pts_{0..#F8abcd-1},pt-> apply(F8abcd, f-> pt f)), baseRing X) ; 
rank Mb

rankList = apply(#F8abcd, i-> rank(Mb_{0..i}));
jumps = {0}|for i to #rankList-1 list if rankList#(i) > rankList#(i-1) then i else continue;
jumps
rank Mb
F8abcd = F8abcd_jumps;
time char8abcd = computeChar(F8abcd)
decomposeCharacter(char8abcd) -- this confirms that you don't get anything new here...






-- point evaluation method takes most of the time, the rest of the steps are essentially free

pts = apply(#F8bigd, i-> rndPt());
time Md = sub(matrix apply(pts_{0..#F8bigd-1},pt-> apply(F8bigd, f-> pt f)), baseRing X) ; 
rankList = apply(#F8bigd+1, i-> rank(Md_{0..i-1}^{0..i-1}))
jumps = for i to #F8bigd-1 list if rankList#(i+1) > rankList#i then i else continue;
jumps
rank(Md_jumps^jumps) -- this might not have full rank?
F8d = F8bigd_jumps;
betti ideal F8bigd
time betti mingens ideal F8bigd
time betti mingens (I4 + ideal F8bigd)
time betti mingens (I4 + ideal F8big)

time betti mingens (ideal F8+ideal F8b + ideal F8c + ideal F8d ) 
-- takes about 70s - get a 36 dimensional space from their span
time betti mingens (I4 + ideal F8+ideal F8b + ideal F8c+ ideal F8d)
 -- takes about 125s -- and it tells us that one new polynomial is in the mingens
-- compute the character anyway:
makeMgd = myPerm -> sub(matrix apply(pts_{0..#F8d-1},pt-> apply(F8d, f-> (pt*myPerm) f)), baseRing X) ; 
Mgsd = for rep in reps list time makeMgd(perm rep);
-- maybe this takes about 20 minutes with YS, or with polynomials 2s for each step (over ZZ/101) or 8s over ZZ/nextPrime(10000) per iteration 
-- the way that the partitions list out, the one corresponding to the identity is last:
M = last Mgsd; mychi = apply(Mgsd, MMg->trace(M^(-1)*MMg ))
apply(7, i-> ip(matrix({mychi}), CT^{i}))
-- character {−1,−1,1,−1,5,7,29}
-- decomposition {1,2,3,0,1,0,0}
netList {toList\(partitions 5),  apply(7, i-> ip(matrix({mychi}), CT^{i}))}

-- value for projection to trivial:
time sum(permutations 5, pm-> ((pts#0)*(perm pm)) inv8b)
-- compute the number of copies of the trivial rep:
trivProj = (pt,f)->  sum(permutations 5, pm-> ((pt)*(perm pm)) f) ;

matrix apply(4, i->  trivProj_(pts#i)\{inv8, inv8b, inv8c, inv8d })
rank oo
betti mingens (I4+  ideal {inv8, inv8b, inv8c, inv8d })

#F8bigd

time tmp = inv8%I4;
mingensto8 = flatten entries mingens(I4 + ideal(fullF8) );
betti ideal mingensto8
-- just take the mingens
mg8= mingensto8_{5..25};
betti ideal mg8-- 21 mininmal generators
-- point evaluation method takes most of the time, the rest of the steps are essentially free
pts = apply(#mg8, i-> rndPt());

-- compute the character 
makeMg = myPerm -> sub(matrix apply(pts_{0..#mg8-1},pt-> apply(mg8, f-> (pt*myPerm) f)), baseRing X) ; 
Mgs = for rep in reps list time makeMg(perm rep);
M = last Mgs; mychi = apply(Mgs, MMg->trace(M^(-1)*MMg ))
apply(7, i-> ip(matrix({mychi}), CT^{i}))
netList {toList\(partitions 5),  apply(7, i-> ip(matrix({mychi}), CT^{i}))}
-- we see that the deg 8 only minimal generators module provided is not G-invariant, so the character doesn't work
mg8t = (perm {1,0,2,3,4})\ mg8;
betti mingens ideal(mg8|mg8t) -- strangely, this has dimension even bigger than the ambient space?

-- try again with the entire module?
mg8= mingensto8;
betti ideal mg8-- 5+21 mininmal generators
-- point evaluation method takes most of the time, the rest of the steps are essentially free
-- compute the character 
pts = apply(#mg8, i-> rndPt());
makeMg = myPerm -> sub(matrix apply(pts_{0..#mg8-1},pt-> apply(mg8, f-> (pt*myPerm) f)), baseRing X) ; 
Mgs = for rep in reps list time makeMg(perm rep);
M = last Mgs; mychi = apply(Mgs, MMg->trace(M^(-1)*MMg ))
apply(7, i-> ip(matrix({mychi}), CT^{i}))
netList {toList\(partitions 5),  apply(7, i-> ip(matrix({mychi}), CT^{i}))}
-- we see that the minimal generators module provided is not G-invariant, so the character doesn't work
mg8t = (perm {1,0,2,3,4})\ mg8;
betti mingens ideal(mg8|mg8t)




time tooFull8 = flatten entries mingens ideal tooFull8;
betti ideal tooFull8 -- 36 mininmal generators
-- point evaluation method takes most of the time, the rest of the steps are essentially free
fullF8 = tooFull8;
pts = apply(#fullF8, i-> rndPt());
betti ideal fullF8
--time betti mingens (I4 + ideal fullF8)

-- compute the character anyway:
-- if you just take the minimal generators, they might take longer to evaluate than the polynomials that come from the young symmetrizers
makeMgf = myPerm -> sub(matrix apply(pts_{0..#fullF8-1},pt-> apply(fullF8, f-> (pt*myPerm) f)), baseRing X) ; 
Mgsf = for rep in reps list time makeMgf(perm rep);
-- maybe this takes about 20 minutes with YS, or with polynomials 2s for each step (over ZZ/101) or 8s over ZZ/nextPrime(10000) per iteration 
-- the way that the partitions list out, the one corresponding to the identity is last:
M = last Mgsf; mychi = apply(Mgsf, MMg->trace(M^(-1)*MMg ))
apply(7, i-> ip(matrix({mychi}), CT^{i}))
-- character {−1,−1,1,−1,5,7,29}
-- decomposition {1,2,3,0,1,0,0}
netList {toList\(partitions 5),  apply(7, i-> ip(matrix({mychi}), CT^{i}))}


-- this computes the character decomposition, in the order that the character table is made in:
netList goodTabList10
first goodTabList10
all10 = goodTabList10|{{0,17,24,37,40}, {0,17,24,40,37}, {0,17,40,24,37}, {0,24,17,37,40}, {17,0,24,37,40},{0,12,30,31,38}}
P = {0,0,4,34,40}
tally (sort\goodTabList10)
elapsedTime inv10 = tab2poly(apply(5, i-> syt10#(P#i))); -- time 757s
tgt = goodTabList10_{0..8} -- just grab the known contributors to the mingens
ct = 0; myPerms = flatten for i to #tgt-1 list if ct == i then for xx in permutations({0,1,2,3,4}) list if ct == i and P_xx == tgt#i then (ct = ct+1; xx ) else continue
time F10 = apply(myPerms, p-> (perm p) inv10); --8s
P = last goodTabList10; elapsedTime inv10b = tab2poly(apply(5, i-> syt10#(P#i))); -- time 428s?
F10 = F10|{inv10b};
F46 = apply(F4, f4 -> f4*inv6);
time char10 = computeChar(F10|F46)
decomposeCharacter(char10) -- space is not invariant unless you include what you get from deg 4 times deg 6
betti ideal F10
time char46 =  computeChar(F46)
decomposeCharacter(char46) -- space is not invariant unless you include what you get from deg 4 times deg 6
betti ideal (F10|F46)
I10 = ideal (F10|F46);
betti I8
-- try to just compute with the first orbit span via point evaluations and not computing the polynomial?
time M10 = sub(matrix apply(pts_{0..#goodTabList10-2},pt-> apply(goodTabList10_{0..#goodTabList10-2}, Q -> tab2polyMap(apply(5, i-> syt10#(Q#i)), pt))), baseRing X) ; -- expect 20s
det M10
makeMg10 = myPerm -> sub(matrix apply(pts_{0..#goodTabList10-2},pt-> apply(goodTabList10_{0..#goodTabList10-2}, Q -> tab2polyMap(apply(5, i-> syt10#(Q#i)), pt*myPerm))), baseRing X); 
Mgs10 = for rep in reps list time makeMg10(perm rep);-- expect 200s?
mychi10 = apply(Mgs10, MMg->trace(M10^(-1)*MMg ))
apply(7, i-> ip(matrix({mychi10}), charMat^{i}))


time M10l = sub(matrix apply(pts_{0},pt-> apply(goodTabList10_{#goodTabList10-1}, Q -> tab2polyMap(apply(5, i-> syt10#(Q#i)), pt))), baseRing X) ; 
 M10l
makeMg10l = myPerm -> sub(matrix apply(pts_{0},pt-> apply(goodTabList10_{#goodTabList10-1}, Q -> tab2polyMap(apply(5, i-> syt10#(Q#i)), pt*myPerm))), baseRing X); 

-- this doesn't work because we don't have a spanning vector for the space, but rather a vector that projects non-zero to that space. 
Mgs10l = for rep in reps list time makeMg10l(perm rep);
Mgs10l
mychi10l = apply(Mgs10l, MMg->trace(M10l^(-1)*MMg ))
apply(7, i-> ip(matrix({mychi10l}), charMat^{i}))

-- for higher degree, you could split up by orbits, since the span of an orbit will be a sum of irreducibles.
-- this means that you don't have to make such large matrices. 
-- but we didn't build our list of permutations as the orbits, but rather the new things that increase the number of independent vectors in the list. 

all10 = goodTabList10|{{0,17,24,40,37}, {0,17,40,24,37}, {0,24,17,37,40}, {17,0,24,37,40},{0,12,30,31,38}}
-- try to just compute with the first orbit span
time M10 = sub(matrix apply(pts_{0..#all10-1},pt-> apply(all10, Q -> tab2polyMap(apply(5, i-> syt10#(Q#i)), pt))), baseRing X) ; 
-- expect 112s
rank M10
makeMg10 = myPerm -> sub(matrix apply(pts_{0..#all10-1},pt-> apply(all10, Q -> tab2polyMap(apply(5, i-> syt10#(Q#i)), pt*myPerm))), baseRing X); 

Mgs10 = for rep in reps list time makeMg10(perm rep);-- expect 7*90s?

mychi10 = apply(Mgs10, MMg->trace(M10^(-1)*MMg ))
-- {15,−7,3,−1,0,3,−1}
apply(7, i-> ip(matrix({mychi10}), charMat^{i}))
-- {0,2,0,2,0,0,1}
-- 2U' + 2V' + W'

goodTabList12={{0,6,35,113,117},{0,6,35,117,113},{0,6,113,35,117},{0,6,113,117,35},{0,6,117,35,113},{0,6,117,113,35},{0,35,6,113,117},{0,35,6,117,113},{0,35,113,6,117},{0,35,113,117,6},{0,35,117,6,113},{0,35,117,113,6},{0,113,6,35,117},{0,113,6,117,35},{0,113,35,6,117},{0,113,35,117,6},{0,113,117,6,35},{0,113,117,35,6},{0,117,6,35,113},{0,117,6,113,35},{0,117,35,6,113},{0,117,35,113,6},{0,117,113,6,35},{0,117,113,35,6},{6,0,35,113,117},{6,0,35,117,113},{6,0,113,35,117},{6,0,113,117,35},{6,0,117,35,113},{6,0,117,113,35},{6,35,0,113,117},{6,35,0,117,113},{6,35,113,0,117},{6,35,113,117,0},{6,35,117,0,113},{6,35,117,113,0},{6,113,0,35,117},{6,113,0,117,35},{6,113,35,0,117},{6,113,35,117,0},{6,113,117,0,35},{6,113,117,35,0},{6,117,0,35,113},{6,117,0,113,35},{6,117,35,0,113},{6,117,35,113,0},{6,117,113,0,35},{6,117,113,35,0},{35,0,6,113,117},{35,0,6,117,113},{35,0,113,6,117},{35,0,113,117,6},{35,0,117,6,113},{35,0,117,113,6},{35,6,0,113,117},{35,6,0,117,113},{35,6,113,0,117},{35,6,113,117,0},{35,6,117,0,113},{35,6,117,113,0},{35,113,0,6,117},{35,113,0,117,6},{35,113,6,0,117},{35,113,6,117,0},{35,113,117,0,6},{35,113,117,6,0},{35,117,0,6,113},{35,117,0,113,6},{35,117,6,0,113},{35,117,6,113,0},{35,117,113,0,6},{113,0,6,35,117},{113,0,6,117,35},{113,0,35,6,117},{113,0,35,117,6},{113,0,117,6,35},{113,0,117,35,6},{113,6,0,35,117},{113,6,0,117,35},{113,6,35,0,117},{113,6,35,117,0},{113,6,117,0,35},{113,6,117,35,0},{113,35,0,6,117},{113,35,0,117,6},{113,35,6,0,117},{113,35,6,117,0},{113,35,117,0,6},{113,35,117,6,0},{113,117,0,6,35},{113,117,0,35,6},{113,117,6,0,35},{113,117,6,35,0},{113,117,35,0,6},{117,0,6,35,113},{117,0,6,113,35},{117,0,35,6,113},{117,0,35,113,6},{117,0,113,6,35},{117,0,113,35,6},{117,6,0,35,113},{117,6,0,113,35},{117,6,35,0,113},{117,6,35,113,0},{117,6,113,0,35},{117,35,0,6,113},{117,35,0,113,6},{117,35,113,0,6},{117,113,0,6,35},{117,113,0,35,6},{117,113,6,0,35},{0,5,32,36,91},{0,5,32,91,36},{0,5,36,32,91},{0,5,36,91,32},{0,5,91,32,36},{0,5,91,36,32},{0,32,5,36,91},{0,32,5,91,36},{0,32,36,5,91},{0,32,36,91,5},{0,32,91,5,36},{0,32,91,36,5},{0,36,5,32,91},{0,36,5,91,32},{0,36,32,5,91},{0,36,32,91,5},{0,36,91,5,32},{0,36,91,32,5},{0,91,5,32,36},{0,91,5,36,32},{0,91,32,5,36},{5,0,32,36,91},{5,0,32,91,36},{5,0,36,32,91},{5,0,36,91,32},{5,0,91,32,36},{5,32,0,36,91},{5,32,0,91,36},{5,32,36,0,91},{5,32,36,91,0},{5,32,91,0,36},{5,36,0,32,91},{5,36,0,91,32},{5,36,91,0,32},{5,91,0,32,36},{32,0,5,36,91},{32,0,5,91,36},{32,0,36,5,91},{32,0,36,91,5},{32,0,91,5,36},{32,5,0,36,91},{32,5,36,0,91},{32,36,0,5,91},{32,36,0,91,5},{32,36,91,0,5},{32,91,0,5,36},{36,0,5,32,91},{36,0,5,91,32},{36,0,32,5,91},{36,0,32,91,5},{36,0,91,5,32},{36,91,0,5,32},{91,0,5,32,36},{0,6,64,65,129},{0,6,64,129,65},{0,6,65,64,129},{0,6,65,129,64},{0,6,129,64,65},{0,64,6,65,129},{0,64,6,129,65},{0,64,65,6,129},{0,64,65,129,6},{0,64,129,6,65},{0,64,129,65,6},{0,65,6,64,129},{0,65,6,129,64},{0,65,129,6,64},{0,129,6,64,65},{0,129,64,6,65},{0,129,64,65,6},{6,0,64,65,129},{6,0,64,129,65},{6,0,65,64,129},{6,0,65,129,64},{6,0,129,64,65},{6,64,0,65,129},{6,64,0,129,65},{6,64,65,0,129},{6,64,65,129,0},{6,64,129,0,65},{6,65,0,64,129},{6,65,0,129,64},{6,65,129,0,64},{6,129,0,64,65},{64,0,6,65,129},{64,0,6,129,65},{64,0,65,6,129},{64,0,65,129,6},{64,0,129,6,65},{64,6,0,65,129},{64,6,65,0,129},{64,65,0,6,129},{64,129,0,6,65},{65,0,6,64,129},{129,0,6,64,65},{0,14,87,102,110},{0,14,87,110,102},{0,14,102,87,110},{0,14,102,110,87},{0,14,110,87,102},{0,87,14,102,110},{0,87,102,14,110},{0,87,102,110,14},{0,102,14,87,110},{0,110,14,87,102},{14,0,87,102,110},{87,0,14,102,110},{102,0,14,87,110},{0,41,78,115,116},{0,43,55,121,124},{0,17,34,72,119},{0,2,36,89,110},{0,3,42,48,116},{0,3,13,15,110},{0,1,12,44,46},{0,5,39,40,42},{0,27,58,80,93}};


viewHelp SchurRings
restart
loadPackage "SchurRings"
A1 =schurRing(a,2) 
A2 =schurRing(A1,b,2) 
A3 =schurRing(A2,c,2) 
A4 =schurRing(A3,d,2) 
A5 =schurRing(A4,e,2) 
r1 = a_1*b_1*c_1*d_1*e_1
symmetricPower_2 r1
symmetricPower_4 r1
last symmetricPower_4 r1
last symmetricPower_6 r1
--
Q = schurRing(q,5,GroupActing => "Sn") -- this is the Sn rep
T = schurRing(Q,t,2) -- this is the GL_2 rep
symmetricPower_2 (t_1)^5

S = schurRing(T,s,2) -- this is the GL_1 rep
r1 = q_5 + q_(4,1) + q_(3,2)
r2 = symmetricPower_2 r1
r1*t_1
symmetricPower_2 t_1
-- there's not any function that I know of yet that expands the potential Hilbert series.
-- this function eats a rep r, and spits out the hilbert series expansion of (1/(1-r)) to order ord:
hilbExpand = (rep,ord) -> sum(ord+1, i-> (symmetricPower_i s_1)*symmetricPower(i,rep))
hilbExpand(q_5, 4)
hilbExpand(q_(4,1), 4)
rep = q_5*t_5 + q_(4,1)*t_(4,1) + q_(3,2)*t_(3,2) 
tmp = hilbExpand(q_5*t_5, 2)* hilbExpand(q_(4,1)*t_(4,1), 2) * hilbExpand(q_(3,2)*t_(3,2) , 2);
diff(s_1, tmp)

hilbExpand(rep, 2)
symmetricPower(2,rep)

symmetricPower(2, q_{2,1,1,1})
symmetricPower(2, q_{3,2})
u = q_{5}
r4 = q_{5} + q_{4,1}
symmetricPower(2,q_{5} + q_{4,1})
symmetricPower(3,q_{5} + q_{4,1})
dim oo
r6 = q_{1,1,1,1,1}
r4*r6
r8 = 4*q_{5} + 3*q_{4,1} + 3*q_{3,2} + q_{2,2,1}
mg8 = 2*u + q_{4,1} + 2*q_{3,2} + q_{2,2,1}
symmetricPower(2, r4)
mg8 + symmetricPower(2, r4) -r8

mg10 = q_{2,2,1} + q_{2,1,1,1} + q_{1,1,1,1,1}
r10 = r4*r6 + mg10
-- non-mingens for deg 12:
r6*r6
symmetricPower(3,r4)
r4*mg8
toString(r6*r6 + symmetricPower(3,r4) + r4*mg8)
nmg12 = r6*r6 + symmetricPower(3,r4) + r4*mg8
dim(r6*r6 + symmetricPower(3,r4) + r4*mg8)
--- notes 
rep = r4 + r6 + q_{3,2} + q_{1,1,1,1,1} + q_{2,2,1}
dim rep
-- declare the ring in even degrees, and don't forget the 0 in degree 2
sr2 = schurResolution(rep, {q_5,0,r4,r6,r8,r10,nmg12 + q_{2,2,1}}, DegreeLimit=>10)
netList sr2
apply(sr2, tm-> dim tm#0#1)
sr =  schurResolution(rep, {q_5}, DegreeLimit=>17)
apply(sr, tm-> dim tm#0#1)


sr =  schurResolution(5*q_{5} + 3*q_{4,1}, {q_5}, DegreeLimit=>17)
apply(sr, tm-> dim tm#0#1)


apply(17, i-> binomial(17,i))
netList schurResolution(rep, {q_5}, DegreeLimit=>5)
sr =  schurResolution(rep, {q_5}, DegreeLimit=>5)
apply(sr, tm-> dim tm#0#1)


--time betti mingens (I4 + ideal F8) -- takes about 30s
-- note that an 8-dimensional subspace is not mingens.
-- time betti mingens ideal F8big -- takes about 29s
-- rankList
-- rankList_jumps

-- these iterative methods take 2-4x longer, maybe because of all the re-writes to the matrix?
-- -- doing one part of the matrix at a time can lead to mistakes because the rank over ZZ/p might fail when the det is non-zero over ZZ but reduces to 0 mod p. 
-- M = new MutableMatrix from (0*id_(KK^(#myPerms)))
-- F8 ={inv8};rk =1;
-- M_(0,0) = pts#0 inv8;
-- rankList = {1};
-- elapsedTime for P in myPerms_{1..#myPerms-2} do (
-- rankList = rankList|{rk};
-- tmp = (perm P) inv8; tmpF8 = F8|{tmp};
-- for i to rk do M_(i,rk) = (pts#i ) (tmpF8#rk ); -- just do the outer band?
-- for j to rk do M_(rk,j) = (pts#rk ) (tmpF8#j ); -- just do the outer band?
-- if rank (matrix M)>rk then (F8 = tmpF8; rk=rk + 1;) else continue);


-- M = new MutableMatrix from (0*id_(KK^(#myPerms)))
-- F8 ={inv8};rk =1;
-- pmList = {myPerms#0};
-- M_(0,0) = pts#0 inv8;
-- rankList = {1};
-- elapsedTime for P in myPerms_{1..#myPerms-2} do (
-- rankList = rankList|{rk};
-- tmp = (perm P) inv8; tmppmList = pmList|{P};
-- for i to rk do M_(i,rk) = ((pts#i )*(perm(tmppmList#rk))) inv8 ; -- just do the outer band?
-- for j to rk do M_(rk,j) = ((pts#rk )*(perm(tmppmList#j))) inv8 ; -- just do the outer band?
-- if rank (matrix M)>rk then (pmList = tmppmList; rk=rk + 1;) else continue);
-- matrix M
-- ((pts#1)*(perm(tmppmList#1)) )inv8
-- rankList 
-- rank matrix M
-- Mb_{0..7}^{0..7} - matrix(M)_{0..7}^{0..7}
-- matrix M
-- M
-- #F8
-- takes 2s instead of 125s each iteration to apply permutations to polynomials rather than use the YS again and again
--makeMg = myPerm -> sub(matrix apply(pts,pt-> apply(goodTabList8, Q -> tab2polyMap(apply(5, i-> syt8#(Q#i)), pt*myPerm))), baseRing X);
--time Mg = sub(matrix apply(pts,pt-> apply(goodTabList8, Q -> tab2polyMap(apply(5, i-> syt8#(Q#i)), pt*myPerm))), baseRing X) ; -- expect 125s
-- M = makeMg (perm {0,1,2,3,4});
-- Mp = makeMg (perm {1,0,2,3,4});
-- M^(-1)*Mp
-- trace(M^(-1)*Mp)

-- first way to find conjiugacy classes
-- G = permutations(5);
-- con = h-> sort unique apply(G, g-> ((inversePermutation(g))_h)_g)
-- conjClasses = (unique (con\G));
-- conjClasses = conjClasses_{0,1,2,4,6,3,5}; -- this is the order in FH p28
-- orders = apply(conjClasses, c-> #c) -- sizes of conjugacy classes.
-- reps = apply(conjClasses, cl -> first cl)
