-- Young symmetrizers for Specht modules

--KK = ZZ/nextPrime(100); --KK = QQ;
-- KK = QQ;
--X = KK[x_(0,0,0,0,0)..x_(1,1,1,1,1)];
--perm = P -> map(X,X, apply(toList((0,0,0,0,0)..(1,1,1,1,1)), ll -> x_(toSequence apply(P, p-> ll#p))));
-- try the 3-factor case:
-- KK = QQ;
-- X = KK[x_(0,0,0)..x_(1,1,1)];
-- perm = P -> map(X,X, apply(toList((0,0,0)..(1,1,1)), ll -> x_(toSequence apply(P, p-> ll#p))));
-- need a way to make a sub-permutation, only using some of the letters, and apply that to a list

-- this function takes P as a permutation on its own content, and applies that permutation to the content of L
subSym = (P, L) -> (Lm := new MutableList from L; 
for i to #P-1 do Lm#((sort P)#i) = L#(P#i);
toList Lm)


-- given a filling, produce the Young symmetrizer for the Specht module
--- input tableau, output the map that symmetrizes the content of the rows:
-- want to perform a sequence of operations to avoid this map getting to large?
-- however, the sum of two ring maps is not the ring map of the sum of the matrix of the ring maps... 
-- you can see this by acting on anything other than degree <2
-- f = sum(10, j-> product(2, i-> (gens X)_(random(32))))
-- m1 = perm((permutations(5))_(random(120)))
-- m2 = perm((permutations(5))_(random(120)))
-- f12 = (m1 f) + (m2 f)
-- f12p = (map(X,X, matrix m1 + matrix m2) )f
-- f12p-f12

-- function to transpose a tableau
tabTranspose = tab -> (
tmpTab := new MutableList from tab;
toList for j to #first tmpTab-1 list( 
    tmpTab = delete({},tmpTab);
    for i to #tmpTab-1 list( c := first tmpTab#i; tmpTab#i = delete(first tmpTab#i, tmpTab#i ); c)
)
)

-- this way of computing the sign of a permutation works for any list of numbers, even if it's not 0..d:
sgn = p-> det((id_(KK^(max(p)+1)))_p^(sort p))

-- alambda needs to be the sum of the result of the maps applied to the input
-- if you sum the maps, you can't just sum their matrices that define the map in degree 1
alambda = (tab, pol)->(
    tmp := pol;ftab := flatten tab;
    for i to #tab-1 do( 
        plist := apply(permutations(tab#i), x-> subSym(x, sort ftab));
        tmp = sum(plist, p -> (perm p) tmp); 
    );
tmp)
-- you could put in a factor like 1/((#(tab#i))!)* this so that you project polynomials to themselves
--  however, that can cause problems in other characteristics, and it's not the true definition of the map
--- input tableau, output the map that skew-symmetrizes the content of the columns:
-- blambda on the tableau (need to transpose)
blambda = (tmptab, pol)->(
    tmp := pol;
    tab := tabTranspose tmptab;
    ftab := flatten tab;
    for i to #tab -1 do( 
        plist := apply(permutations(tab#i), x-> subSym(x, sort ftab));
        tmp = sum(plist, p -> sgn(p)*((perm p) tmp)) 
    );
tmp)

clambda = (tab, pol) -> alambda_tab blambda_tab pol
end;
restart
input"/Users/lao0004/Documents/GitHub/YoungSymmetrizers/YoungSpecht.m2"
-- try the 2-factor case:
KK = QQ;
X = KK[x_(0,0)..x_(1,1)];
perm = P -> map(X,X, apply(toList((0,0)..(1,1)), ll -> x_(toSequence apply(P, p-> ll#p))));
gx = gens X
myTab = {{0,1}}; rank  diff(transpose matrix {gx}, matrix({clambda_myTab \ gens X}))
myTab = {{1,0}}; rank  diff(transpose matrix {gx}, matrix({clambda_myTab \ gens X}))
myTab = {{0},{1}}; rank diff(transpose matrix {gx}, matrix({clambda_myTab \ gens X}))


myTab = {{0,1}}
I = ideal(apply(5, i->  random(4, X)));
betti(ideal (clambda_myTab\ flatten entries gens I))

betti(alambda_myTab I)
betti(blambda_myTab I)
-- 3 factor examples:
KK = QQ;
X = KK[x_(0,0,0)..x_(1,1,1)];
perm = P -> map(X,X, apply(toList((0,0,0)..(1,1,1)), ll -> x_(toSequence apply(P, p-> ll#p))));
gx = gens X
myTab = {{0,1,2}}; rank diff(transpose matrix {gx}, matrix({clambda_myTab \ gens X}))
myTab = {{0,1},{2}}; rank diff(transpose matrix {gx}, matrix({alambda_myTab \ gens X}))
myTab = {{0,1},{2}}; rank diff(transpose matrix {gx}, matrix({blambda_myTab \ gens X}))
myTab = {{0,1},{2}}; rank diff(transpose matrix {gx}, matrix({clambda_myTab \ gens X}))
myTab = {{0,2},{1}}; rank diff(transpose matrix {gx}, matrix({clambda_myTab \ gens X}))
myTab = {{0},{1},{2}}; rank diff(transpose matrix {gx}, matrix({clambda_myTab \ gens X}))


-- 
-- try the 4-factor case:
KK = QQ;
X = KK[x_(0,0,0,0)..x_(1,1,1,1)];
perm = P -> map(X,X, apply(toList((0,0,0,0)..(1,1,1,1)), ll -> x_(toSequence apply(P, p-> ll#p))));
gx = gens X
alambda_{{0,1},{2,3}} x_(0,1,0,1)
myTab = {{0,1,2,3}}; rank diff(transpose matrix {gx}, matrix({clambda_myTab \ gens X}))
myTab = {{0,1,2},{3}}; rank diff(transpose matrix {gx}, matrix({clambda_myTab \ gens X}))

myTab = {{0,1},{2,3}}; rank diff(transpose matrix {gx}, matrix({alambda_myTab \ gens X}))
myTab = {{0,1},{2,3}}; rank diff(transpose matrix {gx}, matrix({blambda_myTab \ gens X}))
myTab = {{0,1},{2,3}}; rank diff(transpose matrix {gx}, matrix({clambda_myTab \ gens X}))


myTab = {{0,2},{1,3}}; rank diff(transpose matrix {gx}, matrix({alambda_myTab \ gens X})) 
myTab = {{0,2},{1,3}}; rank diff(transpose matrix {gx}, matrix({blambda_myTab \ gens X}))
myTab = {{0,2},{1,3}}; rank diff(transpose matrix {gx}, matrix({clambda_myTab \ gens X}))

myTab = {{0,1},{2},{3}}; rank diff(transpose matrix {gx}, matrix({clambda_myTab \ gens X}))
myTab = {{0},{1},{2},{3}}; rank diff(transpose matrix {gx}, matrix({clambda_myTab \ gens X}))


-- 5-factor examples:
KK = QQ;
X = KK[x_(0,0,0,0,0)..x_(1,1,1,1,1)];
perm = P -> map(X,X, apply(toList((0,0,0,0,0)..(1,1,1,1,1)), ll -> x_(toSequence apply(P, p-> ll#p))));
gx = gens X
myTab = {{0,1,2},{3,4}}

f = sum(10, j-> product(2, i-> (gens X)_(random(32))))

fp = alambda_myTab f
-- check that we actually have a projection:
fp - alambda_myTab fp -- might not be zero because you don't know the multiplication factor
rank matrix{{fp , alambda_myTab fp}} -- has rank 1 because it is a projection onto the linear space:

fp = blambda_myTab f
fp - blambda_myTab fp -- might not be zero because you don't know the multiplication factor
rank matrix{{fp , blambda_myTab fp}} -- has rank 1 because it is a projection onto the linear space:

-- see what happens to some of the variables
alambda_myTab x_(0,0,0,0,1)
alambda_myTab x_(0,1,0,0,0)
alambda_myTab x_(0,1,0,0,1)


blambda_myTab (f)
blambda_myTab (x_(0,1,0,0,1)*x_(0,0,1,1,1))
blambda_myTab x_(0,1,0,0,1)
-- notice that if the tableau has more than two rows, then blambda sends all variables to zero:
blambda_myTab \ gens X
-- notice that if the tableau has more than two rows, then blambda sends all variables to zero:
blambda_(tabTranspose myTab) \ gens X

blambda_(tabTranspose myTab)  (f)

-- compute the dimensions of the projections:
myTab = {{0,1,2,3,4}}; rank diff(transpose matrix {gx}, matrix({clambda_myTab \ gx})) -- 6
myTab = {{0,1,2,3},{4}}; rank diff(transpose matrix {gx}, matrix({clambda_myTab \ gx})) -- 4
myTab = {{0,1,2},{3,4}}; rank diff(transpose matrix {gx}, matrix({clambda_myTab \ gx})) --2
myTab = {{0,1,2},{3,4}}; rank diff(transpose matrix {gx}, matrix({alambda_myTab \ gx})) -- 4*3 = 12
myTab = {{0,1,2},{3,4}}; rank diff(transpose matrix {gx}, matrix({blambda_myTab \ gx})) -- 1*1*2 = 2
myTab = {{0,1},{2,3},{4}}; rank diff(transpose matrix {gx}, matrix({clambda_myTab \ gx})) 

-- this computes the decomposition of S^2 (C^2)^5 as a sum of images of these Young symmetrizers, but we only get one per standard filling, so to get 528, you have to multiply by the dimensions of the S_5 modules
-- this is the decomposition as an S_5 module.
-- 34+45*4+33*5+15*6+11*5+1*4 = 528
-- the irreducible decomposition for SL_2^5 goes like this:
-- 243+10*27+5*3 = 528

gx2 = flatten entries basis(2, X);
myTab = {{0,1,2,3,4}}; time betti mingens ideal(clambda_myTab \gx2) -- 34
myTab = {{0,1,2,3},{4}}; time rank diff(transpose matrix {gx2}, matrix({clambda_myTab \ gx2})) -- 45
myTab = {{0,1,2},{3,4}}; time rank diff(transpose matrix {gx2}, matrix({clambda_myTab \ gx2})) -- 33
myTab = {{0,1,2},{3},{4}}; time rank diff(transpose matrix {gx2}, matrix({clambda_myTab \ gx2})) -- 15
myTab = {{0,1},{2,3},{4}}; time rank diff(transpose matrix {gx2}, matrix({clambda_myTab \ gx2})) -- 11
myTab = {{0,1},{2},{3},{4}}; time rank diff(transpose matrix {gx2}, matrix({clambda_myTab \ gx2})) -- 1
myTab = {{0},{1},{2},{3},{4}}; time rank diff(transpose matrix {gx2}, matrix({clambda_myTab \ gx2})) -- 0


myTab = {{0,1},{2},{3},{4}};  tmp = ideal mingens ideal (clambda_myTab\ gx2); betti mingens tmp
myTab = {{0,2},{1},{3},{4}};  tmp = tmp + ideal mingens ideal (clambda_myTab\ gx2); betti mingens tmp
myTab = {{0,3},{2},{1},{4}};  tmp = tmp + ideal mingens ideal (clambda_myTab\ gx2); betti mingens tmp
myTab = {{0,4},{2},{3},{1}};  tmp = tmp + ideal mingens ideal (clambda_myTab\ gx2); betti mingens tmp
betti mingens tmp
34+45+33+11+1


myTab = {{0,1,2,3},{4}}; M =  diff(transpose matrix {gx}, matrix({clambda_myTab \ gx}))
mingens image M

bx = {x_(1,0,0,0,0), x_(0,1,0,0,0), x_(0,0,1,0,0), x_(0,0,0,1,0), x_(0,0,0,0,1)}
myTab = {{0,3,1,4},{2}}; M =  diff(transpose matrix {gx}, matrix({clambda_myTab \ bx}))
codim ideal(clambda_myTab \gx)
netList transpose {clambda_myTab \bx}
mingens image M

myTab = {{0,1,2,3,4}}; M =  diff(transpose matrix {gx}, matrix({clambda_myTab \ bx}))
M =  diff(transpose matrix {gx}, matrix({clambda_myTab \ gx}))
ker M
{{0,1,2,3,4},{0,1,3,2,4},{0,2,1,3,4},{0,2,3,1,4},{0,3,1,2,4},{0,3,2,1,4},{1,0,2,3,4},{1,0,3,2,4},{1,2,0,3,4},{1,2,3,0,4},{1,3,0,2,4},{1,3,2,0,4},{2,0,1,3,4},{2,0,3,1,4},{2,1,0,3,4},{2,1,3,0,4},{2,3,0,1,4},{2,3,1,0,4},{3,0,1,2,4},{3,0,2,1,4},{3,1,0,2,4},{3,1,2,0,4},{3,2,0,1,4},{3,2,1,0,4}}  
{{0,1,2,3,4}}  

-- 6 factor examples:
KK = QQ;
X = KK[x_(0,0,0,0,0,0)..x_(1,1,1,1,1,1)];
perm = P -> map(X,X, apply(toList((0,0,0,0,0,0)..(1,1,1,1,1,1)), ll -> x_(toSequence apply(P, p-> ll#p))));
gx = gens X;

myTab = {{0,1,2,3,4,5}}; rank diff(transpose matrix {gx}, matrix({clambda_myTab \ gx})) -- 7
myTab = {{0,1,2,3,4},{5}}; rank diff(transpose matrix {gx}, matrix({clambda_myTab \ gx})) -- 5
myTab = {{0,1,2,3},{4,5}}; rank diff(transpose matrix {gx}, matrix({clambda_myTab \ gx})) --3
myTab = {{0,1,2},{3,4,5}}; rank diff(transpose matrix {gx}, matrix({clambda_myTab \ gx})) --1
myTab = {{0,1},{2,3},{4,5}}; rank diff(transpose matrix {gx}, matrix({clambda_myTab \ gx})) --0


subSym = (P, L) -> (Lm := new MutableList from L; 
for i to #P-1 do Lm#((sort P)#i) = L#(P#i);
toList Lm)

subSym_{4,0}\ apply(permutations({0,1,2,3}), xx-> {4}|xx)
subSym_{4,0} {0,1,2,3,4}
subSym_{4,0} {1,2,3,4,0}

apply(permutations({0,1,2,3}), xx-> {4}|xx)
