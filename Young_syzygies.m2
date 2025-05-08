---Young Symmetrizers  -- trying to compute syzygies, try the Koszul complex
restart
aa = 2; bb = 2;
list2mats = (L,a) -> matrix apply(#L, i-> apply( L, j-> a_(i,j-1))) -- make matrices of the right size... 
X = QQ[x_(0,0)..x_(aa,bb)] -- ring on A \otimes B
tab2polyList = (T,L)->(    
    ta:=T#0;    tb:=T#1; dg:=length flatten ta;
    R := X[a_(0,0)..a_(#ta#0-1,dg-1),b_(0,0)..b_(#tb#0-1,dg-1)];
    A := apply(ta, Lm -> list2mats(Lm, a));
    B := apply(tb, Lm -> list2mats(Lm, b));
    usedA := {};usedB := {};
    F :=1_R; T :=1_R;
    for z in L do time ( T=product(
        (for i to sub(dg/2-1,ZZ) list if isMember(z+1,ta#i) and not isMember(i,usedA) then (usedA = usedA|{i}; det A_i) else continue)|
        (for i to sub(dg/2-1,ZZ) list if isMember(z+1,tb#i) and not isMember(i,usedB) then (usedB = usedB|{i}; det B_i) else continue));
        F = F*T;
        F= sum(#tb#0-1, j-> sum(#ta#0-1,i->valmap(x_(i,j))*diff(a_(i,z)*b_(j,z),F)));
    ); 
    F
)

tab2polyListEval = (T,L,valmap)->(    
    --T = {{{1,2,3}},{{1,2,3}}}; L = {0}; valmap= val;
    ta:=T#0;    tb:=T#1; dg:=length flatten ta; 
    R := X[a_(0,0)..a_(#ta#0-1,dg-1),b_(0,0)..b_(#tb#0-1,dg-1)]; 
    A := apply(ta, Lm -> list2mats(Lm, a)); 
    B := apply(tb, Lm -> list2mats(Lm, b)); 
    usedA := {};usedB := {};
    F :=1_R; T :=1_R;
    for z in L do time ( T=product(
        (for i to sub(dg/2-1,ZZ) list if isMember(z+1,ta#i) and not isMember(i,usedA) then (usedA = usedA|{i}; det A_i) else continue)|
        (for i to sub(dg/2-1,ZZ) list if isMember(z+1,tb#i) and not isMember(i,usedB) then (usedB = usedB|{i}; det B_i) else continue));
        F = F*T; 
        F= sum(#tb#0-1, j-> sum(#ta#0-1,i->valmap(x_(i,j))*diff(a_(i,z)*b_(j,z),F)));
    ); 
    F
)

tab2polyListEval2 = (T,L,M,rules)->(    
    ta:=T#0;    tb:=T#1; dg:=length flatten ta;
    R := X[a_(0,0)..a_(#ta#0-1,dg-1),b_(0,0)..b_(#tb#0-1,dg-1)];
    A := apply(ta, Lm -> list2mats(Lm, a)); 
    B := apply(tb, Lm -> list2mats(Lm, b)); 
    usedA := {};usedB := {};
    F :=1_R; T :=1_R;
    for z in L do time ( T=product(
        (for i to sub(dg/2-1,ZZ) list if isMember(z+1,ta#i) and not isMember(i,usedA) then (usedA = usedA|{i}; det A_i) else continue)|
        (for i to sub(dg/2-1,ZZ) list if isMember(z+1,tb#i) and not isMember(i,usedB) then (usedB = usedB|{i}; det B_i) else continue));
        F = F*T;
        F= sum(#tb#0-1, j-> sum(#ta#0-1,i->sub(x_(i,j),rules)*diff(a_(i,z)*b_(j,z),F)));
    );
    for z in M do time ( T=product(
        (for i to sub(dg/2-1,ZZ) list if isMember(z+1,ta#i) and not isMember(i,usedA) then (usedA = usedA|{i}; det A_i) else continue)|
        (for i to sub(dg/2-1,ZZ) list if isMember(z+1,tb#i) and not isMember(i,usedB) then (usedB = usedB|{i}; det B_i) else continue));
        F = F*T;
        F= sum(#tb#0-1, j-> sum(#ta#0-1,i-> x_(i,j)*diff(a_(i,z)*b_(j,z),F)));
    ); 
    F
)

-- suppose we already know that the bi-partition {{1},{1}} is the first module M, so we compute the second one by computing all the possible syzygy modules in the tensor product M ** R
point2eval = pt -> map(X,X, apply(gens X, xx-> xx=> contract(xx,pt)))
point2rules = pt -> apply(gens X, xx-> xx=> contract(xx,pt))
val = point2eval x_(0,0)
myRules = point2rules x_(0,0)
-- two modules show up:
tab2polyListEval({{{1, 2}},{{1,2}}}, {0,1}, val)
tab2polyListEval({{{1, 2}},{{1,2}}}, {0}, val)

tab2polyListEval({{{1},{2}},{{1},{2}}}, {0,1}, val)
tab2polyListEval({{{1, 2}},{{1,2}}}, {1}, val)

-- the idea here is to evaluate part of the tableau and then put variables in for the rest. 
tab2polyListEval2({{{1},{1}},{{1},{1}}}, {0,1},{}, val)
tab2polyListEval2({{{1, 2}},{{1,2}}}, {0},{1}, val)


tab2polyListEval({{{1,2,3}},{{1,2,3}}}, {0}, val)
tab2polyListEval2({{{1,2,3}},{{1,2,3}}}, {0},{1}, myRules)