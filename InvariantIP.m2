-- invariantIP
-- make an invariant inner product via point evaluations that works with few points by 
KK = QQ;
R = QQ[x_0..x_4]
perm = P-> map(R,R, apply(#P, i-> x_i => x_(P#i)));
perm({0,1,2,4,3})
f = random(2,R)
g = random(2,R)

rndPt = () -> map(KK, R, random(KK^1, KK^(numgens R))) 
rp = rndPt();

points = apply(2, i-> rndPt())
pointIP = (myPoints, a,b) -> sum(myPoints, pt -> pt(a)*pt(b))
pointIP (f,g)
fp = f - pointIP (f,g)/pointIP (g,g)*g

-- if you use a global variable points:
-- invariantPointIP = (a,b) -> sum(points, pt -> sum(permutations 5, P-> (pt*(perm P))(a) * (pt*(perm P))(b) ))
--if you want to specify what points to use:
invariantPointIP = (myPoints, a,b) -> sum(myPoints, pt -> sum(permutations 5, P-> (pt*(perm P))(a) * (pt*(perm P))(b) ))
invariantPointIP_points(f,g)
--invariantProj = (u,v) -> invariantPointIP_myPoints(u,v)/ invariantPointIP_myPoints(u,u)*u
invariantProj = (myPoints, u,v) -> invariantPointIP_myPoints(u,v)/ invariantPointIP_myPoints(u,u)*u
gp = g - invariantProj_points(f,g) -- the claim is that this polynomial is in the span of f and g, but it is orthogonal to f
invariantPointIP_points(f,gp)
-- this shows that the polynomial we constructed is indeed still in the span of the previous two:
betti mingens ideal {f,g,gp}

-- this works, but it is very slow. It probably makes more sense to make a matrix of values, then GS the matrix, then apply that to the forms rather than compute all these inner products over and over again....
GSprocess = polys-> (
    mypoints :=  apply(#polys, i-> rndPt());
    ulist := new MutableList from apply(#polys, i-> 0_(ring polys#0));
    elist := new MutableList from apply(#polys, i-> 0_(ring polys#0));
    for k from 0 to #polys-1 do(
        ulist#k =  polys#k - sum(k, j-> invariantProj_mypoints(ulist#j, polys#k));
        elist#k = sub(ulist#k/invariantPointIP_mypoints(ulist#k,ulist#k),ring ulist#k);
    );
    toList elist
)

-- each polynomial will get sent to a matrix of points x permutations:
-- ordering is polynomial, point, vector of values of permuted inputs
GSprocessVecs = polys-> (
    mypoints :=  apply(#polys, i-> rndPt());
    vlist := new HashTable from for i to #polys-1 list  i=> new HashTable from apply(#mypoints, j-> j=> matrix apply(permutations 5, P-> {(mypoints#j*(perm P)) polys#i}) );
    vlist
)
 v4v =  GSprocessVecs V4;
 (transpose v4v#0#0)*v4v#0#1
time f448a =  GSprocessVecs F448a; -- approximately 12 s per matrix for the first batch, 25s each for the second... takes about 675s total to make this array.
-- but the inner products should be able to be computed quickly now.

V4onb = GSprocess(V4);
mypts= apply(5, i-> rndPt())
f4onb = U4|V4onb;
-- this shows that even if you change your points for evaluations, the orthogonality still holds:
matrix apply(f4onb, f-> apply(f4onb, g-> invariantPointIP_mypts(f,g)))
