-- SL_2(Z_3) is a group, what is its character table?
restart
-- get all the elements by brute force:
R = ZZ/3
SL2Z3 = flatten flatten flatten for a to 2 list for b to 2 list for c to 2 list for d to 2 list (tmpA = sub(matrix{{a,b},{c,d}},R); if det tmpA ==1 then tmpA else continue)
#SL2Z3 -- it's a group with 24 elements
-- check closure:
# unique(SL2Z3*SL2Z3)
-- build the conjugacy classes:
remaining = SL2Z3;
CL = for z in SL2Z3 list (if member(z,remaining) then (
    tmp = unique for x in SL2Z3 list x*z*x^(-1); 
    for xx in tmp do remaining = delete(xx,remaining);
    tmp) else continue )
-- this is all the conjugacy classes, but there's some redundancy, and M2 doesn't know how to sort matrices
netList CL
-- # unique apply(CL, L -> sort apply(L, l -> flatten entries l))
-- FL =  unique apply(CL, L -> sort apply(L, l -> flatten entries l))
-- conjClasses = for L in FL list for x in L list reshape(R^2, R^2, matrix{x})
-- netList conjClasses

-- so there are 7 conjugacy classes - I put them in order by # elements and by order of elements
CL = CL_{5,6,4,2,3,1,0}
-- these are the sizes of conjugacy classes:
apply(CL, l-> #l)
netList CL
-- this is the verification of the dimension formula - there's only one solution.
24 - (1+1+1+4+4+4+9)
apply(6, i->(first CL#0)^(i+1))
apply(6, i->(first CL#1)^(i+1))
apply(6, i->(first CL#2)^(i+1))
apply(6, i->(first CL#3)^(i+1))
apply(6, i->(first CL#4)^(i+1))
apply(6, i->(first CL#5)^(i+1))
apply(6, i->(first CL#6)^(i+1))
CL#6
 unique flatten ((CL#2)*(CL#3))
 netList unique for x in CL#2 list x*CL#2
 CL#3
 unique flatten ((CL#1)*(CL#2))
 CL#5
 findOrder = a -> (tmp:=1; min for i from 1 to 6 list if a^i == a^0 then i else continue)
 findOrder\ CL#2
 findOrder\ CL#5
 findOrder\ CL#6
-- this points out that the conjugacy classes don't multiply
findOrder\ unique flatten (CL#2 * CL#3)
apply(unique flatten (CL#2 * CL#3), xx-> member(xx, CL#6))
apply(unique flatten (CL#2 * CL#3), xx-> member(xx, CL#0))
-- however, if y is in class [x] then y^2 is in class [x^2] 
-- this is because if y = hxh^(-1), then y^2 = h x^2 h^(-1) by direct computation
-- we can check if any equivalence class has an element and its inverse in the same class
for x in CL#6 do for y in CL#6 do if x*y ==(x*y)^0 then print(x,y)
for x in CL#5 do for y in CL#5 do if x*y ==(x*y)^0 then print(x,y)
for x in CL#4 do for y in CL#4 do if x*y ==(x*y)^0 then print(x,y)
for x in CL#3 do for y in CL#3 do if x*y ==(x*y)^0 then print(x,y)
for x in CL#2 do for y in CL#2 do if x*y ==(x*y)^0 then print(x,y)
-- checks if the equivalence class has an element of order 2
for i to 6 do for x in CL#i do if x*x ==(x*x)^0 then print(i,x)

-- character table so far:
--restart
S = QQ[a,b,c,d,e,f, w]/(w^2 + w +1 )
conj = map(S,S, apply({a,b,c,d,e,f}, xx-> xx=>xx)| {w=> w^2})
A = matrix {apply(7, i->1), {1,1,w,w^2,w^2, w, 1}, {1,1,w^2,w,w, w^2, 1},
{2, a,b,c,d,e,f}
}
hadamardProduct = (u,v) -> apply(#u, i-> u_i *v_i)
wts = {1,1,4,4,4,4,6}
ip = (x,y) -> 1/24*sum(7, i-> conj(x_i)*y_i*wts_i)
A = A||matrix {hadamardProduct(flatten(entries A^{1}), flatten entries( A^{3}))}||matrix {hadamardProduct(flatten(entries A^{2}), flatten entries( A^{3}))}
-- the orthogonality conditions
ip(flatten entries A^{0}, flatten entries A^{0}) -- = 1
ip(flatten entries A^{0}, flatten entries A^{1}) -- = 0 
ip(flatten entries A^{0}, flatten entries A^{2}) -- = 0 
ip(flatten entries A^{0}, flatten entries A^{3}) --?
ip(flatten entries A^{0}, flatten entries A^{4}) --?
ip(flatten entries A^{0}, flatten entries A^{5}) --?
ip(flatten entries A^{1}, flatten entries A^{1}) -- = 1 
ip(flatten entries A^{1}, flatten entries A^{2}) -- = 0 
ip(flatten entries A^{1}, flatten entries A^{3}) --?
ip(flatten entries A^{1}, flatten entries A^{4}) --?
ip(flatten entries A^{1}, flatten entries A^{5}) --?
ip(flatten entries A^{2}, flatten entries A^{2}) -- = 1
ip(flatten entries A^{2}, flatten entries A^{3}) --?
ip(flatten entries A^{2}, flatten entries A^{4}) --?
ip(flatten entries A^{2}, flatten entries A^{5}) --?
ip(flatten entries A^{3}, flatten entries A^{3}) --?
ip(flatten entries A^{3}, flatten entries A^{4}) --?
ip(flatten entries A^{3}, flatten entries A^{5}) --?
ip(flatten entries A^{4}, flatten entries A^{4}) --?
ip(flatten entries A^{4}, flatten entries A^{5}) --?
ip(flatten entries A^{5}, flatten entries A^{5}) --?
kdelta = (i,j) -> (if i==j then 1 else 0) 
I = ideal matrix apply(6, i -> apply(6, j-> ip(flatten entries A^{i}, flatten entries A^{j}) - kdelta(i,j) ) ) 
netList flatten entries mingens I
-- 
chiV2 = matrix{ hadamardProduct(flatten entries (A^{3}), flatten entries (A^{3}))}
-- this list computes the landing place for g^2 for the first g in the class
--g2List =  for j to 6 list flatten for mem in CL#j list for i to 6 list if member((mem)^2, CL#i) then i else 
g2List =  flatten for j to 6 list for i to 6 list if member((first CL#j)^2, CL#i) then i else continue
for xx in CL#6 list(member(xx^2,CL#6 ))
for xx in unique flatten((CL#6)*(CL#6)) list(member(xx,CL#6 ))
for xx in CL#5 list(member(xx^2,CL#5 ))
chiVg2 = (A^{3})_g2List
-- notes S^2 V is 3-dimensional, so it's either irreducible (there's only one irrep of that dimension), or it's a sum either of type 1+1+1 or of type 2+1
chiS2V = 1/2*(chiV2 + chiVg2)
-- ip(chiS2V , chiS2V) =? if 1, then it's irreducible
ip(flatten entries chiS2V ,flatten entries  chiS2V)
factor oo 
-- note wedge^2 V is 1-dimensional, so it's irreducible. 
chiW2V = 1/2*(chiV2 - chiVg2)
ip(flatten entries chiW2V ,flatten entries  chiW2V)
-- if this is zero, then they're both irreducible:
ip(flatten entries chiS2V ,flatten entries  chiS2V) - ip(flatten entries chiW2V ,flatten entries  chiW2V)

mingens ideal (chiW2V - A^{0})
mingens ideal (chiW2V - A^{1})
mingens ideal (chiW2V - A^{2})
sub(sub(chiW2V, S/(a^2-4)), S)
sub(sub(chiS2V, S/(a^2-4)), S)
-- we don't know which 1-d irrep it is, but they all say that a^2 -4 = 0, so a^2 = 4, and a = +/- 2
sub(chiS2V, a=>-2)
sub(chiS2V, a=>2)

sub(chiW2V, a=>-2)
sub(chiW2V, a=>2)

-- this isn't enough unless you know more, like the 3-dimensional rep lifted from a quotient G/N


-- try the regular representation?

-- Reg Rep is permutation rep on the 24 elements

apply(CL#2, xx-> xx^2)
CL#3
-- is it the case that the group is generated by 2 elements?
myGens = {CL#1#0, last CL#2}
myGens = CL#1|CL#2
soFar = myGens; for i to 10 do (soFar = unique flatten (soFar*myGens); print(#soFar))
soFar
apply(CL#2,x ->  member((x)^2, CL#3))
apply(CL#3,x ->  member((x)^2, CL#2))
apply(CL# 2,x ->  member((x)^3, CL#0))
A
isq = {0,0,3,2, 3,2,1}
chiw2v = i -> 1/2*(A_(3,i)^2 - A_(3, isq_i))
chiw2v\{0,1,2,3,4,5,6}
mingens (ideal(chiw2v\{0,1,2,3,4,5,6} - apply(7, i-> 1)) )
