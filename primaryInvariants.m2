-- finding primary invariants for 5 qubits
--restart
input"5factor.m2"
---- try to find the minimal degree algebraically independent generators of the ring of invariants:
gx = gens X
nnz = 22; -- what happens when nothing is zeroed out?
vl =  unique apply(32-nnz, i-> random(32)); -- pick some variables to zero
while #vl !=32-nnz do vl = unique(vl|apply(32-nnz-#vl, i-> random(32)))
randomLinear = apply(vl, i-> gx_i=>0); --random dense points takes longer than just evaluation. Try random sparse...
bx = sub(basis(1,X),randomLinear);
jc = I-> diff(transpose bx, gens I );
jacEval = (fns, rnp) -> (
    bxn1 := reverse subsets( gx, #gx-1);
    matrix for i to 31 list (
    rules = apply(bxn1#i, xx-> xx=>rnp xx);
    fns0 = matrix {apply(fns,fn -> sub(fn, rules ))};
    flatten entries rnp diff(gx#i,  fns0 ))
)

L4 = {syt4#0,syt4#0,syt4#0,syt4#0,syt4#1};
inv4p1 = tab2polyEval(L4,randomLinear); 
inv4p2 = tab2polyEval(L4_{0,1,2,4,3},randomLinear); 
inv4p3 = tab2polyEval(L4_{0,1,4,2,3},randomLinear); 
inv4p4 = tab2polyEval(L4_{0,4,2,2,3},randomLinear); 
inv4p5 = tab2polyEval(L4_{4,0,1,2,3},randomLinear); 
rp = rndPt();
li = ideal(inv4p1,inv4p2, inv4p3, inv4p4, inv4p5 );
rank rp jacobian li, betti mingens li

inv6 = tab2polyEval(syt6,randomLinear);
li = li = ideal(inv4p1,inv4p2, inv4p3, inv4p4, inv4p5, inv6);
rank rp jc li
betti mingens li

time I8 = ideal apply(mgTabList8_{0..4}, pp-> tab2polyEval(apply(5, i-> syt8#(pp#i)), randomLinear));
I = li + I8;
betti I
rank rp jc I
time I10 = ideal apply(goodTabList10_{0}, pp-> tab2polyEval(apply(5, i-> syt10#(pp#i)), randomLinear));
I = I+ I10;
Ito10 = I;
betti I
rank rp jc I
I12 = {};
for L in minList12_{0..4} do (elapsedTime I12 =  I12 |{ tab2polyEval(apply(L, l-> syt12#l),randomLinear) };
I = Ito10+ ideal I12;
print(betti I);
--time print(rank jacEval(flatten entries gens I, rp));
time print(rank rp jc I);
)
I12 = ideal I12;
I = I+ I12;
betti I
rank jacEval(flatten entries gens I, rp)
rank rp jc I


end;
restart
input"primaryInvariants.m2"
I;
--I = ideal {inv4p1,inv4p2, inv4p3, inv4p4, inv4p5, inv6}+  I8 +  I10;
betti I
rank rp jc I
-- try to find the 5,1,5,1,5 pattern:
gns  = (flatten entries gens I)_({0,1,2,3,4}|{5}|{6,7,8,9,10}});
degree\ gns
for i from 1 to #gns -1 do print (i+1, rank jacEval(gns_{0..i}, rp)) 
rank rp jc I
betti I

