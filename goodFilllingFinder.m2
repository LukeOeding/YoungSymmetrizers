-- find all the linearly independent good fillings


indFill8 = {};tmp = 0; tmp2 = 0;
L8a = {syt8#0,syt8#0,syt8#2,syt8#9,syt8#12}
L8b = {syt8#0,syt8#2,syt8#3,syt8#12,syt8#13}
L8c = {syt8#0,syt8#0,syt8#12,syt8#12,syt8#12}
L8d = {syt8#0,syt8#3,syt8#6,syt8#9,syt8#13}

plist = unique permutations(L8a);
fl8 = for pp in plist list (tmpInv = tab2polyEval(pp, randomLinear); tmp2 = numgens ideal mingens ideal (indFill8|{tmpInv}); print(tmp2);
    if tmp2>tmp then (tmp = tmp2; indFill8 = indFill8|{tmpInv}; pp) else continue); -- 28 invariants
plist = unique permutations(L8b);
fl8 = fl8 | for pp in plist list (tmpInv = tab2polyEval(pp, randomLinear); tmp2 = numgens ideal mingens ideal (indFill8|{tmpInv}); print(tmp2);
    if tmp2>tmp then (tmp = tmp2; indFill8 = indFill8|{tmpInv}; pp) else continue); -- bumps up to 34
plist = {L8c}; -- I know that permutations aren't necessary from prior computation
fl8 = fl8 |for pp in plist list (tmpInv = tab2polyEval(pp, randomLinear); tmp2 = numgens ideal mingens ideal (indFill8|{tmpInv}); print(tmp2);
    if tmp2>tmp then (tmp = tmp2; indFill8 = indFill8|{tmpInv}; pp) else continue);
plist = {L8d};-- I know that permutations aren't necessary from prior computation
fl8 = fl8 |for pp in plist list (tmpInv = tab2polyEval(pp, randomLinear); tmp2 = numgens ideal mingens ideal (indFill8|{tmpInv}); print(tmp2);    if tmp2>tmp then (tmp = tmp2; indFill8 = indFill8|{tmpInv}; pp) else continue);

ff8 = for xx in fl8 list flatten for j to 4 list (for i to #syt8-1 list if syt8#i ==xx#j then i else continue)
toString ff8

-- L8a = {syt8#0,syt8#0,syt8#2,syt8#9,syt8#12}
-- elapsedTime inv8p0 = tab2polyEval(L8a_{0,0,2,3,4}, randomLinear); 
-- elapsedTime inv8p1 = tab2polyEval(L8a_{0,0,3,2,4}, randomLinear); 
-- elapsedTime inv8p2 = tab2polyEval(L8a_{0,0,3,4,2}, randomLinear); 
-- elapsedTime inv8p3 = tab2polyEval(L8a_{0,2,0,3,4}, randomLinear); 
-- elapsedTime inv8p4 = tab2polyEval(L8a_{0,2,3,0,4}, randomLinear); 
-- elapsedTime inv8p5 = tab2polyEval(L8a_{0,2,3,4,0}, randomLinear); 
-- elapsedTime inv8p6 = tab2polyEval(L8a_{0,3,0,2,4}, randomLinear); 
-- elapsedTime inv8p7 = tab2polyEval(L8a_{0,3,0,4,2}, randomLinear); 
-- elapsedTime inv8p8 = tab2polyEval(L8a_{0,3,2,0,4}, randomLinear); 
-- elapsedTime inv8p9 = tab2polyEval(L8a_{0,3,2,4,0}, randomLinear); 
-- elapsedTime inv8p10 = tab2polyEval(L8a_{0,3,4,0,2}, randomLinear); 

-- indFill8 = {};tmp = 0; tmp2 = 0;
L8a = {syt8#0,syt8#0,syt8#2,syt8#9,syt8#12}
L8b = {syt8#0,syt8#2,syt8#3,syt8#12,syt8#13}
L8c = {syt8#0,syt8#0,syt8#12,syt8#12,syt8#12}
L8d = {syt8#0,syt8#3,syt8#6,syt8#9,syt8#13}


-- plist = unique permutations(L8a);
-- fl8 = for pp in plist list (tmpInv = tab2polyEval(pp, randomLinear); tmp2 = numgens ideal mingens ideal (indFill8|{tmpInv}); print(tmp2);
--     if tmp2>tmp then (tmp = tmp2; indFill8 = indFill8|{tmpInv}; pp) else continue); -- 28 invariants
-- plist = unique permutations(L8b);
-- fl8 = fl8 | for pp in plist list (tmpInv = tab2polyEval(pp, randomLinear); tmp2 = numgens ideal mingens ideal (indFill8|{tmpInv}); print(tmp2);
--     if tmp2>tmp then (tmp = tmp2; indFill8 = indFill8|{tmpInv}; pp) else continue); -- bumps up to 34
-- plist = {L8c}; -- I know that permutations aren't necessary from prior computation
-- fl8 = fl8 |for pp in plist list (tmpInv = tab2polyEval(pp, randomLinear); tmp2 = numgens ideal mingens ideal (indFill8|{tmpInv}); print(tmp2);
--     if tmp2>tmp then (tmp = tmp2; indFill8 = indFill8|{tmpInv}; pp) else continue);
-- plist = {L8d};-- I know that permutations aren't necessary from prior computation
-- fl8 = fl8 |for pp in plist list (tmpInv = tab2polyEval(pp, randomLinear); tmp2 = numgens ideal mingens ideal (indFill8|{tmpInv}); print(tmp2);    if tmp2>tmp then (tmp = tmp2; indFill8 = indFill8|{tmpInv}; pp) else continue);

-- ff8 = for xx in fl8 list flatten for j to 4 list (for i to #syt8-1 list if syt8#i ==xx#j then i else continue)
-- toString ff8
inv8a = tab2polyEval(L8a, randomLinear); -- 
inv8b = tab2polyEval(L8b, randomLinear); -- 
inv8c = tab2polyEval(L8c, randomLinear); -- 
inv8d = tab2polyEval(L8d, randomLinear); -- 
testRelns = fn -> (
myList = {inv4p1,inv4p2, inv4p3, inv4p4, inv4p5};
qd = flatten for i to 4 list for j from i to 4 list myList#i*myList#j;
Mtmp = rp matrix {qd|{fn}};
for i to 17 do (rp = rndPt(); Mtmp = Mtmp || rp matrix {qd|{fn}});
ker Mtmp)

testRelns(inv8a)
testRelns(inv8b)
testRelns(inv8c)
testRelns(inv8d)

time I8 = apply(goodTabList8, pp-> tab2polyEval(apply(5, i-> syt8#(pp#i)), randomLinear));
myList = {inv4p1,inv4p2, inv4p3, inv4p4, inv4p5};
qd = flatten for i to 4 list for j from i to 4 list myList#i*myList#j;
#I8
Mtmp = rp matrix {qd|I8};
for j to 53 do (rp = rndPt(); Mtmp = Mtmp || rp matrix {qd|I8});
ker Mtmp
rank Mtmp
-- find minimal generators:
li = ideal(inv4p1,inv4p2, inv4p3, inv4p4, inv4p5, inv6);
indFill8 = {};tmp = 0; tmp2 = 0;
plist = goodTabList8;
mg8 = for pp in plist list (tmpInv = tab2polyEval(apply(5, i-> syt8#(pp#i)), randomLinear); tmp2 = numgens ideal mingens (li + ideal (indFill8|{tmpInv})); print(tmp2);
    if tmp2>tmp then (tmp = tmp2; indFill8 = indFill8|{tmpInv}; pp) else continue); -- 
toString mg8
I8 = ideal indFill8;
rank rp jc (li + ideal I8) -- rank 15 (previously we had 17...?)
betti (li + ideal I8), betti mingens (li + ideal I8)
betti mingens (li + ideal indFill8)
time betti mingens (li + I8)
rank jacEval(flatten entries gens (li+I8), rp)
rank rp jc (li+I8)
-- primaryInvariants = {inv4p1,inv4p2, inv4p3, inv4p4, inv4p5, inv6, inv8p0,inv8p1,inv8p2, inv8p3, inv8p4, inv8p5, inv8p6, inv8p7, inv8p8, inv8p9, inv8p10};
-- li = I8+ ideal primaryInvariants;
-- betti mingens li, betti mingens ideal primaryInvariants
-- this would say that there are more than 5 degree 8 invariants that are algebraically independent, which is not what we expect.

-- plist = unique apply(permutations({1,2,3,4}), L-> ({L8a#0})|(L8a_L));
-- tmp = rank rp jc li
-- rp = rndPt();
-- tmp = rank rp jc (li + I8) -- the jacobian predicts something lower rank over Z/101 with a 
-- betti mingens (li + I8)

-- for L in plist do (
--     time invtmp = tab2polyEval(L, randomLinear); 
--     li = li+ ideal invtmp;
--     tmp2 = rank rp jc (li); print(tmp2);
--     if tmp2 >tmp then (tmp = tmp2; primaryInvariants = primaryInvariants |{invtmp}; li = li + ideal invtmp; print (L); print(betti li););
-- )
-- betti mingens li, betti mingens ideal primaryInvariants

-- L8ar = L8a_{2,0,1,3,4};
-- plist = unique apply(permutations({1,2,3,4}), L-> ({L8ar#0})|(L8ar_L));
-- tmp = rank rp jc li
-- for L in plist do (
--     time invtmp = tab2polyEval(L, randomLinear); 
--     li = li+ ideal invtmp;
--     tmp2 = rank rp jc (li); print(tmp2);
--     if tmp2 >tmp then (tmp = tmp2; primaryInvariants = primaryInvariants |{invtmp}; li = li + ideal invtmp; print (L); print(betti li););
-- )
-- betti mingens li, betti mingens ideal primaryInvariants
-- L8ar = L8a_{3,0,1,2,4};
-- plist = unique apply(permutations({1,2,3,4}), L-> ({L8ar#0})|(L8ar_L));
-- tmp = rank rp jc li
-- for L in plist do (
--     time invtmp = tab2polyEval(L, randomLinear); 
--     li = li+ ideal invtmp;
--     tmp2 = rank rp jc (li); print(tmp2);
--     if tmp2 >tmp then (tmp = tmp2; primaryInvariants = primaryInvariants |{invtmp}; li = li + ideal invtmp; print (L); print(betti li););
-- )
-- betti mingens li, betti mingens ideal primaryInvariants

-- L8ar = L8a_{4,0,1,2,3};
-- plist = unique apply(permutations({1,2,3,4}), L-> ({L8ar#0})|(L8ar_L));
-- tmp = rank rp jc li
-- for L in plist do (
--     time invtmp = tab2polyEval(L, randomLinear); 
--     li = li+ ideal invtmp;
--     tmp2 = rank rp jc (li); print(tmp2);
--     if tmp2 >tmp then (tmp = tmp2; primaryInvariants = primaryInvariants |{invtmp}; li = li + ideal invtmp; print (L); print(betti li););
-- )
-- betti mingens li, betti mingens ideal primaryInvariants
--tl = {{{1, 2}, {3, 4}, {5, 6}, {7, 8}}, {{1, 3}, {2, 6}, {4, 7}, {5, 8}}, {{1, 2}, {3, 4}, {5, 6}, {7, 8}}, {{1, 4}, {2, 6}, {3, 7}, {5, 8}}, {{1, 2}, {3, 5}, {4, 6}, {7, 8}}};
--for j to 4 list first for i to 4 list (if L8a#i == tl#j then i else continue)

L8b = {syt8#0,syt8#3,syt8#6,syt8#9,syt8#13}
plist = unique permutations(L8b);
fl8 = fl8|for pp in plist list (tmpInv = tab2polyEval(pp, randomLinear); tmp2 = numgens ideal mingens ideal (indFill8|{tmpInv}); print(tmp2);
    if tmp2>tmp then (tmp = tmp2; indFill8 = indFill8|{tmpInv}; pp) else continue)
ff8 = for xx in fl8 list flatten for j to 4 list (for i to #syt8-1 list if syt8#i ==xx#j then i else continue)
toString ff8
betti mingens ideal indFill8

L8c = {syt8#0,syt8#0,syt8#12,syt8#12,syt8#12}
plist = unique permutations(L8c);
fl8 = fl8|for pp in plist list (tmpInv = tab2polyEval(pp, randomLinear); tmp2 = numgens ideal mingens ideal (indFill8|{tmpInv}); print(tmp2);
    if tmp2>tmp then (tmp = tmp2; indFill8 = indFill8|{tmpInv}; pp) else continue)


plist = unique permutations(L8d);
fl8 = fl8|for pp in plist list (tmpInv = tab2polyEval(pp, randomLinear); tmp2 = numgens ideal mingens ideal (indFill8|{tmpInv}); print(tmp2);
    if tmp2>tmp then (tmp = tmp2; indFill8 = indFill8|{tmpInv}; pp) else continue)

ff8 = for xx in fl8 list flatten for j to 4 list (for i to #syt8-1 list if syt8#i ==xx#j then i else continue)
toString ff8
betti mingens ideal indFill8


I8 = I8+ ideal inv8b;
betti mingens (li + I8)
rank rp jc (li + I8)
betti (li + I8)
time betti mingens (li + I8)
rank jacEval(flatten entries gens (li+I8), rp)

-- betti mingens li, betti mingens ideal primaryInvariants
--toString (transpose\ syt8)

-- L8c = {syt8#0,syt8#0,syt8#12,syt8#12,syt8#12}
--  elapsedTime inv8c = tab2polyEval(L8c, randomLinear); 

-- myList = {inv4p1,inv4p2, inv4p3, inv4p4, inv4p5, inv8c};
-- qd = flatten for i to 4 list for j from i to 4 list myList#i*myList#j;
-- flatten for i to 4 list for j from i to 4 list {i,j}
-- rndPt = () -> map(KK, X, random(KK^1, KK^(numgens X))) 
-- rp = rndPt();
-- Mtmp = rp matrix {qd|{inv8c}}
-- for i to 17 do (rp = rndPt(); Mtmp = Mtmp || rp matrix {qd|{inv8c}};)
-- ker Mtmp
-- Tm = QQ[a,b]
-- factor (1/16*a^2 + 3/8*a*b + 1/16*b^2)

--
-- fl = {{1,1,5,35,41},{1,1,35,5,41},{1,1,35,41,5},{1,5,1,35,41},{1,5,35,1,41},{1,5,35,41,1},{1,35,1,5,41},{1,35,1,41,5},{1,35,41,1,5},{1,18,25,38,41}}
-- apply(fl, ff -> ff-{1,1,1,1,1})
-- #fl
-- -- L10b ={syt10#0,syt10#17,syt10#24,syt10#37,syt10#40}
-- -- elapsedTime inv10b = tab2polyEval(L10b,randomLinear); -- approx 448s -> 6s
-- I10 = ideal for ff in fl list time tab2polyEval(apply(ff, i-> syt10#i),randomLinear);

-- betti I10
-- betti mingens I10
---

-- L10 = {syt10#0,syt10#0,syt10#4,syt10#34,syt10#40};
-- L10b = {syt10#0,syt10#17,syt10#24,syt10#37,syt10#40};
-- tmp = rank rp jc li
-- tmp2 = tmp;
-- plist = unique apply(permutations({1,2,3,4}), L-> ({L10#0})|(L10_L));

-- for L in plist do (
--     time invtmp = tab2polyEval(L, randomLinear); 
--     tmp2 = rank rp jc (li+ ideal invtmp); print(tmp2);
--     if tmp2 >tmp then (tmp = tmp2; li = li + ideal invtmp; print (L); print(betti li););
-- )
-- betti li
-- time betti mingens li

-- this would say that there are more than 5 degree 8 invariants that are algebraically independent, which is not what we expect.
fl = {{1,1,5,35,41},{1,1,35,5,41},{1,1,35,41,5},{1,5,1,35,41},{1,5,35,1,41},{1,5,35,41,1},{1,35,1,5,41},{1,35,1,41,5},{1,35,41,1,5},{1,18,25,38,41}}
fl = apply(fl, ff -> ff-{1,1,1,1,1})
# fl
--L10b ={syt10#0,syt10#17,syt10#24,syt10#37,syt10#40}
--elapsedTime inv10b = tab2polyEval(L10b,randomLinear); -- approx 448s -> 6s
I10 = ideal for ff in fl list time tab2polyEval(apply(ff, i-> syt10#i),randomLinear);
betti mingens I10

I = ideal {inv4p1,inv4p2, inv4p3, inv4p4, inv4p5, inv6}+  I8 +  I10;
betti I
rank rp jc I
gns  = flatten entries gens I;
for i from 1 to #gns -1 do print (i, rank jacEval(gns_{0..i}, rp)) 
rank rp jc I
betti I
-- li = li + I10;
-- betti li
-- time betti mingens li
-- the whole list for degree 10:
--{1,1,5,35,41},{1,1,35,5,41},{1,1,35,41,5},{1,5,1,35,41},{1,5,35,1,41},{1,5,35,41,1},{1,35,1,5,41},{1,35,1,41,5},{1,35,41,1,5},{1,18,25,38,41},{1,18,25,41,38},{1,18,41,25,38},{1,25,18,38,41},{18,1,25,38,41},{1,13,31,32,39}
--time print(betti mingens I12)
I12 = {};
for L in minList12 do (elapsedTime I12 =  I12 |{ tab2polyEval(apply(L, l-> syt12#l),randomLinear) };)
betti ideal I12
I12 = ideal I12;
--I12 = ideal mingens I12;

betti mingens I12
li = li + I12;
time betti mingens li
{flatten entries gens I12}
primaryInvariants = {inv4p1,inv4p2, inv4p3, inv4p4, inv4p5, inv6, inv8p0,inv8p1,inv8p2, inv8p3, inv8p4}|{last flatten entries gens I10}|(flatten entries gens I12);
betti ideal primaryInvariants
--time rank rp jc ideal primaryInvariants
-- this takes a long time - so instead we should evaluate each form on all but one variable and then take that derivative...
rank jacEval(primaryInvariants, rp) -- this happens in an instant!
I12;
I = ideal {inv4p1,inv4p2, inv4p3, inv4p4, inv4p5, inv6}+  I8+  I10 + I12;
rank jacEval(gens I, rp) -- this happens in an instant!
gns  = flatten entries gens I;
for i from 1 to #gns -1 do print (i, rank jacEval(gns_{0..i}, rp)) 

-- jc = I-> diff(transpose bx, gens I ); -- this can take a long time...


-- L12 = {syt12#0, syt12#43, syt12#54, syt12#121, syt12#124}
-- elapsedTime inv12 = tab2polyEval(L12, randomLinear); --
-- tmp = rank rp jc (li + ideal inv12)
-- tmp2 = tmp;
-- plist = unique apply(permutations({1,2,3,4}), L-> ({L12#0})|(L12_L));
-- li = li + ideal inv12;
-- time betti mingens li
-- for L in plist do (
--     time invtmp = tab2polyEval(L, randomLinear); 
--     tmp2 = rank rp jc (li+ ideal invtmp); print(tmp2);
--     if tmp2 >tmp then (tmp = tmp2; li = li + ideal invtmp; print (L); print(betti li););
-- )
-- betti li
-- time betti mingens li


-- elapsedTime inv12 = tab2polyEval(L12,randomLinear); -- 52.77s
-- #terms inv12
-- primaryInvariants = {inv4a,inv4b, inv4c, inv4d, inv4e , inv6, inv8,inv813,inv814,inv815, inv8b,inv10b, inv12};
-- li = ideal(primaryInvariants);
-- print(rank rp jacobian li, betti li, betti mingens li);
-- -- this suggests that the deg 12 invariant we chose is not algebraically independent. 

-- tmp = 12;
-- for L12 in fill12 do if tmp==12 then(
--     print(L12);
-- elapsedTime inv12 = tab2polyEval(L12,randomLinear); -- 52.77s
-- primaryInvariants = {inv4a,inv4b, inv4c, inv4d, inv4e , inv6, inv8,inv813,inv814,inv815, inv8b,inv10b, inv12};
-- li = ideal(primaryInvariants);
-- tmp= rank rp jc li; print(tmp);
-- )
-- print(tmp, betti li, betti mingens li);



-- random

-- p0 = x_(1,0,0,0,0)+x_(0,1,1,1,1)
-- p1 = x_(0,1,0,0,0)+x_(1,0,1,1,1)
-- p2 = x_(0,0,1,0,0)+x_(1,1,0,1,1)
-- p3 = x_(0,0,0,1,0)+x_(1,1,1,0,1)
-- p4 = x_(0,0,0,0,1)+x_(1,1,1,1,0)
-- --  p_{0,\pm} = e_{1}e_{2}e_{4}e_{6}e_{8} \pm e_{0}e_{3}e_{5}e_{7}e_{9},&
-- --  p_{1,\pm} = e_{0}e_{3}e_{4}e_{6}e_{8} \pm e_{1}e_{2}e_{5}e_{7}e_{9},\\
-- --  p_{2,\pm} = e_{0}e_{2}e_{5}e_{6}e_{8} \pm e_{1}e_{3}e_{4}e_{7}e_{9},&
-- --  p_{3,\pm} = e_{0}e_{2}e_{4}e_{7}e_{8} \pm e_{1}e_{3}e_{5}e_{6}e_{9},\\
-- --  p_{4,\pm} = e_{1}e_{3}e_{5}e_{7}e_{8} \pm e_{0}e_{2}e_{4}e_{6}e_{9}.
-- usedList = (terms p0)|(terms p1)|terms p2|terms p3|terms p4
-- r0 = random(KK);
-- r1 = random(KK);
-- r2 = random(KK);
-- r3 = random(KK);
-- r4 = random(KK);

-- pseudoCartan = {x_(1,0,0,0,0)=>r0*x_(0,1,1,1,1),x_(0,1,1,1,1 )=> -r0*x_(0,1,1,1,1),
--                 x_(0,1,0,0,0)=>r1*x_(1,0,1,1,1),x_(1,0,1,1,1) => -r1*x_(1,0,1,1,1),
--                 x_(0,0,1,0,0)=>r2*x_(1,1,0,1,1),x_(1,1,0,1,1) => -r2*x_(1,1,0,1,1),
--                 x_(0,0,0,1,0)=>r3*x_(1,1,1,0,1),x_(1,1,1,0,1) => -r3*x_(1,1,1,0,1),
--                 x_(0,0,0,0,1)=>r4*x_(1,1,1,1,0),x_(1,1,1,1,0) => -r4*x_(1,1,1,1,0)}|  for xx in gens X list if not member(xx, usedList) then xx=>0 else continue
-- ---- random evaluations
vl =  unique apply(15, i-> random(31)) -- pick 15 variables to zero
while #vl !=15 do vl = unique(vl|apply(15-#vl, i-> random(31)))
#vl

randomLinear = apply(vl, i-> gx_i=>0); --random dense points takes longer than just evaluation. Try random sparse...
bx = sub(basis(1,X),randomLinear)
jc = I-> diff(transpose bx, gens I )
L4 = {syt4#0,syt4#0,syt4#0,syt4#0,syt4#1}
inv4e = tab2polyEval(L4,randomLinear); 
elapsedTime inv4d = tab2polyEval(L4_{0,1,2,4,3},randomLinear); 
inv4c = tab2polyEval(L4_{0,1,4,2,3},randomLinear); 
inv4b = tab2polyEval(L4_{0,4,2,2,3},randomLinear); 
inv4a = tab2polyEval(L4_{4,2,1,2,3},randomLinear); 
rp = map(KK, X, random(KK^1, KK^(numgens X))) 

li = ideal(inv4a,inv4b, inv4c, inv4d, inv4e );
rank rp jacobian li, betti mingens li

inv6 = tab2polyEval(syt6,randomLinear);
li = li + ideal(inv6);
rank rp jacobian li, betti li, betti mingens li
rank rp jc li

L8a = {syt8#0,syt8#0,syt8#2,syt8#9,syt8#12}
elapsedTime inv8 = tab2polyEval(L8a, randomLinear); 
inv813 = tab2polyEval(L8a_{0,2,1,3,4}, randomLinear); 
inv814 = tab2polyEval(L8a_{0,2,3,1,4}, randomLinear); 
inv815 = tab2polyEval(L8a_{0,2,3,4,1}, randomLinear); 

L8b = {syt8#0,syt8#2,syt8#3,syt8#12,syt8#13}
inv8b = tab2polyEval(L8b, randomLinear); -- 3.49s
inv8b13 = tab2polyEval(L8b_{0,2,1,3,4}, randomLinear); -- 3.49s
--inv8c = tab2polyEval({syt8#0,syt8#0,syt8#12,syt8#12,syt8#12}, randomLinear); -- 4.88
L8d = {syt8#0,syt8#3,syt8#6,syt8#9,syt8#13}
inv8d = tab2polyEval(L8d, randomLinear); -- 4.3s

li = li + ideal(inv8, inv8b,inv8d); -- inv8c,
betti li, betti mingens li
primaryInvariants = {inv4a,inv4b, inv4c, inv4d, inv4e , inv6, inv8,inv813,inv814,inv815, inv8b,inv8d,inv8b13};
for i from 1 to #primaryInvariants do(li = ideal(primaryInvariants_{0..i-1});
print(rank rp jacobian li, betti li, betti mingens li);
)

--- symbolic evaluations
elapsedTime f4 = tab2poly(syt4#0,syt4#0,syt4#0,syt4#0,syt4#1);
elapsedTime f4 = tab2polyList(syt4#0,syt4#0,syt4#0,syt4#0,syt4#1,{0,3,1,2});
I4 = ideal(f4);
betti mingens I4 -- these are the first 5 primary invariants - but let's do a randomized jacobian to make sure we think they're algebraically independent.
time I4 = I4 + cycle I4 + cycle2 I4 + cycle3 I4 + cycle4 I4;
J4 = jacobian I4;
phi = rndPt();
rank phi J4 -- random rank 5 gives us some confidence that it is rank 5. -- I think: lower bounds are certain, upper bounds are probabilistic
--time ind4 = det J4^{0..4}; -- if this is not the zero polynomial, then we'd have confirmation that they're actually independent. However it is likely to take foreover. On the other hand, the zero polynomnial doesn't have non-zero values. So we could just evaluate it at random points and if it's not zero, then the rank is 5. That's what the randomized test does. 
--hilbertFunction(4,I4)
primaryInvariants = flatten entries gens I4;
elapsedTime inv6 = tab2poly(syt6#0, syt6#1, syt6#2, syt6#3, syt6#4); -- 1.09s (tensor product ring) vs .934s (flat ring)
elapsedTime inv6pc = tab2polyEval(syt6#0, syt6#1, syt6#2, syt6#3, syt6#4,randomLinear); 
elapsedTime inv6p = tab2polyEval(syt6#0, syt6#1, syt6#2, syt6#3, syt6#4,randomPoint)
factor inv6l
--for L in permutations(6) do (print L; elapsedTime test = tab2polyList(syt6#0, syt6#1, syt6#2, syt6#3, syt6#4, L)); -- 1.09s (tensor product ring) vs .934s (flat ring)
--test == inv6

I6 = ideal inv6;
I = I4+I6;
J46 = jacobian(I); -- to show that these polynomials are algebraically independent, you could have a general point in their intersection, then compute their Jacobian. 
rp = rndPt();
rank rp J46
primaryInvariants = primaryInvariants|{inv6};


elapsedTime inv8 = tab2poly(syt8#0,syt8#0,syt8#2,syt8#9,syt8#12); -- 20s
elapsedTime inv8b = tab2poly(syt8#0,syt8#2,syt8#3,syt8#12,syt8#13); -- 30s
elapsedTime inv8c = tab2poly(syt8#0,syt8#0,syt8#12,syt8#12,syt8#12); -- 24s -- not independent of deg 4
elapsedTime inv8d = tab2poly(syt8#0,syt8#3,syt8#6,syt8#9,syt8#13); -- 39s


--I8 = ideal{ inv8, inv8c, inv8b,inv8d};  -- leave out inv8c from the non-minimal generators. 
--primaryInvariants = primaryInvariants|{inv8, inv8c, inv8d};
rank rp jacobian ideal primaryInvariants
PI8 = {inv8b , cycle inv8b , cycle2 inv8b , cycle3 inv8b , cycle4 inv8b};
rp = rndPt();
rank rp jacobian( ideal(PI8) + ideal primaryInvariants)
PI = ideal(PI8) + ideal primaryInvariants;
time betti mingens PI
primaryInvariants = primaryInvariants |PI8;

I = I4+I6+I8;
betti mingens (I)
betti I
J468 = jacobian(I); -- to show that these polynomials are algebraically independent, you could have a general point in their intersection, then compute their Jacobian. 
rp = rndPt();
rank rp J468
--
--J = jacobian(I4 + I6 + ideal{inv8, cycle inv8, cycle2 inv8, cycle3 inv8, cycle4 inv8}); 
--rp = rndPt();
--rank rp J

--time I8 = ideal(inv8, cycle inv8, cycle2 inv8, cycle3 inv8, cycle4 inv8); 
--time I8s = swap I8; betti I8
--time I8 = (I8+I8s); betti I8
time I8s = swap I8;
print betti (I8+I8s)
time I8 = I8 + I8s + cycle I8 + cycle2 I8 + cycle3 I8 + cycle4 I8  + cycle I8s + cycle2 I8s + cycle3 I8s + cycle4 I8s; 
time I8 = ideal mingens I8; betti I8 -- 25 polynomials when you only take inv8,inv8b,inv8d. 30 if you include inv8c
time I8s = swap I8; betti I8s
time I8 = ideal mingens (I8+I8s); betti I8 -- 36 poly's when you start from inv8, inv8b, inv8c, inv8d

I = I4+I6+I8;
betti I
tex betti mingens (I)
betti I
J468 = jacobian(I); -- to show that these polynomials are algebraically independent, you could have a general point in their intersection, then compute their Jacobian. 
rp = rndPt();
rank rp J468
--apply(9, i-> hilbertFunction(i, I))
netList\ transpose\ syt10
elapsedTime inv10b = tab2poly(syt10#0,syt10#17,syt10#24,syt10#37,syt10#40); -- approx 448s

time rank rp jacobian( ideal(inv10b) + ideal primaryInvariants)
PI = ideal(inv10b) + ideal primaryInvariants;
betti  PI

I10 = ideal inv10b;
time I = ideal mingens (I+I10);
betti I
elapsedTime inv10 = tab2poly(syt10#0,syt10#0,syt10#4,syt10#34,syt10#40); -- approx 724s
I10 = I10 +  ideal inv10b;
I = I+I10;
time I = ideal mingens I;
betti I

-- try to get the rest of the generators:
I10 = I10 + I6*I4;
I10s = swap I10;
I10 = I10 + I10s;
time I10 = I10 + cycle I10 + cycle2 I10 + cycle3 I10 + cycle4 I10;
time I10 = ideal mingens I10; 
time betti I10
I10s = swap I10;
I10 = I10 + I10s;
--time hilbertFunction(10,I10);
time I10 = ideal mingens I10; 
time betti I10  

I = I+I10;
betti I
I = ideal mingens (I)
tex oo
--apply(11, i-> hilbertFunction(i, I))

m12to2row \last fill12


g12 = {{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}, {1, 3, 2, 4, 5, 6, 7, 8, 9, 11, 10, 12}, {1, 3, 2, 4, 5, 9, 6, 10, 7, 11, 8, 12}, {1, 5, 2, 7, 3, 8, 4, 9, 6, 10, 11, 12}, {1, 5, 2, 7, 3, 9, 4, 10, 6, 11, 8, 12}}
for gg in g12 list for ss to #syt12 -1 list if gg==(syt12#ss) then ss else continue 
g = "g"
elapsedTime inv12 = tab2poly(syt12#0, syt12#43, syt12#54, syt12#121, syt12#124); -- approx ?
time rank rp jacobian( ideal(inv12) + PI)
PI = PI + ideal(inv12) 
betti  PI

time rank rp jacobian( ideal(inv12) + ideal primaryInvariants)
PI = ideal(inv10b) + ideal primaryInvariants;
betti  PI
