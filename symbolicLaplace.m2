--symbolic (block) Laplace expansion of determinants
--------------------
-- expand the determinant as a sum of determinants, but the determinants are symoblic until requested
--------------------
--restart
--------------------
-- Basic functions
--------------------

--------------------
-- Set Complement
--------------------
setComplement:=method();
setComplement (List, ZZ) := (I,n) ->  toList(set(0..(n-1)) - set I);
setComplement (List, List):= (I,P) -> if isSubset(set I, set P) then sort(toList(set P - set I)) else print("error, not a subset");
setComplement({1,2},{1,2,3})
setComplement({1,2},3)
--------------------
-- subsets containing a fixed element
--------------------
subsetsContaining = (i,j,k) -> for ss in subsets(i,j) list if member(k,ss) then ss else continue
subsetsContaining(4,2,1)
--------------------
-- locate in list
--------------------
-- assume I is a list with no repeats, find the position of i in I:
locateInList:=method();
locateInList (ZZ,List):= (i,I) -> for s to length I-1 do if I_s ==i then return s else continue;

----------------------------
-- tableaux to determinant (symbolic)
----------------------------
tab2det :=method();
tab2det (List, IndexedVariableTable) := (L, aa) -> product(L, l-> aa_(toList(0..(length(l) -1)), l ))
----------------------------
-- take the derivative of an expression that has symbolic determinants in it, but noting which determinants might be there.
----------------------------
Ldiff = (ij,tab,aa, F) -> (i:=ij#0, j:=ij#1;  diff(aa_(i,j),F) + sum flatten for J in tab list if (member(j,J) and #J>1) then (for I in subsetsContaining(3,#J,i) list(
        (-1)^(locateInList(i,I)- locateInList(j,J))*( aa_(setComplement({i},I), setComplement({j},J) ))* diff(aa_(I,J), F))) else continue
        )
----------------------------


----------------------------
-- try the 3 factor case for Young symmetrizers
----------------------------

R0 = QQ; 
R1 = R0[a_(0,0)..a_(3,3)]
R2 = R1[flatten apply(subsets(4,2), j-> apply(subsets(4,2), i-> a_(i,j)))]
R3 = R2[flatten apply(subsets(4,3), j-> apply(subsets(4,3), i-> a_(i,j)))]
R4 = R3[flatten apply(subsets(4,4), j-> apply(subsets(4,4), i-> a_(i,j)))]

-- Strassen's Degree 4 invariant;
--restart
R0= QQ[x_(0,0,0)..x_(2,2,2)];
RA = R0[a_(0,0)..a_(3,3),flatten flatten apply({2,3,4}, k-> apply(subsets(4,k), L->apply(subsets(4,k), M-> a_(L,M)  )))]
for i to 3 do for j to 3 do a_({i},{j}) = a_(i,j)
RB = RA[b_(0,0)..b_(3,3),flatten flatten apply({2,3,4}, k-> apply(subsets(4,k), L->apply(subsets(4,k), M-> b_(L,M)  )))]
for i to 3 do for j to 3 do b_({i},{j}) = b_(i,j)
RC = RB[c_(0,0)..c_(3,3),flatten flatten apply({2,3,4}, k-> apply(subsets(4,k), L->apply(subsets(4,k), M-> c_(L,M)  )))]
for i to 3 do for j to 3 do c_({i},{j}) = c_(i,j)

--R = RA**RB**RC
-- the list of lists -- inner lists are the column labels
tab1 = {{0,1,2},{3}}; tab1ss = tab1|flatten apply(tab1, tb -> subsets(tb,2));
tab2 = {{0,1,3},{2}}; tab2ss = tab2|flatten apply(tab2, tb -> subsets(tb,2));
tab3 = {{0,2,3},{1}}; tab3ss = tab3|flatten apply(tab3, tb -> subsets(tb,2))


f1 = tab2det(tab1, a)
f2 = tab2det(tab2, b)
f3 = tab2det(tab3, c)
F4 = f1*f2*f3;

F = F4;for dg to 3 do( time F = sum(3, i-> sum(3, j-> sum(3, k-> x_(i,j,k)* Ldiff({k,dg},tab3ss,c,Ldiff({j,dg},tab2ss,b,Ldiff({i,dg},tab1ss,a,F)))))); print( #terms F))

-- Strassen's Degree 9 invariant;
R0= QQ;
RA = R0[a_(0,0)..a_(3,8),flatten flatten apply({2,3}, k-> apply(subsets(3,k), L->apply(subsets(9,k), M-> a_(L,M)  )))]
for i to 3 do for j to 8 do a_({i},{j}) = a_(i,j)
RB = RA[b_(0,0)..b_(3,8),flatten flatten apply({2,3}, k-> apply(subsets(3,k), L->apply(subsets(9,k), M-> b_(L,M)  )))]
for i to 3 do for j to 8 do b_({i},{j}) = b_(i,j)
RC = RB[c_(0,0)..c_(3,8),flatten flatten apply({2,3}, k-> apply(subsets(3,k), L->apply(subsets(9,k), M-> c_(L,M)  )))]
for i to 3 do for j to 8 do c_({i},{j}) = c_(i,j)
R = RC[x_(0,0,0)..x_(2,2,2)]
tab1 = {{0,1,2},{3,4,5},{6,7,8}}; tab1ss = tab1|flatten apply(tab1, tb -> subsets(tb,2))
tab2 = {{0,3,6},{1,4,7},{2,5,8}}; tab2ss = tab2|flatten apply(tab2, tb -> subsets(tb,2))
tab3 = {{0,3,5},{1,4,6},{2,7,8}}; tab3ss = tab3|flatten apply(tab3, tb -> subsets(tb,2))

f1 = tab2det(tab1, a)
f2 = tab2det(tab2, b)
f3 = tab2det(tab3, c)
F9 = f1*f2*f3;

Ldiff({0,0},tab1ss,a,f1)
F9;#terms oo
F=F9;
for dg to 8 do( time F = sum(3, i-> sum(3, j-> sum(3, k-> x_(i,j,k)* Ldiff({k,dg},tab3ss,c,Ldiff({j,dg},tab2ss,b,Ldiff({i,dg},tab1ss,a,F)))))); print(#terms F))


R0= QQ[x_(0,0,0)..x_(2,2,2)];--- this way seems the fastest, but it doesn't beat just using the determinants...
RA = R0[a_(0,0)..a_(3,8),flatten flatten apply({2,3}, k-> apply(subsets(3,k), L->apply(subsets(9,k), M-> a_(L,M)  )))];
for i to 3 do for j to 8 do a_({i},{j}) = a_(i,j);
RB = RA[b_(0,0)..b_(3,8),flatten flatten apply({2,3}, k-> apply(subsets(3,k), L->apply(subsets(9,k), M-> b_(L,M)  )))];
for i to 3 do for j to 8 do b_({i},{j}) = b_(i,j);
RC = RB[c_(0,0)..c_(3,8),flatten flatten apply({2,3}, k-> apply(subsets(3,k), L->apply(subsets(9,k), M-> c_(L,M)  )))];
for i to 3 do for j to 8 do c_({i},{j}) = c_(i,j);
tab1 = {{0,1,2},{3,4,5},{6,7,8}}; tab1ss = tab1|flatten apply(tab1, tb -> subsets(tb,2));
tab2 = {{0,3,6},{1,4,7},{2,5,8}}; tab2ss = tab2|flatten apply(tab2, tb -> subsets(tb,2));
tab3 = {{0,3,5},{1,4,6},{2,7,8}}; tab3ss = tab3|flatten apply(tab3, tb -> subsets(tb,2));

f1 = tab2det(tab1, a);
f2 = tab2det(tab2, b);
f3 = tab2det(tab3, c);
F9 = f1*f2*f3;


F9;#terms oo
F=F9;
for dg to 8 do( time F = sum(3, i-> sum(3, j-> sum(3, k-> x_(i,j,k)* Ldiff({k,dg},tab3ss,c,Ldiff({j,dg},tab2ss,b,Ldiff({i,dg},tab1ss,a,F)))))); print( #terms F))


R0= QQ[x_(0,0,0)..x_(2,2,2),a_(0,0)..a_(3,8),flatten flatten apply({2,3}, k-> apply(subsets(3,k), L->apply(subsets(9,k), M-> a_(L,M)  ))),
b_(0,0)..b_(3,8),flatten flatten apply({2,3}, k-> apply(subsets(3,k), L->apply(subsets(9,k), M-> b_(L,M)  ))),
c_(0,0)..c_(3,8),flatten flatten apply({2,3}, k-> apply(subsets(3,k), L->apply(subsets(9,k), M-> c_(L,M)  )))];

for i to 3 do for j to 8 do a_({i},{j}) = a_(i,j);for i to 3 do for j to 8 do b_({i},{j}) = b_(i,j);for i to 3 do for j to 8 do c_({i},{j}) = c_(i,j);
Ldiff = (ij,tab,aa, F) -> (i:=ij#0, j:=ij#1;  diff(aa_(i,j),F) + sum flatten for J in tab list if member(j,J) then (for I in subsetsContaining(3,#J,i) list(
        (-1)^(locateInList(i,I)- locateInList(j,J))*( aa_(setComplement({i},I), setComplement({j},J) ))* diff(aa_(I,J), F))) else continue
        )
tab1 = {{0,1,2},{3,4,5},{6,7,8}}; tab1ss = tab1|flatten apply(tab1, tb -> subsets(tb,2));
tab2 = {{0,3,6},{1,4,7},{2,5,8}}; tab2ss = tab2|flatten apply(tab2, tb -> subsets(tb,2));
tab3 = {{0,3,5},{1,4,6},{2,7,8}}; tab3ss = tab3|flatten apply(tab3, tb -> subsets(tb,2));

f1 = tab2det(tab1, a);
f2 = tab2det(tab2, b);
f3 = tab2det(tab3, c);
F9 = f1*f2*f3;

F9;#terms oo
F=F9;
for dg to 8 do( time F = sum(3, i-> sum(3, j-> sum(3, k-> x_(i,j,k)* Ldiff({k,dg},tab3ss,c,Ldiff({j,dg},tab2ss,b,Ldiff({i,dg},tab1ss,a,F)))))); print( #terms F))
F
RX = QQ[x_(0,0,0)..x_(2,2,2)];
#terms sub(F,RX)