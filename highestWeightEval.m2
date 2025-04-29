restart
highestWeightEval = (tab1,tab2,tab3, rules) -> (
    -- inputs 3 tableaux and some rules for the values of x_i
(dg, aa, bb, cc) = (length flatten tab1, length tab1_0 ,length tab2_0 ,length tab3_0);
-- print(dg, aa, bb, cc);
Rx := QQ[x_(1,1,1)..x_(aa,bb,cc)];
--for i from 1 to aa do for j from 1 to bb do for k from 1 to cc do X_(i,j,k) = sub(x_(i,j,k), rules);
 -- print("---dg---");
 -- print(dg);
--Rx := QQ;--QQ[x_(1,1,1)..x_(aa,bb,cc)];
 -- print(Rx);
genMat := (Rg, a,b) -> value toString transpose genericMatrix(Rg,Rg_0, a,b);
listRa = {};
for i from 1 to length tab1 list (
if i > 1 then listRa = listRa | {(last listRa)[flatten apply(length tab1_(i-1), k-> apply(tab1_(i-1), j-> a_(k+1,j)))]}
else listRa = listRa | {Rx[flatten apply(length tab1_(i-1), k-> apply(tab1_(i-1), j->(a_(k+1,j))))]}
);
 -- print("---space---");
 -- print(listRa);
 -- print(genMat(listRa_0, length tab1_0,length tab1_0));
--  -- print("---space---");
--  -- print(genMat(listRa_1, length tab1_1, length tab1_1));
detListA = {};
for i from 0 to ((length listRa) -1) do (
detListA = detListA | {det(genMat(listRa_i, length tab1_i,length tab1_i))}
);
 -- print("---det---");
 -- print(detListA);
prodDetA = product detListA;
 -- print("---product---");
 -- print(prodDetA);
listRb = {};
for i from 1 to length tab2 list (
if i > 1 then (
listRb = listRb | {(last listRb)[flatten apply(length tab2_(i-1), k->
apply(tab2_(i-1), j-> b_(k+1,j)))]}
)
else (
listRb = listRb | {(last listRa)[flatten apply(length tab2_(i-1), k->
apply(tab2_(i-1), j-> b_(k+1,j)))]}
)
);
 -- print("---space---");
 -- print(listRb);
 -- print(genMat(listRb_0, length tab2_0,length tab2_0));
--  -- print("---space---");
--  -- print(genMat(listRb_1, length tab2_1, length tab2_1));
detListB = {};
for i from 0 to ((length listRb) -1) do (
detListB = detListB | {det(genMat(listRb_i, length tab2_i,length tab2_i))}
);
 -- print("---det---");
 -- print(detListB);
prodDetB = product detListB;
 -- print("---product---");
 -- print(prodDetB);
listRc = {};
for i from 1 to length tab3 list (
if i > 1 then (
listRc = listRc | {(last listRc)[flatten apply(length tab3_(i-1), k->
apply(tab3_(i-1), j-> c_(k+1,j)))]}
)
else (
listRc = listRc | {(last listRb)[flatten apply(length tab3_(i-1), k->
apply(tab3_(i-1), j-> c_(k+1,j)))]}
)
);
 -- print("---space---");
 -- print(genMat(listRc_0, length tab3_0,length tab3_0));
--  -- print("---space---");
detListC = {};
for i from 0 to ((length listRc) -1) do (
detListC = detListC | {det(genMat(listRc_i, length tab3_i,length tab3_i))}
);
 -- print("---det---");
 -- print(detListC);
prodDetC = product detListC;
 -- print("---product---");
 -- print(prodDetC);
limitTab := (d,tab) -> (
for item in tab do (
if (isMember(d,item)) then (
return {length item, position(tab,i->i==item)}
)
)
);
F1 = 1;
templistA={};
templistB={};
templistC={};
for d from 1 to dg do (
Lc:=limitTab(d,tab3);
Lb:=limitTab(d,tab2);
La:=limitTab(d,tab1);
La1:=La_1;
Lb1:=Lb_1;
Lc1:=Lc_1;
if (isMember(La1,templistA)==false) then (
F1=F1*detListA_La1;
templistA=templistA|{La1});
if(isMember(Lb1,templistB)==false) then (
F1=F1*detListB_Lb1;
templistB=templistB|{Lb1} );
if (isMember(Lc1,templistC)==false) then (
F1=F1*detListC_Lc1;
templistC=templistC|{Lc1} );
F1 = sum(Lc_0,k-> sum(Lb_0, j-> sum(La_0,i->sub(x_(i+1,j+1,k+1),rules)*diff(a_(i+1,d)*b_(j+1,d)*c_(k+1,d),F1))));
);
 -- print("---F1---");
sub(F1,QQ)
)
tab1 = {{1,2}, {3,4}}
tab2 = {{1,2}, {3,4}}
tab3 = {{1,3}, {2,4}}
u1 = flatten entries random(QQ^3);u2 = flatten entries random(QQ^3);u3 = flatten entries random(QQ^3);u4 = flatten entries random(QQ^3);u5 = flatten entries random(QQ^3);
v1 = flatten entries random(QQ^3);v2 = flatten entries random(QQ^3);v3 = flatten entries random(QQ^3);v4 = flatten entries random(QQ^3);
v5 = flatten entries random(QQ^3);
w1 = flatten entries random(QQ^3);w2 = flatten entries random(QQ^3);w3 = flatten entries random(QQ^3);w4 = flatten entries random(QQ^3);
w5 = flatten entries random(QQ^3);

rules4 = flatten flatten  apply(3, i-> apply(3, j-> apply(3, k->  (x_(i+1,j+1,k+1) => 
u1_i*v1_j*w1_k+
u2_i*v2_j*w2_k+
u3_i*v3_j*w3_k+
u4_i*v4_j*w4_k
) )))
elapsedTime highestWeightEval(tab1,tab2,tab3, rules4)

rules5 = flatten flatten  apply(3, i-> apply(3, j-> apply(3, k->  (x_(i+1,j+1,k+1) => 
u1_i*v1_j*w1_k +
u2_i*v2_j*w2_k +
u3_i*v3_j*w3_k +
u4_i*v4_j*w4_k + 
u5_i*v5_j*w5_k
) )))

elapsedTime highestWeightEval(tab1,tab2,tab3, rules4)