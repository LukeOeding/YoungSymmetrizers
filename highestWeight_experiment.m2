restart
highestWeight = (tab1,tab2,tab3) -> (
(dg, aa, bb, cc) = (length flatten tab1, length tab1_0 ,length tab2_0 ,length tab3_0);
Rx := QQ[x_(1,1,1)..x_(aa,bb,cc)];
genMat := (Rg, a,b) -> value toString transpose genericMatrix(Rg,Rg_0, a,b);
listRa = {};
for i from 1 to length tab1 list (
if i > 1 then listRa = listRa | {(last listRa)[flatten apply(length tab1_(i-1), k-> apply(tab1_(i-1), j-> a_(k+1,j)))]}
else listRa = listRa | {Rx[flatten apply(length tab1_(i-1), k-> apply(tab1_(i-1), j->(a_(k+1,j))))]}
);
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
R = Rx[(unique flatten (gens\ listRa))|(unique flatten (gens\ listRb))|(unique flatten (gens\ listRc))];
detListA = {};
for i from 0 to ((length listRa) -1) do (
detListA = detListA | {sub(det(genMat(listRa_i, length tab1_i,length tab1_i)), R)}
);
detListB = {};
for i from 0 to ((length listRb) -1) do (
detListB = detListB | {sub(det(genMat(listRb_i, length tab2_i,length tab2_i)),R)}
);
detListC = {};
for i from 0 to ((length listRc) -1) do (
detListC = detListC | {sub(det(genMat(listRc_i, length tab3_i,length tab3_i)),R)}
);
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
F1 = sum(Lc_0,k-> sum(Lb_0, j-> sum(La_0,i->x_(i+1,j+1,k+1)*diff(a_(i+1,d)*b_(j+1,d)*c_(k+1,d),F1))));
);
 print("---F1---");
sub(F1, Rx)
)
tab1 = {{1,2}, {3,4}}
tab2 = {{1,2}, {3,4}}
tab3 = {{1,3}, {2,4}}
highestWeight(tab1,tab2,tab3)
tab1 = {{1,2,3},{4,5,6},{7,8,9}}
tab2 = {{1,4,7},{2,5,8},{3,6,9}}
tab3 = {{1,4,6},{2,5,7},{3,8,9}}

elapsedTime strassen9 = highestWeight(tab1,tab2,tab3);
#terms strassen9


tab1 = {{1,2,3},{4}}
tab2 = {{1,2,4},{3}}
tab3 = {{1,3,4},{2}}

highestWeight(tab1,tab2,tab3)