tmp = (a) -> (x := local x; R := QQ[x_0..x_2]; factor(random(10^10)))
parallelApply(1..10, i-> tmp(i))

tmp = (a) -> (x := local x; R := QQ[x_0..x_2]; factor(random(10^10)))
parallelApply(1..10, i-> tmp(i))


tmp = (a) -> (x := local x; x_0; x_1; R := QQ[x_0,x_1]; 1); parallelApply(1..10, i-> tmp(i))

QQ(new Array from apply(2, i-> value( "local x_"|toString(i))))
tmp = (a) -> (
   for i to 1 do value( "local x"|toString(i));
    QQ(new Array from apply(2, i-> value( "x"|toString(i)))); 1); parallelApply(1..10, i-> tmp(i))

tmp = (a) -> (x := local x; x_0; x_1; R := QQ[x_0,x_1]; 1); parallelApply(1..10, i-> tmp(i))
QQ[Variables => 2]
tmp = (a) -> (QQ["x"|toString(a)]; 1);  parallelApply(1..10, i-> tmp(i))
?Variables
methods lock
tmp = (a) -> (x := local x; R := QQ[x]; factor(random(10^10)))
parallelApply(1..10, i-> tmp(i))

tmp = (a) -> (local x; local y; R := QQ[x,y]; factor(random(10^10)))
parallelApply(1..10, i-> tmp(i))

tmp = (a) -> (x:=local x; R := QQ[x_0..x_3]; factor(random(10^10)))
parallelApply(1..10, i-> tmp(i))

tmp = (a) -> (x:=local x; R := QQ[x_0,x_1];)
parallelApply(1..10, i-> tmp(i))

tmp = (a) -> (local x; local y; R := QQ[x,y]; 1)
parallelApply(1..10, i-> tmp(i))


restart

M = random(QQ^100,QQ^90);
N = random(QQ^90,QQ^100);
A = M*N;
v = random(QQ^90,QQ^1);
(A_{0..89})*v

rank((A_{0..89})*v |A_{90})