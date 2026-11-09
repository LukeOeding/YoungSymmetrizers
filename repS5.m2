-- representations of S_5
S5 = {e, s12, s123, s1234, s12345, s12s34, s12s345}
wts = {1,10,20,30,24,15,20}
ip = (x,y) -> 1/120*sum(7, i-> x_i*y_i*wts_i)
R = QQ[a_1..a_6]
A = matrix{{1, 1, 1, 1, 1, 1, 1}, 
           {1,-1, 1,-1, 1, 1,-1},
           {4, 2, 1, 0,-1, 0,-1},
           {4,-2, 1, 0,-1, 0, 1},
           {6, 0, 0, 0, 1,-2, 0},
           {5,a_1,a_2,a_3,a_4,a_5,a_6},
           {5,-a_1,a_2,-a_3,a_4,a_5,-a_6}}

ip(flatten entries A^{1}, flatten entries A^{5})

M = (transpose matrix {apply(5,i->  ip(flatten entries A^{i}, flatten entries A^{5}))})||
(transpose matrix {apply(5,i->  ip(flatten entries A^{i}, flatten entries A^{6}))})
-- Either solve the linear system on the a_i given by the inner prodcuts, or do this:
mingens ideal M -- these are the minimal expressions that have to be zero.

-- use this and the fact that the a_i are all +/-1 to conclude once you pick a_6=1, a_1 = a_6, and a_3=-a_6
X = sub(A, {a_5=> 1, a_4=>0, a_2=>-1, a_6=>1, a_1=>1, a_3=>-1})

matrix apply(7, i-> apply(7,j->  ip(flatten entries X^{i}, flatten entries X^{j})))