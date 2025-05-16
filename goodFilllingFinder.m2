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