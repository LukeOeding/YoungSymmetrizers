# YoungSymmetrizers on Julia, using the modelKit variables from HomotopyContinuation, and matrices from LinearAlgebra
#
using HomotopyContinuation
using LinearAlgebra
using Base.Threads
using Combinatorics
syt2 = [[[1,2]]];
syt4 = [[[1,2],[3,4]],[[1,3],[2,4]]];
syt6 = [[[1,2],[3,4],[5,6]], [[1,2],[3,5],[4,6]],[[1,3],[2,4],[5,6]], [[1,3],[2,5],[4,6]], [[1,4],[2,5],[3,6]]];
syt8 = [[[1,2],[3,4],[5,6],[7,8]],[[1,2],[3,4],[5,7],[6,8]],[[1,2],[3,5],[4,6],[7,8]],[[1,2],[3,5],[4,7],[6,8]],[[1,2],[3,6],[4,7],[5,8]],[[1,3],[2,4],[5,6],[7,8]],[[1,3],[2,4],[5,7],[6,8]],[[1,3],[2,5],[4,6],[7,8]],[[1,3],[2,5],[4,7],[6,8]],[[1,3],[2,6],[4,7],[5,8]],[[1,4],[2,5],[3,6],[7,8]],[[1,4],[2,5],[3,7],[6,8]],[[1,4],[2,6],[3,7],[5,8]],[[1,5],[2,6],[3,7],[4,8]]];
gf81 = [[1, 2, 3, 4, 5, 6, 7, 8], [1, 2, 3, 4, 5, 6, 7, 8], [1, 2, 3, 5, 4, 6, 7, 8], [1, 3, 2, 6, 4, 7, 5, 8], [1, 4, 2, 6, 3, 7, 5, 8]]
@var x[1:2, 1:2, 1:2, 1:2, 1:2, 1:2]
@var a[1:2,1:8]
@var b[1:2,1:8]
@var c[1:2,1:8]
@var d[1:2,1:8]
@var e[1:2,1:8]
@var f[1:2,1:8]

list2mats = function(L,atmp)
     [stack([[atmp[i,j] for i=1:2] for j in ell]) for ell in L];
end
diffList = function(f,L)
local tmp;
    tmp = f;
    for xx in L tmp = differentiate(tmp,xx);end
tmp        
end
varSymbols = [a,b,c,d,e,f]
# this function is for 2^6 tensors
unfactor = function()
P = [1,1,1,1,1,1]
F = prod([prod(det.(list2mats(syt2[P[i]],varSymbols[i]))) for i=1:6])


P = [1,1,1,1,2,2]
F = prod([prod(det.(list2mats(syt4[P[i]],varSymbols[i]))) for i=1:6])
tmp = F; @time for dg=1:4
    @time tmp = sum([diffList(tmp, [a[i1,dg],b[i2,dg],c[i3,dg],d[i4,dg],e[i5,dg],f[i6,dg]])*x[i1,i2,i3,i4,i5,i6] for i1=1:2, i2=1:2, i3=1:2, i4=1:2, i5=1:2, i6=1:2]);
end

J = [differentiate(tmp,xx) for xx in x[1:2, 1:2, 1:2, 1:2, 1:2, 1:2]];

# compare to Maple:
P = [1,1,1,13,13]
@var x[1:2, 1:2, 1:2, 1:2, 1:2]
xvars = reduce(vcat, variables(:x, 1:2, 1:2, 1:2, 1:2, 1:2))
F = prod([prod(det.(list2mats(syt8[P[i]],varSymbols[i]))) for i=1:5])
tmp = F; @time for dg=1:8
    tmp = sum([diffList(tmp, [a[i1,dg],b[i2,dg],c[i3,dg],d[i4,dg],e[i5,dg]])*x[i1,i2,i3,i4,i5] for i1=1:2, i2=1:2, i3=1:2, i4=1:2, i5=1:2]);
end
xvars = reduce(vcat, variables(:x, 1:2, 1:2, 1:2, 1:2, 1:2))
F85qubit = tmp;

ctmp = @time coefficients(F85qubit, xvars);
size(ctmp)
# can you make a permutation ring map?

# can you solve for the lines on one of these hypersurfaces? That could explain the null cone?

P = [1,1,1,1,4,4]
F = prod([prod(det.(list2mats(syt6[P[i]],varSymbols[i]))) for i=1:6])
tmp = F; @time for dg=1:6
    tmp = sum([diffList(tmp, [a[i1,dg],b[i2,dg],c[i3,dg],d[i4,dg],e[i5,dg],f[i6,dg]])*x[i1,i2,i3,i4,i5,i6] for i1=1:2, i2=1:2, i3=1:2, i4=1:2, i5=1:2, i6=1:2]);
end
# 15s on server, 4.8s on mbp

P = [1,1,1,1,2,2];
F = prod([prod(det.(list2mats(syt8[P[i]],varSymbols[i]))) for i=1:6]);
tmp = F; @time for dg=1:8
    tmp = sum([diffList(tmp, [a[i1,dg],b[i2,dg],c[i3,dg],d[i4,dg],e[i5,dg],f[i6,dg]])*x[i1,i2,i3,i4,i5,i6] for i1=1:2, i2=1:2, i3=1:2, i4=1:2, i5=1:2, i6=1:2]);
end
# ?s on server, 9.6s on mbp
using Combinatorics
permList = unique(collect(Combinatorics.permutations(P)))
I8q=[0*x[1,1,1,1,1,1] for i=1:15 ]
@time for q=1:15
    I8q[q] = prod([prod(det.(list2mats(syt8[permList[q][i]],varSymbols[i]))) for i=1:6]);
    for dg=1:8
    I8q[q] = sum([diffList(I8q[q], [a[i1,dg],b[i2,dg],c[i3,dg],d[i4,dg],e[i5,dg],f[i6,dg]])*x[i1,i2,i3,i4,i5,i6] for i1=1:2, i2=1:2, i3=1:2, i4=1:2, i5=1:2, i6=1:2]);
    end
end
# on mbp 169.533696 seconds (63.96 M allocations: 981.951 MiB, 0.40% gc time, 0.03% compilation time)

J8q=[0*x[1,1,1,1,1,1] for i=1:15 ]
@time @threads for q=1:15
    J8q[q] = prod([prod(det.(list2mats(syt8[permList[q][i]],varSymbols[i]))) for i=1:6]);
    for dg=1:8
    J8q[q] = sum([diffList(J8q[q], [a[i1,dg],b[i2,dg],c[i3,dg],d[i4,dg],e[i5,dg],f[i6,dg]])*x[i1,i2,i3,i4,i5,i6] for i1=1:2, i2=1:2, i3=1:2, i4=1:2, i5=1:2, i6=1:2]);
    end
end
# on server: used 14-15 threads most of the time
# 266.413723 seconds (140.79 M allocations: 3.173 GiB, 0.14% gc time, 8.34% compilation time)
# on mbp: used 10 threads most of the time 
# 207.850514 seconds (64.02 M allocations: 985.855 MiB, 0.30% compilation time)

permList = unique(collect(Combinatorics.permutations(P)))
K8q=[0*x[1,1,1,1,1,1] for i=1:15 ]
holder = [[0*x[1,1,1,1,1,1] for i=1:64] for q=1:15]
sumList = reshape([[i1,i2,i3,i4,i5,i6] for i1=1:2, i2=1:2, i3=1:2, i4=1:2, i5=1:2, i6=1:2],1,64);
@time @threads for q=1:15
    K8q[q] = prod([prod(det.(list2mats(syt8[permList[q][i]],varSymbols[i]))) for i=1:6]);
    for dg=1:8
        @threads for i=1:64
            holder[q][i] = diffList(K8q[q], [a[sumList[i][1],dg],b[sumList[i][2],dg],c[sumList[i][3],dg],d[sumList[i][4],dg],e[sumList[i][5],dg],f[sumList[i][6],dg]])*x[CartesianIndex(sumList[i]...)];
        end
        K8q[q] = sum(holder[q]);
    end
end
# on server 225.461205 seconds (140.48 M allocations: 3.176 GiB, 0.33% gc time, 5.92% compilation time)
# server used 90 cores for most of the computation
# on MBP 337.558598 seconds (64.04 M allocations: 987.182 MiB, 0.06% gc time, 0.17% compilation time)
# used all cores most of the time


P = [1,1,1,1,4,4]
F = prod([prod(det.(list2mats(syt6[P[i]],varSymbols[i]))) for i=1:6])
sumList = reshape([[i1,i2,i3,i4,i5,i6] for i1=1:2, i2=1:2, i3=1:2, i4=1:2, i5=1:2, i6=1:2],1,64);
tmp = F; @time for dg=1:6
    holder = [0*x[1,1,1,1,1,1] for i=1:64]
    for i=1:64
    holder[i] = diffList(tmp, [a[sumList[i][1],dg],b[sumList[i][2],dg],c[sumList[i][3],dg],d[sumList[i][4],dg],e[sumList[i][5],dg],f[sumList[i][6],dg]])*x[CartesianIndex(sumList[i]...)];
    end
    tmp = sum(holder);
end

# 9.4s on mbp
P = [1,1,1,1,4,4]
F = prod([prod(det.(list2mats(syt6[P[i]],varSymbols[i]))) for i=1:6])
sumList = reshape([[i1,i2,i3,i4,i5,i6] for i1=1:2, i2=1:2, i3=1:2, i4=1:2, i5=1:2, i6=1:2],1,64);
tmp = F; @time for dg=1:6
    holder = [0*x[1,1,1,1,1,1] for i=1:64]
    @threads for i=1:64
    holder[i] = diffList(tmp, [a[sumList[i][1],dg],b[sumList[i][2],dg],c[sumList[i][3],dg],d[sumList[i][4],dg],e[sumList[i][5],dg],f[sumList[i][6],dg]])*x[CartesianIndex(sumList[i]...)];
    end
    tmp = sum(holder);
end
 #10.537976 s on server in parallel, 16.861424 seconds in serial

P = [1,1,1,1,2,2]
F = prod([prod(det.(list2mats(syt8[P[i]],varSymbols[i]))) for i=1:6])
sumList = reshape([[i1,i2,i3,i4,i5,i6] for i1=1:2, i2=1:2, i3=1:2, i4=1:2, i5=1:2, i6=1:2],1,64);
tmp = F; @time for dg=1:8
    tmp = sum([diffList(tmp, [a[sumList[i][1],dg],b[sumList[i][2],dg],c[sumList[i][3],dg],d[sumList[i][4],dg],e[sumList[i][5],dg],f[sumList[i][6],dg]])*x[CartesianIndex(sumList[i]...)] for i=1:64]);
end
#8.53s on mbp 
#27.5s on server
# on server
@time F8 = expand(tmp);
#6.3s on mbp

P = [1,1,1,1,2,2]
F = prod([prod(det.(list2mats(syt8[P[i]],varSymbols[i]))) for i=1:6])
sumList = reshape([[i1,i2,i3,i4,i5,i6] for i1=1:2, i2=1:2, i3=1:2, i4=1:2, i5=1:2, i6=1:2],1,64);
tmp = F; @time for dg=1:8
    holder = [0*x[1,1,1,1,1,1] for i=1:64]
    for i=1:64
    holder[i] = diffList(tmp, [a[sumList[i][1],dg],b[sumList[i][2],dg],c[sumList[i][3],dg],d[sumList[i][4],dg],e[sumList[i][5],dg],f[sumList[i][6],dg]])*x[CartesianIndex(sumList[i]...)];
    end
    tmp = sum(holder);
end
#29s on server serial
# 8s on mbp serial

tmp = F; @time for dg=1:8
    holder = [0*x[1,1,1,1,1,1] for i=1:64]
    @threads for i=1:64
    holder[i] = diffList(tmp, [a[sumList[i][1],dg],b[sumList[i][2],dg],c[sumList[i][3],dg],d[sumList[i][4],dg],e[sumList[i][5],dg],f[sumList[i][6],dg]])*x[CartesianIndex(sumList[i]...)];
    end
    tmp = sum(holder);
end;

# 8.161510 seconds on server parallel
# 14s on mbp parallel

# using Oscar
# R, x = polynomial_ring(QQ, :x => (1:2, 1:2, 1:2, 1:2, 1:2, 1:2))
# x[1,1,1,1,1,1]
# I = ideal(R, I8q[1]);
# # note that OSCAR automatically expands expressions... 
# (x[1,1,1,1,1,1]*x[1,1,1,1,1,2] - x[1,1,1,2,1,1])*x[1,1,2,1,1,1]
# # using Oscar
# # R = ExpressionRing(QQ)