2023.09.20

%%%%%%%
第 8 章
%%%%%%%

図8.1

cd1d_std(5, 20)

cd1d_std(5, 40)


図8.2

cd1d_uwd(5,20)

cd1d_uwd(5,40)


図8.3

meshSquare2(1,20)

square2.dat のファイル名をsquare2020.dat に変更

p1femStab1("square2020.dat", @femFunc3, 0.24, 0.01, 0, 0)

p1femStab1("square2020.dat", @femFunc3, 0.12, 0.01, 0, 0)

p1femStab1("square2020.dat", @femFunc3, 0.06, 0.01, 0, 0)

p1femStab1("square2020.dat", @femFunc3, 0.03, 0.01, 0,0)

Remark: dirichlet1.m, getmesh.m, matrixStab1.m, P1grad.m, plotP1fem.m が必要



図8.3

p1femStab1("square2020.dat", @femFunc3, 0.12, 0.01, 0.8, 0)

p1femStab1("square2020.dat", @femFunc3, 0.06, 0.01, 0.8, 0)

p1femStab1("square2020.dat", @femFunc3, 0.03, 0.01, 0.8, 0)

p1femStab1("square2020.dat", @femFunc3, 0.015, 0.01, 0.8, 0)

Remark: dirichlet1.m, getmesh.m, matrixStab1.m, P1grad.m, plotP1fem.m が必要



図8.8

meshSquare4(1,40)

p1femBT("square4.dat",@femFunc4,1)

Remark: dirichlet1.m, getmesh.m, matrixBT.m, P1grad.m, plotP1fem.m が必要


