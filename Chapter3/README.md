2023.09.20

%%%%%%%
第 3 章
%%%%%%%


図3.5
meshSquare2(1,10)
square2.dat のファイル名を square2010.dat に変更
plotTriangulation("square2010.dat")

meshSquare2(1,20)
square2.dat のファイル名をsquare2020.dat に変更
plotTriangulation("square2020.dat")



図3.6
はじめに、Freefem++を利用して、
% freefem++ mk_domain2.edp 5
mk_domain2.dat のファイル名を mk_domain2_5.dat に変更
% freefem++ mk_domain2.edp 10   
mk_domain2.dat のファイル名を mk_domain2_510dat に変更

その後MTALABで、以下を行う
plotTriangulation("mk_domain2_5.dat")
plotTriangulation("mk_domain2_10.dat")


図3.7
はじめに、Freefem++を利用して、
% freefem++ mk_domain1.edp 10
mk_domain1.dat のファイル名を mk_domain1_10.dat に変更
% freefem++ mk_domain1.edp 20
mk_domain1.dat のファイル名を mk_domain1_20.dat に変更

その後MTALABで、以下を行う
plotTriangulation("mk_domain1_10.dat")
plotTriangulation("mk_domain1_20.dat")


図3.9
p1fem1("square2010.dat",@femFunc2)
p1fem1("square2020.dat",@femFunc2)

Remark: dirichlet1.m, getmesh.m, matrix1.m, P1grad.m, plotP1fem.m, robin1.m が必要


図3.11
meshSquare2(1,40)
square2.dat のファイル名をsquare2040.dat に変更
meshSquare2(1,80)
square2.dat のファイル名をsquare2080.dat に変更
meshSquare2(1,160)
square2.dat のファイル名をsquare2160.dat に変更

p1fem2error

Remark: dirichlet1.m, getmesh.m, matrix1.m, P1grad.m, plotP1fem.m, robin1.m, gauss2D.m, p1error.m, p1fem2.m が必要


図3.14
はじめに、Freefem++を利用して、
% freefem++ nonconvex.edp 0.7
nonconvex.dat のファイル名を nonconvex1.dat に変更
% freefem++ nonconvex.edp 1.5
nonconvex.dat のファイル名を nonconvex2.dat に変更
% freefem++ nonconvex.edp 1.8
nonconvex.dat のファイル名を nonconvex3.dat に変更

その後MATLABで、以下を行う
convex_error


%%%%%%
付録 C
%%%%%%


図C.4
meshSquare1(1, 5)
square1.dat のファイル名を square1005.dat に変更
plotTriangulation("square1005.dat")

meshSquare1(1, 20)
square1.dat のファイル名を square1020.dat に変更
plotTriangulation("square1020.dat")

Remark: getmesh.m が必要


図C.5
p1fem1("square1005.dat", @femFunc1)
p1fem1("square1020.dat", @femFunc1)

Remark: dirichlet1.m, getmesh.m, matrix1.m, P1grad.m, plotP1fem.m, robin1.m が必要
