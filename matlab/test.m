libname = 'anode';
[notfound,warnings]=loadlibrary(libname, 'anode_wrap.h');

z = 64
e = 10
v = 100
i = 6
a = 26
tic
out = calllib(libname, 'c_anode_int', z, e, v, i, a);
toc
