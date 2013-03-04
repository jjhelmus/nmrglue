import nmrglue as ng

# process the direct dimension of the sum data set
sdic, sdata = ng.pipe.read('test_sum.fid')
sdic, sdata = ng.pipe_proc.sp(sdic, sdata, off=0.45, end=0.95, pow=1, c=1.0)
sdic, sdata = ng.pipe_proc.zf(sdic, sdata, size=8192)
sdic, sdata = ng.pipe_proc.ft(sdic, sdata)
uc = ng.pipe.make_uc(sdic, sdata, dim=1)
pts = uc.f('27.5 Hz') - uc.f('0 Hz')
sdic, sdata = ng.pipe_proc.fsh(sdic, sdata, dir='ls', pts=pts)
sdic, sdata = ng.pipe_proc.ps(sdic, sdata, p0=-79.0, p1=0.0)
sdic, sdata = ng.pipe_proc.di(sdic, sdata)

# process the direct dimension of the difference data set
ddic, ddata = ng.pipe.read('test_dif.fid')
ddic, ddata = ng.pipe_proc.sp(ddic, ddata, off=0.45, end=0.95, pow=1, c=1.0)
ddic, ddata = ng.pipe_proc.zf(ddic, ddata, size=8192)
ddic, ddata = ng.pipe_proc.ft(ddic, ddata)
ddic, ddata = ng.pipe_proc.ps(ddic, ddata, p0=-90.0, p1=0.0)
uc = ng.pipe.make_uc(ddic, ddata, dim=1)
pts = uc.f('27.5 Hz') - uc.f('0 Hz')
ddic, ddata = ng.pipe_proc.fsh(ddic, ddata, dir='rs', pts=pts)
ddic, ddata = ng.pipe_proc.ps(ddic, ddata, p0=-79.0, p1=0.0)
ddic, ddata = ng.pipe_proc.di(ddic, ddata)

# sum the different and sum data sets
data = sdata + ddata
dic = ddic

# process the indirect dimension
dic, data = ng.pipe_proc.tp(dic, data)
dic, data = ng.pipe_proc.sp(dic, data, off=0.45, end=0.95, pow=1, c=1.0)
dic, data = ng.pipe_proc.zf(dic, data, size=2048)
dic, data = ng.pipe_proc.ft(dic, data, neg=True)
dic, data = ng.pipe_proc.ps(dic, data, p0=0.0, p1=0.0)
dic, data = ng.pipe_proc.di(dic, data)
dic, data = ng.pipe_proc.tp(dic, data)

# write out the results
ng.pipe.write('test.ft2', dic, data, overwrite=True)
