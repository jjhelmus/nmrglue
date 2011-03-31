#! /usr/bin/env python
# List of tests.


suite = {}

#########################
# Apodization functions #
#########################
apod_list = ["apod1.dat","apod1.glue","apod2.dat","apod2.glue",
             "apod3.dat","apod3.glue","apod4.dat","apod4.glue",
             "apod5.dat","apod5.glue","apod6.dat","apod6.glue",
             "apod7.dat","apod7.glue"]
apod = {"py":"apod.py","sh":"apod.com","f":apod_list}
suite["apod"] = apod

em_list = ["em.dat","em.glue","em2.dat","em2.glue"]
em  = {"py":"em.py","sh":"em.com","f":em_list}
suite["em"] = em

gm_list = ["gm.dat","gm.glue","gm2.dat","gm2.glue","gm3.dat","gm3.glue"]
gm = {"py":"gm.py","sh":"gm.com","f":gm_list}
suite["gm"] = gm

gmb_list = ["gmb.dat","gmb.glue","gmb2.dat","gmb2.glue"]
gmb = {"py":"gmb.py","sh":"gmb.com","f":gmb_list}
suite["gmb"] = gmb

jmod_list = ["jmod.dat","jmod.glue","jmod2.dat","jmod2.glue"]
jmod= {"py":"jmod.py","sh":"jmod.com","f":jmod_list}
suite["jmod"] = jmod

sp = {"py":"sp.py","sh":"sp.com","f":["sp.dat","sp.glue","sp2.dat","sp2.glue"]}
suite["sp"] = sp

sine_list = ["sine.dat","sine.glue","sine2.dat","sine2.glue"]
sine = {"py":"sine.py","sh":"sine.com","f":sine_list}
suite["sine"] = sine

tm_list = ["tm.dat","tm.glue","tm2.dat","tm2.glue"]
tm  = {"py":"tm.py","sh":"tm.com","f":tm_list}
suite["tm"] = tm

tri_list = ["tri.dat","tri.glue","tri2.dat","tri2.glue"]
tri = {"py":"tri.py","sh":"tri.com","f":tri_list}
suite["tri"] = tri

###################
# Shift functions #
###################
rs_list = ["rs1.dat","rs1.glue","rs2.dat","rs2.glue","rs3.dat","rs3.glue",
           "rs4.dat","rs4.glue","rs5.dat","rs5.glue"]
rs = {"py":"rs.py","sh":"rs.com","f":rs_list}
suite["rs"] = rs

ls_list = ["ls1.dat","ls1.glue","ls2.dat","ls2.glue","ls3.dat","ls3.glue",
           "ls4.dat","ls4.glue","ls5.dat","ls5.glue"]
ls = {"py":"ls.py","sh":"ls.com","f":ls_list}
suite["ls"] = ls

cs_list = ["cs1.dat","cs1.glue","cs2.dat","cs2.glue","cs3.dat","cs3.glue",
           "cs4.dat","cs4.glue","cs5.dat","cs5.glue","cs6.dat","cs6.glue",
           "cs7.dat","cs7.glue","cs8.dat","cs8.glue"]
cs = {"py":"cs.py","sh":"cs.com","f":cs_list}
suite["cs"] = cs

# fsh does not have a test which passes, see broken for examples why.
# fsh_list= []
# fsh = {"py":"fsh.py","sh":"fsh.com","f":fsh_list}
# suite["fsh"] = fsh

##############
# Transforms #
##############
ft_list = ["ft1.dat","ft1.glue","ft2.dat","ft2.glue","ft3.dat","ft3.glue",
           "ft4.dat","ft4.glue","ft5.dat","ft5.glue","ft6.dat","ft6.glue",
           "ft7.dat","ft7.glue","ft8.dat","ft8.glue"]
ft = {"py":"ft.py","sh":"ft.com","f":ft_list}
suite["ft"] = ft

# a few of these fail...see broken
rft_list = ["rft1.dat","rft1.glue","rft2.dat","rft2.glue",
            "rft3.dat","rft3.glue","rft4.dat","rft4.glue",
            "rft5.dat","rft5.glue","rft6.dat","rft6.glue",
            "rft7.dat","rft7.glue","rft8.dat","rft8.glue",
                                   "rft12.dat","rft12.glue",
            "rft13.dat","rft13.glue","rft14.dat","rft14.glue"]
rft = {"py":"rft.py","sh":"rft.com","f":rft_list}
suite["rft"] = rft

ha_list = ["ha1.dat","ha1.glue","ha2.dat","ha2.glue"]
ha = {"py":"ha.py","sh":"ha.com","f":ha_list}
suite["ha"] = ha

# ps90-180 mode fails..see broken
ht_list = ["ht1.dat","ht1.glue","ht2.dat","ht2.glue","ht3.dat","ht3.glue",
           "ht5.dat","ht5.glue","ht6.dat","ht6.glue",
           "ht7.dat","ht7.glue","ht8.dat","ht8.glue"]
ht = {"py":"ht.py","sh":"ht.com","f":ht_list}
suite["ht"] = ht


##########################
# Standard NMR Functions #
##########################

ps_list = ["ps1.dat","ps1.glue","ps2.dat","ps2.glue","ps3.dat","ps3.glue",
           "ps4.dat","ps4.glue","ps5.dat","ps5.glue","ps6.dat","ps6.glue"]
ps = {"py":"ps.py","sh":"ps.com","f":ps_list}
suite["ps"] = ps

tp_list = ["tp.dat","tp.glue","tp2.dat","tp2.glue","tp3.dat","tp3.glue"]
tp = {"py":"tp.py","sh":"tp.com","f":tp_list}
suite["tp"] = tp

ytp_list = ["ytp.dat","ytp.glue","ytp2.dat","ytp2.glue","ytp3.dat","ytp3.glue"]
ytp = {"py":"ytp.py","sh":"ytp.com","f":ytp_list}
suite["ytp"] = ytp

xy2yx_list = ["xy.dat","xy.glue","xy2.dat","xy2.glue","xy3.dat","xy3.glue"]
xy2yx = {"py":"xy2yx.py","sh":"xy2yx.com","f":xy2yx_list}
suite["xy2yx"] = xy2yx

zf_list = ["zf.dat","zf.glue","zf2.dat","zf2.glue","zf3.dat","zf3.glue",
           "zf4.dat","zf4.glue","zf5.dat","zf5.glue"]
zf = {"py":"zf.py","sh":"zf.com","f":zf_list}
suite["zf"] = zf


###################
# Basic Utilities #
###################

# when multiple parameter provided produces different results...see broken
add_list = ["add1.dat","add1.glue","add2.dat","add2.glue",
            "add3.dat","add3.glue","add4.dat","add4.glue"]
add = {"py":"add.py","sh":"add.com","f":add_list}
suite["add"] = add

dx_list = ["dx.dat","dx.glue"]
dx = {"py":"dx.py","sh":"dx.com","f":dx_list}
suite["dx"] = dx

ext_list = ["ext1.dat","ext1.glue","ext2.dat","ext2.glue",
            "ext3.dat","ext3.glue","ext4.dat","ext4.glue",
            "ext5.dat","ext5.glue","ext6.dat","ext6.glue",
            "ext7.dat","ext7.glue"]
ext = {"py":"ext.py","sh":"ext.com","f":ext_list}
suite["ext"] = ext

integ_list = ["integ.dat","integ.glue"]
integ = {"py":"integ.py","sh":"integ.com","f":integ_list}
suite["integ"] = integ

mc_list = ["mc.dat","mc.glue","mc2.dat","mc2.glue"]
mc = {"py":"mc.py","sh":"mc.com","f":mc_list}
suite["mc"] = mc

mir_list = ["mir1.dat","mir1.glue","mir2.dat","mir2.glue",
            "mir3.dat","mir3.glue","mir4.dat","mir4.glue",
            "mir5.dat","mir5.glue","mir6.dat","mir6.glue",
            "mir7.dat","mir7.glue","mir8.dat","mir8.glue",
            "mir9.dat","mir9.glue","mir10.dat","mir10.glue",
            "mir11.dat","mir11.glue","mir12.dat","mir12.glue",
            "mir13.dat","mir13.glue","mir14.dat","mir14.glue"]
mir = {"py":"mir.py","sh":"mir.com","f":mir_list}
suite["mir"] = mir

# when multiple parameter provided produces different results...see broken
mult_list = ["mult1.dat","mult1.glue","mult2.dat","mult2.glue",
             "mult3.dat","mult3.glue"]
mult = {"py":"mult.py","sh":"mult.com","f":mult_list}
suite["mult"] = mult

rev_list = ["rev1.dat","rev1.glue","rev2.dat","rev2.glue",
            "rev3.dat","rev3.glue"]
rev = {"py":"rev.py","sh":"rev.com","f":rev_list}
suite["rev"] = rev

set_list = ["set1.dat","set1.glue","set2.dat","set2.glue",
            "set3.dat","set3.glue","set4.dat","set4.glue"]
set = {"py":"set.py","sh":"set.com","f":set_list}
suite["set"] = set

# insignificant differences and no int mode...see broken
shuf_list = ["s1.dat","s1.glue","s2.dat","s2.glue","s3.dat","s3.glue",
             "s4.dat","s4.glue","s5.dat","s5.glue","s6.dat","s6.glue",   
             "s7.dat","s7.glue"]
shuf = {"py":"shuf.py","sh":"shuf.com","f":shuf_list}
suite["shuf"] = shuf

sign_list = ["s1.dat","s1.glue","s2.dat","s2.glue","s3.dat","s3.glue",
             "s4.dat","s4.glue","s5.dat","s5.glue","s6.dat","s6.glue",
             "s7.dat","s7.glue","s8.dat","s8.glue"]
sign = {"py":"sign.py","sh":"sign.com","f":sign_list}
suite["sign"] = sign

########
# Misc #
########

coadd_list = ["ca1.dat","ca1.glue","ca2.dat","ca2.glue"]
coadd = {"py":"coadd.py","sh":"coadd.com","f":coadd_list}
suite["coadd"] = coadd

coad_list = ["co1.dat","co1.glue","co2.dat","co2.glue"]
coad = {"py":"coad.py","sh":"coad.com","f":coad_list}
suite["coad"] = coad

dev_list = []
dev = {"py":"dev.py","sh":"dev.com","f":dev_list}
suite["dev"] = dev

# XXX make tests
#img_list = []
#img = {"py":"","sh":"","f":img_list}
#suite["img"] = img

null_list = ["null.dat","null.glue","null2.dat","null2.glue"]
null = {"py":"null.py","sh":"null.com","f":null_list}
suite["null"] = null

qart_list = ["qart.dat","qart.glue","qart2.dat","qart2.glue"]
qart = {"py":"qart.py","sh":"qart.com","f":qart_list}
suite["qart"] = qart

qmix_list = ["qmix.dat","qmix.glue","qmix2.dat","qmix2.glue"]
qmix = {"py":"qmix.py","sh":"qmix.com","f":qmix_list}
suite["qmix"] = qmix

#save_list = ["save.dat","save.glue","save2.dat","save2.glue"]
#save = {"py":"save.py","sh":"save.com","f":save_list}
#suite["save"] = save

smo_list = ["smo.dat","smo.glue","smo2.dat","smo2.glue","smo3.dat","smo3.glue"]
smo = {"py":"smo.py","sh":"smo.com","f":smo_list}
suite["smo"] = smo

zd_list = ["zd.dat","zd.glue","zd2.dat","zd2.glue","zd3.dat","zd3.glue","zd4.dat","zd4.glue"]
zd = {"py":"zd.py","sh":"zd.com","f":zd_list}
suite["zd"] = zd

######################
# Baseline functions #
######################

base_list = ["base1.dat","base1.glue","base2.dat","base2.glue",
             "base3.dat","base3.glue","base4.dat","base4.glue",
             "base5.dat","base5.glue","base6.dat","base6.glue",
             "base7.dat","base7.glue"]
base = {"py":"base.py","sh":"base.com","f":base_list}
suite["base"] = base

cbf_list = ["cbf1.dat","cbf1.glue","cbf2.dat","cbf2.glue",
            "cbf3.dat","cbf3.glue","cbf4.dat","cbf4.glue"]
cbf = {"py":"cbf.py","sh":"cbf.com","f":cbf_list}
suite["cbf"] = cbf

# XXX make tests
#med_list = []
#med = {"py":"","sh":"","f":med_list}
#suite["med"] = med

# XXX make tests
#sol_list = []
#sol = {"py":"","sh":"","f":sol_list}
#suite["sol"] = sol

# Additional code tests

units_list = ["u1.dat","u1.glue","u2.dat","u2.glue","u3.dat","u3.glue",
              "u4.dat","u4.glue","u5.dat","u5.glue","u6.dat","u6.glue",
              "u7.dat","u7.glue","u8.dat","u8.glue","u9.dat","u9.glue",
              "u10.dat","u10.glue","u11.dat","u11.glue","u12.dat","u12.glue",
              "u13.dat","u13.glue","u14.dat","u14.glue","u15.dat","u15.glue",
              "u16.dat","u16.glue","u17.dat","u17.glue"]

units = {"py":"units.py","sh":"units.com","f":units_list}
suite["units"] = units

# This is a collection of tests which fail for one reason or another.

broken_list = ["fsh1.dat","fsh1.glue","fsh2.dat","fsh2.glue",
               "fsh3.dat","fsh3.glue","fsh4.dat","fsh4.glue" ,
               "rft9.dat","rft9.glue","rft10.dat","rft10.glue",
               "rft11.dat","rft11.glue",
               "ht4.dat","ht4.glue",
               "add5.dat","add5.glue",
               "mult4.dat","mult4.glue",
               "s8.dat","s8.glue","s9.dat","s9.glue","s10.dat","s10.glue"]
broken = {"py":"broken.py","sh":"broken.com","f":broken_list}
