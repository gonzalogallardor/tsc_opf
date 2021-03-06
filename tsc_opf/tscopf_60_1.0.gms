Scalars

Sb      /100.0/
Dt      /0.02/
angle_limit     /60/
w0      /314.1592653589793/


Sets

s       /1 * 100/
sfirst(s)
sf(s)   /2 * 11/
spf(s)  /12 * 100/
sdfig_1(s)      /2 * 19/
sdfig_2(s)      /20 * 100/
spitch_1(s)     /1 * 19/
spitch_2(s)     /19 * 100/
b       /1 * 20/
rb(b)   /1 * 6/
sg(rb)  /1, 2, 3, 4, 5/
dfig(rb)        /6/
ngb(b)  /7 * 20/
i       /1 * 5/;

sfirst(s) = yes$(ord(s) eq 1);
alias (sg,sgp);
alias (dfig,dfigp);
alias (i,j);


Parameters

a1(rb)  /1 70.0, 2 80.0, 3 40.0, 4 100.0, 5 120.0, 6 0.0/
Pl(b)   /set.b 0/
Ql(b)   /set.b 0/
pmax(sg)        /1 600.0, 2 470.0, 3 510.0, 4 275.0, 5 350.0/
pmin(sg)        /1 0.0, 2 0.0, 3 0.0, 4 0.0, 5 0.0/
qmax(sg)        /1 320.0, 2 250.0, 3 270.0, 4 150.0, 5 185.0/
qmin(sg)        /1 -320.0, 2 -250.0, 3 -270.0, 4 -150.0, 5 -185.0/
e_fd_max(sg)    /1 2.0, 2 2.0, 3 2.0, 4 2.0, 5 2.0/
e_fd_min(sg)    /1 0.0, 2 0.0, 3 0.0, 4 0.0, 5 0.0/
ra(sg)  /1 0.0, 2 0.0, 3 0.0, 4 0.0, 5 0.0/
xd(sg)  /1 1.5, 2 1.5, 3 1.5, 4 1.5, 5 1.5/
xd_p(sg)        /1 0.3, 2 0.3, 3 0.3, 4 0.3, 5 0.3/
xq(sg)  /1 1.5, 2 1.5, 3 1.5, 4 1.5, 5 1.5/
xq_p(sg)        /1 0.3, 2 0.3, 3 0.3, 4 0.3, 5 0.3/
H(sg)   /1 3.2, 2 3.0, 3 3.0, 4 2.0, 5 2.0/
D(sg)   /1 2.0, 2 2.0, 3 2.0, 4 2.0, 5 2.0/
Td_p(sg)        /1 6.0, 2 6.0, 3 6.0, 4 6.0, 5 6.0/
Tq_p(sg)        /1 1.0, 2 1.0, 3 1.0, 4 1.0, 5 1.0/
Ngen_w(dfig)    /6 30/
u_wind(dfig)    /6 11.0/
Sb_w(dfig)      /6 1.5/
w0_w(dfig)      /6 377.0/
To_p(dfig)      /6 1.5/
xss(dfig)       /6 4.0/
xs_p(dfig)      /6 0.35/
rs(dfig)        /6 0.01/
xl_s(dfig)      /6 0.15/
H_w(dfig)       /6 3.0/
kc_crowbar(dfig)        /6 10.0/
initial_pitch_angle(dfig)       /6 0.0/
Kb(dfig)        /6 56.6/
Kwind(dfig)     /6 0.00159/
T_pc(dfig)      /6 0.05/
Kptrq(dfig)     /6 3.0/
Kitrq(dfig)     /6 0.6/
wr_ref(dfig)    /6 1.2/
T_p(dfig)       /6 0.3/

Sm(sg)
lm(dfig)
lrr(dfig)
rr(dfig)
r_crowbar(dfig)
xm(dfig)
xrr(dfig)
xss(dfig)
lss(dfig)
rr_1(dfig)
rr_2(dfig)
To_p_1(dfig)
To_p_2(dfig);
Pmax(sg) = Pmax(sg)/Sb;
Pmin(sg) = Pmin(sg)/Sb;
Qmax(sg) = Qmax(sg)/Sb;
Qmin(sg) = Qmin(sg)/Sb;
Sm(sg) = 1.1*Pmax(sg);
Ra(sg) = Ra(sg)/Sm(sg);
Xd(sg) = Xd(sg)/Sm(sg);
Xd_p(sg) = Xd_p(sg)/Sm(sg);
Xq(sg) = Xq(sg)/Sm(sg);
Xq_p(sg) = Xq_p(sg)/Sm(sg);
H(sg) = H(sg)*Sm(sg);
D(sg) = D(sg)*Sm(sg);
xss(dfig) = xss(dfig)*Sb/Sb_w(dfig)*w0/w0_w(dfig);
xs_p(dfig) = xs_p(dfig)*Sb/Sb_w(dfig)*w0/w0_w(dfig);
rs(dfig) = rs(dfig)*Sb/Sb_w(dfig);
xl_s(dfig) = xl_s(dfig)*Sb/Sb_w(dfig)*w0/w0_w(dfig);
lss(dfig) = xss(dfig);
H_w(dfig) = H_w(dfig)*Sb_w(dfig)/Sb;
xm(dfig) = xss(dfig)-xl_s(dfig);
lm(dfig) = xm(dfig);
lrr(dfig) = sqr(lm(dfig))/(xss(dfig)-xs_p(dfig));
xrr(dfig) = lrr(dfig);
rr(dfig) =  lrr(dfig)/(w0*To_p(dfig));
r_crowbar(dfig) = kc_crowbar(dfig)*rr(dfig);
rr_1(dfig) = rr(dfig) + r_crowbar(dfig);
rr_2(dfig) = rr(dfig);
To_p_1(dfig) = To_p(dfig)*rr(dfig)/(rr(dfig) + r_crowbar(dfig));
To_p_2(dfig) = To_p(dfig);
Kptrq(dfig) = Kptrq(dfig)*Sb_w(dfig)/Sb;
Kitrq(dfig) = Kitrq(dfig)*Sb_w(dfig)/Sb;


Table Y(b,b)
                 1          2          3          4          5          6          7          8          9         10         11         12         13         14         15         16         17         18         19         20
         1   85.68980    0.00000    0.00000    0.00000    0.00000    0.00000   85.68980    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000
         2    0.00000   72.83321    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000   72.83321    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000
         3    0.00000    0.00000   72.83321    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000   72.83321    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000
         4    0.00000    0.00000    0.00000   39.29273    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000   39.29273    0.00000    0.00000    0.00000    0.00000
         5    0.00000    0.00000    0.00000    0.00000   39.29273    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000   39.29273    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000
         6    0.00000    0.00000    0.00000    0.00000    0.00000  210.73077  136.32822    0.00000    0.00000    0.00000    0.00000    0.00000   37.49531    0.00000    0.00000    0.00000   37.49531    0.00000    0.00000    0.00000
         7   85.68980    0.00000    0.00000    0.00000    0.00000  136.32822  370.15731   81.75117    0.00000   29.91482    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000   37.49531    0.00000    0.00000
         8    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000   81.75117  190.98423   36.07302   36.07302    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000   37.49531    0.00000
         9    0.00000   72.83321    0.00000    0.00000    0.00000    0.00000    0.00000   36.07302  108.60222    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000
        10    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000   29.91482   36.07302    0.00000  317.61714  122.62676  129.10867    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000
        11    0.00000    0.00000   72.83321    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000  122.62676  195.56923    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000
        12    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000  129.10867    0.00000  194.19042    0.00000    0.00000   37.49531    0.00000    0.00000    0.00000    0.00000   28.12148
        13    0.00000    0.00000    0.00000    0.00000    0.00000   37.49531    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000  170.41117  137.87858    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000
        14    0.00000    0.00000    0.00000    0.00000   39.29273    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000  137.87858  173.27246    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000
        15    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000   37.49531    0.00000    0.00000   90.36594   52.93321    0.00000    0.00000    0.00000    0.00000
        16    0.00000    0.00000    0.00000   39.29273    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000   52.93321   92.41534    0.00000    0.00000    0.00000    0.00000
        17    0.00000    0.00000    0.00000    0.00000    0.00000   37.49531    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000   63.12706   21.41957    5.23601    0.00000
        18    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000   37.49531    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000   21.41957   81.04284   23.56182    0.00000
        19    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000   37.49531    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    5.23601   23.56182   70.98174    5.92761
        20    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000   28.12148    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    5.92761   34.01624 ;


Table theta(b,b)
                 1          2          3          4          5          6          7          8          9         10         11         12         13         14         15         16         17         18         19         20
         1   -1.57080    0.00000    0.00000    0.00000    0.00000    0.00000    1.57080    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000
         2    0.00000   -1.57080    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    1.57080    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000
         3    0.00000    0.00000   -1.57080    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    1.57080    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000
         4    0.00000    0.00000    0.00000   -1.57080    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    1.57080    0.00000    0.00000    0.00000    0.00000
         5    0.00000    0.00000    0.00000    0.00000   -1.57080    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    1.57080    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000
         6    0.00000    0.00000    0.00000    0.00000    0.00000   -1.47230    1.72339    0.00000    0.00000    0.00000    0.00000    0.00000    1.57080    0.00000    0.00000    0.00000    1.57080    0.00000    0.00000    0.00000
         7    1.57080    0.00000    0.00000    0.00000    0.00000    1.72339   -1.46877    1.72345    0.00000    1.72335    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    1.57080    0.00000    0.00000
         8    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    1.72345   -1.44803    1.72325    1.72325    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    1.57080    0.00000
         9    0.00000    1.57080    0.00000    0.00000    0.00000    0.00000    0.00000    1.72325   -1.52033    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000
        10    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    1.72335    1.72325    0.00000   -1.41801    1.72345    1.72374    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000
        11    0.00000    0.00000    1.57080    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    1.72345   -1.46246    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000
        12    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    1.72374    0.00000   -1.46933    0.00000    0.00000    1.57080    0.00000    0.00000    0.00000    0.00000    1.57080
        13    0.00000    0.00000    0.00000    0.00000    0.00000    1.57080    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000   -1.10825    2.15497    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000
        14    0.00000    0.00000    0.00000    0.00000    1.57080    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    2.15497   -1.10294    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000
        15    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    1.57080    0.00000    0.00000   -1.53187    1.63728    0.00000    0.00000    0.00000    0.00000
        16    0.00000    0.00000    0.00000    1.57080    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    1.63728   -1.51974    0.00000    0.00000    0.00000    0.00000
        17    0.00000    0.00000    0.00000    0.00000    0.00000    1.57080    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000   -1.36284    2.01825    2.01829    0.00000
        18    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    1.57080    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    2.01825   -1.30912    2.01843    0.00000
        19    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    1.57080    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    2.01829    2.01843   -1.33578    2.01830
        20    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    1.57080    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000    2.01830   -1.44360 ;


Table Y_f_red(rb,rb)
                 1          2          3          4          5          6
         1   15.05975    0.43041    0.00005    0.00894    0.00000   10.59443
         2    0.43041    8.72996    0.00005    0.01190    0.00000    4.19605
         3    0.00005    0.00005   12.41474    0.00007    0.00000    0.00047
         4    0.00894    0.01190    0.00007    5.94103    0.00000    0.10952
         5    0.00000    0.00000    0.00000    0.00000    7.76923    7.60345
         6   10.59443    4.19605    0.00047    0.10952    7.60345   61.47865 ;


Table theta_f_red(rb,rb)
                 1          2          3          4          5          6
         1   -1.55932    1.45515    1.57320    1.46634    0.00000    1.57354
         2    1.45515   -1.52014    1.60633    1.44208    0.00000    1.61128
         3    1.57320    1.60633   -1.55356    1.52563    0.00000    1.73145
         4    1.46634    1.44208    1.52563   -1.53694    0.00000    1.68476
         5    0.00000    0.00000    0.00000    0.00000   -1.52404    1.55161
         6    1.57354    1.61128    1.73145    1.68476    1.55161   -1.39218 ;


Table Y_pf_red(rb,rb)
                 1          2          3          4          5          6
         1   14.78805    0.65085    0.56735    0.27177    0.00000   13.17002
         2    0.65085    8.46163    0.90953    0.43166    0.00000    6.31925
         3    0.56735    0.90953    9.17746    1.53947    0.00000    5.60238
         4    0.27177    0.43166    1.53947    5.26971    0.00000    2.70273
         5    0.00000    0.00000    0.00000    0.00000    7.76923    7.60345
         6   13.17002    6.31925    5.60238    2.70273    7.60345   37.37895 ;


Table theta_pf_red(rb,rb)
                 1          2          3          4          5          6
         1   -1.55622    1.44841    1.38410    1.34799    0.00000    1.57292
         2    1.44841   -1.51498    1.41636    1.37756    0.00000    1.60582
         3    1.38410    1.41636   -1.47685    1.33505    0.00000    1.54815
         4    1.34799    1.37756    1.33505   -1.49580    0.00000    1.51504
         5    0.00000    0.00000    0.00000    0.00000   -1.52404    1.55161
         6    1.57292    1.60582    1.54815    1.51504    1.55161   -1.37712 ;
Table
alpha_ij(i,j)
               1              2               3              4              5
1   -4.1909e-001    2.1808e-001    -1.2406e-002   -1.3365e-004    1.1524e-005
2   -6.7606e-002    6.0405e-002    -1.3934e-002    1.0683e-003   -2.3895e-005
3    1.5727e-002   -1.0996e-002     2.1495e-003   -1.4855e-004    2.7937e-006
4   -8.6018e-004    5.7051e-004    -1.0479e-004    5.9924e-006   -8.9194e-008
5    1.4787e-005   -9.4839e-006     1.6167e-006   -7.1535e-008    4.9686e-010;

Variables

Pg(rb)
z
Qg(rb)
V(b)
alpha(b)
Ps_w(dfig,s)
Qs_w(dfig,s)
P_gsc(dfig,s)
P_mech(dfig,s)
VDr(dfig,s)
VQr(dfig,s)
IDr(dfig,s)
IQr(dfig,s)
iDs_w(dfig,s)
iQs_w(dfig,s)
eD_p_w(dfig,s)
eQ_p_w(dfig,s)
Ig(sg)
phi(sg)
delta(sg,s)
ed_p(sg,s)
eq_p(sg,s)
id(sg,s)
iq(sg,s)
e_fd(sg)
Domega(sg,s)
Pe(sg,s)
delta_COI(s)
vD_w(dfig,s)
vQ_w(dfig,s)
rr_var(dfig,s)
To_p_var(dfig,s)
T_el(dfig,s)
T_mech(dfig,s)
wr(dfig,s)
delta_wr(dfig,s)
lambda(dfig,s)
Cp(dfig,s)
v_w(dfig,s)
alpha_w(dfig,s)
IDr_c(dfig,sdfig_2)
IQr_c(dfig,sdfig_2)
iD_w(dfig,s)
iQ_w(dfig,s)
wr_ref_var(dfig,s)
error_wr(dfig,s)
T_el_ref(dfig,s)
P_ord_pre(dfig,s)
P_ord(dfig,s)
pitch_angle(dfig,s)
pll_angle(dfig,s)
Ps_w_ref(dfig,sdfig_2)
iDr_c_ref(dfig,sdfig_2)
iQr_c_ref(dfig,sdfig_2)
Pg_w(dfig,s)
Qg_w(dfig,s);

Pg.l('1') = 300.0/Sb;
Pg.l('2') = 100.0/Sb;
Pg.l('3') = 300.0/Sb;
Pg.l('4') = 50.0/Sb;
Pg.l('5') = 50.0/Sb;
Pg.lo(sg) = Pmin(sg);
Pg.up(sg) = Pmax(sg);
Qg.lo(sg) = Qmin(sg);
Qg.up(sg) = Qmax(sg);
V.l(b) = 0.95;
V.lo(b) = 0.95;
V.up(b) = 1.05;
alpha.l(b) = 0;
alpha.lo(b) = -pi;
alpha.up(b) = pi;
Ig.l(sg) = 1;
Ig.lo(sg) = 0.001;
Ig.up(sg) = Sm(sg);
phi.l(sg) = 0;
phi.lo(sg) = -pi/2;
phi.up(sg) = pi/2;
Qg.l('1') = 0.0/Sb;
Qg.l('2') = 0.0/Sb;
Qg.l('3') = 0.0/Sb;
Qg.l('4') = 0.0/Sb;
Qg.l('5') = 0.0/Sb;
delta.l(sg,s)= 0;
delta.lo(sg,s) = -9999;
delta.up(sg,s) = 9999;
ed_p.l(sg,s) = 0.2;
ed_p.lo(sg,s) = 0;
ed_p.up(sg,s) = 1.5;
eq_p.l(sg,s) = 1.0;
eq_p.lo(sg,s) = 0;
eq_p.up(sg,s) = 1.5;
id.l(sg,s) = 1;
id.lo(sg,s) = -Sm(sg);
id.up(sg,s) = 3*Sm(sg);
iq.l(sg,s) = 1;
iq.lo(sg,s) = -Sm(sg);
iq.up(sg,s) = 3*Sm(sg);
e_fd.l(sg) = 1;
Domega.l(sg,s) = 0;
Domega.lo(sg,s) = -1;
Domega.up(sg,s) = 1;
Pe.l(sg,s) = 3.0;
Pe.lo(sg,s) = -99;
Pe.up(sg,s) = 99;
delta_COI.l(s) = 0;
delta_COI.lo(s) = -9999;
delta_COI.up(s) = 9999;
v_w.l(dfig,s) = 1;
v_w.lo(dfig,s) = 0.1;
rr_var.fx(dfig,'1') = rr(dfig);
rr_var.fx(dfig,sdfig_1) = rr_1(dfig);
rr_var.fx(dfig,sdfig_2) = rr_2(dfig);
To_p_var.fx(dfig,'1') = To_p(dfig);
To_p_var.fx(dfig,sdfig_1) = To_p_1(dfig);
To_p_var.fx(dfig,sdfig_2) = To_p_2(dfig);
wr_ref_var.fx(dfig,s) = wr_ref(dfig);
delta_wr.fx(dfig,'1') = wr_ref(dfig) - 1;
pitch_angle.fx(dfig,'1') = initial_pitch_angle(dfig);
vDr.fx(dfig,sdfig_1) = 0;
vQr.fx(dfig,sdfig_1) = 0;


Equations

total_cost
p_balance_gen
p_balance_nongen
q_balance_gen
q_balance_nongen
ref_bus
dfig_load_flow_f1
dfig_load_flow_f2
dfig_load_flow_f3
dfig_load_flow_f4
dfig_load_flow_f5
dfig_load_flow_f6
P_mech_ini
I_brach_lim_inf_1
I_brach_lim_inf_2
I_brach_lim_sup_1
I_brach_lim_sup_2
Generator_currents
Power_factor
ed_p_initialization
eq_p_initialization
Vd_initialization
Vq_initialization
Domega_initialization
id_initialization
iq_initialization
e_fd_lim_inf
e_fd_lim_sup
Internal_voltaje_d
Internal_voltaje_q
oscilation_omega
oscilation_delta
electric_power
center_of_inertia
angular_deviation_min
angular_deviation_max
id_stator_fault
iq_stator_fault
id_stator_postfault
iq_stator_postfault
iD_w_fault
iQ_w_fault
iD_w_postfault
iQ_w_postfault
eD_p_w_equation
eQ_p_w_equation
rotor_flux_D
rotor_flux_Q
stator_equation_D
stator_equation_Q
T_el_equation
delta_wr_equation
wr_equation
lambda_equation
Cp_equation
P_mech_equation
T_mech_equation
error_wr_equation
T_el_ref_ini
T_el_ref_equation
P_ord_pre_equation
P_ord_ini
P_ord_equation
pitch_angle_1
pitch_angle_2
vD_w_equation
vQ_w_equation
pll_angle_ini
Ps_w_ref_equation
iQr_c_ref_equation
IDr_c_ref_2
IDr_equation
IQr_equation
IDr_c_2
IQr_c_2
P_gsc_equation
Ps_w_equation
Qs_w_equation
Pg_w_equation
Qg_w_equation
P_balance
Q_balance;

total_cost .. z =e= sum(sg,a1(sg)*(Pg(sg)*Sb)) + sum(dfig,a1(dfig)*(Pg(dfig)*Sb));
p_balance_gen(rb) .. Pg(rb) - Pl(rb) - V(rb)*sum(b,V(b)*Y(rb,b)*cos(alpha(rb) - alpha(b) - theta(rb,b))) =e= 0;
p_balance_nongen(ngb) ..  - Pl(ngb) - V(ngb)*sum(b,V(b)*Y(ngb,b)*cos(alpha(ngb) - alpha(b) - theta(ngb,b))) =e= 0;
q_balance_gen(rb) .. Qg(rb) - Ql(rb) - V(rb)*sum(b,V(b)*Y(rb,b)*sin(alpha(rb) - alpha(b) - theta(rb,b))) =e= 0;
q_balance_nongen(ngb) ..  - Ql(ngb) - V(ngb)*sum(b,V(b)*Y(ngb,b)*sin(alpha(ngb) - alpha(b) - theta(ngb,b))) =e= 0;
ref_bus .. alpha('1') =e= 0;
dfig_load_flow_f1(dfig) .. Pg(dfig)/Ngen_w(dfig) =e= Ps_w(dfig,'1') - P_gsc(dfig,'1');
dfig_load_flow_f2(dfig) .. Qg(dfig)/Ngen_w(dfig) =e= Qs_w(dfig,'1');
dfig_load_flow_f3(dfig) .. VDr(dfig,'1') =e= IDr(dfig,'1') + w0*To_p(dfig)*delta_wr(dfig,'1')*eD_p_w(dfig,'1');
dfig_load_flow_f4(dfig) .. VQr(dfig,'1') =e= IQr(dfig,'1') - w0*To_p(dfig)*delta_wr(dfig,'1')*eQ_p_w(dfig,'1');
dfig_load_flow_f5(dfig) .. V(dfig)*cos(alpha(dfig)) =e= xs_p(dfig)*iQs_w(dfig,'1') - eD_p_w(dfig,'1');
dfig_load_flow_f6(dfig) .. V(dfig)*sin(alpha(dfig)) =e= - xs_p(dfig)*iDs_w(dfig,'1') + eQ_p_w(dfig,'1');
P_mech_ini(dfig) .. Pg(dfig)/Ngen_w(dfig) =e= P_mech(dfig,'1') - rr(dfig)/sqr(lm(dfig))*(sqr(IDr(dfig,'1')) + sqr(IQr(dfig,'1')));
I_brach_lim_inf_1('1') .. 0 =l= (sqr(V('13')*cos(alpha('13')) - V('14')*cos(alpha('14'))) + sqr(V('13')*sin(alpha('13')) - V('14')*sin(alpha('14'))))*sqr(Y('13','14'));
I_brach_lim_inf_2('1') .. 0 =l= (sqr(V('15')*cos(alpha('15')) - V('16')*cos(alpha('16'))) + sqr(V('15')*sin(alpha('15')) - V('16')*sin(alpha('16'))))*sqr(Y('15','16'));
I_brach_lim_sup_1('1') .. (sqr(V('13')*cos(alpha('13')) - V('14')*cos(alpha('14'))) + sqr(V('13')*sin(alpha('13')) - V('14')*sin(alpha('14'))))*sqr(Y('13','14')) =l= sqr(2.0);
I_brach_lim_sup_2('1') .. (sqr(V('15')*cos(alpha('15')) - V('16')*cos(alpha('16'))) + sqr(V('15')*sin(alpha('15')) - V('16')*sin(alpha('16'))))*sqr(Y('15','16')) =l= sqr(1.8);
Generator_currents(sg) .. sqr(Ig(sg)*V(sg)) - sqr(Pg(sg)) - sqr(Qg(sg)) =e= 0;
Power_factor(sg) .. sin(phi(sg)) - Qg(sg)/(V(sg)*Ig(sg)) =e= 0;
ed_p_initialization(sg) .. ed_p(sg,'1') - (Xq(sg) - Xq_p(sg))*Ig(sg)*cos(delta(sg,'1') - alpha(sg) + phi(sg)) =e= 0;
eq_p_initialization(sg) .. eq_p(sg,'1') + (Xd(sg) - Xd_p(sg))*Ig(sg)*sin(delta(sg,'1') - alpha(sg) + phi(sg)) - e_fd(sg) =e= 0;
Vd_initialization(sg) .. V(sg)*sin(delta(sg,'1') - alpha(sg)) - ed_p(sg,'1') + (Ra(sg)*sin(delta(sg,'1') - alpha(sg) +
                   phi(sg)) - Xq_p(sg)*cos(delta(sg,'1') - alpha(sg) + phi(sg)))*Ig(sg) =e= 0;
Vq_initialization(sg) .. V(sg)*cos(delta(sg,'1') - alpha(sg)) - eq_p(sg,'1') + (Ra(sg)*cos(delta(sg,'1') - alpha(sg) +
                   phi(sg)) + Xd_p(sg)*sin(delta(sg,'1') - alpha(sg) + phi(sg)))*Ig(sg) =e= 0;
Domega_initialization(sg) .. Domega(sg,'1') =e= 0;
id_initialization(sg) .. id(sg,'1') - Ig(sg)*sin(delta(sg,'1') - alpha(sg) + phi(sg)) =e= 0;
iq_initialization(sg) .. iq(sg,'1') - Ig(sg)*cos(delta(sg,'1') - alpha(sg) + phi(sg)) =e= 0;
e_fd_lim_inf(sg) .. e_fd_min(sg) =l= e_fd(sg);
e_fd_lim_sup(sg) .. e_fd(sg) =l= e_fd_max(sg);
Internal_voltaje_d(sg,s)$(not sfirst(s)) .. ed_p(sg,s)*(1 + Dt/(2*Tq_p(sg))) -
                   ed_p(sg,s-1)*(1 - Dt/(2*Tq_p(sg))) - (Dt/(2*Tq_p(sg)))*(Xq(sg) - Xq_p(sg))*(iq(sg,s) + iq(sg,s-1)) =e= 0;
Internal_voltaje_q(sg,s)$(not sfirst(s)) .. eq_p(sg,s)*(1 + Dt/(2*Td_p(sg))) -
                   eq_p(sg,s-1)*(1 - Dt/(2*Td_p(sg))) - (Dt/(2*Td_p(sg)))*(2*e_fd(sg) - (Xd(sg) - Xd_p(sg))*(id(sg,s) + id(sg,s-1))) =e= 0;
oscilation_omega(sg,s)$(not sfirst(s)) .. Domega(sg,s)*(1 + Dt*D(sg)/(4*H(sg))) -
                   Domega(sg,s-1)*(1 - Dt*D(sg)/(4*H(sg))) - (Dt/(4*H(sg)))*(2*Pg(sg) - Pe(sg,s) - Pe(sg,s-1)) =e= 0;
oscilation_delta(sg,s)$(not sfirst(s)) .. delta(sg,s) - delta(sg,s-1) -
                   (Dt*100*pi/2)*(Domega(sg,s) + Domega(sg,s-1)) =e= 0;
electric_power(sg,s) .. Pe(sg,s) - ed_p(sg,s)*id(sg,s) - eq_p(sg,s)*iq(sg,s) =e= 0;
center_of_inertia(s) .. delta_COI(s) - sum(sg,H(sg)*delta(sg,s)) / sum(sg,H(sg)) =e= 0;
angular_deviation_min(sg,s) .. - (angle_limit*pi/180) =l= delta(sg,s) - delta_COI(s);
angular_deviation_max(sg,s) .. delta(sg,s) - delta_COI(s) =l= (angle_limit*pi/180);
id_stator_fault(sg,sf) .. id(sg,sf)
                                               - (sum(sgp,Y_f_red(sg,sgp)*(ed_p(sgp,sf)*cos(delta(sg,sf) - delta(sgp,sf) - theta_f_red(sg,sgp))
                                                                         + eq_p(sgp,sf)*sin(delta(sg,sf) - delta(sgp,sf) - theta_f_red(sg,sgp))))
                                               +  sum(dfig,Y_f_red(sg,dfig)*(vD_w(dfig,sf)*sin(delta(sg,sf) - theta_f_red(sg,dfig))
                                                                           - vQ_w(dfig,sf)*cos(delta(sg,sf) - theta_f_red(sg,dfig))))) =e= 0;
iq_stator_fault(sg,sf) .. iq(sg,sf)
                                               - (sum(sgp,Y_f_red(sg,sgp)*(eq_p(sgp,sf)*cos(delta(sg,sf) - delta(sgp,sf) - theta_f_red(sg,sgp))
                                                                         - ed_p(sgp,sf)*sin(delta(sg,sf) - delta(sgp,sf) - theta_f_red(sg,sgp))))
                                               +  sum(dfig,Y_f_red(sg,dfig)*(vD_w(dfig,sf)*cos(delta(sg,sf) - theta_f_red(sg,dfig))
                                                                           + vQ_w(dfig,sf)*sin(delta(sg,sf) - theta_f_red(sg,dfig))))) =e= 0;
id_stator_postfault(sg,spf) .. id(sg,spf)
                                     - (sum(sgp,Y_pf_red(sg,sgp)*(ed_p(sgp,spf)*cos(delta(sg,spf) - delta(sgp,spf) - theta_pf_red(sg,sgp))
                                                                + eq_p(sgp,spf)*sin(delta(sg,spf) - delta(sgp,spf) - theta_pf_red(sg,sgp))))
                                     +  sum(dfig,Y_pf_red(sg,dfig)*(vD_w(dfig,spf)*sin(delta(sg,spf) - theta_pf_red(sg,dfig))
                                                                  - vQ_w(dfig,spf)*cos(delta(sg,spf) - theta_pf_red(sg,dfig))))) =e= 0;
iq_stator_postfault(sg,spf) .. iq(sg,spf)
                                     - (sum(sgp,Y_pf_red(sg,sgp)*(eq_p(sgp,spf)*cos(delta(sg,spf) - delta(sgp,spf) - theta_pf_red(sg,sgp))
                                                                - ed_p(sgp,spf)*sin(delta(sg,spf) - delta(sgp,spf) - theta_pf_red(sg,sgp))))
                                     +  sum(dfig,Y_pf_red(sg,dfig)*(vD_w(dfig,spf)*cos(delta(sg,spf) - theta_pf_red(sg,dfig))
                                                                  + vQ_w(dfig,spf)*sin(delta(sg,spf) - theta_pf_red(sg,dfig))))) =e= 0;
iD_w_fault(dfig,sf) .. iD_w(dfig,sf)*Ngen_w(dfig)
                                        - (sum(sg,Y_f_red(dfig,sg)*(ed_p(sg,sf)*sin(delta(sg,sf) + theta_f_red(dfig,sg))
                                                                  + eq_p(sg,sf)*cos(delta(sg,sf) + theta_f_red(dfig,sg))))
                                        +  sum(dfigp,Y_f_red(dfig,dfigp)*(vD_w(dfigp,sf)*cos(theta_f_red(dfig,dfigp))
                                                                        - vQ_w(dfigp,sf)*sin(theta_f_red(dfig,dfigp))))) =e= 0;
iQ_w_fault(dfig,sf) .. iQ_w(dfig,sf)*Ngen_w(dfig)
                                        - (sum(sg,Y_f_red(dfig,sg)*(- ed_p(sg,sf)*cos(delta(sg,sf) + theta_f_red(dfig,sg))
                                                                    + eq_p(sg,sf)*sin(delta(sg,sf) + theta_f_red(dfig,sg))))
                                        +  sum(dfigp,Y_f_red(dfig,dfigp)*(vD_w(dfigp,sf)*sin(theta_f_red(dfig,dfigp))
                                                                        + vQ_w(dfigp,sf)*cos(theta_f_red(dfig,dfigp))))) =e= 0;
iD_w_postfault(dfig,spf) .. iD_w(dfig,spf)*Ngen_w(dfig)
                                             - (sum(sg,Y_pf_red(dfig,sg)*(ed_p(sg,spf)*sin(delta(sg,spf) + theta_pf_red(dfig,sg))
                                                                        + eq_p(sg,spf)*cos(delta(sg,spf) + theta_pf_red(dfig,sg))))
                                             +  sum(dfigp,Y_pf_red(dfig,dfigp)*(vD_w(dfigp,spf)*cos(theta_pf_red(dfig,dfigp))
                                                                              - vQ_w(dfigp,spf)*sin(theta_pf_red(dfig,dfigp))))) =e= 0;
iQ_w_postfault(dfig,spf) .. iQ_w(dfig,spf)*Ngen_w(dfig)
                                             - (sum(sg,Y_pf_red(dfig,sg)*(- ed_p(sg,spf)*cos(delta(sg,spf) + theta_pf_red(dfig,sg))
                                                                          + eq_p(sg,spf)*sin(delta(sg,spf) + theta_pf_red(dfig,sg))))
                                             +  sum(dfigp,Y_pf_red(dfig,dfigp)*(vD_w(dfigp,spf)*sin(theta_pf_red(dfig,dfigp))
                                                                              + vQ_w(dfigp,spf)*cos(theta_pf_red(dfig,dfigp))))) =e= 0;
eD_p_w_equation(dfig,s)$(ord(s) gt 1) .. eD_p_w(dfig,s)*To_p_var(dfig,s) =e= eD_p_w(dfig,s-1)*To_p_var(dfig,s-1) + Dt/2*((VQr(dfig,s) + VQr(dfig,s-1)) -
                           (IQr(dfig,s) + IQr(dfig,s-1)) + w0*(To_p_var(dfig,s)*eQ_p_w(dfig,s)*delta_wr(dfig,s) + To_p_var(dfig,s-1)*eQ_p_w(dfig,s-1)*delta_wr(dfig,s-1)));
eQ_p_w_equation(dfig,s)$(ord(s) gt 1) .. eQ_p_w(dfig,s)*To_p_var(dfig,s) =e= eQ_p_w(dfig,s-1)*To_p_var(dfig,s-1) + Dt/2*((VDr(dfig,s) + VDr(dfig,s-1)) -
                           (IDr(dfig,s) + IDr(dfig,s-1)) - w0*(To_p_var(dfig,s)*eD_p_w(dfig,s)*delta_wr(dfig,s) + To_p_var(dfig,s-1)*eD_p_w(dfig,s-1)*delta_wr(dfig,s-1)));
rotor_flux_D(dfig,s) .. IDr(dfig,s) =e= eQ_p_w(dfig,s) + (xss(dfig) - xs_p(dfig))*iDs_w(dfig,s);
rotor_flux_Q(dfig,s) .. IQr(dfig,s) =e= eD_p_w(dfig,s) + (xss(dfig) - xs_p(dfig))*iQs_w(dfig,s);
stator_equation_D(dfig,s) .. eD_p_w(dfig,s) =e= xs_p(dfig)*iQs_w(dfig,s) - vD_w(dfig,s);
stator_equation_Q(dfig,s) .. eQ_p_w(dfig,s) =e= vQ_w(dfig,s) + xs_p(dfig)*iDs_w(dfig,s);
vD_w_equation(dfig,s) .. vD_w(dfig,s) =e= v_w(dfig,s)*cos(alpha_w(dfig,s));
vQ_w_equation(dfig,s) .. vQ_w(dfig,s) =e= v_w(dfig,s)*sin(alpha_w(dfig,s));
T_el_equation(dfig,s) .. T_el(dfig,s) =e= iQs_w(dfig,s)*IDr(dfig,s) - iDs_w(dfig,s)*IQr(dfig,s);
delta_wr_equation(dfig,s)$(ord(s) gt 1) .. delta_wr(dfig,s) =e= delta_wr(dfig,s-1) + Dt/(4*H_w(dfig))*((T_mech(dfig,s-1) - T_el(dfig,s-1)) + (T_mech(dfig,s) - T_el(dfig,s)));
wr_equation(dfig,s) .. wr(dfig,s) =e= 1 + delta_wr(dfig,s);
lambda_equation(dfig,s) .. lambda(dfig,s) =e= Kb(dfig)*(wr(dfig,s)/u_wind(dfig));
Cp_equation(dfig,s) .. Cp(dfig,s) =e= sum(i,sum(j,alpha_ij(i,j)*power(pitch_angle(dfig,s),ord(i)-1)*power(lambda(dfig,s),ord(j)-1)));
P_mech_equation(dfig,s) .. P_mech(dfig,s) =e= Cp(dfig,s)*Kwind(dfig)*power(u_wind(dfig),3)*Sb_w(dfig)/Sb;
T_mech_equation(dfig,s) .. T_mech(dfig,s)*wr(dfig,s) =e= P_mech(dfig,s);
error_wr_equation(dfig,s) .. error_wr(dfig,s) =e= wr(dfig,s) - wr_ref(dfig);
T_el_ref_ini(dfig,s)$(ord(s) eq 1) .. T_el_ref(dfig,s) =e= Ps_w(dfig,s);
T_el_ref_equation(dfig,s)$(ord(s) gt 1) .. T_el_ref(dfig,s) =e= T_el_ref(dfig,s-1)
                                                                                      + Kptrq(dfig)*(error_wr(dfig,s) - error_wr(dfig,s-1))
                                                                                      + Kitrq(dfig)*Dt/2*(error_wr(dfig,s) + error_wr(dfig,s-1));
P_ord_pre_equation(dfig,s) .. P_ord_pre(dfig,s) =e= T_el_ref(dfig,s);
P_ord_ini(dfig,s)$(ord(s) eq 1) .. P_ord(dfig,s) =e= P_ord_pre(dfig,s);
P_ord_equation(dfig,s)$(ord(s) gt 1) .. P_ord(dfig,s) =e= P_ord(dfig,s-1)*(2*T_pc(dfig) - Dt)/(2*T_pc(dfig) + Dt)
                                                                            + Dt/(2*T_pc(dfig) + Dt)*(P_ord_pre(dfig,s) + P_ord_pre(dfig,s-1));
pitch_angle_1(dfig,spitch_1)$(ord(spitch_1) gt 1) .. pitch_angle(dfig,spitch_1) =e= pitch_angle(dfig,spitch_1-1) + 10*Dt;
pitch_angle_2(dfig,spitch_2)$(ord(spitch_2) gt 1) .. pitch_angle(dfig,spitch_2) =e= pitch_angle(dfig,spitch_2-1);
pll_angle_ini(dfig,s) .. pll_angle(dfig,s) =e= alpha_w(dfig,s) - pi/2;
IDr_equation(dfig,sdfig_2) .. IDr(dfig,sdfig_2) =e= IDr_c(dfig,sdfig_2)*cos(pll_angle(dfig,sdfig_2)) - IQr_c(dfig,sdfig_2)*sin(pll_angle(dfig,sdfig_2));
IQr_equation(dfig,sdfig_2) .. IQr(dfig,sdfig_2) =e= IDr_c(dfig,sdfig_2)*sin(pll_angle(dfig,sdfig_2)) + IQr_c(dfig,sdfig_2)*cos(pll_angle(dfig,sdfig_2));
Ps_w_ref_equation(dfig,sdfig_2) .. Ps_w_ref(dfig,sdfig_2) =e= P_ord(dfig,sdfig_2);
iQr_c_ref_equation(dfig,sdfig_2) .. IQr_c_ref(dfig,sdfig_2) =e= Ps_w_ref(dfig,sdfig_2)*xss(dfig)/v_w(dfig,sdfig_2);
IDr_c_ref_2(dfig,sdfig_2) .. IDr_c_ref(dfig,sdfig_2) =e= v_w(dfig,sdfig_2);
IDr_c_2(dfig,sdfig_2) .. IDr_c(dfig,sdfig_2) =e= IDr_c_ref(dfig,sdfig_2);
IQr_c_2(dfig,sdfig_2) .. IQr_c(dfig,sdfig_2) =e= IQr_c_ref(dfig,sdfig_2);
P_gsc_equation(dfig,s) .. P_gsc(dfig,s) =e= rr_var(dfig,s)/sqr(lm(dfig))*(VDr(dfig,s)*IDr(dfig,s) + VQr(dfig,s)*IQr(dfig,s));
Ps_w_equation(dfig,s) .. Ps_w(dfig,s) =e= vD_w(dfig,s)*iDs_w(dfig,s) + vQ_w(dfig,s)*iQs_w(dfig,s);
Qs_w_equation(dfig,s) .. Qs_w(dfig,s) =e= vQ_w(dfig,s)*iDs_w(dfig,s) - vD_w(dfig,s)*iQs_w(dfig,s);
Pg_w_equation(dfig,s) .. Pg_w(dfig,s) =e= vD_w(dfig,s)*iD_w(dfig,s) + vQ_w(dfig,s)*iQ_w(dfig,s);
Qg_w_equation(dfig,s) .. Qg_w(dfig,s) =e= vQ_w(dfig,s)*iD_w(dfig,s) - vD_w(dfig,s)*iQ_w(dfig,s);
P_balance(dfig,s) .. Pg_w(dfig,s) - Ps_w(dfig,s) + P_gsc(dfig,s) =e= 0;
Q_balance(dfig,s) .. Qg_w(dfig,s) - Qs_w(dfig,s) =e= 0;


Model tscopf /all/;
        tscopf.workfactor = 100;
        Option nlp = ipopt
        iterlim = 10000;
        Solve tscopf using nlp minimizing z;
        Display Pg.l, Qg.l, V.l, alpha.l, ed_p.l, eq_p.l, delta.l, Domega.l, Pe.l;

file salida /galla_plot_60_1_0.m/;
salida.nd = 8;
put salida
put 'f_obj = [...'/put z.l/put'];'/
put 'Pg = [...'/loop(rb, put Pg.l(rb))put'];'/
put 'Qg = [...'/loop(rb, put Qg.l(rb))put'];'/
put 'V = [...'/loop(b, put V.l(b))put'];'/
put 'alpha = [...'/loop(b, put alpha.l(b))put'];'/
put 't = [...'/loop(s, put (ord(s)*Dt)/)put'];'/
put 'delta = [...'/loop(s,loop(sg, put delta.l(sg,s))put /);put'];'/
put 'delta_COI = [...'/loop(s, put delta_COI.l(s) /)put'];'/
put 'Domega = [...'/loop(s,loop(sg, put Domega.l(sg,s))put /);put'];'/
put 'ed_p = [...'/loop(s,loop(sg, put ed_p.l(sg,s))put /);put'];'/
put 'eq_p = [...'/loop(s,loop(sg, put eq_p.l(sg,s))put /);put'];'/
put 'Pe = [...'/loop(s,loop(sg, put Pe.l(sg,s))put /);put'];'/

put 'figure(1)'/
put 'plot(t, delta(:, 1))'/
put 'hold on'/
put 'plot(t, delta(:, 2))'/
put 'hold on'/
put 'plot(t, delta(:, 3))'/
put 'hold on'/
put 'plot(t, delta(:, 4))'/
put 'hold on'/
put 'plot(t, delta(:, 5))'/

put ' rr_var = [...'/loop(s,    loop(dfig, put rr_var.l(dfig,s))    put /);put'];'/
put ' To_p_var = [...'/loop(s,    loop(dfig, put To_p_var.l(dfig,s))    put /);put'];'/
put ' T_el = [...'/loop(s,    loop(dfig, put T_el.l(dfig,s))    put /);put'];'/
put ' T_mech = [...'/loop(s,    loop(dfig, put T_mech.l(dfig,s))    put /);put'];'/
put ' wr = [...'/loop(s,    loop(dfig, put wr.l(dfig,s))    put /);put'];'/
put ' delta_wr = [...'/loop(s,    loop(dfig, put delta_wr.l(dfig,s))    put /);put'];'/
put ' P_mech = [...'/loop(s,    loop(dfig, put P_mech.l(dfig,s))    put /);put'];'/
put ' lambda = [...'/loop(s,    loop(dfig, put lambda.l(dfig,s))    put /);put'];'/
put ' Cp = [...'/loop(s,    loop(dfig, put Cp.l(dfig,s))    put /);put'];'/
put ' vD_w = [...'/loop(s,    loop(dfig, put vD_w.l(dfig,s))    put /);put'];'/
put ' vQ_w = [...'/loop(s,    loop(dfig, put vQ_w.l(dfig,s))    put /);put'];'/
put ' v_w = [...'/loop(s,    loop(dfig, put v_w.l(dfig,s))    put /);put'];'/
put ' alpha_w = [...'/loop(s,    loop(dfig, put alpha_w.l(dfig,s))    put /);put'];'/
put ' vDr = [...'/loop(s,    loop(dfig, put vDr.l(dfig,s))    put /);put'];'/
put ' vQr = [...'/loop(s,    loop(dfig, put vQr.l(dfig,s))    put /);put'];'/
put ' iDr = [...'/loop(s,    loop(dfig, put iDr.l(dfig,s))    put /);put'];'/
put ' iQr = [...'/loop(s,    loop(dfig, put iQr.l(dfig,s))    put /);put'];'/
put ' iDr_c = [...'/loop(sdfig_2,    loop(dfig, put iDr_c.l(dfig,sdfig_2))    put /);put'];'/
put ' iQr_c = [...'/loop(sdfig_2,    loop(dfig, put iQr_c.l(dfig,sdfig_2))    put /);put'];'/
put ' eD_p_w = [...'/loop(s,    loop(dfig, put eD_p_w.l(dfig,s))    put /);put'];'/
put ' eQ_p_w = [...'/loop(s,    loop(dfig, put eQ_p_w.l(dfig,s))    put /);put'];'/
put ' iDs_w = [...'/loop(s,    loop(dfig, put iDs_w.l(dfig,s))    put /);put'];'/
put ' iQs_w = [...'/loop(s,    loop(dfig, put iQs_w.l(dfig,s))    put /);put'];'/
put ' iD_w = [...'/loop(s,    loop(dfig, put iD_w.l(dfig,s))    put /);put'];'/
put ' iQ_w = [...'/loop(s,    loop(dfig, put iQ_w.l(dfig,s))    put /);put'];'/
put ' wr_ref_var = [...'/loop(s,    loop(dfig, put wr_ref_var.l(dfig,s))    put /);put'];'/
put ' error_wr = [...'/loop(s,    loop(dfig, put error_wr.l(dfig,s))    put /);put'];'/
put ' T_el_ref = [...'/loop(s,    loop(dfig, put T_el_ref.l(dfig,s))    put /);put'];'/
put ' P_ord_pre = [...'/loop(s,    loop(dfig, put P_ord_pre.l(dfig,s))    put /);put'];'/
put ' P_ord = [...'/loop(s,    loop(dfig, put P_ord.l(dfig,s))    put /);put'];'/
put ' pitch_angle = [...'/loop(s,    loop(dfig, put pitch_angle.l(dfig,s))    put /);put'];'/
put ' pll_angle = [...'/loop(s,    loop(dfig, put pll_angle.l(dfig,s))    put /);put'];'/
*put ' u_pll = [...'/loop(s,    loop(dfig, put u_pll.l(dfig,s))    put /);put'];'/
*put ' vD_w_c = [...'/loop(s,    loop(dfig, put vD_w_c.l(dfig,s))    put /);put'];'/
put ' Ps_w_ref = [...'/loop(sdfig_2,    loop(dfig, put Ps_w_ref.l(dfig,sdfig_2))    put /);put'];'/
put ' iDr_c_ref = [...'/loop(sdfig_2,    loop(dfig, put iDr_c_ref.l(dfig,sdfig_2))    put /);put'];'/
put ' iQr_c_ref = [...'/loop(sdfig_2,    loop(dfig, put iQr_c_ref.l(dfig,sdfig_2))    put /);put'];'/
put ' P_gsc = [...'/loop(s,    loop(dfig, put P_gsc.l(dfig,s))    put /);put'];'/
put ' Ps_w = [...'/loop(s,    loop(dfig, put Ps_w.l(dfig,s))    put /);put'];'/
put ' Qs_w = [...'/loop(s,    loop(dfig, put Qs_w.l(dfig,s))    put /);put'];'/
put ' Pg_w = [...'/loop(s,    loop(dfig, put Pg_w.l(dfig,s))    put /);put'];'/
put ' Qg_w = [...'/loop(s,    loop(dfig, put Qg_w.l(dfig,s))    put /);put'];'/



put "figure(2), plot(t(1:end, 1), rr_var(:, 1)), title('rr_var')"/
put "figure(3), plot(t(1:end, 1), To_p_var(:, 1)), title('To_p_var')"/
put "figure(4), plot(t(1:end, 1), T_el(:, 1)), title('T_el')"/
put "figure(5), plot(t(1:end, 1), T_mech(:, 1)), title('T_mech')"/
put "figure(6), plot(t(1:end, 1), wr(:, 1)), title('wr')"/
put "figure(7), plot(t(1:end, 1), delta_wr(:, 1)), title('delta_wr')"/
put "figure(8), plot(t(1:end, 1), P_mech(:, 1)), title('P_mech')"/
put "figure(9), plot(t(1:end, 1), lambda(:, 1)), title('lambda')"/
put "figure(10), plot(t(1:end, 1), Cp(:, 1)), title('Cp')"/
put "figure(11), plot(t(1:end, 1), vD_w(:, 1)), title('vD_w')"/
put "figure(12), plot(t(1:end, 1), vQ_w(:, 1)), title('vQ_w')"/
put "figure(13), plot(t(1:end, 1), v_w(:, 1)), title('v_w')"/
put "figure(14), plot(t(1:end, 1), alpha_w(:, 1)), title('alpha_w')"/
put "figure(15), plot(t(1:end, 1), vDr(:, 1)), title('vDr')"/
put "figure(16), plot(t(1:end, 1), vQr(:, 1)), title('vQr')"/
put "figure(17), plot(t(1:end, 1), iDr(:, 1)), title('iDr')"/
put "figure(18), plot(t(1:end, 1), iQr(:, 1)), title('iQr')"/
put "figure(19), plot(t(20:end, 1), iDr_c(:, 1)), title('iDr_c')"/
put "figure(20), plot(t(20:end, 1), iQr_c(:, 1)), title('iQr_c')"/
put "figure(21), plot(t(1:end, 1), eD_p_w(:, 1)), title('eD_p_w')"/
put "figure(22), plot(t(1:end, 1), eQ_p_w(:, 1)), title('eQ_p_w')"/
put "figure(23), plot(t(1:end, 1), iDs_w(:, 1)), title('iDs_w')"/
put "figure(24), plot(t(1:end, 1), iQs_w(:, 1)), title('iQs_w')"/
put "figure(25), plot(t(1:end, 1), iD_w(:, 1)), title('iD_w')"/
put "figure(26), plot(t(1:end, 1), iQ_w(:, 1)), title('iQ_w')"/
put "figure(27), plot(t(1:end, 1), wr_ref_var(:, 1)), title('wr_ref_var')"/
put "figure(28), plot(t(1:end, 1), error_wr(:, 1)), title('error_wr')"/
put "figure(29), plot(t(1:end, 1), T_el_ref(:, 1)), title('T_el_ref')"/
put "figure(30), plot(t(1:end, 1), P_ord_pre(:, 1)), title('P_ord_pre')"/
put "figure(31), plot(t(1:end, 1), P_ord(:, 1)), title('P_ord')"/
put "figure(32), plot(t(1:end, 1), pitch_angle(:, 1)), title('pitch_angle')"/
put "figure(33), plot(t(1:end, 1), pll_angle(:, 1)), title('pll_angle')"/
*put "figure(34), plot(t(1:end, 1), u_pll(:, 1)), title('u_pll')"/
*put "figure(35), plot(t(1:end, 1), vD_w_c(:, 1)), title('vD_w_c')"/
put "figure(36), plot(t(20:end, 1), Ps_w_ref(:, 1)), title('Ps_w_ref')"/
put "figure(37), plot(t(20:end, 1), iDr_c_ref(:, 1)), title('iDr_c_ref')"/
put "figure(38), plot(t(20:end, 1), iQr_c_ref(:, 1)), title('iQr_c_ref')"/
put "figure(39), plot(t(1:end, 1), P_gsc(:, 1)), title('P_gsc')"/
put "figure(40), plot(t(1:end, 1), Ps_w(:, 1)), title('Ps_w')"/
put "figure(41), plot(t(1:end, 1), Qs_w(:, 1)), title('Qs_w')"/
put "figure(42), plot(t(1:end, 1), Pg_w(:, 1)), title('Pg_w')"/
put "figure(43), plot(t(1:end, 1), Qg_w(:, 1)), title('Qg_w')"/
