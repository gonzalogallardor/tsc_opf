from functionalities.read_module import *
from functionalities.construction_module import *
from functionalities.write_module import *
import numpy as np
import copy

data_path = '../tests/gallardo/raw_data.raw'
gen_parameters_path = '../tests/gallardo/thermal_units_param.csv'
hvdc_param_path = '../tests/gallardo/hvdc_param.csv'
dfig_param_path = '../tests/gallardo/dfig_param.csv'
tsc_param_path = '../tests/gallardo/tsc_param.csv'
opf_param_path = '../tests/gallardo/opf_param.csv'

### DATA READ
s_base, buses, loads, generators, \
    branches, transformers = read_raw_data(data_path)
add_gen_parameters(gen_parameters_path, generators)
#hvdc_param = read_hvdc_param(hvdc_param_path)
dfig_param = read_dfig_param(dfig_param_path)
opf_param = read_opf_param(opf_param_path)
tsc_param = read_tsc_param(tsc_param_path)

### EXTRACTING DATA 
total_time = float(tsc_param['total_time'])
clearing_time = float(tsc_param['clearing_time'])
delta_t = float(tsc_param['delta_t'])
faulted_bus_number = int(tsc_param['faulted_bus'])
[faulted_bus_index] = get_index(buses, [faulted_bus_number])
isolated_branch = get_branch_index(branches,
                                   tsc_param['isolated_branch_bus1'],
                                   tsc_param['isolated_branch_bus2'])
dfig_bus_number = int(dfig_param['bus_number'])  # bus_number
[dfig_bus_index] = get_index(buses, [dfig_bus_number])
f = 50

# postfault_time_analisis / total_postfault_time

load_factors = tuple(np.arange(1.55, 2.5, 0.05))
angle_limits = (60,)
for load_factor in load_factors:
    for angle_limit in angle_limits:
        
        ### Initial Ybus construction considering all buses, branches and transformers,
        ### indexed according to read_data function and saved in the respective system
        ### class of each element

        Y = Ybus_construction(buses, branches, transformers)

        # Adding "load_factor" effect to Y matrix
        for load in loads:
            y_shunt_equivalent =  np.conjugate(complex(load.p, 
                                                       load.q)) * \
                                                       load_factor / s_base
            Y = Ybus_add_shunt_element(Y, 
                                       y_shunt_equivalent,
                                       load.bus_index)

#        Y[hvdc_bus_index, hvdc_bus_index] += complex(0.0,1.23)    #  HVDC link capacitor: Q_consumption = 1/2 P_hvdc


        ###
        ### Creation of extended Ybus matrix (prefault) : [Y_ext]
        ###

        total_buses = len(buses) + len(generators)
        Y_ext = np.zeros((total_buses, total_buses), complex)
        Y_ext[:len(buses), :len(buses)] = Y

        # Adding generator transformer impedances to Y_ext
        for generator in generators:
            
            y_tran = 1 / complex(generator.r_tran, 
                                 generator.x_tran) * \
                        generator.mva_base / s_base

            Y_ext = Ybus_add_series_element(Y_ext, 
                                            y_tran, 
                                            generator.bus_index, 
                                            generator.internal_bus_index)


        ###
        ### Creation of faulted Ybus matrix (during the fault): [Y_fault_ext] and [Y_fault_reduced]
        ###

        Y_fault_ext = copy.deepcopy(Y_ext)
        Y_fault_ext[faulted_bus_index, faulted_bus_index] += complex(0,-1e6)  # 3-f fault representation

        retained_bus_indexes = []
        retained_bus_numbers = []
        for generator in generators:
            retained_bus_indexes.append(generator.internal_bus_index)
            retained_bus_numbers.append(generator.bus_number)

        
        retained_bus_indexes.append(dfig_bus_index)
        retained_bus_numbers.append(dfig_bus_number)
        
        # Creation of a reduced faulted Ybus matrix considering generators internal buses and HVDC-link bus (only)
        Y_fault_reduced = kron_reduction(Y_fault_ext, 
                                         retained_bus_indexes)

        ###
        ### Creation of post fault Ybus matrix: [Y_postfault_ext] and [Y_postfault_reduced]
        ###

        Y_postfault_ext = copy.deepcopy(Y_ext)

        # removing isolated branch from Ybus
        y_series = 1 / complex(isolated_branch.r, isolated_branch.x)
        y_shunt = 1 / 2 * complex(isolated_branch.g, isolated_branch.b)

        Y_postfault_ext = Ybus_add_series_element(Y_postfault_ext,
                                                  - y_series,
                                                  isolated_branch.bus_index_1,
                                                  isolated_branch.bus_index_2)

        Y_postfault_ext = Ybus_add_shunt_element(Y_postfault_ext,
                                                 - y_shunt,
                                                 isolated_branch.bus_index_1)

        Y_postfault_ext = Ybus_add_shunt_element(Y_postfault_ext,
                                                 - y_shunt,
                                                 isolated_branch.bus_index_2)

        Y_postfault_reduced = kron_reduction(Y_postfault_ext, 
                                             retained_bus_indexes)


        ###
        ### Writing the gams file
        ###

        scalars = write_scalars(s_base,       # MVA
                                delta_t,      # in seconds
                                angle_limit,  # degrees
                                f)            # frequency[Hz]

        sets = write_sets(buses,                 # all buses instances in the list
                          generators,            # all generator instances in the list
                          branches,              # all branch instances in the list
                          dfig_param,            # dfig parameters
                          retained_bus_numbers,  # retained bus numbers as list
                          total_time,            # seconds
                          clearing_time,         # seconds
                          delta_t)               # seconds
        
        sets += "\nsfirst(s) = yes$(ord(s) eq 1);\n" + \
                "alias (sg,sgp);\n" + \
                "alias (dfig,dfigp);\n" + \
                "alias (i,j);\n"

        parameters = write_parameters(buses,
                                      generators,
                                      branches,
                                      tsc_param,
                                      dfig_param,  # dfig parameters as structured array
                                      s_base)


        tables = write_tables(Y,                     # Y bus, all buses, pre fault stage  (array)
                              Y_fault_reduced,       # Y bus, retained buses only, during fault (array)
                              Y_postfault_reduced,   # Y bus, retained buses only, post fault stage (array)
                              buses,                 # All buses instances(class) as list
                              retained_bus_numbers)  # All retained in a list

        tables += 'Table\n\
alpha_ij(i,j)\n\
               1              2               3              4              5\n\
1   -4.1909e-001    2.1808e-001    -1.2406e-002   -1.3365e-004    1.1524e-005\n\
2   -6.7606e-002    6.0405e-002    -1.3934e-002    1.0683e-003   -2.3895e-005\n\
3    1.5727e-002   -1.0996e-002     2.1495e-003   -1.4855e-004    2.7937e-006\n\
4   -8.6018e-004    5.7051e-004    -1.0479e-004    5.9924e-006   -8.9194e-008\n\
5    1.4787e-005   -9.4839e-006     1.6167e-006   -7.1535e-008    4.9686e-010;'

        variables, var_restrictions, equation_names, equations = write_equations(generators,
                                                                                branches,
                                                                                dfig_param,
                                                                                opf_param,
                                                                                s_base)
        gams_file = "galla_opf_" + str(load_factor) + ".gms"
        file = open(gams_file,'w')
        file.write("""
Scalars
    Sb "Potencia base en [MVA]" /100/
    w0 /314.16/;


Sets
    b buses /1 * 20/
    rb(b) retained buses /1 * 6/
    sg(rb) synchronous generators /1 * 5/
    dfig(rb) dfig wind-plant bus /6/
    ngb(b) nongenbuses /7 * 20/
    i /1 * 5/;

* añade a sfirst el primer index de s
alias (sg,sgp);
alias (dfig,dfigp)
alias (i,j);

Parameters
    a1(b) /1 7.0e1, 2 8.0e1, 3 4.00e1, 4 1.00e2, 5 1.20e2, 6 0/
    Pl(b) /1 0, 2 0, 3 0, 4 0, 5 0, 6 0, 7 0, 8 0, 9 0, 10 0, 11 0, 12 0, 13 0, 14 0, 15 0, 16 0, 17 0, 18 0, 19 0, 20 0/
    Ql(b) /1 0, 2 0, 3 0, 4 0, 5 0, 6 0, 7 0, 8 0, 9 0, 10 0, 11 0, 12 0, 13 0, 14 0, 15 0, 16 0, 17 0, 18 0, 19 0, 20 0/
    Pmax(sg) /1 600, 2 470, 3 510, 4 275, 5 350/
    Pmax(sg);
    Pmax(sg) = pmax(sg)/Sb;
Parameters
    Pmin(sg) /1 0.000000, 2 0.000000, 3 0.000000, 4 0.000000, 5 0.000000/
    Qmax(sg);
    Qmax(sg) = Pmax(sg)*0.5;
Parameters
    Qmin(sg);
    Qmin(sg) = -Qmax(sg);
Parameters
    Sm(sg);
    Sm(sg) = Pmax(sg)*1.1;
Parameters
    e_fd_max(sg) /1 2.0, 2 2.0, 3 2.0, 4 2.0, 5 2.0/
    e_fd_min(sg) /1 0.0, 2 0.0, 3 0.0, 4 0.0, 5 0.0/
    Ra(sg) /1 0.0, 2 0.0, 3 0.0, 4 0.0, 5 0.0/
    Ra(sg);
    Ra(sg) = ra(sg)/Sm(sg);
Parameters
    Xd(sg) /1 1.5, 2 1.5, 3 1.5, 4 1.5, 5 1.5/
    Xd(sg);
    Xd(sg) = xd(sg)/Sm(sg);
Parameters
    Xd_p(sg) /1 0.3, 2 0.3, 3 0.3, 4 0.3, 5 0.3/
    Xd_p(sg);
    Xd_p(sg) = xd_p(sg)/Sm(sg);
Parameters
    Xq(sg) /1 1.5, 2 1.5, 3 1.5, 4 1.5, 5 1.5/
    Xq(sg);
    Xq(sg) = xq(sg)/Sm(sg);
Parameters
    Xq_p(sg) /1 0.3, 2 0.3, 3 0.3, 4 0.3, 5 0.3/
    Xq_p(sg);
    Xq_p(sg) = xq_p(sg)/Sm(sg);
Parameters
    H(sg) /1 3.2, 2 3.0, 3 3.0, 4 2.0, 5 2.0/
    H(sg);
    H(sg) = h(sg)*Sm(sg);
Parameters
    D(sg) /1 2.0, 2 2.0, 3 2.0, 4 2.0, 5 2.0/
    D(sg);
    D(sg) = d(sg)*Sm(sg);
Parameters
    Td_p(sg) /1 6.0, 2 6.0, 3 6.0, 4 6.0, 5 6.0/
    Tq_p(sg) /1 1.0, 2 1.0, 3 1.0, 4 1.0, 5 1.0/;

Parameters
Ngen_w(dfig) /6 30/
*MACHINE, in self nominal base!
Sb_w(dfig) /6 1.5/
w0_w(dfig) /6 377/
To_p(dfig) /6 1.5/
xss(dfig) /6 4/
xs_p(dfig) /6 0.35/
rs(dfig) /6 0.01/
xl_s(dfig) /6 0.15/
H_w(dfig) /6 3/
Kc_crowbar /6 10/

*PLL
Kpll(dfig) /6 30/
angle_pll_max(dfig) /6 0.1/

*CURRENT CONTROLLER
Kp(dfig) /6 10/
Ki(dfig) /6 1/

*ELECTRICAL CONTROLLER
KQi(dfig) /6 0.1/
KVi(dfig) /6 0.5/

*WIND POWER
pitch_angle(dfig) /6 0/
Pg_w_ini(dfig) /6 1/
Ps_w_max(dfig) /6 1/
Kb(dfig) /6 56.6/
Kwind(dfig) /6 0.00159/
u_wind(dfig) /6 11/

*TORQUE CONTROLLER
T_pc(dfig) /6 0.05/
*Kptrq(dfig) /6 3/
*Kitrq(dfig) /6 0.6/
Kptrq(dfig) /6 0.045/
Kitrq(dfig) /6 0.009/

*PITCH CONTROLLER
wr_ref(dfig) /6 1.2/
Kpp(dfig) /6 150/
Kip(dfig) /6 25/
Kpc(dfig) /6 3/
Kic(dfig) /6 30/
T_p(dfig) /6 0.3/

lm(dfig)
lrr(dfig)
rr(dfig)

* Initial real power inyected to the grid in System bases (p.u), SINGLE MACHINE
Pg_w_ini_SM(dfig)
Ps_w_ini_SM(dfig)
Ps_w_max_SM(dfig)

r_crowbar(dfig)

xm(dfig)
xrr(dfig)
xss(dfig)
lss(dfig)

*1: crowbar connected, 2:crowbar disconnected (control activated)
rr_1(dfig)
rr_2(dfig)
To_p_1(dfig)
To_p_2(dfig)

delta_wr(dfig);

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

*initial (estimated) inyected active power of a single dfig wind turbine, with v_wind given, in Sb[MVA] base of the system
Pg_w_ini_SM(dfig) = Pg_w_ini(dfig)*Sb_w(dfig)/Sb;
Ps_w_ini_SM(dfig) = Pg_w_ini_SM(dfig)*0.8;
Ps_w_max_SM(dfig) = Ps_w_max(dfig)*Sb_w(dfig)/Sb;

rr_1(dfig) = rr(dfig) + r_crowbar(dfig);
rr_2(dfig) = rr(dfig);
To_p_1(dfig) = To_p(dfig)*rr(dfig)/(rr(dfig) + r_crowbar(dfig));
To_p_2(dfig) = To_p(dfig);

delta_wr(dfig) = wr_ref(dfig) - 1;
        """)

        file.write(tables)

        file.write("""
Variables
        Pg(rb)
        Qg(rb)
        Ig(sg)
        phi(sg)
        V(b)
        alpha(b)
        z;

        Pg.l('1') = 3; Pg.l('2') = 1; Pg.l('3') = 3; Pg.l('4') = 0.5; Pg.l('5') = 0.5;
        Pg.lo(dfig) = 0.48*1.25*Ps_w_max_SM(dfig)*Ngen_w(dfig); Pg.up(dfig) = 1.25*Ps_w_max_SM(dfig)*Ngen_w(dfig);
        Pg.lo(sg)=Pmin(sg); Pg.up(sg)=Pmax(sg);
        Qg.lo(sg)=Qmin(sg); Qg.up(sg)=Qmax(sg);
        Qg.fx(dfig) = 0;
        Ig.l(sg) = 1; Ig.lo(sg) = 0.001; Ig.up(sg) = Sm(sg);
        phi.l(sg) = 0; phi.lo(sg) = -pi/2; phi.up(sg) = pi/2;
        V.l(b) = 1; V.lo(b) = 0.95; V.up(b) = 1.05;
        alpha.l(b) = 0; alpha.lo(b) = -pi; alpha.up(b) = pi;

Variables
*DFIG SYSTEM VARIABLES
Cp(dfig)
lambda(dfig)
P_mech(dfig)
VDr(dfig)
VQr(dfig)
IDr(dfig)
IQr(dfig)
iDs_w(dfig)
iQs_w(dfig)
eD_p_w(dfig)
eQ_p_w(dfig)
wr(dfig)
Ps_w(dfig)
Pgsc_w(dfig)
Qs_w(dfig);
*********************************************************************DFIG INITIAL VALUES
lambda.l(dfig) = Kb(dfig)*wr_ref(dfig)/u_wind(dfig); lambda.lo(dfig) = 3; lambda.up(dfig) = 15;
Cp.l(dfig) = 0.5; Cp.lo(dfig) = 0;
P_mech.l(dfig) = 1.2*Ps_w_ini_SM(dfig); P_mech.lo(dfig) = 0;

Equations
total_cost

* Power Flow equations
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

Ps_eq
Pgsc_eq
Qs_eq

* Current limits for the braches
        I_brach_lim_inf_1
        I_brach_lim_inf_2
        I_brach_lim_sup_1
        I_brach_lim_sup_2

* Auxiliary equations
        Generators_current
        Power_factor

*Mechanical power & torque
wr_eq
lambda_equation
Cp_equation
P_mech_equation;

* Objetive function
total_cost .. z =e= sum(sg,a1(sg)*(Pg(sg)*Sb)) + sum(dfig,a1(dfig)*(Pg(dfig)*Sb));

* Power Flow equations
p_balance_gen(rb) .. Pg(rb) - Pl(rb) - V(rb)*sum(b,V(b)*Y(rb,b)*cos(alpha(rb) - alpha(b) - theta(rb,b))) =e= 0;
p_balance_nongen(ngb) ..  - Pl(ngb) - V(ngb)*sum(b,V(b)*Y(ngb,b)*cos(alpha(ngb) - alpha(b) - theta(ngb,b))) =e= 0;
q_balance_gen(rb) .. Qg(rb) - Ql(rb) - V(rb)*sum(b,V(b)*Y(rb,b)*sin(alpha(rb) - alpha(b) - theta(rb,b))) =e= 0;
q_balance_nongen(ngb) ..  - Ql(ngb) - V(ngb)*sum(b,V(b)*Y(ngb,b)*sin(alpha(ngb) - alpha(b) - theta(ngb,b))) =e= 0;
ref_bus .. alpha('1') =e= 0;

dfig_load_flow_f1(dfig) .. Pg(dfig)/Ngen_w(dfig) =e= Ps_w(dfig) - Pgsc_w(dfig);
dfig_load_flow_f2(dfig) .. Qg(dfig)/Ngen_w(dfig) =e= Qs_w(dfig);
dfig_load_flow_f3(dfig) .. VDr(dfig) =e= IDr(dfig) + w0*To_p(dfig)*delta_wr(dfig)*eD_p_w(dfig);
dfig_load_flow_f4(dfig) .. VQr(dfig) =e= IQr(dfig) - w0*To_p(dfig)*delta_wr(dfig)*eQ_p_w(dfig);
dfig_load_flow_f5(dfig) .. V(dfig)*cos(alpha(dfig)) =e= xs_p(dfig)*iQs_w(dfig) - eD_p_w(dfig);
dfig_load_flow_f6(dfig) .. V(dfig)*sin(alpha(dfig)) =e= - xs_p(dfig)*iDs_w(dfig) + eQ_p_w(dfig);
P_mech_ini(dfig) .. Pg(dfig)/Ngen_w(dfig) =e= P_mech(dfig) - rr(dfig)/sqr(lm(dfig))*(sqr(IDr(dfig)) + sqr(IQr(dfig)));

Ps_eq(dfig) .. Ps_w(dfig) =e= V(dfig)*cos(alpha(dfig))*iDs_w(dfig) + V(dfig)*sin(alpha(dfig))*iQs_w(dfig);
Qs_eq(dfig) .. Qs_w(dfig) =e= V(dfig)*sin(alpha(dfig))*iDs_w(dfig) - V(dfig)*cos(alpha(dfig))*iQs_w(dfig);
Pgsc_eq(dfig) .. Pgsc_w(dfig) =e= rr(dfig)/sqr(lm(dfig))*(VDr(dfig)*IDr(dfig) + VQr(dfig)*IQr(dfig));

* Current limits for the braches
        I_brach_lim_inf_1('1') .. 0 =l= (sqr(V('13')*cos(alpha('13')) - V('14')*cos(alpha('14'))) + sqr(V('13')*sin(alpha('13')) - V('14')*sin(alpha('14'))))*sqr(Y('13','14'));
        I_brach_lim_inf_2('1') .. 0 =l= (sqr(V('15')*cos(alpha('15')) - V('16')*cos(alpha('16'))) + sqr(V('15')*sin(alpha('15')) - V('16')*sin(alpha('16'))))*sqr(Y('15','16'));
        I_brach_lim_sup_1('1') .. (sqr(V('13')*cos(alpha('13')) - V('14')*cos(alpha('14'))) + sqr(V('13')*sin(alpha('13')) - V('14')*sin(alpha('14'))))*sqr(Y('13','14')) =l= sqr(2.0);
        I_brach_lim_sup_2('1') .. (sqr(V('15')*cos(alpha('15')) - V('16')*cos(alpha('16'))) + sqr(V('15')*sin(alpha('15')) - V('16')*sin(alpha('16'))))*sqr(Y('15','16')) =l= sqr(1.8);

* Auxiliary equations
        Generators_current(sg) .. sqr(Ig(sg)*V(sg)) - sqr(Pg(sg)) - sqr(Qg(sg)) =e= 0;
        Power_factor(sg) .. sin(phi(sg)) - Qg(sg)/(V(sg)*Ig(sg)) =e=0;

* Initial condition equations
wr_eq(dfig) .. wr(dfig) =e= 1 + delta_wr(dfig);
lambda_equation(dfig) .. lambda(dfig) =e= Kb(dfig)*(1.2/u_wind(dfig));
Cp_equation(dfig) .. Cp(dfig) =e= sum(i,sum(j,alpha_ij(i,j)*power(0,ord(i)-1)*power(lambda(dfig),ord(j)-1)));
P_mech_equation(dfig) .. P_mech(dfig) =e= Cp(dfig)*Kwind(dfig)*power(u_wind(dfig),3)*Sb_w(dfig)/Sb;
        """)

        str_loadfactor = str(load_factor)
        str_loadfactor_split = str_loadfactor.split('.')

        file.write("\n\nModel tscopf /all/;\n\
        Option nlp = ipopt\n\
        iterlim = 100000;\n\
        Solve tscopf using nlp minimizing z;\n\
        Display Pg.l, Qg.l, V.l, alpha.l;")

        file.write("\n\nfile salida /galla_opf_" + str_loadfactor_split[0] + "_" + str_loadfactor_split[1] + ".m/;")
        file.write("""
salida.nd = 8;
put salida

put 'f_obj = [...'/put z.l/put'];'/
put 'Pg = [...'/loop(rb, put Pg.l(rb))put'];'/
put 'Qg = [...'/loop(rb, put Qg.l(rb))put'];'/
put 'V = [...'/loop(b, put V.l(b))put'];'/
put 'alpha = [...'/loop(b, put alpha.l(b))put'];'
        """)


        file.close()