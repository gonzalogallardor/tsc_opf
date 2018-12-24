from functionalities.read_module import *
from functionalities.construction_module import *
from functionalities.write_module import *
import numpy as np
import copy

# Defining data paths
data_path = '../tests/gallardo/raw_data.raw'
gen_parameters_path = '../tests/gallardo/thermal_units_param.csv'
hvdc_param_path = '../tests/gallardo/hvdc_param.csv'
dfig_param_path = '../tests/gallardo/dfig_param.csv'
tsc_param_path = '../tests/gallardo/tsc_param.csv'
opf_param_path = '../tests/gallardo/opf_param.csv'

# Data read
s_base, buses, loads, generators, \
    branches, transformers_2w = read_raw_data(data_path)
add_gen_parameters(gen_parameters_path, generators)
dfig_param = read_dfig_param(dfig_param_path)
opf_param = read_opf_param(opf_param_path)
tsc_param = read_tsc_param(tsc_param_path)

# Extracting data
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
f = 50  # frequency

# load_factors = tuple(np.arange(1.05, 2, 0.05))
load_factors = (1.0,)
angle_limits = (60,)
for load_factor in load_factors:
    for angle_limit in angle_limits:
        
        # [Ybus] construction considering all buses, branches and transformers,
        Y = Ybus_construction(buses, branches, transformers_2w)

        # Adding "load_factor" to [Ybus:
        for load in loads:
            y_shunt_equivalent =  np.conjugate(complex(load.p, 
                                                       load.q)) * \
                                                       load_factor / s_base
            Y = Ybus_add_shunt_element(Y, 
                                       y_shunt_equivalent,
                                       load.bus_index)

        # Creation of extended Ybus matrix (prefault) : [Y_ext]

        total_buses = len(buses) + len(generators)
        Y_ext = np.zeros((total_buses, total_buses), complex)
        Y_ext[:len(buses), :len(buses)] = Y

        # Adding generator transformer impedances to Y_ext
        for generator in generators:
            
            y_transformer = 1 / complex(generator.r_transformer, 
                                        generator.x_transformer) * \
                                        generator.mva_base / s_base

            Y_ext = Ybus_add_series_element(Y_ext, 
                                            y_transformer, 
                                            generator.bus_index, 
                                            generator.internal_bus_index)

        ### Creation of faulted Ybus matrix (during the fault): [Y_fault_ext] and [Y_fault_reduced]
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

        # Creation of post fault Ybus matrix: [Y_postfault_ext] and [Y_postfault_reduced]
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

        # Writing the gams file
        gams_filename = "tscopf_" + str(angle_limit) + "_" + str(load_factor) + ".gms"
        file = open(gams_filename,'w')

        scalars = write_scalars(s_base,       # MVA base
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

        file.write(scalars)
        file.write(sets)
        file.write(parameters)
        file.write(tables)
        file.write(variables)
        file.write(var_restrictions)
        file.write(equation_names)
        file.write(equations)

        str_loadfactor = str(load_factor)
        str_loadfactor_split = str_loadfactor.split('.')

        file.write("\n\nModel tscopf /all/;\n\
        tscopf.workfactor = 100;\n\
        Option nlp = ipopt\n\
        iterlim = 10000;\n\
        Solve tscopf using nlp minimizing z;\n\
        Display Pg.l, Qg.l, V.l, alpha.l, ed_p.l, eq_p.l, delta.l, Domega.l, Pe.l;")

        file.write("\n\nfile salida /galla_plot_" + str(angle_limit) + "_" + str_loadfactor_split[0] + "_" + str_loadfactor_split[1] + ".m/;")
        file.write("""
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
put "figure(43), plot(t(1:end, 1), Qg_w(:, 1)), title('Qg_w')"/""")


        file.close()