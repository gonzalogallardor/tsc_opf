import numpy as np

def write_scalars(s_base, delta_t, angle_limit, f):

    scalars = "Scalars\n\n"

    scalars += "Sb\t"+ "/" + str(s_base) + "/\n"
    scalars += "Dt\t"+ "/" + str(delta_t) + "/\n"
    scalars += "angle_limit\t"+ "/" + str(angle_limit) + "/\n"
    scalars += "w0\t"+ "/" + str(2*np.pi*f) + "/\n"

    return scalars


def write_sets(buses, generators, 
               branches, dfig_param, 
               retained_buses, total_time, 
               clearing_time, delta_t):
    
    sets = "\n\nSets\n\n"

    total_samples = int(total_time/delta_t + 1)
    sets += "s\t"+ "/" + "1 * " + str(total_samples) + "/\n"
    sets += "sfirst(s)\n"

    fault_sample = int(clearing_time/delta_t + 1)  # pesimist -> include fault_sample inside the fault period
    sets += "sf(s)\t"+ "/" + str(2) + " * " + str(fault_sample) + "/\n"
    
    fault_clearing_sample = int(fault_sample + 1)  # first sample inside postfault period = next sample after fault_sample
    sets += "spf(s)\t"+ "/" + str(fault_clearing_sample) + \
            " * " + str(total_samples) + "/\n"

    sets += "sdfig_1(s)\t"+ "/" + str(2) + \
            " * " + str(dfig_param['s_crow_end']) + "/\n"
    
    sets += "sdfig_2(s)\t"+ "/" + str(dfig_param['s_crow_end'] + 1) + \
            " * " + str(total_samples) + "/\n"

    sets += "spitch_1(s)\t"+ "/" + str(1) + \
            " * " + str(int(dfig_param['s_pitch_end'])) + "/\n"

    sets += "spitch_2(s)\t"+ "/" + str(dfig_param['s_pitch_end']) + \
            " * " + str(total_samples) + "/\n"

    bus_numbers = []
    for bus in buses:
        bus_numbers.append(bus.number)
    sets += "b\t" + "/"
    sets = setsFunc_bus_writing(sets, bus_numbers)
    sets += "/\n"
    
    retained_buses.sort()
    sets += "rb(b)\t" + "/"
    sets = setsFunc_bus_writing(sets, retained_buses)
    sets += "/\n"

    generator_buses = []
    for generator in generators:
        generator_buses.append(generator.bus_number)
    generator_buses.sort()
    sets += "sg(rb)\t" + "/"
    sets = setsFunc_bus_writing(sets, generator_buses)
    sets += "/\n"

    sets += "dfig(rb)\t" + "/" 
    sets = setsFunc_bus_writing(sets, [dfig_param['bus_number']])
    sets += "/\n"

    non_gen_buses = list(set(bus_numbers) - set(retained_buses))
    non_gen_buses.sort()
    sets += "ngb(b)\t" + "/"
    sets = setsFunc_bus_writing(sets, non_gen_buses)
    sets += "/\n"
    
    sets += "i\t" + "/"
    sets += "1 * 5"
    sets += "/;\n"

    return sets

def write_parameters(buses, 
                     generators, 
                     branches, 
                     tsc_param, 
                     dfig_param, 
                     s_base):
    
    parameters = "\n\nParameters\n\n"

    parameters += "a1(rb)\t" + "/"
    counter = 0
    for generator in generators:
            parameters += str(generator.bus_number)+ " " + str(generator.a_1) + ", "
            counter += 1

    parameters += str(dfig_param['bus_number']) + " " + str(dfig_param['a_1'])
    parameters += "/\n"

    parameters += "Pl(b)\t" + "/set.b 0/\n"

    parameters += "Ql(b)\t" + "/set.b 0/\n"

    sg_set_names = ('pmax',
             'pmin',
             'qmax',
             'qmin',
             'e_fd_max',
             'e_fd_min',
             'ra',
             'xd',
             'xd_p',
             'xq',
             'xq_p',
             'H',
             'D',
             'Td_p',
             'Tq_p')

    dfig_set_names = ('Ngen_w',
             'u_wind',
             'Sb_w',
             'w0_w',
             'To_p',
             'xss',
             'xs_p',
             'rs',
             'xl_s',
             'H_w',
             'kc_crowbar',
#             'Kpll',
#             'Kp',
#             'Ki',
#             'KQi',
#             'KVi',
             'initial_pitch_angle',
             'Kb',
             'Kwind',
             'T_pc',
             'Kptrq',
             'Kitrq',
             'wr_ref',
#             'Kpp',
#             'Kip',
#             'Kpc',
#             'Kic',
             'T_p')

    for set_name in sg_set_names:
        parameters += set_name + '(sg)\t' + '/'
        parameters = parametersFun_gen_writing(parameters,
                                               generators,
                                               set_name)
        parameters += '/\n'

    for set_name in dfig_set_names:
        parameters += set_name + '(dfig)\t' + '/'
        parameters = parametersFun_dfig_writing(parameters,
                                                dfig_param,
                                                set_name)
        parameters += '/\n'

    parameters += "\n"

    parameters += "Sm(sg)\n"
    parameters += "lm(dfig)\n"
    parameters += "lrr(dfig)\n"
    parameters += "rr(dfig)\n"
    parameters += "r_crowbar(dfig)\n"
    parameters += "xm(dfig)\n"
    parameters += "xrr(dfig)\n"
    parameters += "xss(dfig)\n"
    parameters += "lss(dfig)\n"
    parameters += "rr_1(dfig)\n"
    parameters += "rr_2(dfig)\n"
    parameters += "To_p_1(dfig)\n"
    parameters += "To_p_2(dfig);\n"

    parameters += "Pmax(sg) = Pmax(sg)/Sb;\n"
    parameters += "Pmin(sg) = Pmin(sg)/Sb;\n"
    parameters += "Qmax(sg) = Qmax(sg)/Sb;\n"
    parameters += "Qmin(sg) = Qmin(sg)/Sb;\n"
    parameters += "Sm(sg) = " + str(tsc_param['Sm/Pmax']) + "*Pmax(sg);\n"
    parameters += "Ra(sg) = Ra(sg)/Sm(sg);\n"
    parameters += "Xd(sg) = Xd(sg)/Sm(sg);\n"
    parameters += "Xd_p(sg) = Xd_p(sg)/Sm(sg);\n"
    parameters += "Xq(sg) = Xq(sg)/Sm(sg);\n"
    parameters += "Xq_p(sg) = Xq_p(sg)/Sm(sg);\n"
    parameters += "H(sg) = H(sg)*Sm(sg);\n"
    parameters += "D(sg) = D(sg)*Sm(sg);\n"

    parameters += "xss(dfig) = xss(dfig)*Sb/Sb_w(dfig)*w0/w0_w(dfig);\n"
    parameters += "xs_p(dfig) = xs_p(dfig)*Sb/Sb_w(dfig)*w0/w0_w(dfig);\n"
    parameters += "rs(dfig) = rs(dfig)*Sb/Sb_w(dfig);\n"
    parameters += "xl_s(dfig) = xl_s(dfig)*Sb/Sb_w(dfig)*w0/w0_w(dfig);\n"
    parameters += "lss(dfig) = xss(dfig);\n"
    parameters += "H_w(dfig) = H_w(dfig)*Sb_w(dfig)/Sb;\n"
    parameters += "xm(dfig) = xss(dfig)-xl_s(dfig);\n"
    parameters += "lm(dfig) = xm(dfig);\n"
    parameters += "lrr(dfig) = sqr(lm(dfig))/(xss(dfig)-xs_p(dfig));\n"
    parameters += "xrr(dfig) = lrr(dfig);\n"
    parameters += "rr(dfig) =  lrr(dfig)/(w0*To_p(dfig));\n"
    parameters += "r_crowbar(dfig) = kc_crowbar(dfig)*rr(dfig);\n"
    parameters += "rr_1(dfig) = rr(dfig) + r_crowbar(dfig);\n"
    parameters += "rr_2(dfig) = rr(dfig);\n"
    parameters += "To_p_1(dfig) = To_p(dfig)*rr(dfig)/(rr(dfig) + r_crowbar(dfig));\n"
    parameters += "To_p_2(dfig) = To_p(dfig);\n"
    parameters += "Kptrq(dfig) = Kptrq(dfig)*Sb_w(dfig)/Sb;\n"
    parameters += "Kitrq(dfig) = Kitrq(dfig)*Sb_w(dfig)/Sb;\n"

    return parameters

def write_tables(Y, Y_fault_reduced, Y_postfault_reduced, 
                 buses, retained_bus_numbers):

    tables = "\n\nTable Y(b,b)\n"
    tables = tablesFunc_write_ybus(tables,      # string
                                   Y,           # array to be writen n x n
                                   buses,       # list of numbers or list of instances of the bus class; n-elements
                                   "absolute")  # string; "absolute" or "angle"

    tables += "\n\nTable theta(b,b)\n"
    tables = tablesFunc_write_ybus(tables,
                                   Y,
                                   buses,
                                   "angle")

    tables += "\n\nTable Y_f_red(rb,rb)\n"
    tables = tablesFunc_write_ybus(tables,
                                   Y_fault_reduced,
                                   retained_bus_numbers,
                                   "absolute")

    tables += "\n\nTable theta_f_red(rb,rb)\n"
    tables = tablesFunc_write_ybus(tables,
                                   Y_fault_reduced,
                                   retained_bus_numbers,
                                   "angle")

    tables += "\n\nTable Y_pf_red(rb,rb)\n"
    tables = tablesFunc_write_ybus(tables,
                                   Y_postfault_reduced,
                                   retained_bus_numbers,
                                   "absolute")

    tables += "\n\nTable theta_pf_red(rb,rb)\n"
    tables = tablesFunc_write_ybus(tables,
                                   Y_postfault_reduced,
                                   retained_bus_numbers,
                                   "angle")

    return tables

def write_equations(generators,
                    branches,
                    dfig_param,
                    opf_param,
                    s_base):
    
    cost = def_cost_eq(generators, 
                       dfig_param)
    power_flow = def_power_flow_eq(generators,
                                   dfig_param,
                                   opf_param)
    auxiliary = def_auxiliary_eq(opf_param, 
                                 generators)
    initial_conditions = def_initial_conditions_eq(opf_param,
                                                   generators)
    field_voltage = def_field_voltage_eq()
    internal_gen_voltage = def_internal_gen_voltage_eq()
    mechanical_eq = def_mechanical_eq(generators)
    dq_gen_currents = def_dq_gen_currents_eq()
    stability_criterion = def_stability_criterion_eq()
    dfig_eq = def_dfig_eq(dfig_param)
    power_output_eq = def_power_output_eq()
    branch_restrictions_eq = def_branch_restrictions_eq()

    total_equations = [cost,
                       power_flow,
                       branch_restrictions_eq,
                       auxiliary,
                       initial_conditions,
                       field_voltage,
                       internal_gen_voltage,
                       mechanical_eq,
                       power_output_eq,
                       stability_criterion,
                       dq_gen_currents,
                       dfig_eq]
    
    variables = []
    var_restrictions = []
    equation_names = []
    equations = []

    for equation in total_equations:
        for variable in equation.variables:
            if variable not in variables:
                variables.append(variable)
        for var_restriction in equation.var_restrictions:
            if var_restriction not in var_restrictions:
                var_restrictions.append(var_restriction)
        for equation_name in equation.equation_names:
            if equation_name not in equation_names:
                equation_names.append(equation_name)
        for equation in equation.equations:
            if equation not in equations:
                equations.append(equation)

    string_variables = "\n\nVariables\n\n"
    string_var_restrictions = "\n"
    string_equation_names = "\n\nEquations\n\n"
    string_equations = "\n"
    
    counter = 0
    for variable in variables:
        string_variables += variable
        counter += 1
        if counter == len(variables): string_variables += ";\n"
        else: string_variables += "\n"

    for var_restriction in var_restrictions:
        string_var_restrictions += var_restriction + ";\n"

    counter = 0
    for equation_name in equation_names:
        string_equation_names += equation_name
        counter += 1
        if counter == len(equation_names): string_equation_names += ";\n"
        else: string_equation_names += "\n"

    for equation in equations:
        string_equations += equation + ";\n"

    return string_variables, string_var_restrictions, string_equation_names, string_equations




###
### writing sub_module
###

def setsFunc_bus_writing(sets, bus_list):  # Efficient writing
    
    sets += str(bus_list[0])

    if len(bus_list) <= 5:
        for i in range (1, len(bus_list)):
            sets += ", " + str(bus_list[i])
    else:
        i = 1
        sequence = False
        initiaL_sequence = False
        aux = 0
        while i <= len(bus_list) - 1:
            if i < len(bus_list) - 1:
                while bus_list[i + 1] == bus_list[i] + 1:
                    if sequence == False:
                        aux = bus_list[i]
                        if i == 1: initiaL_sequence = True
                    i += 1
                    sequence = True
                    if i == len(bus_list) - 1: break
            else:
                if bus_list[-2] == bus_list[-1]:
                    if sequence == False:
                        aux = bus_list[i]
                    i += 1
                    sequence = True
                    if i == len(bus_list) - 1: break
            
            if sequence == False:
                sets += ", " + str(bus_list[i])
            else:
                if initiaL_sequence == True:  # if the initial bus_number is part of the first sequence
                    sets += " * " + str(bus_list[i])
                else:
                    sets += ", " + str(aux) + " * " + str(bus_list[i])

            if i == len(bus_list) - 1: break
            i += 1
            sequence = False

    return sets

def parametersFun_gen_writing(parameters, generators, attribute):  # Efficient writing

    counter = 0
    for generator in generators:
            parameters += str(generator.bus_number)+ " " + str(getattr(generator, attribute))
            if counter == len(generators) - 1: break
            else: parameters += ", "
            counter += 1

    return parameters

def parametersFun_dfig_writing(parameters, dfig_param, attribute):  # Efficient writing

    parameters += str(dfig_param['bus_number'])+ " " + str(dfig_param[attribute])

    return parameters

def tablesFunc_write_ybus(tables, Ybus, bus_list, decision):

    Y = np.zeros((Ybus.shape[0] + 1, Ybus.shape[0] + 1))

    if decision == "absolute":
        Y[1:, 1:] = np.absolute(Ybus)
    if decision == "angle":
        Y[1:, 1:] = np.angle(Ybus)

    if isinstance(bus_list[0], int):  # then bus_list is no a list of instances of the bus class
        counter = 1
        for bus_number in bus_list:
            Y[0, counter] = bus_number
            Y[counter, 0] = bus_number
            counter += 1
    else:
        counter = 1
        for bus in bus_list:
            Y[0, counter] = bus.number
            Y[counter, 0] = bus.number
            counter += 1
        
    #  Write the rest of the Ybus elements as floats
    for i in range(0, Y.shape[0]):
        for j in range(0, Y.shape[0]):
            if i == 0 and j == 0:
                tables += "\t" 
            else:
                if i == 0 or j == 0:
                    tables += "%10i " %int(Y[i, j])
                else:
                    tables += "%10.5f " %Y[i, j]

        if (i == Y.shape[0] - 1) and (j == Y.shape[0] - 1):
            tables += ";\n"
            break
        else:
            tables += "\n"

    return tables

class Gams_equation:

    def __init__(self,
                 name,
                 variables,
                 var_restrictions,
                 equation_names,
                 equations):

        self.name = name  # string
        for variable in variables:
            variable.replace(" ", "")
        self.variables = variables
        for var_restriction in var_restrictions:
            var_restriction.replace(" ", "")
            var_restriction.replace("=", " = ")
        self.var_restrictions = var_restrictions
        for equation_name in equation_names:
            equation_name.replace(" ", "")
        self.equation_names = equation_names
        self.equations = equations

def def_power_flow_eq(generators,
                      dfig_param,
                      opf_param):
    
    name = 'Power Flow'  # Name of the equations to be defined
    variables = []       # List of strings; varibles to be added in GAMS format
    restrictions = []    # List of strings; restriction to be added in GAMS format (to the variables already defined)
    equation_names = []  # List of strings; equation names defined before the equations definition
    equations = []       # List of strings; equations to be added in GAMS format (involving all "variables" in variables list)

    #  Adding Pg(rb) variable and restrictions associated to it
    variables.append("Pg(rb)")
    for generator in generators:
        initial_real_power = "Pg.l('" + str(generator.bus_number) + "') = " + str(generator.P_ini) + "/Sb"
        restrictions.append(initial_real_power)

    variables.append("Qg(rb)")
    for generator in generators:
        initial_imag_power = "Qg.l('" + str(generator.bus_number) + "') = " + str(generator.Q_ini) + "/Sb"
        restrictions.append(initial_real_power)

    restrictions.append("Pg.lo(sg) = Pmin(sg)")
    restrictions.append("Pg.up(sg) = Pmax(sg)")
    restrictions.append("Qg.lo(sg) = Qmin(sg)")
    restrictions.append("Qg.up(sg) = Qmax(sg)")
    
    variables.append("V(b)")
    restrictions.append("V.l(b) = " + str(float(opf_param['v_ini'])))
    restrictions.append("V.lo(b) = " + str(float(opf_param['v_min'])))
    restrictions.append("V.up(b) = " + str(float(opf_param['v_max'])))

    variables.append("alpha(b)")                  
    restrictions.append("alpha.l(b) = 0")
    restrictions.append("alpha.lo(b) = -pi")
    restrictions.append("alpha.up(b) = pi")

    variables.append("Ps_w(dfig,s)")
    variables.append("Qs_w(dfig,s)")
    variables.append("P_gsc(dfig,s)")
    variables.append("P_mech(dfig,s)")
    variables.append("VDr(dfig,s)")
    variables.append("VQr(dfig,s)")
    variables.append("IDr(dfig,s)")
    variables.append("IQr(dfig,s)")
    variables.append("iDs_w(dfig,s)")
    variables.append("iQs_w(dfig,s)")
    variables.append("eD_p_w(dfig,s)")
    variables.append("eQ_p_w(dfig,s)")
    
    equation_names = ["p_balance_gen",
                      "p_balance_nongen",
                      "q_balance_gen",
                      "q_balance_nongen",
                      "ref_bus",
                      "dfig_load_flow_f1",
                      "dfig_load_flow_f2",
                      "dfig_load_flow_f3",
                      "dfig_load_flow_f4",
                      "dfig_load_flow_f5",
                      "dfig_load_flow_f6",
                      "P_mech_ini"]

    equations =  ["p_balance_gen(rb) .. Pg(rb) - Pl(rb) - V(rb)*sum(b,V(b)*Y(rb,b)*cos(alpha(rb) - alpha(b) - theta(rb,b))) =e= 0",
                  "p_balance_nongen(ngb) ..  - Pl(ngb) - V(ngb)*sum(b,V(b)*Y(ngb,b)*cos(alpha(ngb) - alpha(b) - theta(ngb,b))) =e= 0",
                  "q_balance_gen(rb) .. Qg(rb) - Ql(rb) - V(rb)*sum(b,V(b)*Y(rb,b)*sin(alpha(rb) - alpha(b) - theta(rb,b))) =e= 0",
                  "q_balance_nongen(ngb) ..  - Ql(ngb) - V(ngb)*sum(b,V(b)*Y(ngb,b)*sin(alpha(ngb) - alpha(b) - theta(ngb,b))) =e= 0",
                  "ref_bus .. alpha('"+ str(int(opf_param['slack_bus'])) + "') =e= 0",
                  "dfig_load_flow_f1(dfig) .. Pg(dfig)/Ngen_w(dfig) =e= Ps_w(dfig,'1') - P_gsc(dfig,'1')",
                  "dfig_load_flow_f2(dfig) .. Qg(dfig)/Ngen_w(dfig) =e= Qs_w(dfig,'1')",
                  "dfig_load_flow_f3(dfig) .. VDr(dfig,'1') =e= IDr(dfig,'1') + w0*To_p(dfig)*delta_wr(dfig,'1')*eD_p_w(dfig,'1')",
                  "dfig_load_flow_f4(dfig) .. VQr(dfig,'1') =e= IQr(dfig,'1') - w0*To_p(dfig)*delta_wr(dfig,'1')*eQ_p_w(dfig,'1')",
                  "dfig_load_flow_f5(dfig) .. V(dfig)*cos(alpha(dfig)) =e= xs_p(dfig)*iQs_w(dfig,'1') - eD_p_w(dfig,'1')",
                  "dfig_load_flow_f6(dfig) .. V(dfig)*sin(alpha(dfig)) =e= - xs_p(dfig)*iDs_w(dfig,'1') + eQ_p_w(dfig,'1')",
                  "P_mech_ini(dfig) .. Pg(dfig)/Ngen_w(dfig) =e= P_mech(dfig,'1') - rr(dfig)/sqr(lm(dfig))*(sqr(IDr(dfig,'1')) + sqr(IQr(dfig,'1')))"]

    gams_equation_instance = Gams_equation(name,
                                           variables,
                                           restrictions,
                                           equation_names,
                                           equations)

    return gams_equation_instance


def def_branch_restrictions_eq():

    name = 'Branch restrictions' 
    variables = []
    restrictions = []
    equation_names = []
    equations = []

    equation_names = ["I_brach_lim_inf_1",
                      "I_brach_lim_inf_2",
                      "I_brach_lim_sup_1",
                      "I_brach_lim_sup_2"]
      
    equations =  ["I_brach_lim_inf_1('1') .. 0 =l= (sqr(V('13')*cos(alpha('13')) - V('14')*cos(alpha('14'))) + sqr(V('13')*sin(alpha('13')) - V('14')*sin(alpha('14'))))*sqr(Y('13','14'))",
                  "I_brach_lim_inf_2('1') .. 0 =l= (sqr(V('15')*cos(alpha('15')) - V('16')*cos(alpha('16'))) + sqr(V('15')*sin(alpha('15')) - V('16')*sin(alpha('16'))))*sqr(Y('15','16'))",
                  "I_brach_lim_sup_1('1') .. (sqr(V('13')*cos(alpha('13')) - V('14')*cos(alpha('14'))) + sqr(V('13')*sin(alpha('13')) - V('14')*sin(alpha('14'))))*sqr(Y('13','14')) =l= sqr(2.0)",
                  "I_brach_lim_sup_2('1') .. (sqr(V('15')*cos(alpha('15')) - V('16')*cos(alpha('16'))) + sqr(V('15')*sin(alpha('15')) - V('16')*sin(alpha('16'))))*sqr(Y('15','16')) =l= sqr(1.8)"]
    
    gams_equation_instance = Gams_equation(name,
                                           variables,
                                           restrictions,
                                           equation_names,
                                           equations)  

    return gams_equation_instance


def def_auxiliary_eq(opf_param, 
                     generators):

    name = 'Auxiliary equations' 
    variables = []
    restrictions = []
    equation_names = []
    equations = []

    variables.append("Ig(sg)")
    restrictions.append("Ig.l(sg) = 1")
    restrictions.append("Ig.lo(sg) = 0.001")
    restrictions.append("Ig.up(sg) = Sm(sg)")

    variables.append("phi(sg)")
    restrictions.append("phi.l(sg) = 0")
    restrictions.append("phi.lo(sg) = -pi/2")
    restrictions.append("phi.up(sg) = pi/2")

    variables.append("V(b)")
    restrictions.append("V.l(b) = " + str(float(opf_param['v_ini'])))
    restrictions.append("V.lo(b) = " + str(float(opf_param['v_min'])))
    restrictions.append("V.up(b) = " + str(float(opf_param['v_max'])))

    variables.append("Pg(rb)")
    for generator in generators:
        initial_real_power = "Pg.l('" + str(generator.bus_number) + "') = " + str(generator.P_ini) + "/Sb"
        restrictions.append(initial_real_power)

    variables.append("Qg(rb)")
    for generator in generators:
        initial_imag_power = "Qg.l('" + str(generator.bus_number) + "') = " + str(generator.Q_ini) + "/Sb"
        restrictions.append(initial_imag_power)

    restrictions.append("Pg.lo(sg) = Pmin(sg)")
    restrictions.append("Pg.up(sg) = Pmax(sg)")
    restrictions.append("Qg.lo(sg) = Qmin(sg)")
    restrictions.append("Qg.up(sg) = Qmax(sg)")                
    
    equation_names = ["Generator_currents",
                      "Power_factor"]

    equations =  ["Generator_currents(sg) .. sqr(Ig(sg)*V(sg)) - sqr(Pg(sg)) - sqr(Qg(sg)) =e= 0",
                  "Power_factor(sg) .. sin(phi(sg)) - Qg(sg)/(V(sg)*Ig(sg)) =e= 0"]


    gams_equation_instance = Gams_equation(name,
                                           variables,
                                           restrictions,
                                           equation_names,
                                           equations)

    return gams_equation_instance

def def_initial_conditions_eq(opf_param,
                              generators):

    name = 'Initial condition equations' 
    variables = []
    restrictions = []
    equation_names = []
    equations = []

    variables.append("Ig(sg)")
    restrictions.append("Ig.l(sg) = 1")
    restrictions.append("Ig.lo(sg) = 0.001")
    restrictions.append("Ig.up(sg) = Sm(sg)")

    variables.append("Pg(rb)")
    for generator in generators:
        initial_real_power = "Pg.l('" + str(generator.bus_number) + "') = " + str(generator.P_ini) + "/Sb"
        restrictions.append(initial_real_power)

    restrictions.append("Pg.lo(sg) = Pmin(sg)")
    restrictions.append("Pg.up(sg) = Pmax(sg)")

    variables.append("phi(sg)")
    restrictions.append("phi.l(sg) = 0")
    restrictions.append("phi.lo(sg) = -pi/2")
    restrictions.append("phi.up(sg) = pi/2")

    variables.append("alpha(b)")                  
    restrictions.append("alpha.l(b) = 0")
    restrictions.append("alpha.lo(b) = -pi")
    restrictions.append("alpha.up(b) = pi")
    
    variables.append("delta(sg,s)")
    restrictions.append("delta.l(sg,s)= 0")
    restrictions.append("delta.lo(sg,s) = -9999")
    restrictions.append("delta.up(sg,s) = 9999")

    variables.append("V(b)")
    restrictions.append("V.l(b) = " + str(float(opf_param['v_ini'])))
    restrictions.append("V.lo(b) = " + str(float(opf_param['v_min'])))
    restrictions.append("V.up(b) = " + str(float(opf_param['v_max'])))

    variables.append("ed_p(sg,s)")
    restrictions.append("ed_p.l(sg,s) = 0.2")
    restrictions.append("ed_p.lo(sg,s) = 0") 
    restrictions.append("ed_p.up(sg,s) = 1.5")

    variables.append("eq_p(sg,s)")
    restrictions.append("eq_p.l(sg,s) = 1.0")
    restrictions.append("eq_p.lo(sg,s) = 0")
    restrictions.append("eq_p.up(sg,s) = 1.5")

    variables.append("id(sg,s)")
    restrictions.append("id.l(sg,s) = 1")
    restrictions.append("id.lo(sg,s) = -Sm(sg)")
    restrictions.append("id.up(sg,s) = 3*Sm(sg)")

    variables.append("iq(sg,s)")
    restrictions.append("iq.l(sg,s) = 1")
    restrictions.append("iq.lo(sg,s) = -Sm(sg)")
    restrictions.append("iq.up(sg,s) = 3*Sm(sg)")

    variables.append("e_fd(sg)")
    restrictions.append("e_fd.l(sg) = 1")
    
    variables.append("Domega(sg,s)")
    restrictions.append("Domega.l(sg,s) = 0")
    restrictions.append("Domega.lo(sg,s) = -1")
    restrictions.append("Domega.up(sg,s) = 1")

    variables.append("Pe(sg,s)")
    restrictions.append("Pe.l(sg,s) = 3.0")
    restrictions.append("Pe.lo(sg,s) = -99")
    restrictions.append("Pe.up(sg,s) = 99")

    
    equation_names = ["ed_p_initialization",
                      "eq_p_initialization",
                      "Vd_initialization",
                      "Vq_initialization",
                      "Domega_initialization",
                      "id_initialization",
                      "iq_initialization"]

    equations =  ["ed_p_initialization(sg) .. ed_p(sg,'1') - (Xq(sg) - Xq_p(sg))*Ig(sg)*cos(delta(sg,'1') - alpha(sg) + phi(sg)) =e= 0",
                  "eq_p_initialization(sg) .. eq_p(sg,'1') + (Xd(sg) - Xd_p(sg))*Ig(sg)*sin(delta(sg,'1') - alpha(sg) + phi(sg)) - e_fd(sg) =e= 0",
                  "Vd_initialization(sg) .. V(sg)*sin(delta(sg,'1') - alpha(sg)) - ed_p(sg,'1') + (Ra(sg)*sin(delta(sg,'1') - alpha(sg) +\n\
                   phi(sg)) - Xq_p(sg)*cos(delta(sg,'1') - alpha(sg) + phi(sg)))*Ig(sg) =e= 0",
                  "Vq_initialization(sg) .. V(sg)*cos(delta(sg,'1') - alpha(sg)) - eq_p(sg,'1') + (Ra(sg)*cos(delta(sg,'1') - alpha(sg) +\n\
                   phi(sg)) + Xd_p(sg)*sin(delta(sg,'1') - alpha(sg) + phi(sg)))*Ig(sg) =e= 0",
                  "Domega_initialization(sg) .. Domega(sg,'1') =e= 0",
                  "id_initialization(sg) .. id(sg,'1') - Ig(sg)*sin(delta(sg,'1') - alpha(sg) + phi(sg)) =e= 0",
                  "iq_initialization(sg) .. iq(sg,'1') - Ig(sg)*cos(delta(sg,'1') - alpha(sg) + phi(sg)) =e= 0"]


    gams_equation_instance = Gams_equation(name,
                                           variables,
                                           restrictions,
                                           equation_names,
                                           equations)

    return gams_equation_instance

def def_field_voltage_eq():

    name = 'Field voltage limit' 
    variables = []
    restrictions = []
    equation_names = []
    equations = []

    variables.append("e_fd(sg)")
    restrictions.append("e_fd.l(sg) = 1")

    equation_names = ["e_fd_lim_inf",
                      "e_fd_lim_sup"]

    equations =  ["e_fd_lim_inf(sg) .. e_fd_min(sg) =l= e_fd(sg)",
                  "e_fd_lim_sup(sg) .. e_fd(sg) =l= e_fd_max(sg)"]

    gams_equation_instance = Gams_equation(name,
                                           variables,
                                           restrictions,
                                           equation_names,
                                           equations)

    return gams_equation_instance


def def_internal_gen_voltage_eq():

    name = 'Internal d-q genrator voltage equations' 
    variables = []
    restrictions = []
    equation_names = []
    equations = []

    variables.append("e_fd(sg)")
    restrictions.append("e_fd.l(sg) = 1")

    variables.append("id(sg,s)")
    restrictions.append("id.l(sg,s) = 1")
    restrictions.append("id.lo(sg,s) = -Sm(sg)")
    restrictions.append("id.up(sg,s) = 3*Sm(sg)")

    variables.append("iq(sg,s)")
    restrictions.append("iq.l(sg,s) = 1")
    restrictions.append("iq.lo(sg,s) = -Sm(sg)")
    restrictions.append("iq.up(sg,s) = 3*Sm(sg)")
    
    equation_names = ["Internal_voltaje_d",
                      "Internal_voltaje_q"]

    equations =  ["Internal_voltaje_d(sg,s)$(not sfirst(s)) .. ed_p(sg,s)*(1 + Dt/(2*Tq_p(sg))) - \n\
                   ed_p(sg,s-1)*(1 - Dt/(2*Tq_p(sg))) - (Dt/(2*Tq_p(sg)))*(Xq(sg) - Xq_p(sg))*(iq(sg,s) + iq(sg,s-1)) =e= 0",
                  "Internal_voltaje_q(sg,s)$(not sfirst(s)) .. eq_p(sg,s)*(1 + Dt/(2*Td_p(sg))) - \n\
                   eq_p(sg,s-1)*(1 - Dt/(2*Td_p(sg))) - (Dt/(2*Td_p(sg)))*(2*e_fd(sg) - (Xd(sg) - Xd_p(sg))*(id(sg,s) + id(sg,s-1))) =e= 0"]


    gams_equation_instance = Gams_equation(name,
                                           variables,
                                           restrictions,
                                           equation_names,
                                           equations)

    return gams_equation_instance


def def_power_output_eq():

    name = 'Power output equations' 
    variables = []
    restrictions = []
    equation_names = []
    equations = []
    
    variables.append("delta(sg,s)")
    restrictions.append("delta.l(sg,s)= 0")
    restrictions.append("delta.lo(sg,s) = -9999")
    restrictions.append("delta.up(sg,s) = 9999")

    variables.append("ed_p(sg,s)")
    restrictions.append("ed_p.l(sg,s) = 0.2")
    restrictions.append("ed_p.lo(sg,s) = 0")
    restrictions.append("ed_p.up(sg,s) = 1.5")

    variables.append("eq_p(sg,s)")
    restrictions.append("eq_p.l(sg,s) = 1.0")
    restrictions.append("eq_p.lo(sg,s) = 0")
    restrictions.append("eq_p.up(sg,s) = 1.5")

    variables.append("id(sg,s)")
    restrictions.append("id.l(sg,s) = 1")
    restrictions.append("id.lo(sg,s) = -Sm(sg)")
    restrictions.append("id.up(sg,s) = 3*Sm(sg)")

    variables.append("iq(sg,s)")
    restrictions.append("iq.l(sg,s) = 1")
    restrictions.append("iq.lo(sg,s) = -Sm(sg)")
    restrictions.append("iq.up(sg,s) = 3*Sm(sg)")

    variables.append("Pe(sg,s)")
    restrictions.append("Pe.l(sg,s) = 3.0")
    restrictions.append("Pe.lo(sg,s) = -99")
    restrictions.append("Pe.up(sg,s) = 99")

    equation_names = ["electric_power"]

    equations =  ["electric_power(sg,s) .. Pe(sg,s) - ed_p(sg,s)*id(sg,s) - eq_p(sg,s)*iq(sg,s) =e= 0"]


    gams_equation_instance = Gams_equation(name,
                                           variables,
                                           restrictions,
                                           equation_names,
                                           equations)

    return gams_equation_instance


def def_mechanical_eq(generators):

    name = 'Mechanical equations' 
    variables = []
    restrictions = []
    equation_names = []
    equations = []

    variables.append("Domega(sg,s)")
    restrictions.append("Domega.l(sg,s) = 0")
    restrictions.append("Domega.lo(sg,s) = -1")
    restrictions.append("Domega.up(sg,s) = 1")

    variables.append("Pg(rb)")
    for generator in generators:
        initial_real_power = "Pg.l('" + str(generator.bus_number) + "') = " + str(generator.P_ini) + "/Sb"
        restrictions.append(initial_real_power)

    restrictions.append("Pg.lo(sg) = Pmin(sg)")
    restrictions.append("Pg.up(sg) = Pmax(sg)")

    variables.append("delta(sg,s)")
    restrictions.append("delta.l(sg,s)= 0")
    restrictions.append("delta.lo(sg,s) = -9999")
    restrictions.append("delta.up(sg,s) = 9999")

    variables.append("Pe(sg,s)")
    restrictions.append("Pe.l(sg,s) = 3.0")
    restrictions.append("Pe.lo(sg,s) = -99")
    restrictions.append("Pe.up(sg,s) = 99")
    
    equation_names = ["oscilation_omega",
                      "oscilation_delta"]

    equations =  ["oscilation_omega(sg,s)$(not sfirst(s)) .. Domega(sg,s)*(1 + Dt*D(sg)/(4*H(sg))) - \n\
                   Domega(sg,s-1)*(1 - Dt*D(sg)/(4*H(sg))) - (Dt/(4*H(sg)))*(2*Pg(sg) - Pe(sg,s) - Pe(sg,s-1)) =e= 0",
                  "oscilation_delta(sg,s)$(not sfirst(s)) .. delta(sg,s) - delta(sg,s-1) - \n\
                   (Dt*100*pi/2)*(Domega(sg,s) + Domega(sg,s-1)) =e= 0"]


    gams_equation_instance = Gams_equation(name,
                                           variables,
                                           restrictions,
                                           equation_names,
                                           equations)

    return gams_equation_instance

def def_stability_criterion_eq():

    name = 'Stability criterion equations' 
    variables = []
    restrictions = []
    equation_names = []
    equations = []

    variables.append("delta(sg,s)")
    restrictions.append("delta.l(sg,s)= 0")
    restrictions.append("delta.lo(sg,s) = -9999")
    restrictions.append("delta.up(sg,s) = 9999")

    variables.append("delta_COI(s)")
    restrictions.append("delta_COI.l(s) = 0")
    restrictions.append("delta_COI.lo(s) = -9999")
    restrictions.append("delta_COI.up(s) = 9999")
    
    equation_names = ["center_of_inertia",
                      "angular_deviation_min",
                      "angular_deviation_max"]

    equations =  ["center_of_inertia(s) .. delta_COI(s) - sum(sg,H(sg)*delta(sg,s)) / sum(sg,H(sg)) =e= 0",
                  "angular_deviation_min(sg,s) .. - (angle_limit*pi/180) =l= delta(sg,s) - delta_COI(s)",
                  "angular_deviation_max(sg,s) .. delta(sg,s) - delta_COI(s) =l= (angle_limit*pi/180)"]

    gams_equation_instance = Gams_equation(name,
                                           variables,
                                           restrictions,
                                           equation_names,
                                           equations)

    return gams_equation_instance


def def_cost_eq(generators, 
                dfig_param):

    name = 'Cost equation' 
    variables = []
    restrictions = []
    equation_names = []
    equations = []

    variables.append("Pg(rb)")
    for generator in generators:
        initial_real_power = "Pg.l('" + str(generator.bus_number) + "') = " + str(generator.P_ini) + "/Sb"
        restrictions.append(initial_real_power)

    variables.append("z")
    
    equation_names = ["total_cost"]

    equations =  ["total_cost .. z =e= sum(sg,a1(sg)*(Pg(sg)*Sb)) + sum(dfig,a1(dfig)*(Pg(dfig)*Sb))"]

    gams_equation_instance = Gams_equation(name,
                                           variables,
                                           restrictions,
                                           equation_names,
                                           equations)

    return gams_equation_instance


def def_dfig_eq(dfig_param):

    name = 'DFIG equations' 
    variables = []
    restrictions = []
    equation_names = []
    equations = []

    variables.append('rr_var(dfig,s)')
    variables.append('To_p_var(dfig,s)')
    variables.append('T_el(dfig,s)')
    variables.append('T_mech(dfig,s)')
    variables.append('wr(dfig,s)')
    variables.append('delta_wr(dfig,s)')
    variables.append('P_mech(dfig,s)')
    variables.append('lambda(dfig,s)')
    variables.append('Cp(dfig,s)')
    variables.append('vD_w(dfig,s)')
    variables.append('vQ_w(dfig,s)')
    variables.append('v_w(dfig,s)')
    variables.append('alpha_w(dfig,s)')
    variables.append('VDr(dfig,s)')
    variables.append('VQr(dfig,s)')
    variables.append('IDr(dfig,s)')
    variables.append('IQr(dfig,s)')
    variables.append('IDr_c(dfig,sdfig_2)')
    variables.append('IQr_c(dfig,sdfig_2)')
    variables.append('eD_p_w(dfig,s)')
    variables.append('eQ_p_w(dfig,s)')
    variables.append('iDs_w(dfig,s)')
    variables.append('iQs_w(dfig,s)')
    variables.append('iD_w(dfig,s)')
    variables.append('iQ_w(dfig,s)')
    variables.append('wr_ref_var(dfig,s)')
    variables.append('error_wr(dfig,s)')
    variables.append('T_el_ref(dfig,s)')
    variables.append('P_ord_pre(dfig,s)')
    variables.append('P_ord(dfig,s)')
    variables.append('pitch_angle(dfig,s)')
    variables.append('pll_angle(dfig,s)')
    variables.append('Ps_w_ref(dfig,sdfig_2)')
    variables.append('iDr_c_ref(dfig,sdfig_2)')
    variables.append('iQr_c_ref(dfig,sdfig_2)')
    variables.append('P_gsc(dfig,s)')
    variables.append('Ps_w(dfig,s)')
    variables.append('Qs_w(dfig,s)')
    variables.append('Pg_w(dfig,s)')
    variables.append('Qg_w(dfig,s)')

    restrictions.append("v_w.l(dfig,s) = 1") 
    restrictions.append("v_w.lo(dfig,s) = 0.1")
    restrictions.append("rr_var.fx(dfig,'1') = rr(dfig)")
    restrictions.append("rr_var.fx(dfig,sdfig_1) = rr_1(dfig)")
    restrictions.append("rr_var.fx(dfig,sdfig_2) = rr_2(dfig)")
    restrictions.append("To_p_var.fx(dfig,'1') = To_p(dfig)")
    restrictions.append("To_p_var.fx(dfig,sdfig_1) = To_p_1(dfig)")
    restrictions.append("To_p_var.fx(dfig,sdfig_2) = To_p_2(dfig)")
    restrictions.append("wr_ref_var.fx(dfig,s) = wr_ref(dfig)")
    restrictions.append("delta_wr.fx(dfig,'1') = wr_ref(dfig) - 1")
    restrictions.append("pitch_angle.fx(dfig,'1') = initial_pitch_angle(dfig)")
    restrictions.append("vDr.fx(dfig,sdfig_1) = 0")
    restrictions.append("vQr.fx(dfig,sdfig_1) = 0")
    
    equation_names = ["iD_w_fault",
                      "iQ_w_fault",
                      "iD_w_postfault",
                      "iQ_w_postfault",
                      "P_mech_ini",
                      "eD_p_w_equation",
                      "eQ_p_w_equation",
                      "rotor_flux_D",
                      "rotor_flux_Q",
                      "stator_equation_D",
                      "stator_equation_Q",
                      "T_el_equation",
                      "delta_wr_equation",
                      "wr_equation",
                      "lambda_equation",
                      "Cp_equation",
                      "P_mech_equation",
                      "T_mech_equation",
                      "error_wr_equation",
                      "T_el_ref_ini",
                      "T_el_ref_equation",
                      "P_ord_pre_equation",
                      "P_ord_ini",
                      "P_ord_equation",
                      "pitch_angle_1",
                      "pitch_angle_2",
                      "vD_w_equation",
                      "vQ_w_equation",
                      "pll_angle_ini",
                      "Ps_w_ref_equation",
                      "iQr_c_ref_equation",
                      "IDr_c_ref_2",
                      "IDr_equation",
                      "IQr_equation",
                      "IDr_c_2",
                      "IQr_c_2",
                      "P_gsc_equation",
                      "Ps_w_equation",
                      "Qs_w_equation",
                      "Pg_w_equation",
                      "Qg_w_equation",
                      "P_balance",
                      "Q_balance"]

    equations =  ["iD_w_fault(dfig,sf) .. iD_w(dfig,sf)*Ngen_w(dfig)\n\
                                        - (sum(sg,Y_f_red(dfig,sg)*(ed_p(sg,sf)*sin(delta(sg,sf) + theta_f_red(dfig,sg))\n\
                                                                  + eq_p(sg,sf)*cos(delta(sg,sf) + theta_f_red(dfig,sg))))\n\
                                        +  sum(dfigp,Y_f_red(dfig,dfigp)*(vD_w(dfigp,sf)*cos(theta_f_red(dfig,dfigp))\n\
                                                                        - vQ_w(dfigp,sf)*sin(theta_f_red(dfig,dfigp))))) =e= 0",
                  "iQ_w_fault(dfig,sf) .. iQ_w(dfig,sf)*Ngen_w(dfig)\n\
                                        - (sum(sg,Y_f_red(dfig,sg)*(- ed_p(sg,sf)*cos(delta(sg,sf) + theta_f_red(dfig,sg))\n\
                                                                    + eq_p(sg,sf)*sin(delta(sg,sf) + theta_f_red(dfig,sg))))\n\
                                        +  sum(dfigp,Y_f_red(dfig,dfigp)*(vD_w(dfigp,sf)*sin(theta_f_red(dfig,dfigp))\n\
                                                                        + vQ_w(dfigp,sf)*cos(theta_f_red(dfig,dfigp))))) =e= 0",
                  "iD_w_postfault(dfig,spf) .. iD_w(dfig,spf)*Ngen_w(dfig)\n\
                                             - (sum(sg,Y_pf_red(dfig,sg)*(ed_p(sg,spf)*sin(delta(sg,spf) + theta_pf_red(dfig,sg))\n\
                                                                        + eq_p(sg,spf)*cos(delta(sg,spf) + theta_pf_red(dfig,sg))))\n\
                                             +  sum(dfigp,Y_pf_red(dfig,dfigp)*(vD_w(dfigp,spf)*cos(theta_pf_red(dfig,dfigp))\n\
                                                                              - vQ_w(dfigp,spf)*sin(theta_pf_red(dfig,dfigp))))) =e= 0",
                  "iQ_w_postfault(dfig,spf) .. iQ_w(dfig,spf)*Ngen_w(dfig)\n\
                                             - (sum(sg,Y_pf_red(dfig,sg)*(- ed_p(sg,spf)*cos(delta(sg,spf) + theta_pf_red(dfig,sg))\n\
                                                                          + eq_p(sg,spf)*sin(delta(sg,spf) + theta_pf_red(dfig,sg))))\n\
                                             +  sum(dfigp,Y_pf_red(dfig,dfigp)*(vD_w(dfigp,spf)*sin(theta_pf_red(dfig,dfigp))\n\
                                                                              + vQ_w(dfigp,spf)*cos(theta_pf_red(dfig,dfigp))))) =e= 0",
                  "eD_p_w_equation(dfig,s)$(ord(s) gt 1) .. eD_p_w(dfig,s)*To_p_var(dfig,s) =e= eD_p_w(dfig,s-1)*To_p_var(dfig,s-1) + Dt/2*((VQr(dfig,s) + VQr(dfig,s-1)) -\n\
                           (IQr(dfig,s) + IQr(dfig,s-1)) + w0*(To_p_var(dfig,s)*eQ_p_w(dfig,s)*delta_wr(dfig,s) + To_p_var(dfig,s-1)*eQ_p_w(dfig,s-1)*delta_wr(dfig,s-1)))",
                  "eQ_p_w_equation(dfig,s)$(ord(s) gt 1) .. eQ_p_w(dfig,s)*To_p_var(dfig,s) =e= eQ_p_w(dfig,s-1)*To_p_var(dfig,s-1) + Dt/2*((VDr(dfig,s) + VDr(dfig,s-1)) -\n\
                           (IDr(dfig,s) + IDr(dfig,s-1)) - w0*(To_p_var(dfig,s)*eD_p_w(dfig,s)*delta_wr(dfig,s) + To_p_var(dfig,s-1)*eD_p_w(dfig,s-1)*delta_wr(dfig,s-1)))",
                  "rotor_flux_D(dfig,s) .. IDr(dfig,s) =e= eQ_p_w(dfig,s) + (xss(dfig) - xs_p(dfig))*iDs_w(dfig,s)",
                  "rotor_flux_Q(dfig,s) .. IQr(dfig,s) =e= eD_p_w(dfig,s) + (xss(dfig) - xs_p(dfig))*iQs_w(dfig,s)",
                  "stator_equation_D(dfig,s) .. eD_p_w(dfig,s) =e= xs_p(dfig)*iQs_w(dfig,s) - vD_w(dfig,s)",
                  "stator_equation_Q(dfig,s) .. eQ_p_w(dfig,s) =e= vQ_w(dfig,s) + xs_p(dfig)*iDs_w(dfig,s)",
                  "vD_w_equation(dfig,s) .. vD_w(dfig,s) =e= v_w(dfig,s)*cos(alpha_w(dfig,s))",
                  "vQ_w_equation(dfig,s) .. vQ_w(dfig,s) =e= v_w(dfig,s)*sin(alpha_w(dfig,s))",
                  "T_el_equation(dfig,s) .. T_el(dfig,s) =e= iQs_w(dfig,s)*IDr(dfig,s) - iDs_w(dfig,s)*IQr(dfig,s)",
                  "delta_wr_equation(dfig,s)$(ord(s) gt 1) .. delta_wr(dfig,s) =e= delta_wr(dfig,s-1) + Dt/(4*H_w(dfig))*((T_mech(dfig,s-1) - T_el(dfig,s-1)) + (T_mech(dfig,s) - T_el(dfig,s)))",
                  "wr_equation(dfig,s) .. wr(dfig,s) =e= 1 + delta_wr(dfig,s)",
                  "lambda_equation(dfig,s) .. lambda(dfig,s) =e= Kb(dfig)*(wr(dfig,s)/u_wind(dfig))",
                  "Cp_equation(dfig,s) .. Cp(dfig,s) =e= sum(i,sum(j,alpha_ij(i,j)*power(pitch_angle(dfig,s),ord(i)-1)*power(lambda(dfig,s),ord(j)-1)))",
                  "P_mech_equation(dfig,s) .. P_mech(dfig,s) =e= Cp(dfig,s)*Kwind(dfig)*power(u_wind(dfig),3)*Sb_w(dfig)/Sb",
                  "T_mech_equation(dfig,s) .. T_mech(dfig,s)*wr(dfig,s) =e= P_mech(dfig,s)",
                  "error_wr_equation(dfig,s) .. error_wr(dfig,s) =e= wr(dfig,s) - wr_ref(dfig)",
                  "T_el_ref_ini(dfig,s)$(ord(s) eq 1) .. T_el_ref(dfig,s) =e= Ps_w(dfig,s)",
                  "T_el_ref_equation(dfig,s)$(ord(s) gt 1) .. T_el_ref(dfig,s) =e= T_el_ref(dfig,s-1)\n\
                                                                                      + Kptrq(dfig)*(error_wr(dfig,s) - error_wr(dfig,s-1))\n\
                                                                                      + Kitrq(dfig)*Dt/2*(error_wr(dfig,s) + error_wr(dfig,s-1))",
                  "P_ord_pre_equation(dfig,s) .. P_ord_pre(dfig,s) =e= T_el_ref(dfig,s)",
                  "P_ord_ini(dfig,s)$(ord(s) eq 1) .. P_ord(dfig,s) =e= P_ord_pre(dfig,s)",
                  "P_ord_equation(dfig,s)$(ord(s) gt 1) .. P_ord(dfig,s) =e= P_ord(dfig,s-1)*(2*T_pc(dfig) - Dt)/(2*T_pc(dfig) + Dt)\n\
                                                                            + Dt/(2*T_pc(dfig) + Dt)*(P_ord_pre(dfig,s) + P_ord_pre(dfig,s-1))",
                  "pitch_angle_1(dfig,spitch_1)$(ord(spitch_1) gt 1) .. pitch_angle(dfig,spitch_1) =e= pitch_angle(dfig,spitch_1-1) + 10*Dt",
                  "pitch_angle_2(dfig,spitch_2)$(ord(spitch_2) gt 1) .. pitch_angle(dfig,spitch_2) =e= pitch_angle(dfig,spitch_2-1)",
                  "pll_angle_ini(dfig,s) .. pll_angle(dfig,s) =e= alpha_w(dfig,s) - pi/2",
                  "IDr_equation(dfig,sdfig_2) .. IDr(dfig,sdfig_2) =e= IDr_c(dfig,sdfig_2)*cos(pll_angle(dfig,sdfig_2)) - IQr_c(dfig,sdfig_2)*sin(pll_angle(dfig,sdfig_2))",
                  "IQr_equation(dfig,sdfig_2) .. IQr(dfig,sdfig_2) =e= IDr_c(dfig,sdfig_2)*sin(pll_angle(dfig,sdfig_2)) + IQr_c(dfig,sdfig_2)*cos(pll_angle(dfig,sdfig_2))",
                  "Ps_w_ref_equation(dfig,sdfig_2) .. Ps_w_ref(dfig,sdfig_2) =e= P_ord(dfig,sdfig_2)",
                  "iQr_c_ref_equation(dfig,sdfig_2) .. IQr_c_ref(dfig,sdfig_2) =e= Ps_w_ref(dfig,sdfig_2)*xss(dfig)/v_w(dfig,sdfig_2)",
                  "IDr_c_ref_2(dfig,sdfig_2) .. IDr_c_ref(dfig,sdfig_2) =e= v_w(dfig,sdfig_2)",
                  "IDr_c_2(dfig,sdfig_2) .. IDr_c(dfig,sdfig_2) =e= IDr_c_ref(dfig,sdfig_2)",
                  "IQr_c_2(dfig,sdfig_2) .. IQr_c(dfig,sdfig_2) =e= IQr_c_ref(dfig,sdfig_2)",
                  "P_gsc_equation(dfig,s) .. P_gsc(dfig,s) =e= rr_var(dfig,s)/sqr(lm(dfig))*(VDr(dfig,s)*IDr(dfig,s) + VQr(dfig,s)*IQr(dfig,s))",
                  "Ps_w_equation(dfig,s) .. Ps_w(dfig,s) =e= vD_w(dfig,s)*iDs_w(dfig,s) + vQ_w(dfig,s)*iQs_w(dfig,s)",
                  "Qs_w_equation(dfig,s) .. Qs_w(dfig,s) =e= vQ_w(dfig,s)*iDs_w(dfig,s) - vD_w(dfig,s)*iQs_w(dfig,s)",
                  "Pg_w_equation(dfig,s) .. Pg_w(dfig,s) =e= vD_w(dfig,s)*iD_w(dfig,s) + vQ_w(dfig,s)*iQ_w(dfig,s)",
                  "Qg_w_equation(dfig,s) .. Qg_w(dfig,s) =e= vQ_w(dfig,s)*iD_w(dfig,s) - vD_w(dfig,s)*iQ_w(dfig,s)",
                  "P_balance(dfig,s) .. Pg_w(dfig,s) - Ps_w(dfig,s) + P_gsc(dfig,s) =e= 0",
                  "Q_balance(dfig,s) .. Qg_w(dfig,s) - Qs_w(dfig,s) =e= 0"]


    gams_equation_instance = Gams_equation(name,
                                           variables,
                                           restrictions,
                                           equation_names,
                                           equations)

    return gams_equation_instance


def def_dq_gen_currents_eq():

    name = 'Generator d-q current equations' 
    variables = []
    restrictions = []
    equation_names = []
    equations = []

    variables.append("delta(sg,s)")
    restrictions.append("delta.l(sg,s)= 0")
    restrictions.append("delta.lo(sg,s) = -9999")
    restrictions.append("delta.up(sg,s) = 9999")

    variables.append("ed_p(sg,s)")
    restrictions.append("ed_p.l(sg,s) = 0.2")
    restrictions.append("ed_p.lo(sg,s) = 0")
    restrictions.append("ed_p.up(sg,s) = 1.5")

    variables.append("eq_p(sg,s)")
    restrictions.append("eq_p.l(sg,s) = 1.0")
    restrictions.append("eq_p.lo(sg,s) = 0")
    restrictions.append("eq_p.up(sg,s) = 1.5")

    variables.append("id(sg,s)")
    restrictions.append("id.l(sg,s) = 1")
    restrictions.append("id.lo(sg,s) = -Sm(sg)")
    restrictions.append("id.up(sg,s) = 3*Sm(sg)")

    variables.append("iq(sg,s)")
    restrictions.append("iq.l(sg,s) = 1")
    restrictions.append("iq.lo(sg,s) = -Sm(sg)")
    restrictions.append("iq.up(sg,s) = 3*Sm(sg)")

    variables.append("vD_w(dfig,s)")
    variables.append("vQ_w(dfig,s)")

    
    equation_names = ["id_stator_fault",
                      "iq_stator_fault",
                      "id_stator_postfault",
                      "iq_stator_postfault"]

    equations =  ["id_stator_fault(sg,sf) .. id(sg,sf)\n\
                                               - (sum(sgp,Y_f_red(sg,sgp)*(ed_p(sgp,sf)*cos(delta(sg,sf) - delta(sgp,sf) - theta_f_red(sg,sgp))\n\
                                                                         + eq_p(sgp,sf)*sin(delta(sg,sf) - delta(sgp,sf) - theta_f_red(sg,sgp))))\n\
                                               +  sum(dfig,Y_f_red(sg,dfig)*(vD_w(dfig,sf)*sin(delta(sg,sf) - theta_f_red(sg,dfig))\n\
                                                                           - vQ_w(dfig,sf)*cos(delta(sg,sf) - theta_f_red(sg,dfig))))) =e= 0",
                  "iq_stator_fault(sg,sf) .. iq(sg,sf)\n\
                                               - (sum(sgp,Y_f_red(sg,sgp)*(eq_p(sgp,sf)*cos(delta(sg,sf) - delta(sgp,sf) - theta_f_red(sg,sgp))\n\
                                                                         - ed_p(sgp,sf)*sin(delta(sg,sf) - delta(sgp,sf) - theta_f_red(sg,sgp))))\n\
                                               +  sum(dfig,Y_f_red(sg,dfig)*(vD_w(dfig,sf)*cos(delta(sg,sf) - theta_f_red(sg,dfig))\n\
                                                                           + vQ_w(dfig,sf)*sin(delta(sg,sf) - theta_f_red(sg,dfig))))) =e= 0",
                  "id_stator_postfault(sg,spf) .. id(sg,spf)\n\
                                     - (sum(sgp,Y_pf_red(sg,sgp)*(ed_p(sgp,spf)*cos(delta(sg,spf) - delta(sgp,spf) - theta_pf_red(sg,sgp))\n\
                                                                + eq_p(sgp,spf)*sin(delta(sg,spf) - delta(sgp,spf) - theta_pf_red(sg,sgp))))\n\
                                     +  sum(dfig,Y_pf_red(sg,dfig)*(vD_w(dfig,spf)*sin(delta(sg,spf) - theta_pf_red(sg,dfig))\n\
                                                                  - vQ_w(dfig,spf)*cos(delta(sg,spf) - theta_pf_red(sg,dfig))))) =e= 0",
                  "iq_stator_postfault(sg,spf) .. iq(sg,spf)\n\
                                     - (sum(sgp,Y_pf_red(sg,sgp)*(eq_p(sgp,spf)*cos(delta(sg,spf) - delta(sgp,spf) - theta_pf_red(sg,sgp))\n\
                                                                - ed_p(sgp,spf)*sin(delta(sg,spf) - delta(sgp,spf) - theta_pf_red(sg,sgp))))\n\
                                     +  sum(dfig,Y_pf_red(sg,dfig)*(vD_w(dfig,spf)*cos(delta(sg,spf) - theta_pf_red(sg,dfig))\n\
                                                                  + vQ_w(dfig,spf)*sin(delta(sg,spf) - theta_pf_red(sg,dfig))))) =e= 0"]


    gams_equation_instance = Gams_equation(name,
                                           variables,
                                           restrictions,
                                           equation_names,
                                           equations)

    return gams_equation_instance