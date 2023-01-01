################################################################################
###                                                                          ###
###   GH_Gcode for HW                                                        ###
###                                                                          ###
###       Component : (2) Points to gcode                                    ###
###                                                                          ###
###                                                                          ###
###   Base Script >>> GH_Gcode                                               ###
###                                                                          ###
###       Repository : https://github.com/naysok/GH_Gcode                    ###
###       Coding : naoki yoshioka (naysok)                                   ###
###       License : MIT License                                              ###
###                                                                          ###
###   Update, 210822 / dd Leveling)                                          ###
###   Update, 210826 / Add E_Retract)                                        ###
###   Update, 221128 / Bug_Fix for E_Retract                                 ###
###   Update, 221203 / Update E_Retract                                      ###
###   Update, 221231 / Add F_Retract/ E_Travel                               ###
###                                                                          ###
################################################################################


import datetime
import math

import rhinoscriptsyntax as rs
import Rhino.Geometry as rg


class Util():

    @staticmethod
    def get_current_time():
        return str(datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S"))

    @staticmethod
    def remap_number(src, old_min, old_max, new_min, new_max):
        return ((src - old_min) / (old_max - old_min) * (new_max - new_min) + new_min)

    @staticmethod
    def flatten_runtime_list(list_):
        all_ = []
        
        for i in xrange(len(list_)):
            # print(i)
            sub = list_[i]
            for j in xrange(len(sub)):
                # print(j)
                all_.append(sub[j])

        return all_

    @staticmethod
    def export_gcode(dir_path, now, txt):

        file_path = dir_path + now + ".gcode"

        ### Export
        with open(file_path, 'w') as f:
            f.write(txt)

        print("********************\n*** Export GCode ***\n********************\n{}".format(file_path))

    @staticmethod
    def zip_matrix(mat):
        ### https://note.nkmk.me/python-list-transpose/
        return [list(x) for x in zip(*mat)]

    @staticmethod
    def padding_previous_value(list_):
        
        list_pad = []

        for i in xrange(len(list_)):
            
            item_ = list_[i]

            ### First 
            if i == 0:
                if (item_ == None):
                    list_pad.append(0)
                else:
                    list_pad.append(item_)
            
            ### Not Frist
            else:
                if (item_ == None):
                    tmp = list_pad[i-1]
                    list_pad.append(tmp)
                else:
                    list_pad.append(item_)
        
        return list_pad

    @staticmethod
    def remove_previous_elements(a_list):

        ### Remove Same Element as the Previous One
        new_list = []
        src_length = len(a_list)

        for i in range(src_length):
            tmp = a_list[i]

            ### 
            if i < src_length-1:
                if a_list[i] != a_list[i+1]:
                    new_list.append(tmp)
            ### Last
            else:
                new_list.append(tmp)
                
        return new_list

    @staticmethod
    def bitwise_or_arrays(arrays):

        ### Merge Bool from Brep.isPointInside

        # print(len(arrays[0]))
        # print(len(arrays[0][0]))

        if len(arrays) == 1:
            return arrays[0]

        else:
            arrays_zip = self.zip_matrix(arrays)
            # print(len(arrays_zip))
            # print(len(arrays_zip[0]))
            # print(len(arrays_zip[0][0]))

            bool_inside = []

            for i in xrange(len(arrays_zip)):
                # print(i)

                sub_ = []

                item = arrays_zip[i]
                item_zip = self.zip_matrix(item)

                # print(len(item_zip))
                # print(len(item_zip[0]))

                for j in xrange(len(item_zip)):

                    values = item_zip[j]
                    
                    if True in values:
                        sub_.append(True)
                    else:
                        sub_.append(False)

                bool_inside.append(sub_)

            # print(len(bool_inside))
            # print(len(bool_inside[0]))

            return bool_inside


class Transform():


    def pt_pt_add(self, pt0, pt1):
        pt = [
            float(pt1[0]) + float(pt0[0]),
            float(pt1[1]) + float(pt0[1]),
            float(pt1[2]) + float(pt0[2])]
        return pt


    def pt_pt_subtract(self, pt0, pt1):
        pt = [
            float(pt0[0]) - float(pt1[0]),
            float(pt0[1]) - float(pt1[1]),
            float(pt0[2]) - float(pt1[2])]
        return pt


    def vector_multiplicate(self, vector, value):
        vec = [
            float(vector[0]) * value,
            float(vector[1]) * value,
            float(vector[2]) * value]
        return vec


    def vector_unitize(self, vector):
        length = math.sqrt(
            math.pow(float(vector[0]), 2) + 
            math.pow(float(vector[1]), 2) + 
            math.pow(float(vector[2]), 2))
        new_vector = self.vector_multiplicate(vector, (1.0 / length))
        return new_vector


    def move_pt_vec(self, pt, vec):
        p = [
            float(pt[0]) + float(vec[0]),
            float(pt[1]) + float(vec[1]),
            float(pt[2]) + float(vec[2])]
        return p
        
tr = Transform()


class Curve():


    def polyline_to_points(self, polyline):

        ### Polyline to Points by RhinoCommon
        pl = rs.coercegeometry(polyline)
        new_pl = rg.PolylineCurve.ToPolyline(pl)
        points = new_pl.ToArray()
        
        return points


    def polylines_to_points(self, polylines):

        points_array = []

        for i in xrange(len(polylines)):

            points = self.polyline_to_points(polylines[i])
            points_array.append(points)

        return points_array


    def offset_z_curves(self, curves, z_offset_value):

        move_vec = rs.CreateVector(0, 0, z_offset_value)
        geos_off = rs.MoveObjects(curves, move_vec)
        
        return geos_off

op_c = Curve()


################################################################


class CalcVector():

    @staticmethod
    def calc_distance_2pt(x0, y0, z0, x1, y1, z1):

        d = math.sqrt(
            ((x0 - x1) * (x0 - x1)) + 
            ((y0 - y1) * (y0 - y1)) + 
            ((z0 - z1) * (z0 - z1)))

        return d


################################################################


class ParametersExtrude():

    def __init__(self, extrude_amp=None, extrude_retract=None, extrude_retract_back=None):
        
        self.extrude_amp = extrude_amp
        self.extrude_retract = extrude_retract
        self.extrude_retract_back = extrude_retract_back

    def define_print_parameter(self):
        
        ### For Header
        gcode_extrude_amp =          "; EXTRUDE_AMP, ratio       : {}\n".format(str(self.extrude_amp))
        gcode_extrude_retract =      "; EXTRUDE_RETRACT, mm      : {}\n".format(str(self.extrude_retract))
        gcode_extrude_retract_back = "; EXTRUDE_RETRACT_BACK, mm : {}\n".format(str(self.extrude_retract_back))

        gcode_extrude_parameter = gcode_extrude_amp + gcode_extrude_retract + gcode_extrude_retract_back

        return gcode_extrude_parameter


class ParametersFeed():

    def __init__(self, feed_print=None, feed_retract=None, feed_travel=None):
        
        self.feed_print = feed_print
        self.feed_retract = feed_retract
        self.feed_travel = feed_travel

    def define_print_parameter(self):
            
        ### For Header
        gcode_feed_print =   "; FEED_PRINT   : {}\n".format(str(self.feed_print))
        gcode_feed_retract = "; FEED_RETRACT : {}\n".format(str(self.feed_retract))
        gcode_feed_travel =  "; FEED_TRAVEL  : {}\n".format(str(self.feed_travel))

        gcode_feed_parameter = gcode_feed_print + gcode_feed_retract + gcode_feed_travel

        return gcode_feed_parameter


class ParametersTemperature():

    def __init__(self, temp_nozzle=None, temp_bed=None):
        
        self.temp_nozzle = temp_nozzle
        self.temp_bed = temp_bed

    def define_print_parameter(self):
            
        ### For Header
        gcode_temp_nozzle = "; TEMP_NOZZLE : {}\n".format(str(self.temp_nozzle))
        gcode_temp_bed =    "; TEMP_BED    : {}\n".format(str(self.temp_bed))

        gcode_temp_parameter = gcode_temp_nozzle + gcode_temp_bed

        return gcode_temp_parameter


################################################################


class MarlinGcodeHeader():
    
    def __init__(self, now, component_info, param_e, param_feed, param_temp, fan, z_buffer):
        
        self.now = now
        self.component_info = component_info
        self.param_e = param_e
        self.param_feed = param_feed
        self.param_temp = param_temp
        self.fan = fan
        self.z_buffer = z_buffer

    def define_header_top(self):

        t0 = "; GH_Gcode for HW\n"
        t1 = "; Polyline to Gcode by Grasshopper\n"
        t2 = "; For Mariln 3D Printer\n"

        return t0 + t1 + t2

    def define_print_parameter(self):

        prm_e = self.param_e
        prm_f = self.param_feed
        prm_t = self.param_temp

        ____line____ = "; ----- Print Parameter -----\n"

        ### Time Stamp
        time = "; Export Time : {}\n".format(self.now)
        
        ### Print Parameter
        print_comp = "; Component Info : {}\n".format(self.component_info)

        ### Extrude
        print_extrude = prm_e.define_print_parameter()

        ### Feed
        print_feed = prm_f.define_print_parameter()

        ### Temp
        print_temp = prm_t.define_print_parameter()

        ### Other
        print_fan = "; FAN : {}\n".format(str(self.fan))
        print_z_buffer = "; Z_BUFFER : {}\n".format(str(self.z_buffer))
        
        print_prm = [
            ____line____,
            time,
            print_comp,
            print_extrude,
            print_feed,
            print_temp,
            print_fan,
            print_z_buffer,
            ____line____]

        print_prm_join = "".join(print_prm)

        return print_prm_join

    def define_header(self):

        top = self.define_header_top()
        prms = self.define_print_parameter()

        return top + prms


################################################################


class MarlinGcodeMachineSetup():


    def __init__(self, param_temp, fan):

        self.param_temp = param_temp
        self.fan = fan


    def define_general_settings(self):

        ################################

        ### G90 - Absolute Positioning
        set_g90 = "G90 ; Absolute Positioning\n"

        ### M82 - E Absolute Mode
        set_m82 = "M82 ; E Absolute Mode\n"

        ### Set Tool
        set_t0 = "T0\n"

        ################################

        settings = set_g90 + set_m82 + set_t0

        return settings


    def start_fan(self):
        
        ### M106 - Set Fan Speed
        set_fan = "M106 S{} ; Set - Fan\n".format(str(self.fan))

        return set_fan


    def start_bed(self):
        
        prm_t = self.param_temp

        ### M140 - Set Bed Temperature
        ### M190 - Wait for Bed Temperature

        set_temp_bed = "M140 S{} ; Set - Temp Bed\n".format(str(prm_t.temp_bed))

        if float(prm_t.temp_bed) < 40:
            return set_temp_bed

        else:
            wait_temp_bed = "M190 S{} ; Wait - Temp Bed\n".format(str(prm_t.temp_bed))
            return set_temp_bed + wait_temp_bed


    def start_extruder(self):

        prm_t = self.param_temp

        ### M104 - Set Hotend Temperature
        ### M109 - Wait for Hotend Temperature
        
        set_temp_extruder = "M104 S{} T0 ; Set - Temp Extruder\n".format(str(prm_t.temp_nozzle))
        wait_temp_extruder = "M109 S{} T0 ; Wait - Temp Extruder\n".format(str(prm_t.temp_nozzle))

        return set_temp_extruder + wait_temp_extruder


    def homing_all_axes(self):
        
        ### G28 - Auto Home
        homing = "G28 ; Homing\n"

        return homing


    def homing_x(self):
        
        ### G28 - Auto Home
        homing = "G28 X0 ; Homing X\n"

        return homing

    
    def auto_leveling(self):

        ### G29 - Auto Leveling
        leveling_auto = "G29 ; Auto Leveling\n"

        return leveling_auto


    def reset_extrude_value(self):

        ### G92 - Set Position
        reset_e = "G92 E0 ; Reset Extruder Value, for Start Code\n"

        return reset_e


    def machine_start(self):

        start = []

        ################################

        start.append("; ----- Start Code -----\n")

        ### General Setting
        start.append(self.define_general_settings())
        
        ### Start Fan
        start.append(self.start_fan())

        ### Start Bed
        start.append(self.start_bed())

        ### Start Extruder
        start.append(self.start_extruder())

        ### Homing
        start.append(self.homing_all_axes())

        ### Leveling
            ### Ref
            ### hineri_kazamatsuri.gcode
        # start.append(self.auto_leveling())

        ### Reset Extruder Value
        start.append(self.reset_extrude_value())

        start.append("; ----- Start Code -----\n")

        ################################

        start_join = "".join(start)

        return start_join


    def machine_end(self):

        end = []

        ################################

        end.append("; ----- End Code -----\n")

        ### Homing X
        end.append(self.homing_x())

        ### End Part
        end.append("M106 S0 ; turn off cooling fan\n")
        end.append("M104 S0 ; turn off extruder\n")
        end.append("M140 S0 ; turn off bed\n")
        end.append("M84 ; disable motors\n")

        end.append("; ----- End Code -----\n")

        ################################

        end_join = "".join(end)

        return end_join


################################################################


class MarlinGcodePrinting():


    def travel(self, z_current, z_zuffer):

        ### Travel

        gcode = []

        comment = "; --- Travel ---\n"

        gz = str("{:.4f}".format(z_current + float(z_zuffer)))
        code_travel = "G1 Z{}\n".format(gz)

        ###

        gcode.append(comment)
        gcode.append(code_travel)
        gcode.append(comment)

        ###

        gcode_join = "".join(gcode)

        return gcode_join


    def retract(self, e_retract):

        ### E_Retract

        gcode = []

        comment = "; --- E_Retract ---\n"
        
        ### Reset Extruder Value
        code_reset = self.reset_extrude_value()
        
        v = str("{:.4f}".format(float(e_retract)))
        code_retract = "G0 E-{}\n".format(v)
        
        ###
        
        gcode.append(comment)
        gcode.append(code_reset)
        gcode.append(code_retract)
        gcode.append(code_reset)
        gcode.append(comment)

        ###

        gcode_join = "".join(gcode)

        return gcode_join


    def retract_back(self, e_retract_back):

        ### E_Retract_Back

        gcode = []

        comment = "; --- E_Retract_Back ---\n"

        ### Reset Extruder Value
        code_reset = self.reset_extrude_value()
        
        v = str("{:.4f}".format(float(e_retract_back)))
        code_retract_back = "G0 E{}\n".format(v)
 
        ###

        gcode.append(comment)
        gcode.append(code_reset)
        gcode.append(code_retract_back)
        gcode.append(code_reset)
        gcode.append(comment)

        ###

        gcode_join = "".join(gcode)

        return gcode_join


    def point_to_gcode(self, count, pts, e_amp, e_retract, e_retract_back, feed, z_zuffer):
        
        layer = []

        ### CR-10 / TPU
        PRINTER_PARAMETER = 0.165

        for i in xrange(len(pts)):
            p = pts[i]
            xx, yy, zz = p

            ### Index[0]
            if i == 0:

                ### Initialize
                ee = 0

                ### Layer Info (Comment)
                start_comment = "; ----- Layer : {} / start -----\n".format(count)
                layer.append(start_comment)

                ### Reset Extruder Value
                layer.append(self.reset_extrude_value())

                gx = str("{:.4f}".format(xx))
                gy = str("{:.4f}".format(yy))
                gz = str("{:.4f}".format(zz))
                gf = str("{}".format(feed))

                gcode = "G1 X{} Y{} Z{} E0 F{}\n".format(gx, gy, gz, gf)
                layer.append(gcode)

                if int(count) != 0:
                    ### E_Retract_Back
                    retract_back = self.retract_back(e_retract_back)
                    layer.append(retract_back)


            ### Index[1] - Index[Last]
            else:
                x0, y0, z0 = pts[i - 1]
                x1, y1, z1 = pts[i]
                
                distance = self.calc_distance_2pt(x0, y0, z0, x1, y1, z1)

                ### PRINTER_PARAMETER // 1.85 to 0.4
                ### e_amp // override
                ee += (distance * float(PRINTER_PARAMETER) * float(e_amp))

                gx = str("{:.4f}".format(xx))
                gy = str("{:.4f}".format(yy))
                gz = str("{:.4f}".format(zz))
                ge = str("{:.4f}".format(ee))
                
                gcode = "G1 X{} Y{} Z{} E{}\n".format(gx, gy, gz, ge)
                layer.append(gcode)
            
                ### Index[Last]
                if i == (len(pts) - 1):

                    ### E_Retract
                    retract = self.retract(e_retract)
                    layer.append(retract)

                    ### Travel
                    travel = self.travel(zz, z_zuffer)
                    layer.append(travel)

                    ### Layer Info (Comment)
                    end_comment = "; ----- Layer : {} / end -----\n".format(count)
                    layer.append(end_comment)

        ### Join
        layer_join = "".join(layer)

        return layer_join


    ##################################################################


    ##############################
    ###                        ###
    ###     Generate Gcode     ###
    ###                        ###
    ##############################


    def points_list_to_gcode(self, points_list, now, component_info, params_e, params_feed, params_temp, fan, z_zuffer):

        export = []

        ################################

        header = MarlinGcodeHeader(now, component_info, params_e, params_feed, params_temp, fan, z_zuffer)
        machine_setup = MarlinGcodeMachineSetup(params_temp, fan)


        ### (1) Print Header, Parameters
        export.append(header.define_header())

        ### (2) Machine Start
        export.append(machine_setup.machine_start())


        """
        ### (3) Printing
        for i in xrange(len(points_list)):

            pts = points_list[i]
            layer_count = str(i)
            export.append(self.point_to_gcode(layer_count, pts, e_amp, e_retract, e_retract_back, feed, z_zuffer))

        """

        ### (4) Machine End
        export.append(machine_setup.machine_end())

        ################################


        ### JOIN
        export_join = "".join(export)

        return export_join


op_ml = MarlinGcodePrinting()



##########

##########


BUILD_DAY = "221231"


ghenv.Component.Message = "2) Points to gcode / {}".format(BUILD_DAY)
COMPONENT_INFO = "ver_{}".format(BUILD_DAY)


NOW = Util.get_current_time()

Z_OFFSET_VALUE = INFO
FAN = 0

PARAMS_EXTRUDE = ParametersExtrude(E_AMP, E_RETRACT, E_RETRACT_BACK)
PARAMS_FEED = ParametersFeed(FEED_PRINT, FEED_RETRACT, FEED_TRAVEL)
PARAMS_TEMP = ParametersTemperature(TEMP_NOZZLE, TEMP_BED)


### Points to Gcode (Not Go Through Machine Origin)
if RUN_AND_EXPORT:
    gcode = op_ml.points_list_to_gcode(POINTS, NOW, COMPONENT_INFO, PARAMS_EXTRUDE, PARAMS_FEED, PARAMS_TEMP, FAN, Z_BUFFER)
    Util.export_gcode(EXPORT_DIR, NOW, gcode)

    ### For DEBUG
    print(NOW)