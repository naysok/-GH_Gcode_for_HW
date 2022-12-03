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
###                                                                          ###
################################################################################


import datetime
import math

import rhinoscriptsyntax as rs
import Rhino.Geometry as rg


class Util():


    def get_current_time(self):
        return str(datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S"))


    def remap_number(self, src, old_min, old_max, new_min, new_max):
        return ((src - old_min) / (old_max - old_min) * (new_max - new_min) + new_min)


    def flatten_runtime_list(self, list):
        all = []
        
        for i in xrange(len(list)):
            # print(i)
            sub = list[i]
            for j in xrange(len(sub)):
                # print(j)
                all.append(sub[j])

        return all


    def export_gcode(self, dir_path, txt):

        now = self.get_current_time()
        file_path = dir_path + now + ".gcode"

        ### Export
        with open(file_path, 'w') as f:
            f.write(txt)

        print("Export GCode : {}".format(file_path))


    def zip_matrix(self, mat):
        ### https://note.nkmk.me/python-list-transpose/
        return [list(x) for x in zip(*mat)]


    def padding_previous_value(self, list_):
        
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


    def remove_previous_elements(self, a_list):

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


    def bitwise_or_arrays(self, arrays):

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

ut = Util()


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

        ########## NEW CODE ##########

        ### Polyline to Points by RhinoCommon
        pl = rs.coercegeometry(polyline)
        new_pl = rg.PolylineCurve.ToPolyline(pl)
        points = new_pl.ToArray()

        ########## NEW CODE ##########

        """
        ########## OLD CODE 1 ##########
        ### Start-Point + CurveDiscontinuity + End-Point
        points = []
        ### Start Point
        start_pt = rs.CurveStartPoint(polyline)
        points.append(start_pt)
        ### Segment Points
        ### Style : 3 = C2 - Continuous first and second derivative
        segments = rs.CurveDiscontinuity(polyline, 3)
        for j in range(len(segments)):
            points.append(segments[j])
        ### End Point
        end_pt =  rs.CurveEndPoint(polyline)
        points.append(end_pt)
        ########## OLD CODE 1 ##########
        """

        """
        ########## OLD CODE 2 ##########
        
        ### Polyline to Points by rhinoscriptsuntax
        ### https://developer.rhino3d.com/api/RhinoScriptSyntax/#collapse-PolylineVertices
        
        points = rs.PolylineVertices(polyline)
        
        ########## OLD CODE 2 ##########
        """
        
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


class MarlinGcode():


    def define_msg(self):
        t0 = "; GH_Gcode for HW\n"
        t1 = "; Polyline to Gcode by Grasshopper\n"
        t2 = "; For Mariln 3D Printer\n"
        t = t0 + t1 + t2
        return t

    def define_print_parameter(self, component, e_amp, feed, temp_nozzle, temp_bed, fan, z_zuffer):
        
        line = "; ----- Print Parameter -----\n"

        ### Time Stamp
        now = ut.get_current_time()
        time = "; Export Time : {}\n".format(now)
        
        ### Print Parameter
        print_comp = "; Component Info : {}\n".format(component)
        print_e_amp = "; E_AMP : {}\n".format(str(e_amp))
        print_feed = "; FEED : {}\n".format(str(feed))
        print_temp_nozzle = "; TEMP_NOZZLE : {}\n".format(str(temp_nozzle))
        print_temp_bed = "; TEMP_BED : {}\n".format(str(temp_bed))
        print_fan = "; FAN : {}\n".format(str(fan))
        print_z_buffer = "; Z_BUFFER : {}\n".format(str(z_zuffer))
        
        print_prm = [line, time, print_comp, print_e_amp, print_feed, 
            print_temp_nozzle, print_temp_bed, print_fan, print_z_buffer, line]

        print_prm_join = "".join(print_prm)

        return print_prm_join


    def define_general_settings(self):

        ### G90 - Absolute Positioning
        set_g90 = "G90 ; Absolute Positioning\n"

        ### M82 - E Absolute Mode
        set_m82 = "M82 ; E Absolute Mode\n"

        ### Set Tool
        set_t0 = "T0\n"

        settings = [set_g90, set_m82, set_t0]
        
        settings_join = "".join(settings)

        return settings_join


    def start_fan(self, fan):
        
        ### M106 - Set Fan Speed
        set_fan = "M106 S{} ; Set - Fan\n".format(str(fan))
        return set_fan


    def start_bed(self, temp_bed):
        
        ### M140 - Set Bed Temperature
        ### M190 - Wait for Bed Temperature

        set_temp_bed = "M140 S{} ; Set - Temp Bed\n".format(str(temp_bed))

        if float(temp_bed) < 40:
            return set_temp_bed

        else:
            wait_temp_bed = "M190 S{} ; Wait - Temp Bed\n".format(str(temp_bed))
            return set_temp_bed + wait_temp_bed


    def start_extruder(self, temp_nozzle):

        ### M104 - Set Hotend Temperature
        ### M109 - Wait for Hotend Temperature
        
        set_temp_extruder = "M104 S{} T0 ; Set - Temp Extruder\n".format(str(temp_nozzle))
        wait_temp_extruder = "M109 S{} T0 ; Wait - Temp Extruder\n".format(str(temp_nozzle))

        return set_temp_extruder + wait_temp_extruder


    def homing_all_axes(self):
        
        ### G28 - Auto Home
        homing = "G28 ; Homing\n"
        return homing


    def homing_x(self):
        
        ### G28 - Auto Home
        homing = "G28 X0 ; Homing X\n"
        return homing

    
    def leveling(self):

        ### G29 - Auto Leveling
        leveling = "G29 ; Auto Leveling\n"
        return leveling


    def reset_extrude_value(self):

        ### G92 - Set Position
        reset_e = "G92 E0 ; Reset Extruder Value\n"
        return reset_e


    def calc_distance_2pt(self, x0, y0, z0, x1, y1, z1):

        ### Calc Distance
        xx = x0 - x1
        yy = y0 - y1
        zz = z0 - z1
        dist = math.sqrt((xx * xx) + (yy * yy) + (zz * zz))

        return dist


    ##################################################################


    #################################
    ###                           ###
    ###     Print Start / End     ###
    ###                           ###
    #################################


    def print_start_old(self, fan, temp_bed, temp_nozzle):

        start = []

        start.append("; ----- Start Code -----\n")

        ### General Setting
        start.append(self.define_general_settings())
        
        ### Start Fan
        start.append(self.start_fan(fan))

        ### Start Bed
        start.append(self.start_bed(temp_bed))

        ### Start Extruder
        start.append(self.start_extruder(temp_nozzle))

        ### Homing
        start.append(self.homing_all_axes())

        ### Reset Extruder Value
        start.append(self.reset_extrude_value())

        start.append("; ----- Start Code -----\n")

        start_join = "".join(start)

        return start_join


    def print_start(self, fan, temp_bed, temp_nozzle):
        
        ### Ref
        ### hineri_kazamatsuri.gcode

        start = []

        start.append("; ----- Start Code -----\n")

        ### General Setting
        start.append(self.define_general_settings())
        
        ### Start Fan
        start.append(self.start_fan(fan))

        ### Start Bed
        start.append(self.start_bed(temp_bed))

        ### Start Extruder
        start.append(self.start_extruder(temp_nozzle))

        ### Homing
        start.append(self.homing_all_axes())

        ### Leveling
        start.append(self.leveling())

        ### Reset Extruder Value
        start.append(self.reset_extrude_value())

        start.append("; ----- Start Code -----\n")

        ###

        start_join = "".join(start)

        return start_join


    def print_end(self):

        end = []

        end.append("; ----- End Code -----\n")

        ### Homing X
        end.append(self.homing_x())

        ### End Part
        end.append("M106 S0 ; turn off cooling fan\n")
        end.append("M104 S0 ; turn off extruder\n")
        end.append("M140 S0 ; turn off bed\n")
        end.append("M84 ; disable motors\n")

        end.append("; ----- End Code -----\n")

        end_join = "".join(end)

        return end_join
    

    ##################################################################


    ########################
    ###                  ###
    ###     Printing     ###
    ###                  ###
    ########################


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


    def points_list_to_gcode(self, points_list, component, e_amp, e_retract, e_retract_back, feed, temp_nozzle, temp_bed, fan, z_zuffer):
        
        ### RUN ALL
        export = []

        ### Msg, Parameters
        export.append(self.define_msg())
        export.append(self.define_print_parameter(component, e_amp, feed, temp_nozzle, temp_bed, fan, z_zuffer))

        ### Print Start
        # export.append(self.print_start(fan, temp_bed, temp_nozzle))
        export.append(self.print_start(fan, temp_bed, temp_nozzle))

        
        ### Printing
        for i in xrange(len(points_list)):

            pts = points_list[i]
            layer_count = str(i)
            export.append(self.point_to_gcode(layer_count, pts, e_amp, e_retract, e_retract_back, feed, z_zuffer))

        ### Print End
        export.append(self.print_end())


        ### JOIN
        export_join = "".join(export)

        return export_join


op_ml = MarlinGcode()



##########

##########


BUILD_DAY = "221203"


ghenv.Component.Message = "2) Points to gcode / {}".format(BUILD_DAY)
comp_info = "ver_{}".format(BUILD_DAY)


Z_OFFSET_VALUE = INFO
FAN = 0


### Points to Gcode (Not Go Through Machine Origin)
if RUN_AND_EXPORT:
    gcode = op_ml.points_list_to_gcode(POINTS, comp_info, E_AMP, E_RETRACT, E_RETRACT_BACK, FEED, TEMP_NOZZLE, TEMP_BED, FAN, Z_BUFFER)
    ut.export_gcode(EXPORT_DIR, gcode)
