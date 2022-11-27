################################################################################
###                                                                          ###
###   GH_Gcode for TL                                                        ###
###                                                                          ###
###       Component : (1) Polyline to Points / 201123                        ###
###                                                                          ###
###                                                                          ###
###   Base Script >>> GH_Gcode                                               ###
###                                                                          ###
###       Repository : https://github.com/naysok/GH_Gcode                    ###
###       Coding : naoki yoshioka (naysok)                                   ###
###       License : MIT License                                              ###
###                                                                          ###
################################################################################



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


##########

##########


ghenv.Component.Message = '(1) Polyline to Points / 201123'


### Z OFFSET
offset_curves = op_c.offset_z_curves(POLYLINES, Z_OFFSET_VALUE)
#DEV = offset_curves
INFO = Z_OFFSET_VALUE


### Polylines to Points
pts = op_c.polylines_to_points(offset_curves)
POINTS = pts


def print_preview(points_array):

    pls = []
    lens_ = 0.0

    for i in xrange(len(points_array)):
        pl = rg.Polyline(points_array[i])
        pls.append(pl)
        
        len_ = pl.Length
        lens_ += float(len_)
        
    lens_int = int(lens_)

    return pls, lens_int

def draw_machine(w, h):
    rc = rs.AddRectangle(rs.WorldXYPlane(), w, h)
    return rc


### Print Preview

print_crvs, print_len = print_preview(pts)
__PRINT_CURVES = print_crvs
__PRINT_LENGTH = print_len


### Draw Machine
MACHINE = draw_machine(400, 400)