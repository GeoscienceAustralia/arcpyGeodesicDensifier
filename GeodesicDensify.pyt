import os
import sys

import arcpy

scripts_dir = os.path.join(os.path.dirname(__file__), 'scripts')
sys.path.append(scripts_dir)
# Do not compile .pyc files for the tool modules.
sys.dont_write_bytecode = True

class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Toolbox"
        self.alias = "Geodesic Densification using arcpy"

        # List of tool classes associated with this toolbox
        self.tools = [GeodesicDensification_arcpy]

                      
class GeodesicDensification_arcpy(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Geodesic Densification using arcpy"
        self.description = "Densifies geometries along geodesic lines"
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = []
        in_layer = arcpy.Parameter(
            displayName="Input Layer",
            name="in_layer",
            datatype="GPLayer",
            parameterType="Required",
            direction="Input")
        params.append(in_layer)
                                    
        in_spacing = arcpy.Parameter(
            displayName="Maximum Point Spacing (m)",
            name="in_spacing",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")
        in_spacing.value = 900
        params.append(in_spacing)
        
        dens_type_field = arcpy.Parameter(
            displayName="Field with densification type",
            name="dens_type_field",
            datatype="Field",
            parameterType="Required",
            direction="Input")
        dens_type_field.parameterDependencies = [in_layer.name]
        params.append(dens_type_field)
        
        dens_type_geodesic = arcpy.Parameter(
            displayName="Field value for Geodesic densification",
            name="dens_type_geodesic",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
        dens_type_geodesic.filter.type = "ValueList"
        dens_type_geodesic.filter.list = []
        params.append(dens_type_geodesic)
        
        dens_type_loxodrome = arcpy.Parameter(
            displayName="Field value for Loxodrome densification",
            name="dens_type_loxodrome",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
        dens_type_loxodrome.filter.type = "ValueList"
        dens_type_loxodrome.filter.list = []
        params.append(dens_type_loxodrome)
        
        out_layer = arcpy.Parameter(
            displayName="output layer",
            name="out_path",
            datatype="GPLayer",
            parameterType="Required",
            direction="Output")
        params.append(out_layer)
        
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        if parameters[0].altered or parameters[1].altered:
            in_layer_value = parameters[0].valueAsText
            in_spacing_value = parameters[1].valueAsText
            if in_layer_value is not None and in_spacing_value is not None:
                parameters[5].value = in_layer_value + "_densified_" + str(int(in_spacing_value)) + "m"
            
        if parameters[2].altered:
            with arcpy.da.SearchCursor(parameters[0].valueAsText, parameters[2].valueAsText) as g_rows:
                parameters[3].filter.list = sorted(list(set([row[0] for row in g_rows])))
            with arcpy.da.SearchCursor(parameters[0].valueAsText, parameters[2].valueAsText) as l_rows:
                parameters[4].filter.list = sorted(list(set([row[0] for row in l_rows])))
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        import math
        import os

        # set default values
        out_layer_path = arcpy.env.scratchWorkspace
        if not out_layer_path:
            out_layer_path = arcpy.env.scratchGDB
        in_layer = parameters[0].valueAsText
        spacing = parameters[1].valueAsText
        out_layer_name = parameters[5].valueAsText
        out_layer_abspath = arcpy.os.path.join(out_layer_path, out_layer_name)
        dens_field = parameters[2].valueAsText
        geodesic_field_value = parameters[3].valueAsText
        loxodrome_field_value = parameters[4].valueAsText
        
        def densify_points(in_layer, out_layer, spacing, dens_field):
            counter = 0
            current_point = arcpy.Point()
            in_field_names = ['SHAPE@' if f.type == 'Geometry' else f.name for f in arcpy.ListFields(in_layer)]
            geom_field_index = in_field_names.index('SHAPE@')
            dens_field_index = in_field_names.index(dens_field)
            current_dens_type = None
            # cursor to write output layer
            cur = arcpy.da.InsertCursor(out_layer, "*")
            # loop through the features densify
            with arcpy.da.SearchCursor(in_layer, in_field_names) as cursor:
                for row in cursor:

                    # get the Densification type
                    if row[dens_field_index] == geodesic_field_value:
                        dens_type = "GEODESIC"
                    elif row[dens_field_index] == loxodrome_field_value:
                        dens_type = "LOXODROME"
                    else:
                        dens_type = None

                    # this is to write the first point
                    if counter == 0:
                        # get the point geometry
                        pt_geom = row[geom_field_index]
                        first_point = pt_geom.getPart()
                        # save point geometry for the next point
                        current_point.X = first_point.X
                        current_point.Y = first_point.Y
                        # save the densification type for the next point
                        current_dens_type = dens_type
                        # save the row for later
                        row_list = list(row)
                        row_list.append('Original')
                        row = tuple(row_list)
                        cur.insertRow(row)

                    # This is for subsequent points
                    elif counter > 0:
                        start_pt = current_point
                        start_ptgeom = arcpy.PointGeometry(start_pt, in_srs)
                        next_pt_geom = row[geom_field_index]
                        next_pt = next_pt_geom.getPart()
                        next_ptgeom = arcpy.PointGeometry(next_pt, in_srs)
                        if current_dens_type is not None:  # densify if loxodrome or geodesic
                            angle, distance = start_ptgeom.angleAndDistanceTo(next_ptgeom, current_dens_type)
                            # determine how many densified segments there will be
                            segment_count = int(math.ceil(distance / float(spacing)))
                            # adjust the spacing distance
                            seglen = distance / segment_count
                            # find every waypoint along segment
                            for i in range(1, segment_count):
                                waypoint = start_ptgeom.pointFromAngleAndDistance(angle, seglen * i, current_dens_type)
                                point = arcpy.Point(waypoint.extent.XMax, waypoint.extent.YMax)
                                current_point.X = point.X
                                current_point.Y = point.Y
                                # write to output layer
                                row_list = list(row)
                                row_list.append("Densified")
                                row_list[geom_field_index] = (point.X, point.Y)
                                out_row = tuple(row_list)
                                cur.insertRow(out_row)
                            # save point geometry for the next point
                            current_point.X = next_pt.X
                            current_point.Y = next_pt.Y
                            # save the densification type for the next point
                            current_dens_type = dens_type
                        else:  # don't densify if not loxodrome or geodesic, just copy the point
                            # write to output layer
                            row_list = list(row)
                            row_list.append("not densified")
                            row_list[geom_field_index] = next_pt_geom
                            out_row = tuple(row_list)
                            cur.insertRow(out_row)
                            # save point geometry for the next point
                            current_point.X = next_pt.X
                            current_point.Y = next_pt.Y
                            # save the densification type for the next point
                            current_dens_type = dens_type
                        # write the last point 
                        row_list = list(row)
                        row_list.append('Original')
                        row = tuple(row_list)
                        cur.insertRow(row)
                    counter += 1
            if cur:
                del cur

        # get input geometry type
        desc = arcpy.Describe(in_layer)
        in_srs = desc.spatialReference
        arcpy.env.overwriteOutput = True
        if desc.featureType == "Simple" and desc.datasetType == "FeatureClass":
            # create output layer
            arcpy.CreateFeatureclass_management(out_layer_path, out_layer_name, "POINT", in_layer, spatial_reference=in_srs)

            # add field for pointType (original or densified)
            fields = arcpy.ListFields(in_layer)
            field_name_list = [field.name for field in fields]
            point_type_field = ''
            for fieldName in ["pointType", "pntType", "pntTyp"]:
                if fieldName not in field_name_list:
                    point_type_field = fieldName
                    break
            arcpy.AddField_management(out_layer_abspath, point_type_field, "TEXT", field_length=50)

            # run the densification
            densify_points(in_layer, out_layer_abspath, spacing, dens_field)

        # add output layer to map TOC
        mxd = arcpy.mapping.MapDocument("CURRENT")
        data_frame = arcpy.mapping.ListDataFrames(mxd)[0]
        layer = arcpy.mapping.Layer(os.path.join(out_layer_path, out_layer_name))
        arcpy.mapping.AddLayer(data_frame, layer, "AUTO_ARRANGE")
        arcpy.GetMessages()
