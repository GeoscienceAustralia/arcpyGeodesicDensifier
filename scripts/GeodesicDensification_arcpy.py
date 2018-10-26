import arcpy

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
        
        out_layer = arcpy.Parameter(
            displayName="output path",
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
                parameters[2].value = in_layer_value + "_densified_" + str(int(in_spacing_value)) + "m"
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
        in_layer = parameters[0].valueAsText
        spacing = parameters[1].valueAsText
        out_layer_name = parameters[2].valueAsText
        out_layer_path = arcpy.env.scratchGDB
        out_layer_abspath = arcpy.os.path.join(out_layer_path, out_layer_name)
        
        def densify_point(in_layer, out_layer, spacing, point_type_field):
            counter = 0
            current_point = arcpy.Point()
            in_field_names = ['SHAPE@' if f.type == 'Geometry' else f.name for f in arcpy.ListFields(in_layer)]
            out_field_names = in_field_names[:]
            out_field_names.append(point_type_field)
            geom_field_index = in_field_names.index('SHAPE@')
            # cursor to write output layer
            cur = arcpy.da.InsertCursor(out_layer, "*")
            # loop through the features and retrieve [X,Y]
            with arcpy.da.SearchCursor(in_layer, in_field_names) as cursor:
                for row in cursor:
                
                    # this is to write the first point
                    if counter == 0:
                        # get the point geometry
                        pt_geom = row[geom_field_index]
                        first_point = pt_geom.getPart()
                        # save point geometry for later
                        current_point.X = first_point.X
                        current_point.Y = first_point.Y
                        # write to output layer
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
                        angle, distance = start_ptgeom.angleAndDistanceTo(next_ptgeom, "GEODESIC")
                        arcpy.AddMessage(" ".join([str(counter), str(angle), str(distance)]))
                        # determine how many densified segments there will be
                        n = int(math.ceil(distance / float(spacing)))
                        # adjust the spacing distance
                        seglen = distance / n
                        for i in range(1, n):
                            if i > 0:
                                newWaypoint = start_ptgeom.pointFromAngleAndDistance(angle, seglen * i, "GEODESIC")
                                point = arcpy.Point(newWaypoint.extent.XMax, newWaypoint.extent.YMax)
                                current_point.X = point.X
                                current_point.Y = point.Y
                                # write to output layer
                                row_list = list(row)
                                row_list.append("Densified")
                                row_list[geom_field_index] = point
                                out_row = tuple(row_list)
                                cur.insertRow(out_row)
                        pt_geom = row[geom_field_index].getPart()
                        current_point.X = pt_geom.X
                        current_point.Y = pt_geom.Y
                        
                        # write to output layer
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
            fieldNameList = [field.name for field in fields]
            pointTypeField = ''
            for fieldName in ["pointType", "pntType", "pntTyp"]:
                if fieldName not in fieldNameList:
                    pointTypeField = fieldName
                    break
            arcpy.AddField_management(arcpy.os.path.join(out_layer_path, out_layer_name), pointTypeField, "TEXT", field_length=50)

            # run the densification
            densify_point(in_layer, out_layer_abspath, spacing, pointTypeField)

        # add output layer to map TOC
        mxd = arcpy.mapping.MapDocument("CURRENT")
        data_frame = arcpy.mapping.ListDataFrames(mxd)[0]
        layer = arcpy.mapping.Layer(os.path.join(out_layer_path, out_layer_name))
        arcpy.mapping.AddLayer(data_frame, layer, "AUTO_ARRANGE")
        arcpy.AddMessage("finished processing")
        arcpy.GetMessages()
