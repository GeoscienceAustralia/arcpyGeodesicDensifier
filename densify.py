import math
import os
import sys
print("importing arcpy")
import arcpy
try:
    from geographiclib.geodesic import Geodesic
except ImportError:
    print("attempting to install geographiclib python package")
    import pip
    pip.main(['install', 'geographiclib'])
    try:
        from geographiclib.geodesic import Geodesic
    except ImportError:
        print("Couldn't find geographiclib python package")

# set default values
in_layer = arcpy.GetParameterAsText(0)
out_path = r'in_memory'
#inLayer = r"C:\Users\Jonah\PycharmProjects\arcpy_densifier\testData\testMultiPoly.shp"
#outLayer = r"C:\Users\Jonah\PycharmProjects\arcpy_densifier\testData\outtestMultiPoly.shp"
#outPath = os.path.dirname(outLayer)
#outLayerName = os.path.basename(outLayer)
ellipsoid_a = 6378137.0
ellipsoid_f = 298.2572236
ellipsoid_name = 'WGS84'

# delete output if already exists
if arcpy.Exists(outLayer):
    try:
        arcpy.Delete_management(outLayer)
    except Exception:
        print("couldn't delete old output layer")
        e = sys.exc_info()[1]
        print(e.args[0])


print("inLayer: " + in_layer)

# default point spacing is 900
spacing = 900

ellipsoid_dict = {'165': [6378165.000, 298.3],
                  'ANS': [6378160, 298.25],
                  'CLARKE 1858': [6378293.645, 294.26],
                  'GRS80': [6378137, 298.2572221],
                  'WGS84': [6378137, 298.2572236],
                  'WGS72': [6378135, 298.26],
                  'International 1924': [6378388, 297]}

# get spatial reference
inSr = arcpy.Describe(in_layer).spatialReference
wgs84Sr = arcpy.SpatialReference(4326)

# Create a geographiclib Geodesic object
geod = Geodesic(ellipsoid_a, 1 / ellipsoid_f)


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
            # this is the first point
            if counter == 0:
                # get the point geometry
                pt_geom = row[geom_field_index]
                point = pt_geom.getPart()
                # save point geometry for later
                current_point.X = point.X
                current_point.Y = point.Y
                # write to output layer
                row_list = list(row)
                row_list.append('Original')
                row = tuple(row_list)
                cur.insertRow(row)

            # This is for subsequent points
            elif counter > 0:
                start_pt = current_point
                end_pt_geom = row[geom_field_index]
                end_pt = end_pt_geom.getPart()
                # project to WGS84 if needed
                if inSr.factoryCode != wgs84Sr.factoryCode:
                    start_pt_geom = arcpy.PointGeometry(start_pt, inSr)
                    start_pt = start_pt_geom.projectAs(wgs84Sr).getPart()
                    end_pt_geom = end_pt_geom.projectAs(wgs84Sr)
                    end_pt = end_pt_geom.getPart()
                # create a geographiclib line object
                line_object = geod.InverseLine(start_pt.Y, start_pt.X, end_pt.Y, end_pt.X)
                # determine how many densified segments there will be
                n = int(math.ceil(line_object.s13 / spacing))
                # adjust the spacing distance
                seglen = line_object.s13 / n
                for i in range(1, n):
                    if i > 0:
                        s = seglen * i
                        g = line_object.Position(s,
                                                 Geodesic.LATITUDE |
                                                 Geodesic.LONGITUDE |
                                                 Geodesic.LONG_UNROLL)
                        point = arcpy.Point(g['lon2'], g['lat2'])
                        current_point.X = point.X
                        current_point.Y = point.Y
                        # Convert each point back to the input CRS
                        if inSr.factoryCode != wgs84Sr.factoryCode:
                            pt_geom = arcpy.PointGeometry(point, wgs84Sr)
                            point = pt_geom.projectAs(inSr).getPart()
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


def densify_poly(in_layer, out_layer, spacing):
    # loop through the features
    field_names = ['SHAPE@' if f.type == 'Geometry' else f.name for f in arcpy.ListFields(in_layer)]
    geom_field_index = field_names.index('SHAPE@')
    with arcpy.da.SearchCursor(in_layer, field_names) as cursor:
        for row in cursor:
            part_count = row[geom_field_index].partCount
            for i in range(0, part_count):
                # create array to hold densified points
                point_array = arcpy.Array()
                # get the geometry
                part = row[geom_field_index].getPart(i)
                point_count = len(part)
                start_pt = part[0]
                # project to WGS84 if needed
                if inSr.factoryCode != wgs84Sr.factoryCode:
                    start_pt_geom = arcpy.PointGeometry(start_pt, inSr)
                    start_pt = start_pt_geom.projectAs(wgs84Sr).getPart()
                point_array.add(start_pt)
                # loop through the line segments
                for j in range(1, point_count):
                    end_pt = part[j]
                    # project to WGS84 if needed
                    if inSr.factoryCode != wgs84Sr.factoryCode:
                        end_pt_geom = arcpy.PointGeometry(end_pt, inSr)
                        end_pt = end_pt_geom.projectAs(wgs84Sr).getPart()
                    # create a geographiclib line object
                    line_object = geod.InverseLine(start_pt.Y, start_pt.X, end_pt.Y, end_pt.X)
                    # determine how many densified segments there will be
                    n = int(math.ceil(line_object.s13 / spacing))
                    if line_object.s13 > spacing:
                        seglen = line_object.s13 / n
                        for k in range(1, n):
                            s = seglen * k
                            g = line_object.Position(s,
                                                     Geodesic.LATITUDE |
                                                     Geodesic.LONGITUDE |
                                                     Geodesic.LONG_UNROLL)
                            point = arcpy.Point(g['lon2'], g['lat2'])
                            # add densified points to output array
                            point_array.add(point)

                    # Convert last point back to the input CRS if needed
                    if inSr.factoryCode != wgs84Sr.factoryCode:
                        end_pt_geom = arcpy.PointGeometry(end_pt, wgs84Sr)
                        end_pt = end_pt_geom.projectAs(inSr).getPart()
                    # add the last point to the output array
                    point_array.add(end_pt)
                    start_pt = end_pt
                # Convert each point back to the input CRS if needed
                if inSr.factoryCode != wgs84Sr.factoryCode:
                    for k in range(len(point_array)):
                        pt_geom = arcpy.PointGeometry(point_array[k], wgs84Sr)
                        point_array[k] = pt_geom.projectAs(inSr).getPart()
                # cursor to write output layer
                cur = arcpy.da.InsertCursor(out_layer, field_names)
                if len(point_array) > 0:
                    # write output point array to output layer
                    row_list = list(row)
                    if create_polygon:
                        row_list[geom_field_index] = arcpy.Polygon(point_array)
                    elif create_polyline:
                        row_list[geom_field_index] = arcpy.Polyline(point_array)
                    out_row = tuple(row_list)
                    cur.insertRow(out_row)


# get input geometry type
create_polyline = False
create_polygon = False
desc = arcpy.Describe(in_layer)
if desc.featureType == "Simple" and desc.datasetType == "FeatureClass":

    if desc.shapeType == 'Point':
        print("Input is Point Type")
        # create a point memory layer
        out_layer_name = "Densified_Point_" + str(ellipsoid_name) + "_" + str(spacing) + "m"
        # create output layer
        arcpy.CreateFeatureclass_management(out_path, out_layer_name, "POINT", in_layer, spatial_reference=inSr)

        # add field for pointType (original or densified)
        fields = arcpy.ListFields(in_layer)
        fieldNameList = [field.name for field in fields]
        pointTypeField = ''
        for fieldName in ["pointType", "pntType", "pntTyp"]:
            if fieldName not in fieldNameList:
                pointTypeField = fieldName
                break
        arcpy.AddField_management(outLayer, pointTypeField, "TEXT", field_length=50)

        # run the densification
        densify_point(in_layer, outLayer, spacing, pointTypeField)

    elif desc.shapeType == 'Polyline':
        print("Input is Polyline Type: processing")
        create_polyline = True
        # create a polyline memory layer
        out_layer_name = "Densified_Polyline_" + str(ellipsoid_name) + "_" + str(spacing) + "m"
        arcpy.CreateFeatureclass_management(out_path, out_layer_name, "POLYLINE", in_layer, spatial_reference=inSr)
        print("created output: " + outLayer)
        densify_poly(in_layer, outLayer, spacing)

    elif desc.shapeType == 'Polygon':
        print("Input is Polygon Type: processing")
        create_polygon = True
        # create a polygon memory layer
        out_layer_name = "Densified_Polygon_" + str(ellipsoid_name) + "_" + str(spacing) + "m"
        arcpy.CreateFeatureclass_management(out_path, out_layer_name, "POLYGON", in_layer, spatial_reference=inSr)
        print("created output: " + outLayer)
        densify_poly(in_layer, outLayer, spacing)
else:
    print("Tool only works with simple features (point/polyline/polygon)")

# add memory layer to map TOC
mxd = arcpy.mapping.MapDocument("CURRENT")
data_frame = arcpy.mapping.ListDataFrames(mxd)[0]
layer = arcpy.mapping.Layer(os.path.join(out_path, out_layer_name))
arcpy.mapping.AddLayer(data_frame, layer, "AUTO_ARRANGE")
print("finished processing")
