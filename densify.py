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
inLayer = r"C:\Users\Jonah\PycharmProjects\arcpy_densifier\testData\testLineAu.shp"
outLayer = r"C:\Users\Jonah\PycharmProjects\arcpy_densifier\testData\outtestLineAu.shp"
outPath = os.path.dirname(outLayer)
outLayerName = os.path.basename(outLayer)
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


print("inLayer: " + inLayer)

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
inSr = arcpy.Describe(inLayer).spatialReference
wgs84Sr = arcpy.SpatialReference(4326)

# Create a geographiclib Geodesic object
geod = Geodesic(ellipsoid_a, 1 / ellipsoid_f)


def densifyPoint(inLayer, outLayer, spacing, pointTypeField):
    counter = 0
    currentPoint = arcpy.Point()
    inFieldNames = ['SHAPE@' if f.type == 'Geometry' else f.name for f in arcpy.ListFields(inLayer)]
    outFieldNames = inFieldNames[:]
    outFieldNames.append(pointTypeField)
    geomFieldIndex = inFieldNames.index('SHAPE@')
    # cursor to write output layer
    cur = arcpy.da.InsertCursor(outLayer, "*")
    # loop through the features and retrieve [X,Y]
    with arcpy.da.SearchCursor(inLayer, inFieldNames) as cursor:
        for row in cursor:
            # this is the first point
            if counter == 0:
                # get the point geometry
                ptGeom = row[geomFieldIndex]
                point = ptGeom.getPart()
                # save point geometry for later
                currentPoint.X = point.X
                currentPoint.Y = point.Y
                # write to output layer
                rowList = list(row)
                rowList.append('Original')
                row = tuple(rowList)
                cur.insertRow(row)

            # This is for subsequent points
            elif counter > 0:
                startPt = currentPoint
                endPtGeom = row[geomFieldIndex]
                endPt = endPtGeom.getPart()
                # project to WGS84 if needed
                if inSr.factoryCode != wgs84Sr.factoryCode:
                    startPtGeom = arcpy.PointGeometry(startPt, inSr)
                    startPt = startPtGeom.projectAs(wgs84Sr).getPart()
                    endPtGeom = endPtGeom.projectAs(wgs84Sr)
                    endPt = endPtGeom.getPart()
                # create a geographiclib line object
                lineObject = geod.InverseLine(startPt.Y, startPt.X, endPt.Y, endPt.X)
                # determine how many densified segments there will be
                n = int(math.ceil(lineObject.s13 / spacing))
                # adjust the spacing distance
                seglen = lineObject.s13 / n
                for i in range(1, n):
                    if i > 0:
                        s = seglen * i
                        g = lineObject.Position(s,
                                                Geodesic.LATITUDE |
                                                Geodesic.LONGITUDE |
                                                Geodesic.LONG_UNROLL)
                        point = arcpy.Point(g['lon2'], g['lat2'])
                        currentPoint.X = point.X
                        currentPoint.Y = point.Y
                        # Convert each point back to the input CRS
                        if inSr.factoryCode != wgs84Sr.factoryCode:
                            ptGeom = arcpy.PointGeometry(point, wgs84Sr)
                            point = ptGeom.projectAs(inSr).getPart()
                        # write to output layer
                        rowList = list(row)
                        rowList.append("Densified")
                        rowList[geomFieldIndex] = point
                        outRow = tuple(rowList)
                        cur.insertRow(outRow)
                ptGeom = row[geomFieldIndex].getPart()
                currentPoint.X = ptGeom.X
                currentPoint.Y = ptGeom.Y
                # write to output layer
                rowList = list(row)
                rowList.append('Original')
                row = tuple(rowList)
                cur.insertRow(row)
            counter += 1
    if cur:
        del cur


def densifyLine(inLayer, outLayer, spacing):
    # loop through the features
    fieldNames = ['SHAPE@' if f.type == 'Geometry' else f.name for f in arcpy.ListFields(inLayer)]
    geomFieldIndex = fieldNames.index('SHAPE@')
    with arcpy.da.SearchCursor(inLayer, fieldNames) as cursor:
        for row in cursor:
            # create array to hold densified points
            pointArray = arcpy.Array()
            partCount = row[geomFieldIndex].partCount
            for i in range(0,partCount):
                # get the geometry
                part = row[geomFieldIndex].getPart(i)
                pointCount = len(part)
                startPt = part[0]
                # project to WGS84 if needed
                if inSr.factoryCode != wgs84Sr.factoryCode:
                    startPtGeom = arcpy.PointGeometry(startPt, inSr)
                    startPt = startPtGeom.projectAs(wgs84Sr).getPart()
                pointArray.add(startPt)
                print("startPt:")
                print(startPt)
                # loop through the line segments
                for j in range(1, pointCount):
                    endPt = part[j]
                    # project to WGS84 if needed
                    if inSr.factoryCode != wgs84Sr.factoryCode:
                        endPtGeom = arcpy.PointGeometry(endPt, inSr)
                        endPt = endPtGeom.projectAs(wgs84Sr).getPart()
                    print("endPt:")
                    print(endPt)
                    # create a geographiclib line object
                    lineObject = geod.InverseLine(startPt.Y, startPt.X, endPt.Y, endPt.X)
                    # determine how many densified segments there will be
                    n = int(math.ceil(lineObject.s13 / spacing))
                    print("n:")
                    print(n)
                    if lineObject.s13 > spacing:
                        seglen = lineObject.s13 / n
                        print("seglen:")
                        print(seglen)
                        for k in range(1, n):
                            s = seglen * k
                            g = lineObject.Position(s,
                                                    Geodesic.LATITUDE |
                                                    Geodesic.LONGITUDE |
                                                    Geodesic.LONG_UNROLL)
                            point = arcpy.Point(g['lon2'], g['lat2'])
                            # add densified points to output array
                            pointArray.add(point)

                    # Convert last point back to the input CRS if needed
                    if inSr.factoryCode != wgs84Sr.factoryCode:
                        endPtGeom = arcpy.PointGeometry(endPt, wgs84Sr)
                        endPt = endPtGeom.projectAs(inSr).getPart()
                    # add the last point to the output array
                    pointArray.add(endPt)
                    startPt = endPt
            # Convert each point back to the input CRS if needed
            if inSr.factoryCode != wgs84Sr.factoryCode:
                for i in range(len(pointArray)):
                    ptGeom = arcpy.PointGeometry(pointArray[i], wgs84Sr)
                    pointArray[i] = ptGeom.projectAs(inSr).getPart()
            # cursor to write output layer
            cur = arcpy.da.InsertCursor(outLayer, fieldNames)
            if len(pointArray) > 0:
                # write output point array to output layer
                rowList = list(row)
                rowList[geomFieldIndex] = arcpy.Polyline(pointArray)
                outRow = tuple(rowList)
                cur.insertRow(outRow)


def densifyPolygon(inLayer, outLayer, spacing):
    # loop through the features
    with arcpy.da.SearchCursor(inLayer, ["SHAPE@"]) as cursor:
        for row in cursor:
            polyPointArray = arcpy.Array()
            for part in row:
                # get the geometry
                arr = part.getPart(0)
                pointCount = len(arr)
                startPt = arr[0]
                if inSr != wgs84Sr:
                    startPt = startPt.projectAs(wgs84Sr)
                polyPointArray.add(startPt)
                for i in range(1, pointCount):
                    endPt = arr[i]
                    if inSr != wgs84Sr:
                        endPt = endPt.projectAs(wgs84Sr)
                    lineObject = geod.InverseLine(startPt.Y, startPt.X, endPt.Y, endPt.X)
                    n = int(math.ceil(lineObject.s13 / spacing))
                    seglen = lineObject.s13 / n
                    for j in range(1, n):
                        s = seglen * j
                        g = lineObject.Position(s,
                                                Geodesic.LATITUDE |
                                                Geodesic.LONGITUDE |
                                                Geodesic.LONG_UNROLL)
                        polyPointArray.add(arcpy.Point(g['lon2'], g['lat2']))
                    polyPointArray.add(endPt)
                    startPt = endPt

                    # Convert each point back to the output CRS
                    if inSr != wgs84Sr:
                        for x, pt in enumerate(polyPointArray):
                            polyPointArray[x] = pt.projectAs(wgs84Sr)

            if len(polyPointArray) > 0:
                # cursor to write output layer
                with arcpy.da.InsertCursor(outLayer, ["SHAPE@"]) as writeCursor:
                    # write output point array to output layer
                    writeCursor.insertRow([arcpy.Polygon(polyPointArray)])


# get input geometry type
desc = arcpy.Describe(inLayer)
if desc.shapeType == 'Point':
    print("Input is Point Type")
    # create and add to map canvas a point memory layer
    # layer_name = "Densified Point " + str(ellipsoid_name) + " " + str(spacing) + "m"
    # create output layer
    arcpy.CreateFeatureclass_management(outPath, outLayerName, "POINT", inLayer, spatial_reference=inSr)

    # add field for pointType (original or densified)
    # get the field list
    fields = arcpy.ListFields(inLayer)
    fieldNameList = [field.name for field in fields]
    pointTypeField = ''
    for fieldName in ["pointType", "pntType", "pntTyp"]:
        if fieldName not in fieldNameList:
            pointTypeField = fieldName
            break
    arcpy.AddField_management(outLayer, pointTypeField, "TEXT", field_length=50)

    # new field list
    fields = arcpy.ListFields(inLayer)

    # run the densification
    densifyPoint(inLayer, outLayer, spacing, pointTypeField)

elif desc.shapeType == 'MultiPoint':
    print("Input is MultiPoint Type: not implemented")

elif desc.shapeType == 'Polyline':
    print("Input is Polyline Type: processing")
    create_polyline = True
    # layer_name = "Densified Line " + str(ellipsoid_name) + " " + str(spacing) + "m"
    # out_line_layer = arcpy.CopyFeatures_management(inLayer, r"in_memory\\" + layer_name)
    arcpy.CreateFeatureclass_management(outPath, outLayerName, "POLYLINE", inLayer, spatial_reference=inSr)
    print("created output: " + outLayer)
    densifyLine(inLayer, outLayer, spacing)

elif desc.shapeType == 'Polygon':
    print("Input is Polygon Type: still in development")
    create_polygon = True
    # create and add to map canvas a polyline memory layer
    layer_name = "Densified Polygon " + str(ellipsoid_name) + " " + str(spacing) + "m"
    # out_point_layer = arcpy.CopyFeatures_management(inLayer, r"in_memory\\" + layer_name)
    arcpy.CreateFeatureclass_management(outPath, outLayerName, "POLYGON", inLayer, spatial_reference=inSr)

elif desc.shapeType == 'MultiPatch':
    print("Input is MultiPatch Type: not implemented")

print("finished processing")
