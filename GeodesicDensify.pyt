import sys
import os
import arcpy

scripts_dir = os.path.join(os.path.dirname(__file__), 'scripts')
sys.path.append(scripts_dir)
# Do not compile .pyc files for the tool modules.
sys.dont_write_bytecode = True

from GeodesicDensification_geographiclib import GeodesicDensification_geographiclib
from GeodesicDensification_arcpy import GeodesicDensification_arcpy

class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Toolbox"
        self.alias = "Geodesic Densification using arcpy or geographiclib"

        # List of tool classes associated with this toolbox
        self.tools = [GeodesicDensification_geographiclib,
                      GeodesicDensification_arcpy]
