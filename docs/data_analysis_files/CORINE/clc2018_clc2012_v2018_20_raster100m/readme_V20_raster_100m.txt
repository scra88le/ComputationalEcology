README.TXT

Corine Land Cover European seamless 100m raster database (Version 20 - final delivery)
-------------------------------------------------------------------------------

This directory contains the Corine Land Cover European seamless raster database 
version 20 (dated 05/2019) in GeoTiff format. 

Vector databases was provided by National Teams within I&CLC2000, CARDS, CLC2006, CLC2012 and CLC2018 project. 
All features in original vector database were classified and digitised based on satellite images 
with 100 m positional accuracy (according to CLC specifications) and 25 ha minimum mapping unit (5ha MMU 
for change layer) into the standardized CLC nomenclature (44 CLC classes). EU CORINE Land Cover 
seamless DB represents the final product of European CLC2000 data integration. 
Harmonization on country border done for CLC2000 and CHA9000 (2km strip) for other themes only small artefact 
along borders (0,1ha) has been removed.
All seamless vector layers were rasterized to the 100m resolution grid. CELL_CENTER method was used for the rasterizing.  


100 m raster part of the V20 database delivery contains:

readme_V20_raster_100m.txt  	.... readme file
CLC_country_coverage_v20.pdf  .... country coverage for all CLC themes


CLC2000_CLC1990_V2018.20.tif  .... representing CORINE Land Cover data for reference year 1990 in 100m resolution
CLC2006_CLC2000_V2018.20.tif  .... representing merged CORINE Land Cover data for reference year 2000 and for reference year 2000 - revised in 100m resolution
CLC2012_CLC2006_V2018.20.tif  .... representing merged CORINE Land Cover data for reference year 2006 and for reference year 2006 - revised in 100m resolution
clc2018_clc2012_V2018.20.tif  .... representing CORINE Land Cover data for reference year 2012-revised in 100m resolution
clc2018_clc2018_V2018.20.tif  .... representing CORINE Land Cover data for reference year 2018 in 100m resolution
CLC2000_CHA9000_V2018.20.tif  .... representing CORINE Land Cover Change data (change 1990-2000) in 100m resolution
CLC2006_CHA0006_V2018.20.tif  .... representing CORINE Land Cover Change data (change 2000-2006) in 100m resolution
CLC2012_CHA0612_V2018.20.tif  .... representing CORINE Land Cover Change data (change 2006-2012) in 100m resolution
clc2018_cha1218_V2018.20.tif	.... representing CORINE Land Cover Change data (change 2012-2018) in 100m resolution

clc2018_cha1218_V2018.20_FR_GLP.tif	.... representing CORINE Land Cover Change data (change 2012-2018) in 100m resolution for Guadeloupe
clc2018_cha1218_V2018.20_FR_GUF.tif	.... representing CORINE Land Cover Change data (change 2012-2018) in 100m resolution for French Guiana
clc2018_cha1218_V2018.20_FR_MTQ.tif	.... representing CORINE Land Cover Change data (change 2012-2018) in 100m resolution for Martinique
clc2018_cha1218_V2018.20_FR_MYT.tif	.... representing CORINE Land Cover Change data (change 2012-2018) in 100m resolution for Mayotte
clc2018_cha1218_V2018.20_FR_REU.tif	.... representing CORINE Land Cover Change data (change 2012-2018) in 100m resolution for Reunion

clc2018_clc2012_V2018.20_FR_GLP.tif  .... representing CORINE Land Cover data for reference year 2012-revised in 100m resolution for Guadeloupe
clc2018_clc2012_V2018.20_FR_GUF.tif  .... representing CORINE Land Cover data for reference year 2012-revised in 100m resolution for French Guiana
clc2018_clc2012_V2018.20_FR_MTQ.tif  .... representing CORINE Land Cover data for reference year 2012-revised in 100m resolution for Martinique
clc2018_clc2012_V2018.20_FR_MYT.tif  .... representing CORINE Land Cover data for reference year 2012-revised in 100m resolution for Mayotte
clc2018_clc2012_V2018.20_FR_REU.tif  .... representing CORINE Land Cover data for reference year 2012-revised in 100m resolution for Reunion


clc2018_clc2018_V2018.20_FR_GLP.tif  .... representing CORINE Land Cover data for reference year 2018 in 100m resolution for Guadeloupe
clc2018_clc2018_V2018.20_FR_GUF.tif  .... representing CORINE Land Cover data for reference year 2018 in 100m resolution for French Guiana
clc2018_clc2018_V2018.20_FR_MTQ.tif  .... representing CORINE Land Cover data for reference year 2018 in 100m resolution for Martinique
clc2018_clc2018_V2018.20_FR_MYT.tif  .... representing CORINE Land Cover data for reference year 2018 in 100m resolution for Mayotte
clc2018_clc2018_V2018.20_FR_REU.tif  .... representing CORINE Land Cover data for reference year 2018 in 100m resolution for Reunion

	
The Coordinate Reference Systems (CRS):

For European seamless layers:

CRS Name: EUR_ETRS89/LAEA1052 
 
Projection: Lambert Azimuthal - Equal Area projection
        longitude of origin             10d00'00.0000"E
        latitude of origin              52d00'00.0000"N
        false easting                   4321000.000
        false northing                  3210000.000
Datum: ETRS89 (European Terrestrial Reference System 1989)
        type 				geodetic 
	valid area 			Europe / EUREF
Ellipsoid: GRS 80 (New International) 
	semi major axis 		6 378 137 m 
	inverse flattening 		298.257222101 

For Guadeloupe:

CRS name: WGS_1984_UTM_Zone_20N

Projection: Transverse Mercator
	central meridian	-63� E     
	latitude of origin 	0� N
	false northing 	0 m 
	false easting 		500 000 m 
	scale factor		0.9996

Datum: WGS_1984
	type  			geodetic  
	valid area  		France - Guadeloupe
Ellipsoid: WGS 1984
    semi major axis 		6 378 137 m
    inverse flattening 		298.257223563
   
For French Guiana:

CRS name: WGS_1984_UTM_Zone_22N

Projection: Transverse Mercator
	central meridian	-51� E     
	latitude of origin 	0� N
	false northing 	0 m 
	false easting 		500 000 m 
	scale factor		0.9996
Datum: WGS_1984
	type  			geodetic  
	valid area  		France - French Guiana
Ellipsoid: WGS 1984
    semi major axis 		6 378 137 m
    inverse flattening 		298.257223563
    
For Martinique:

CRS name: WGS_1984_UTM_Zone_20N

Projection: Transverse Mercator
	central meridian	-63� E     
	latitude of origin 	0� N
	false northing 	0 m 
	false easting 		500 000 m 
	scale factor		0.9996
Datum: WGS_1984
	type  			geodetic  
	valid area  		France - Martinique
Ellipsoid: WGS 1984
    semi major axis 		6 378 137 m
    inverse flattening 		298.257223563
    
For Mayotte:

CRS name: WGS_1984_UTM_Zone_38S

Projection: Transverse Mercator
	central meridian	45� E     
	latitude of origin 	0� N
	false northing 	10 000 000 m 
	false easting 		500 000 m 
	scale factor		0.9996

Datum: WGS_1984
	type  			geodetic  
	valid area  		France - Mayotte
Ellipsoid: WGS 1984
    semi major axis 		6 378 137 m
    inverse flattening 		298.257223563
  
For Reunion:

CRS name : WGS_1984_UTM_Zone_40S

Projection: Transverse Mercator
	central meridian	57� E     
	latitude of origin 	0� N
	false northing 	10 000 000 m 
	false easting 		500 000 m 
	scale factor		0.9996
Datum: WGS_1984
	type  			geodetic  
	valid area  		France - R�union
Ellipsoid: WGS 1984
    semi major axis 		6 378 137 m
    inverse flattening 		298.257223563
    


Changes from previous release - change log:

Version 20 (final version)
TR was added to CLC2012, CLC2018, CHA1218 layers. For actual coverage see CLC_country_coverage_V20.pdf
Repaired FR data were uploaded to CLC2012, CLC2018, CHA1218 layers.
Data for FR DOMs were added (CLC2012, CLC2018, CHA1218). For actual coverage see CLC_country_coverage_V20.pdf
Faroe Islands were added (CLC2018). For actual coverage see CLC_country_coverage_V20.pdf

Version 20b2 (prefinal version)
Layer CLC2012 has been revised. For actual coverage see CLC_country_coverage_V20b2.pdf
New layer CLC2018 has been created. For actual coverage see CLC_country_coverage_V20b2.pdf
New layer CHA1218 has been created. For actual coverage see CLC_country_coverage_V20b2.pdf

Version 18_5
Spain and Turkey were completed. For actual coverage see CLC_coverage_v18_5.pdf
CLC00 and CLC00-revised layers were merged to the one CLC00 layer. For actual coverage see CLC_coverage_v18_5.pdf.
CLC06 and CLC06-revised layers were merged to the one CLC06 layer. For actual coverage see CLC_coverage_v18_5.pdf.

Version 18_4
New parts were included to clc90, clc00, cha00 - Guernsey and Jersey islands and to clc06, cha06 - Guernsey and Jersey islands and Greece. For actual coverage see CLC_coverage_v18_4.pdf
New layer clc12 has been created. For actual coverage see CLC_coverage_v18_4.pdf.
New layer cha12 has been created. For actual coverage see CLC_coverage_v18_4.pdf.
New layer clc06 revised has been created. For actual coverage see CLC_coverage_v18_4.pdf.

Version 17
New part was included to clc90, clc00, clc06, cha00 and cha06 - Azores. For actual coverage see CLC_coverage_v17.pdf

version 16
New country was included in CHA00 and CLC90 - Turkey.For actual coverage see CLC_coverage_v16.pdf.
New layer clc00_revised has been created. For actual coverage see CLC_coverage_v16.pdf.
Shift in Malta position has been repaired - all CLC layers for Malta has been re-projected.
Few polygons inconsistences in CLC layers has been repaired - see details in readme for particular CLC theme included.


version 15 (V5)
New country was included in CHA06 and CLC06 as result of CLC2006 mapping - Great Britain.For actual coverage see CLC_coverage_v14.pdf.
A simplified border matching is applied on CLC06 layer-see details in readme for particular CLC theme included.

Version 14 (V4)
New countries were included in CLC00, CHA06 and CLC06 as result of CLC2006 mapping - Switzerland, Northern Ireland, Madeira island. 
New countries were included in CLC90 and CHA00 as available. For actual coverage see CLC_coverage_v14.pdf.
Seamless database has been further improved addressing action points/feedback from EEA on version 13 (V3).
See details in readme for particular CLC theme included.

Version 3
First real seamless release in ESRI Geodatabase format.
New layers (CLC06 and CHA06) included as result of CLC2006 mapping (36 countries) 
New countries were included in CLC90, CLC00 and CHA00 as available from CLC2006 related mapping.
Sea buffers has been introduced (15km as proxy to 12 nautical miles sea zone).

version 2
All issues in tiles or tile borders, which has been reported by users was corrected.
Original idea was to deliver one big seamless layer, but as it has showed to be not feasible in 
ESRI environment yet it was decided to produce seamless tiles again (but free of any tiling artefacts).
New country deliveries available integrated: Iceland, Montenegro, Serbia, Kosovo from CARDS and CLC2006 projects 
were included. Limited country harmonization done for these add-ons (small artefacts cleaned only).


version 1
Initial version includes all I&CLC2000 project countries
Harmonization on country border done (2km strip)


Prepared by
GISAT 08/2019
For more information contact tomas.soukup@gisat.cz, miroslav.kopecky@gisat.cz