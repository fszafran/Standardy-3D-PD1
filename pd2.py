import geopandas as gpd
import rasterio
from rasterio.mask import mask
from owslib.wfs import WebFeatureService
import requests
from pathlib import Path
import os
from typing import Union
from osgeo import gdal
import zipfile 
import open3d
import numpy as np
import shapely

def download_and_save_file(download_url: str, save_path: Union[Path, str]) -> None:
    response = requests.get(download_url)

    if response.status_code != 200:
        raise Exception(f"Failed to download {download_url}. Status code: {response.status_code}")

    with open(save_path, "wb") as file:
        file.write(response.content)

def download(aoiBounds, source):
    version = "1.0.0"
    if source == "NMT":
        wfsServiceUrl = "https://mapy.geoportal.gov.pl/wss/service/PZGIK/NumerycznyModelTerenuEVRF2007/WFS/Skorowidze"
        typeName = "gugik:SkorowidzNMT2023"
    elif source == "NMPT":
        wfsServiceUrl = "https://mapy.geoportal.gov.pl/wss/service/PZGIK/NumerycznyModelPokryciaTerenuEVRF2007/WFS/Skorowidze"
        typeName = "gugik:SkorowidzNMPT2023"
    elif source == "BDOT":
        wfsServiceUrl = "https://mapy.geoportal.gov.pl/wss/service/PZGIK/BDOT/WFS/PobieranieBDOT10k"
        typeName = "ms:BDOT10k_powiaty"
        version= "2.0.0"
    wfsService = WebFeatureService(
        url = wfsServiceUrl,
        version=version
    )
    response = wfsService.getfeature(
        bbox = tuple(aoiBounds),
        typename = [typeName]
    )
    return response

def rasterize(filePath, tile, destinationFolder, source):
    with rasterio.open(filePath, 'r') as src:
        out_image, out_transform = mask(src, [tile], crop=True)
        out_meta = src.meta
    out_meta.update({"driver": "GTiff",
                     "height": out_image.shape[1],
                     "width": out_image.shape[2],
                     "transform": out_transform})
    output_file = f"{destinationFolder}/{source}.tif"
    with rasterio.open(output_file, "w", **out_meta) as dest:
        dest.write(out_image)

def read_raster(filename):
    raster = gdal.Open(filename)
    return raster

def create_vertex_array(raster):
    transform = raster.GetGeoTransform()
    width = raster.RasterXSize
    height = raster.RasterYSize
    x = np.arange(0, width) * transform[1] + transform[0]
    y = np.arange(0, height) * transform[5] + transform[3]
    xx, yy = np.meshgrid(x, y)
    zz = raster.ReadAsArray()
    vertices = np.vstack((xx, yy, zz)).reshape([3, -1]).transpose()
    return vertices

def extrude_polygon(geometry, minHeight, maxHeight):
    polygon = shapely.Polygon.from_bounds(geometry.bounds[0], geometry.bounds[1], geometry.bounds[2], geometry.bounds[3])
    min_height = minHeight
    max_height = maxHeight

    xy = polygon.exterior.coords
    xy = list(xy[:-1])
    xy = np.float64(xy)
    centroid = np.float64(polygon.representative_point().xy).flatten()

    base_num_vertices = xy.shape[0]
    xyz = np.vstack([
        np.hstack([
            xy, np.full((xy.shape[0], 1), min_height)
        ]),
        np.hstack([
            xy, np.full((xy.shape[0], 1), max_height)
        ]),
        np.hstack([centroid, min_height]),
        np.hstack([centroid, max_height])
    ])

    walls_triangles = []
    for i in range(base_num_vertices - 1):
        walls_triangles.append(
            [i, i + 1, i + 1 + base_num_vertices][::-1]
        )
        walls_triangles.append(
            [i + 1 + base_num_vertices, i + base_num_vertices, i][::-1]
        )
    walls_triangles.append(
        [base_num_vertices - 1, 0, base_num_vertices][::-1]
    )
    walls_triangles.append(
        [base_num_vertices, 2 * base_num_vertices - 1, base_num_vertices - 1][::-1]
    )

    base_triangles = []
    for i in range(base_num_vertices - 1):
        base_triangles.append(
            [i, 2 * base_num_vertices, i + 1][::-1]
        )
        base_triangles.append(
            [base_num_vertices + i, 2 * base_num_vertices + 1, base_num_vertices + i + 1]
        )
    base_triangles.append(
        [base_num_vertices - 1, 2 * base_num_vertices, 0][::-1]
    )
    base_triangles.append(
        [2 * base_num_vertices - 1, 2 * base_num_vertices + 1, base_num_vertices]
    )

    mesh = open3d.geometry.TriangleMesh(
        open3d.utility.Vector3dVector(xyz),
        open3d.utility.Vector3iVector(base_triangles + walls_triangles)
    )
    mesh.compute_vertex_normals()
    mesh.paint_uniform_color([0.0, 0.6, 0.0])
    return mesh
    
def create_index_array(raster):
    width = raster.RasterXSize
    height = raster.RasterYSize

    ai = np.arange(0, width - 1)
    aj = np.arange(0, height - 1)
    aii, ajj = np.meshgrid(ai, aj)
    a = aii + ajj * width
    a = a.flatten()

    tria = np.vstack((a, a + width, a + width + 1, a, a + width + 1, a + 1))
    tria = np.transpose(tria).reshape([-1, 3])
    return tria

def makeMosaic(rasterList,source, destinationFolder):
    g = gdal.Warp(f"{destinationFolder}/merged{source}.tif", rasterList, format="GTiff", options=["COMPRESS=LZW", "TILED=YES"])
    g = None
def getBUBD_A(zip_file, file_path, destinationFolder):
    with zipfile.ZipFile(zip_file, 'r') as zip_ref:
        if file_path in zip_ref.namelist():
            zip_ref.extract(file_path, f"{destinationFolder}")
            bubda = gpd.read_file(f"{destinationFolder}/{file_path}")
        else:
            print("File not found in the zip file.")

        gdal.SetConfigOption('SHAPE_RESTORE_SHX', 'YES')
        gdal.VectorTranslate(f"{destinationFolder}/bubda.shp", f"{destinationFolder}/PL.PZGiK.336.2209/BDOT10k/PL.PZGiK.336.2209__OT_BUBD_A.xml", options='-f "ESRI Shapefile"')
        os.remove(f"{destinationFolder}/{file_path}")

def clipSHP(aoiBounds, destinationFolder):
    features=gpd.read_file(f"{destinationFolder}/bubda.shp")
    features=features.clip(aoiBounds, keep_geom_type=True)
    features.to_file(f"{destinationFolder}/bubda.shp")
def save_mesh(plyFile, mesh):
    mesh = open3d.geometry.TriangleMesh.create_sphere()
    mesh.compute_vertex_normals()
    open3d.io.write_triangle_mesh("mesh.ply", mesh)

def ConstructLOD(hextilesPath, tileIndex, destinationFolder):
    tiles = gpd.read_file(hextilesPath)
    tiles = tiles.to_crs("EPSG:2180")
    myTile = tiles.iloc[tileIndex]
    myTile= myTile["geometry"]
    aoiBounds=myTile.bounds
    nmt = download(aoiBounds, "NMT")
    nmpt = download(aoiBounds, "NMPT")
    bdot = download(aoiBounds, "BDOT")
    os.makedirs(f"{destinationFolder}", exist_ok=True)
    with open(f"{destinationFolder}/nmt.xml", 'wb') as file:
        file.write(nmt.read())
    with open(f"{destinationFolder}/nmpt.xml", 'wb') as file:
        file.write(nmpt.read())
    with open(f"{destinationFolder}/bdot.xml", 'wb') as file:
        file.write(bdot.read())
    XMLnmt = gpd.read_file(f"{destinationFolder}/nmt.xml")
    os.remove(f"{destinationFolder}/nmt.xml")
    XMLnmtp = gpd.read_file(f"{destinationFolder}/nmpt.xml")
    os.remove(f"{destinationFolder}/nmpt.xml")
    XMLbdot = gpd.read_file(f"{destinationFolder}/bdot.xml")
    os.remove(f"{destinationFolder}/bdot.xml")
    nmtRasters = []
    nmtpRasters = []
    for i in range (len(XMLnmt['url_do_pobrania'])):
        download_and_save_file(XMLnmt['url_do_pobrania'][i], f"{destinationFolder}/nmt{i}.asc")
        rasterize(f"{destinationFolder}/nmt{i}.asc", myTile, destinationFolder, f"nmt{i}")
        os.remove(f"{destinationFolder}/nmt{i}.asc")
        nmtRasters.append(f"{destinationFolder}/nmt{i}.tif")
    for i in range (len(XMLnmtp['url_do_pobrania'])):
        download_and_save_file(XMLnmtp['url_do_pobrania'][i], f"{destinationFolder}/nmpt{i}.asc")
        rasterize(f"{destinationFolder}/nmpt{i}.asc", myTile, destinationFolder, f"nmpt{i}")
        os.remove(f"{destinationFolder}/nmpt{i}.asc")
        nmtpRasters.append(f"{destinationFolder}/nmpt{i}.tif")
    for i in range (len(XMLbdot['URL_GML'])):
        download_and_save_file(XMLbdot['URL_GML'][i], f"{destinationFolder}/bdot{i}.zip")
        getBUBD_A(f"{destinationFolder}/bdot{i}.zip", "PL.PZGiK.336.2209/BDOT10k/PL.PZGiK.336.2209__OT_BUBD_A.xml", destinationFolder)
    
    makeMosaic(nmtRasters, "Nmt", destinationFolder)
    makeMosaic(nmtpRasters, "Nmpt", destinationFolder)
    for rasters in nmtRasters:
        os.remove(rasters)
    for rasters in nmtpRasters:
        os.remove(rasters)

    clipSHP(myTile, destinationFolder)
    raster = read_raster(f"{destinationFolder}/mergedNmt.tif")
    vertices = create_vertex_array(raster)
    triangles = create_index_array(raster)
    mesh = open3d.geometry.TriangleMesh(
        open3d.utility.Vector3dVector(vertices),
        open3d.utility.Vector3iVector(triangles)
    )
    mesh.compute_vertex_normals()
    open3d.io.write_triangle_mesh(f"{destinationFolder}/mergedNmt_mesh.ply", mesh)
    centerPointsNmp=[]
    nmtHeights=[]
    nmptHeights=[]
    extruded=[]
    buildings =[]
    BubdaObjects = gpd.read_file(f"{destinationFolder}/bubda.shp")
    for geometry in BubdaObjects['geometry']:
        centerPointsNmp.append(geometry.centroid)
        buildings.append(geometry)
    with rasterio.open(f"{destinationFolder}/mergedNmt.tif") as src:
        Inmt = src.read(1)
        for point in centerPointsNmp:
            row, col = src.index(point.x, point.y)
            height = Inmt[row, col]
            nmtHeights.append(height)
    with rasterio.open(f"{destinationFolder}/mergedNmpt.tif") as src:
        Inmpt = src.read(1)
        for point in centerPointsNmp:
            row, col = src.index(point.x, point.y)
            height = Inmpt[row, col]
            nmptHeights.append(height)
    for i in range(len(centerPointsNmp)):
        extruded.append(extrude_polygon(buildings[i], nmtHeights[i], nmptHeights[i]))
    mergedMesh = open3d.geometry.TriangleMesh()
    for mesh in extruded:
        mergedMesh += mesh
    open3d.io.write_triangle_mesh(f"{destinationFolder}/budyneczki.ply", mergedMesh)
    mesh1 = open3d.io.read_triangle_mesh(f"{destinationFolder}/mergedNmt_mesh.ply")
    mesh2 = open3d.io.read_triangle_mesh(f"{destinationFolder}/budyneczki.ply")
    vertices1 = np.asarray(mesh1.vertices)
    faces1 = np.asarray(mesh1.triangles)
    vertices2 = np.asarray(mesh2.vertices)
    faces2 = np.asarray(mesh2.triangles)
    combined_vertices = np.vstack((vertices1, vertices2))
    combined_faces = np.vstack((faces1, faces2 + len(vertices1))) 
    combined_mesh = open3d.geometry.TriangleMesh(
        open3d.utility.Vector3dVector(combined_vertices),
        open3d.utility.Vector3iVector(combined_faces)
    )
    open3d.visualization.draw_geometries([combined_mesh])
    open3d.io.write_triangle_mesh(f"{destinationFolder}/final.ply", combined_mesh)

#main function
if __name__ == '__main__':
    ConstructLOD("C:/Users/filo1/Desktop/szkola_sem4/Standardy_3D/PD2/hextiles1.fgb", 21, "C:/Users/filo1/Desktop/szkola_sem4/Standardy_3D/PD2/destinationFolder")

