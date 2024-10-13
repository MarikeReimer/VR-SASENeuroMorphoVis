
import bpy
import bmesh
import numpy as np 
from datetime import datetime
from pynwb import NWBFile, NWBHDF5IO
from pynwb.ophys import TwoPhotonSeries, OpticalChannel, ImageSegmentation, ImagingPlane
from pynwb.file import Subject
from pynwb.base import Images
#from pynwb.image import RGBImage

root_collection = bpy.data.collections["CollectionToNWB"]

def NWBImageSegmentation():
    #Create NWB objects and metadata
    #Create pynwb subject
    subject = Subject(
        age = 'age',
        description = 'subject_description',
        genotype = 'genotype',
        sex = 'sex',
        species = 'species',
        strain = 'strain',
        subject_id = 'subject_id'
        )

    #Create pywnb File
    nwbfile = NWBFile(
        experimenter = 'experimenter',
        experiment_description = 'experiment_description',
        file_create_date = datetime.now(),
        identifier = 'identifier',
        institution = 'institution',
        lab = 'lab',
        notes = 'notes', 
        pharmacology = 'pharmacology',
        protocol = 'protocol',
        session_description = 'session_description',
        session_start_time = datetime.now(),
        slices = 'slices',
        surgery = 'surgery', 
        subject = subject
        )  

    #Extract Strings (from Blender's Scene)
    plane_name = 'bpy.context.scene.plane_name'
    device = 'bpy.context.scene.device'
    optical_channel_name = 'bpy.context.scene.optical_channel_name'
    optical_channel_description = 'bpy.context.scene.optical_channel_description'
    emission_lambda = 42.1

    #Add selected device        
    device = nwbfile.create_device(name = device)

    #Create optical channel
    optical_channel = OpticalChannel(
    optical_channel_name,
    optical_channel_description,
    emission_lambda)

    #Extract Imaging Plane Strings
    plane_name = 'bpy.context.scene.plane_name'
    plane_description = 'bpy.context.scene.plane_description'
    excitation_lambda = 42.1
    #external_file = 'bpy.context.scene.external_file'
    grid_spacing = 42.1
    imaging_rate = 42.1
    indicator = 'bpy.context.scene.indicator'
    location = 'bpy.context.scene.location'
    grid_spacing_unit = 'grid_spacing_unit'

    #Create imaging plane
    imaging_plane = nwbfile.create_imaging_plane(
        plane_name, 
        optical_channel, 
        plane_description, 
        device, 
        excitation_lambda,
        indicator, 
        location, 
        imaging_rate, 
        grid_spacing = [grid_spacing, grid_spacing, grid_spacing], 
        grid_spacing_unit = grid_spacing_unit)


    module = nwbfile.create_processing_module("MorphologyData", 'Contains processed morphology data from Blender.')   
    image_segmentation = ImageSegmentation()
    module.add(image_segmentation)
    return image_segmentation, nwbfile, imaging_plane

NWBObjects = NWBImageSegmentation()
image_segmentation = NWBObjects[0]
nwbfile = NWBObjects[1]
imaging_plane = NWBObjects[2]


#Extract location, volume, surface area of meshes 
def find_mesh_attributes(self):                  
    #Create and find center of mass
    self.select_set(True) #The origin_set operator only works on selected object
    bpy.ops.object.origin_set(type='ORIGIN_CENTER_OF_VOLUME')
    center_of_mass = self.location
    center_of_mass = [center_of_mass[0], center_of_mass[1], center_of_mass[2]]

    #Get mesh data from the active object in the Scene Collection
    mesh = self.data
    #Create an empty BMesh
    bm = bmesh.new()
    #Fill Bmesh with the mesh data from the object  
    bm.from_mesh(mesh)

    #Find Volume
    volume = bm.calc_volume(signed=False)
    #Find Surface Area
    surface_area = sum(i.calc_area() for i in bm.faces)
    bm.free()

    return center_of_mass, volume, surface_area


#TODO: The first plane segmentation adds for all objects

def BlenderMorphologyDataToNWB(obj, image_segmentation, imaging_plane):
    
    print('Object in NWB', obj)
    #Make dummy data to enable plane segmentation    
    pixel_mask = np.array([[0, 0, 0]] * 1)
    segmentation_name = obj.name + '_plane_segmentation'
    plane_segmentation = image_segmentation.create_plane_segmentation(
        name = segmentation_name,
        description = 'output from segmenting a mesh in Blender',
        imaging_plane = imaging_plane     
        #reference_images = ['images']
            )     
    mesh_attributes = find_mesh_attributes(obj)
    center_of_mass = mesh_attributes[0]
    volume = mesh_attributes[1]
    surface_area = mesh_attributes[2]

    #Create unique name
    segmentation_name = obj.name + '_plane_segmentaton'
    #Add columns to ROI to hold extracted variables
    plane_segmentation.add_column('center_of_mass', 'center_of_mass')
    plane_segmentation.add_column('volume', 'volume')
    plane_segmentation.add_column('surface_area', 'surface_area')        

    plane_segmentation.add_roi(                
        pixel_mask = pixel_mask,
        center_of_mass = center_of_mass,
        volume=volume,
        surface_area = surface_area,
        )    
    
    return plane_segmentation   
     
# # Start from the master collection or any specific collection
# def iterate_collections(collection, image_segmentation, imaging_plane):
#     print("Collection:", collection.name)
    
#     # Iterate through all objects in the collection
#     for obj in collection.objects:
#         print(f"  Object: {obj.name}")
#         BlenderMorphologyDataToNWB(obj, image_segmentation, imaging_plane)
    
#     # Recursively iterate through child collections
#     for child_collection in collection.children:
#         iterate_collections(child_collection, image_segmentation, imaging_plane)

# iterate_collections(root_collection, image_segmentation, imaging_plane)


def iterate_collections(collection, image_segmentation, imaging_plane):
    print("Collection:", collection.name)
    
    # Iterate through all objects in the collection
    for obj in collection.objects:
        print(f"  Object in iterator: {obj.name}")
        
        # Debugging prints for image_segmentation and imaging_plane
        print(f"image_segmentation: {image_segmentation}")
        print(f"imaging_plane: {imaging_plane}")
        
        # Call BlenderMorphologyDataToNWB
        if obj.type == 'MESH':
            BlenderMorphologyDataToNWB(obj, image_segmentation, imaging_plane)
    
    # Recursively iterate through child collections
    for child_collection in collection.children:
        iterate_collections(child_collection, image_segmentation, imaging_plane)

# Ensure image_segmentation and imaging_plane are properly defined
print(f"root_collection: {root_collection}")
print(f"image_segmentation: {image_segmentation}")

# Call the function
iterate_collections(root_collection, image_segmentation, imaging_plane)

#Use DANDI Archives convention:prepend sub-, insert _ses-
#nwbfile_name = 'sub-' + 'subject_id' + '_ses-' + 'identifier' + '.nwb'
print(nwbfile.processing['MorphologyData'].data_interfaces)

with NWBHDF5IO('SimpleAllenSegmentation.nwb', 'w') as io:
    io.write(nwbfile)