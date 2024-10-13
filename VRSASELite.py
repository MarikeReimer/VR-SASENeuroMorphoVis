
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

def BlenderMorphologyDataToNWB(obj, plane_segmentation):
    print('Object in NWB:', obj.name)
    
    # Dummy data for plane segmentation    
    pixel_mask = np.array([[0, 0, 0]] * 1)
    
    # Extract mesh attributes
    mesh_attributes = find_mesh_attributes(obj)
    center_of_mass = mesh_attributes[0]
    volume = mesh_attributes[1]
    surface_area = mesh_attributes[2]

    # Add columns to ROI if not already added (check to avoid duplication)
    if 'center_of_mass' not in plane_segmentation.colnames:
        plane_segmentation.add_column('center_of_mass', 'center_of_mass')
        plane_segmentation.add_column('volume', 'volume')
        plane_segmentation.add_column('surface_area', 'surface_area')        

    # Add ROI for this object
    plane_segmentation.add_roi(
        pixel_mask=pixel_mask,
        center_of_mass=center_of_mass,
        volume=volume,
        surface_area=surface_area
    )
    
    return plane_segmentation

def BlenderMorphologyDataToNWB(obj, plane_segmentation):
    print(f'Processing Object in NWB: {obj.name}')
    
    # Dummy data for plane segmentation    
    pixel_mask = np.array([[0, 0, 0]] * 1)
    
    # Extract mesh attributes
    mesh_attributes = find_mesh_attributes(obj)
    center_of_mass = mesh_attributes[0]
    volume = mesh_attributes[1]
    surface_area = mesh_attributes[2]

    # Add columns to ROI if not already added (check to avoid duplication)
    if 'center_of_mass' not in plane_segmentation.colnames:
        plane_segmentation.add_column('center_of_mass', 'center_of_mass')
        plane_segmentation.add_column('volume', 'volume')
        plane_segmentation.add_column('surface_area', 'surface_area')

    # Add ROI for this object
    plane_segmentation.add_roi(
        pixel_mask=pixel_mask,
        center_of_mass=center_of_mass,
        volume=volume,
        surface_area=surface_area
    )
    
    return plane_segmentation

#OPTION 1: Creates plane_segmentations for only mesh objects (no curves).

def iterate_collections(collection, image_segmentation, imaging_plane):
    print(f"Iterating Collection: {collection.name}")
    
    # Check if the collection contains objects and avoid creating a segmentation for the collection itself
    for obj in collection.objects:
        print(f"Object in iterator: {obj.name}")
        
        # Process only MESH objects
        if obj.type == 'MESH':
            segmentation_name = obj.name + '_plane_segmentation'
            
            # Check if this plane segmentation already exists or create a new one
            if segmentation_name not in [ps.name for ps in image_segmentation.plane_segmentations.values()]:
                plane_segmentation = image_segmentation.create_plane_segmentation(
                    name=segmentation_name,
                    description='Output from segmenting a mesh in Blender',
                    imaging_plane=imaging_plane
                )
                print(f"Created new plane segmentation: {segmentation_name}")
            else:
                plane_segmentation = image_segmentation.plane_segmentations[segmentation_name]
                print(f"Using existing plane segmentation: {segmentation_name}")
            
            # Process the MESH object
            BlenderMorphologyDataToNWB(obj, plane_segmentation)
        else:
            print(f"Skipping non-MESH object: {obj.name}")
    
    # Recursively iterate through child collections, avoiding segmentation for the collection itself
    for child_collection in collection.children:
        if isinstance(child_collection, bpy.types.Collection):
            iterate_collections(child_collection, image_segmentation, imaging_plane)
        else:
            print(f"Skipping child collection that is not valid: {child_collection.name}")


# Debugging root_collection and image_segmentation
print(f"image_segmentation: {image_segmentation}")

# Call the iteration function, ensuring CollectionToNWB is not treated as a plane segmentation
iterate_collections(root_collection, image_segmentation, imaging_plane)


#OPTION 2: Creates plane_segmentations for non-mesh objects like curves and skips missing data values.

def iterate_collections(collection, image_segmentation, imaging_plane):
    print("Iterating Collection:", collection.name)
    
    # Iterate through all objects in the collection
    for obj in collection.objects:
        print(f"Object in iterator: {obj.name}")
        
        # Process only MESH objects
        if obj.type == 'MESH':
            # Create a single PlaneSegmentation for the collection
            segmentation_name = collection.name + '_plane_segmentation'
            plane_segmentation = image_segmentation.create_plane_segmentation(
            name=segmentation_name,
            description='Output from segmenting a mesh in Blender',
            imaging_plane=imaging_plane
            )   
        
            print(f"Processing MESH object: {obj.name}")
            BlenderMorphologyDataToNWB(obj, plane_segmentation)
        else:
            print(f"Skipping non-MESH object: {obj.name}")
    
    # Recursively iterate through child collections
    for child_collection in collection.children:
        iterate_collections(child_collection, image_segmentation, imaging_plane)

#BEST PRACTICE: Use DANDI Archives convention:prepend sub-, insert _ses-
#nwbfile_name = 'sub-' + 'subject_id' + '_ses-' + 'identifier' + '.nwb'
print(nwbfile.processing['MorphologyData'].data_interfaces)

with NWBHDF5IO('SimpleAllenSegmentation.nwb', 'w') as io:
    io.write(nwbfile)