import bpy

# Get the active object in the scene (assuming it is a curve)
obj = bpy.context.object

# Ensure the object is a curve
if obj and obj.type == 'CURVE':
    curve_data = obj.data
    
    # General curve properties
    print(f"Curve Name: {obj.name}")
    print(f"Dimensions: {'2D' if curve_data.dimensions == '2D' else '3D'}")
    print(f"Resolution: {curve_data.resolution_u}")
    
    # Iterate through each spline in the curve
    for spline_index, spline in enumerate(curve_data.splines):
        print(f"\n--- Spline {spline_index} ---")
        print(f"Spline Type: {spline.type}")
        
        if spline.type == 'BEZIER':
            print(f"Number of control points: {len(spline.bezier_points)}")
            
            # Iterate through control points (Bezier Points)
            for i, point in enumerate(spline.bezier_points):
                print(f"\nControl Point {i}:")
                
                # Control point position
                print(f"  Coordinates: {point.co}")
                
                # Handle positions
                print(f"  Left Handle: {point.handle_left}")
                print(f"  Right Handle: {point.handle_right}")
                
                # Handle types (free, aligned, vector, automatic)
                print(f"  Left Handle Type: {point.handle_left_type}")
                print(f"  Right Handle Type: {point.handle_right_type}")
                
                # Tilt and weight
                print(f"  Tilt: {point.tilt}")
                print(f"  Weight: {point.weight_softbody}")
                
        elif spline.type == 'NURBS':
            print(f"Number of control points: {len(spline.points)} (NURBS)")
            # For NURBS splines, handle the point data differently.
            for i, point in enumerate(spline.points):
                print(f"\nNURBS Control Point {i}:")
                print(f"  Coordinates: {point.co}")
                print(f"  Weight: {point.weight}")
        
        elif spline.type == 'POLY':
            print(f"Number of control points: {len(spline.points)} (Poly)")
            # For Poly splines, only coordinates are relevant.
            for i, point in enumerate(spline.points):
                print(f"\nPoly Control Point {i}:")
                print(f"  Coordinates: {point.co}")
                
        else:
            print(f"Spline {spline_index} is of unsupported type: {spline.type}. Skipping.")
    
    # If the curve has a bevel, print bevel details
    if curve_data.bevel_depth > 0 or curve_data.extrude > 0:
        print("\nBevel and Extrude Settings:")
        print(f"  Bevel Depth: {curve_data.bevel_depth}")
        print(f"  Extrude: {curve_data.extrude}")
        
else:
    print("Please select a curve object.")
