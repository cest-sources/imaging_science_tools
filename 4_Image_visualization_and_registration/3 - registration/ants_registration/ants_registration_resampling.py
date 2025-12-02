import numpy as np
import ants

def register_resample_figures(source_path, target_path, output_path):
    """ Coregister and resample a 4D/3D image (source) to match a target 4D/3D image using ANTsPy.    
    Parameters:
    ----------
    source_path : str
        Path to the source 4D image to be resampled.
    target_path : str
        Path to the target 3D image.
    output_path : str
        Path to save the resampled 4D image.
    """

    # -------------------------------------
    #  1. Ensure paths are strings
    # -------------------------------------
    source_path = str(source_path)
    target_path = str(target_path)
    output_path = str(output_path)

    # -------------------------------------
    #  2. Load images as ants objects
    # -------------------------------------
    img4d_ref = ants.image_read(target_path)   # reference
    img4d_mov = ants.image_read(source_path)   # to be resampled

    print("Ref shape:", img4d_ref.shape)
    print("Mov shape:", img4d_mov.shape)

    # -------------------------------------
    #  3. Choose 3D reference volume
    # -------------------------------------
    if len(img4d_ref.shape) == 4:
        ref3d = img4d_ref[:, :, :, 0]   # ANTsImage 3D
        ref3d = ants.from_numpy(ref3d.numpy(), spacing=img4d_ref.spacing[:3], origin=img4d_ref.origin[:3], direction=img4d_ref.direction[:3,:3])
    else:
        ref3d = img4d_ref

    # -------------------------------------
    # 4. Number of volumes to resample
    # -------------------------------------
    resampled_vols = []
    img4d_mov_np = img4d_mov.numpy()
    if len(img4d_mov.shape) == 4:
        T_mov = img4d_mov.shape[3]
    else:
        T_mov = 1
        img4d_mov_np = img4d_mov_np[..., np.newaxis]

        new_spacing = img4d_mov.spacing + (1.0,)
        new_origin = img4d_mov.origin + (0.0,)

        new_direction = np.eye(4)
        new_direction[0:3,0:3] = img4d_mov.direction

        img4d_mov = ants.from_numpy(img4d_mov_np, spacing=new_spacing, origin=new_origin, direction=new_direction)

    # -------------------------------------
    # 5. Resample/register each volume
    # -------------------------------------
    for t in range(T_mov):
        vol3d = img4d_mov[:, :, :, t]   # 3D ANTsImage
        vol3d = ants.from_numpy(vol3d.numpy(), spacing=img4d_mov.spacing[:3], origin=img4d_mov.origin[:3], direction=img4d_mov.direction[0:3,0:3])
        vol_res = ants.resample_image_to_target(vol3d, ref3d, interp_type='linear')
        resampled_vols.append(vol_res)

    # -------------------------------------
    # 6. 3D stack of registered volumes
    #       to 4D image
    # -------------------------------------
    np_vols = [v.numpy() for v in resampled_vols]   # each is (x,y,z)
    stacked = np.stack(np_vols, axis=3)              # result shape (x,y,z,t)
    print("Stacked numpy shape:", stacked.shape)

    # -------------------------------------
    # 7. Prepare metadata 
    # -------------------------------------

    # Spacing
    sx, sy, sz = ref3d.spacing
    st = img4d_mov.spacing[3] if len(img4d_mov.spacing) > 3 else 1.0
    spacing4 = (sx, sy, sz, st)

    # Origin
    def ensure_4d_origin(origin):
        origin = list(origin)
        if len(origin) == 3:
            origin.append(0.0)
        return tuple(origin)

    origin4 = ensure_4d_origin(ref3d.origin)

    # Direction
    dir3 = np.array(ref3d.direction).reshape((3, 3))
    dir4 = np.eye(4, dtype=float)
    dir4[:3, :3] = dir3
    direction4 = tuple(dir4.flatten().tolist())   # antspy expects a flat tuple

    # -------------------------------------
    # 8. Create 4D ANTsImage 
    # -------------------------------------
    stacked = stacked.astype(img4d_mov.numpy().dtype)

    if len(img4d_ref.spacing) == 3:
        spacing4 = img4d_ref.spacing + (1.0, )
        origin4 = img4d_ref.origin + (0.0,)
        direction4 = np.eye(4)
        direction4[0:3,0:3] = img4d_ref.direction
    else:
        spacing4 = img4d_ref.spacing
        origin4 = img4d_ref.origin
        direction4 = img4d_ref.direction

    img4d_resampled = ants.from_numpy(
        stacked,
        spacing=spacing4,
        origin=origin4,
        direction=direction4
    )

    # -------------------------------------
    # 9. Save image
    # -------------------------------------
    print(img4d_resampled)   # should show Dimensions: (x,y,z,t) and Components: 1
    ants.image_write(img4d_resampled, output_path)
    print(f"Saved {output_path}")

    print(f"Coregistered image saved to {output_path}")
