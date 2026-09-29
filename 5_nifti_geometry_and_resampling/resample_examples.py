"""Small NIfTI resampling examples for the geometry tutorial.

Builds a synthetic ball, then:

- resamples a coarse axial 3D grid onto a finer tilted 3D grid, and back
- resamples that volume onto one tilted 2D slice
- paints the slice back into a 3D grid (a slab, not the original ball)

Run from anywhere:

    python resample_examples.py
"""

from pathlib import Path

import nibabel as nib
import numpy as np
from nibabel.affines import apply_affine
from nibabel.processing import resample_from_to
from scipy.ndimage import map_coordinates

HERE = Path(__file__).resolve().parent
OUT = HERE / "output"


def rotation_about_x(degrees):
    """3x3 rotation. Columns are the new directions of the voxel axes."""
    theta = np.deg2rad(degrees)
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[1.0, 0.0, 0.0],
                      [0.0, c, -s],
                      [0.0, s, c]])


def affine_centered(shape, spacing, rotation=None):
    """4x4 voxel-to-world matrix. The center of the grid sits at world (0, 0, 0).

    spacing: millimetres along voxel axes i, j, k
    rotation: 3x3, columns are the world directions of i, j, k
    """
    if rotation is None:
        rotation = np.eye(3)
    shape = np.asarray(shape, dtype=float)
    spacing = np.asarray(spacing, dtype=float)
    center = (shape - 1) / 2.0
    affine = np.eye(4)
    affine[:3, :3] = rotation @ np.diag(spacing)
    affine[:3, 3] = -(rotation @ (spacing * center))
    return affine


def voxel_indices(shape):
    """Index of every voxel, shape (*shape, 3), in i, j, k order."""
    axes = [np.arange(n) for n in shape]
    grids = np.meshgrid(*axes, indexing="ij")
    return np.stack(grids, axis=-1)


def ball_on_grid(shape, affine, radius_mm=20.0):
    """Ball of value 1, value 2 where z > 0, so a tilt is visible."""
    xyz = apply_affine(affine, voxel_indices(shape))
    inside = np.linalg.norm(xyz, axis=-1) <= radius_mm
    data = inside.astype(np.float32)
    data[inside & (xyz[..., 2] > 0)] = 2.0
    return data


def resample_to_grid(source, source_affine, target_shape, target_affine, order=1):
    """Sample a 3D array onto another grid.

    For each target voxel: index -> world mm -> source index, then interpolate.
    order=1 is linear. order=0 is nearest neighbour (labels).
    Voxels that fall outside the source are 0.
    """
    world = apply_affine(target_affine, voxel_indices(target_shape))
    source_ijk = apply_affine(np.linalg.inv(source_affine), world)
    coords = np.moveaxis(source_ijk, -1, 0)
    sampled = map_coordinates(source, coords, order=order, mode="constant", cval=0.0)
    return sampled.astype(np.float32)


def place_slice_in_volume(slice_2d, slice_affine, volume_shape, volume_affine, order=1):
    """Paint one oriented slice into a 3D grid. Voxels outside the slab stay 0.

    slice_2d is (nx, ny). slice_affine is the 4x4 of the slice stored as
    (nx, ny, 1): columns 0 and 1 are the in-plane steps, column 2 is the
    normal, and the length of column 2 is the slice thickness in mm.
    """
    world = apply_affine(volume_affine, voxel_indices(volume_shape))
    on_slice = apply_affine(np.linalg.inv(slice_affine), world)
    coords = np.moveaxis(on_slice[..., :2], -1, 0)
    painted = map_coordinates(slice_2d, coords, order=order, mode="constant", cval=0.0)
    thickness = float(np.linalg.norm(slice_affine[:3, 2]))
    distance = np.abs(on_slice[..., 2]) * thickness
    painted[distance > (thickness / 2.0)] = 0.0
    return painted.astype(np.float32)


def save_nifti(data, affine, name):
    path = OUT / name
    nib.save(nib.Nifti1Image(np.asarray(data, np.float32), affine), path)
    print(f"  wrote {path.name}  shape={tuple(data.shape)}")


def main():
    OUT.mkdir(exist_ok=True)

    # Coarse axial volume, 2 mm voxels. Fine volume, 1.5 mm, tilted 25 degrees.
    source_shape = (32, 32, 32)
    source_affine = affine_centered(source_shape, spacing=(2.0, 2.0, 2.0))
    tilt = rotation_about_x(25)
    target_shape = (48, 48, 48)
    target_affine = affine_centered(target_shape, spacing=(1.5, 1.5, 1.5), rotation=tilt)

    source = ball_on_grid(source_shape, source_affine)
    print("source spacing mm", np.linalg.norm(source_affine[:3, :3], axis=0))
    print("target spacing mm", np.linalg.norm(target_affine[:3, :3], axis=0))

    on_tilted = resample_to_grid(source, source_affine, target_shape, target_affine)
    back = resample_to_grid(on_tilted, target_affine, source_shape, source_affine)

    # Same 3D mapping, via nibabel, to show the function matches the library.
    source_img = nib.Nifti1Image(source, source_affine)
    nib_tilted = resample_from_to(source_img, (target_shape, target_affine), order=1)
    nib_data = np.asanyarray(nib_tilted.dataobj)
    max_diff = float(np.max(np.abs(on_tilted - nib_data)))
    roundtrip = float(np.mean(np.abs(back - source)))
    print(f"3D -> tilted 3D  shape {on_tilted.shape}  max diff vs nibabel {max_diff:.2e}")
    print(f"tilted 3D -> original grid  shape {back.shape}  mean abs diff {roundtrip:.4f}")

    # One tilted slice: 1 mm in plane, 4 mm thick, same 25 degree tilt.
    plane_shape = (80, 80, 1)
    plane_affine = affine_centered(plane_shape, spacing=(1.0, 1.0, 4.0), rotation=tilt)
    plane = resample_to_grid(source, source_affine, plane_shape, plane_affine)
    slice_2d = plane[:, :, 0]
    slab = place_slice_in_volume(slice_2d, plane_affine, source_shape, source_affine)

    thickness = float(np.linalg.norm(plane_affine[:3, 2]))
    world = apply_affine(source_affine, voxel_indices(source_shape))
    normal = plane_affine[:3, 2] / thickness
    distance = np.abs(world @ normal)
    filled = slab > 0
    print(f"3D -> 2D  slice shape {slice_2d.shape}  nonzero {int(np.count_nonzero(slice_2d))}")
    print(f"2D -> 3D  slab nonzero {int(filled.sum())} / {slab.size}")
    print(f"  thickness {thickness:.1f} mm  max distance of filled voxels {distance[filled].max():.2f} mm")

    center = tuple(n // 2 for n in source_shape)
    center_mm = apply_affine(source_affine, np.array(center, dtype=float))
    print(f"  voxel {center} at world mm {np.round(center_mm, 2)}  source {source[center]:.1f}  slab {slab[center]:.1f}")

    save_nifti(source, source_affine, "source.nii.gz")
    save_nifti(on_tilted, target_affine, "tilted.nii.gz")
    save_nifti(back, source_affine, "back.nii.gz")
    save_nifti(plane, plane_affine, "plane.nii.gz")
    save_nifti(slab, source_affine, "slice_in_volume.nii.gz")


if __name__ == "__main__":
    main()
