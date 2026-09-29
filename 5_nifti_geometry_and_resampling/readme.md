# NIfTI geometry and resampling

This is the next tutorial after the tool list in [4. Image visualization and registration](../4_Image_visualization_and_registration/readme.md). Those pages name the programs (dcm2niix, napari, ANTs). This one walks the short path those programs share: get a NIfTI, see how its geometry is stored, and resample it onto another grid.

The examples use a synthetic ball, so they run without a DICOM download. The same functions work on a real file. A runnable copy is [`resample_examples.py`](resample_examples.py).

```
pip install numpy scipy nibabel
```

Pictures need `matplotlib` as well. Napari is optional.

## Plan

1. **DICOM to NIfTI.** One series becomes one `.nii.gz`.
2. **View the NIfTI.** The array on screen is voxel indices. Millimetres come from the affine.
3. **Data array and voxel-to-world matrix.** Intensities live in the array. Position, voxel size, and tilt live in a 4×4 affine.
4. **One resample function.** For each voxel of the grid you want: index → millimetres → index in the image you have, then interpolate.
5. **3D onto a tilted 3D, and back.** Same anatomy, different voxel size, rotated grid. Swap the two grids to come back.
6. **3D onto a 2D plane, and that plane back into 3D.** A slice is a grid with one voxel along its normal. Putting it back fills only that slab.

Left out on purpose: searching for an unknown alignment (that is [registration](../4_Image_visualization_and_registration/3%20-%20registration/readme.md)), and 4D series.

## 1. DICOM to NIfTI

A DICOM folder is many files. A NIfTI is one file with the pixel stack and its geometry. Two routes already noted in this repo:

- [dcm2niix](https://github.com/rordenlab/dcm2niix)
- SimpleITK, in [`sitk_dcm2nii.py`](../4_Image_visualization_and_registration/1%20-%20dcm2nii/sitk_dcm2nii/sitk_dcm2nii.py)

The short form, one series, no file dialog:

```python
import SimpleITK as sitk

def dicom_series_to_nifti(dicom_dir, nifti_path):
    reader = sitk.ImageSeriesReader()
    files = reader.GetGDCMSeriesFileNames(dicom_dir)
    reader.SetFileNames(files)
    sitk.WriteImage(reader.Execute(), nifti_path)
```

SimpleITK stores origin, spacing, and direction. Written as NIfTI, that geometry is nibabel's affine. Load that file with nibabel for the rest of this tutorial. The array inside SimpleITK is ordered z, y, x; the NIfTI array follows the header.

## 2. Viewing

```python
import nibabel as nib
import matplotlib.pyplot as plt

img = nib.load("volume.nii.gz")
data = img.get_fdata()   # intensities, indexed by voxel
affine = img.affine      # 4x4, voxel index -> world millimetres

mid = [n // 2 for n in data.shape[:3]]
fig, axes = plt.subplots(1, 3, figsize=(9, 3))
axes[0].imshow(data[mid[0], :, :].T, cmap="gray", origin="lower")
axes[1].imshow(data[:, mid[1], :].T, cmap="gray", origin="lower")
axes[2].imshow(data[:, :, mid[2]].T, cmap="gray", origin="lower")
for ax in axes:
    ax.set_axis_off()
plt.show()
```

`data[i, j, k]` selects a voxel by index. Millimetres come from the affine, applied to that index. A slice through the middle index follows one voxel axis; whether that axis is axial in the scanner is what the affine says. Napari and the other programs in the [viewers list](../4_Image_visualization_and_registration/2%20-%20viewers/readme.md) draw this array. Give the viewer the affine too, so the tilt is kept.

## 3. Data array and voxel-to-world matrix

A NIfTI is two arrays.

**Data array.** `data` with shape `(ni, nj, nk)`. `data[i, j, k]` is the intensity at that voxel. Voxel size and orientation live in the affine, described next.

**Affine.** A 4×4 matrix from voxel index to a point in world millimetres (usually scanner RAS):

```
[x]   [a11 a12 a13 tx] [i]
[y] = [a21 a22 a23 ty] [j]
[z]   [a31 a32 a33 tz] [k]
[1]   [ 0   0   0   1] [1]
```

In code: `xyz = affine @ [i, j, k, 1]`.

| column | meaning |
| --- | --- |
| 0 | step, in mm, from voxel `(i, j, k)` to `(i+1, j, k)` |
| 1 | step to `(i, j+1, k)` |
| 2 | step to `(i, j, k+1)` |
| 3 | world position of voxel `(0, 0, 0)` |

Voxel size along an axis is the length of that column. The axis direction is the column divided by that length. Columns that are not parallel to x, y, z mean the grid is tilted.

```python
import numpy as np

spacing = np.linalg.norm(affine[:3, :3], axis=0)
direction = affine[:3, :3] / spacing
origin = affine[:3, 3]
```

`img.header.get_zooms()` stores the spacing. The tilt is only in the affine.

`img.affine` is the matrix nibabel will use: the sform if the header has one, otherwise the qform. Index `0` is the center of the first voxel.

For the examples we build the affine ourselves, with the center of the grid at the world origin (for an even shape that point lies between voxels):

```python
def rotation_about_x(degrees):
    theta = np.deg2rad(degrees)
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[1.0, 0.0, 0.0],
                      [0.0, c, -s],
                      [0.0, s, c]])

def affine_centered(shape, spacing, rotation=None):
    """4x4 affine. The center of the grid sits at world (0, 0, 0)."""
    if rotation is None:
        rotation = np.eye(3)
    shape = np.asarray(shape, dtype=float)
    spacing = np.asarray(spacing, dtype=float)
    center = (shape - 1) / 2.0
    affine = np.eye(4)
    affine[:3, :3] = rotation @ np.diag(spacing)
    affine[:3, 3] = -(rotation @ (spacing * center))
    return affine
```

`spacing` is millimetres per voxel along i, j, k. Columns of `rotation` are the world directions of those axes. `rotation_about_x(25)` is the tilt used below. On a real file, keep the affine already stored in `img.affine`.

## 4. Resampling

Resampling writes a new data array on a new grid. The anatomy can be the same. The grid can change shape, voxel size, and tilt.

For every index on the **target** grid:

1. Index → world millimetres, with the target affine.
2. World millimetres → index in the **source** array, with the inverse source affine.
3. Read the source array there. The index is usually fractional, so interpolate.

```python
from nibabel.affines import apply_affine
from scipy.ndimage import map_coordinates

def voxel_indices(shape):
    axes = [np.arange(n) for n in shape]
    grids = np.meshgrid(*axes, indexing="ij")
    return np.stack(grids, axis=-1)

def resample_to_grid(source, source_affine, target_shape, target_affine, order=1):
    """Sample `source` on the grid given by `target_shape` and `target_affine`.

    order=1 is linear interpolation. order=0 is nearest neighbour (use that for labels).
    """
    world = apply_affine(target_affine, voxel_indices(target_shape))
    source_ijk = apply_affine(np.linalg.inv(source_affine), world)
    coords = np.moveaxis(source_ijk, -1, 0)
    sampled = map_coordinates(source, coords, order=order, mode="constant", cval=0.0)
    return sampled.astype(np.float32)
```

Samples that miss the source volume are 0. These examples are about 50 voxels on a side, so the index grid fits in memory. For a full scan, do the same math one slice at a time, or call nibabel, which does this mapping:

```python
from nibabel.processing import resample_from_to

out = resample_from_to(source_img, (target_shape, target_affine), order=1)
```

## 5. 3D onto a tilted 3D, and back

Source: 32×32×32, 2 mm voxels, axes parallel to the world. Target: 48×48×48, 1.5 mm voxels, tilted 25° about x. The phantom is a ball of radius 20 mm (value 1, and 2 in the upper half) so the tilt shows up.

```python
def ball_on_grid(shape, affine, radius_mm=20.0):
    xyz = apply_affine(affine, voxel_indices(shape))
    inside = np.linalg.norm(xyz, axis=-1) <= radius_mm
    data = inside.astype(np.float32)
    data[inside & (xyz[..., 2] > 0)] = 2.0
    return data

source_shape = (32, 32, 32)
source_affine = affine_centered(source_shape, spacing=(2.0, 2.0, 2.0))
tilt = rotation_about_x(25)
target_shape = (48, 48, 48)
target_affine = affine_centered(target_shape, spacing=(1.5, 1.5, 1.5), rotation=tilt)

source = ball_on_grid(source_shape, source_affine)
on_tilted = resample_to_grid(source, source_affine, target_shape, target_affine)
back = resample_to_grid(on_tilted, target_affine, source_shape, source_affine)
```

`on_tilted` is the finer, tilted volume. `back` is that volume sampled onto the original grid. Coming back is the same function with the two grids swapped. Linear interpolation blurs edges, so a round trip changes the values slightly. The script prints the mean absolute difference, and checks `on_tilted` against `resample_from_to`.

## 6. 3D to 2D, and 2D to 3D

Store one slice as a NIfTI of shape `(nx, ny, 1)`.

- Columns 0 and 1 of the affine place the pixels, in millimetres.
- Column 2 is the slice normal. Its length is the slice thickness in millimetres.
- The only index along that axis is `k = 0`, so the pixels lie on the plane.

**3D → 2D.** Same function. The target shape ends in 1. The affine carries the in-plane size and the tilt.

```python
plane_shape = (80, 80, 1)
plane_affine = affine_centered(plane_shape, spacing=(1.0, 1.0, 4.0), rotation=tilt)
plane = resample_to_grid(source, source_affine, plane_shape, plane_affine)
slice_2d = plane[:, :, 0]
```

In-plane pixels are 1 mm. The plane is tilted 25°. The thickness stored on column 2 is 4 mm.

**2D → 3D.** Painting the slice into a volume fills the slab around the plane. Every other voxel stays 0, because the slice only has intensities on that plane.

Use `place_slice_in_volume` for this direction. On an axis of length 1, linear interpolation has a single sample, so `resample_to_grid` would write 0 almost everywhere. Interpolate inside the plane, and keep a voxel when it lies within half the thickness of the plane:

```python
def place_slice_in_volume(slice_2d, slice_affine, volume_shape, volume_affine, order=1):
    """Paint one oriented slice into a 3D grid. Voxels outside the slab stay 0."""
    world = apply_affine(volume_affine, voxel_indices(volume_shape))
    on_slice = apply_affine(np.linalg.inv(slice_affine), world)
    coords = np.moveaxis(on_slice[..., :2], -1, 0)
    painted = map_coordinates(slice_2d, coords, order=order, mode="constant", cval=0.0)
    thickness = float(np.linalg.norm(slice_affine[:3, 2]))
    distance = np.abs(on_slice[..., 2]) * thickness
    painted[distance > (thickness / 2.0)] = 0.0
    return painted.astype(np.float32)

slab = place_slice_in_volume(slice_2d, plane_affine, source_shape, source_affine)
```

`slab` has the source shape. Nonzero voxels sit in a band about 4 mm thick around the tilted plane.

## Run it

```
python resample_examples.py
```

The script writes NIfTIs under `output/` (that folder is local, not part of the repo):

| file | what it is |
| --- | --- |
| `source.nii.gz` | coarse axial ball, 2 mm |
| `tilted.nii.gz` | same ball, 1.5 mm, tilted 25° |
| `back.nii.gz` | tilted volume sampled back onto the coarse grid |
| `plane.nii.gz` | one tilted slice, 1 mm in plane, 4 mm thick |
| `slice_in_volume.nii.gz` | that slice painted into the coarse grid |

Open `source.nii.gz` and `tilted.nii.gz` in napari or MITK. The ball occupies the same place in the world. The voxel grids that sample it differ.

## Still to add

- A real DICOM series in place of the ball. Keep the affine already stored in the NIfTI.
- Labels resampled with `order=0`.
- A napari view that shows the axial grid and the tilted grid together.
