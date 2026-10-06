# NIfTI, resampling, and registration

DICOM series are converted to NIfTI. A NIfTI keeps a data array (intensities at voxel indices) and a voxel-to-world matrix (a 4×4 affine, index to millimetres). The notebook reslices `MR_Gd` onto a tilted slab, compares a footprint average with nibabel, and registers that high-resolution scan onto the slab with a known affine, with ANTs, and with DIPY.

The worked example is [`nifti_and_registration.ipynb`](nifti_and_registration.ipynb). Display is NiiVue.

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/cest-sources/imaging_science_tools/blob/main/5_Image_visualization_and_registration/nifti_and_registration.ipynb)

The first cell installs the packages with `!pip`.

## Here you will learn

- **From DICOM to NIfTI.** A DICOM series is many files, one per slice. Tools such as `dcm2niix` or SimpleITK turn it into one NIfTI: a data array plus an affine.
- **Reading a NIfTI.** Load `MR_Gd`, look at its shape, voxel size, and affine, and see how the affine places every voxel in world millimetres.
- **Building a new grid.** Define a tilted slab of 3 slices at 2.5 × 2.5 × 7.5 mm, rotated by 25°, that covers the full in-plane field of view, just like a clinical low-resolution acquisition.
- **Resampling.** Move the high-resolution image onto the slab and back. Compare a footprint average, which mimics thick slices, with nibabel's `resample_from_to` interpolation. Put both volumes on one 0.5 mm grid in the slab frame to compare them voxel by voxel.
- **Registration.** First with a known affine (a shift plus an 8° rotation, then its inverse), so you know the right answer. Then rigid registration with ANTs and with DIPY using mutual information.
- **Viewing.** Every result opens in NiiVue with buttons for radiologic or neurologic display, world millimetres or voxel indices, and a checkbox for each loaded volume. Slice edges are labelled R, L, A, P, S, and I (S is head, I is foot). For the meaning of the radiologic and neurologic views, see [`4_nii_orientation`](../4_nii_orientation).
- **Where to go next.** Other viewers (NiiVue VS Code extension, 3D Slicer, MITK, FSLeyes, SPM, napari) and other resampling tools (nilearn, ANTs, SimpleITK).

NiiVue panels stay blank on GitHub until you run the cells.
