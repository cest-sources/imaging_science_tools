# NIfTI, resampling, and registration

DICOM series are converted to NIfTI. A NIfTI keeps a data array (intensities at voxel indices) and a voxel-to-world matrix (a 4×4 affine, index to millimetres). The notebook reslices a volume onto a tilted slab, compares a footprint average with nibabel, and registers a phantom with a known affine, with ANTs, and with DIPY.

The worked example is [`nifti_and_registration.ipynb`](nifti_and_registration.ipynb). Display is NiiVue. The phantom used for registration is [`phantom_for_tutorials/dice_test_phantom.nii.gz`](phantom_for_tutorials/dice_test_phantom.nii.gz).

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/cest-sources/imaging_science_tools/blob/cursor/merge-nifti-registration-7d5a/4_Image_visualization_and_registration/nifti_and_registration.ipynb)

[Open in Colab](https://colab.research.google.com/github/cest-sources/imaging_science_tools/blob/cursor/merge-nifti-registration-7d5a/4_Image_visualization_and_registration/nifti_and_registration.ipynb). The first cell installs the packages with `!pip`.

NiiVue panels stay blank on GitHub until you run the cells. The same buttons are on every viewer: radiologic or neurologic, world millimetres or voxel indices, and, when two images are loaded, a checkbox for each file. Slice edges are labeled R, L, A, P, S, and I (S is head, I is foot).
