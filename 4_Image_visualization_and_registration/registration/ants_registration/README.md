# **ANTs registration**

The registration function presented here uses ANTsPy and is importable to python.

The documentation for ANTsPy can be found under this link: [ANTsPyx - https://github.com/ANTsX/ANTsPy](https://github.com/ANTsX/ANTsPy)

ANTsPy is a python wrapper for a well-established C++ biomedical image processing framework ANTs.

In this folder, you can find:
- ants_registration_resampling.py: a registration/reslicing function that works for 3D and 4d data, uses ANTsPy and can be easily imported into your main script.
- tutorial/ants_registration_resampling_example.ipynb: a jupyter notebook where you can test this function.
- tutorial/dice_test_phantom.nii.gz: nifti phantom used in the tutorial (can be found in MR-zero adddon eventually?).
