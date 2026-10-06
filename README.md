# imaging_science_tools
Tools for doing imaging science and slides that explain it.

Full videos on youtube:

ST1 - Matlab Figures: https://www.youtube.com/watch?v=ur0cTY8sQ0M

ST2 - Quick intro to git, tortoisegit, gitlab, github: https://youtu.be/eYHlW0sh8Is 
	- clone public or private repos, branches
	- local changes, diff, commit, revert
	- push to server, merge/pull request

ST3 - folders_data_codes: https://youtu.be/I48tuQhZ34A 
	- separate raw data,  code  and processed data
	- Projects 001
	- raw data - 2022_09_18
	- link to raw data and copy it with code
	- best practice: git clone, readme, get local or web data, run code

ST4 - NIfTI orientation: the right-hand rule phantom
	- NIfTI = data array plus affine; RAS+ world, R x A = S like thumb x index = middle finger
	- a hand in the right-hand rule pose (Hegdé hand template, NIH 3D 17237): thumb +X, index finger +Y, middle finger +Z
	- neurological versus radiological view: what flips, why both are used, read the R/L letters
	- folder: 4_nii_orientation
	- colab: https://colab.research.google.com/github/cest-sources/imaging_science_tools/blob/main/4_nii_orientation/right_hand_rule_phantom.ipynb

ST5 - NIfTI, resampling, and registration
	- DICOM series to one NIfTI: data array plus voxel-to-world matrix
	- resample MR_Gd onto a tilted 3-slice slab (2.5 x 2.5 x 7.5 mm, full in-plane field of view), and back
	- footprint average versus nibabel; register the high-resolution scan onto the low-resolution slab with ANTs and DIPY
	- folder: 5_Image_visualization_and_registration
	- colab: https://colab.research.google.com/github/cest-sources/imaging_science_tools/blob/main/5_Image_visualization_and_registration/nifti_and_registration.ipynb

TODO:

STx - minimal project management
	- getting things done - inbox, next step!!, two minutes rule







----


Used resources:

- Screencast done with OBS (https://obsproject.com/welcome) 
  - screen @ 2560 x 1440,  output @ 1920 x 1080
  - speex noise suppression -30dB, Kompressor  -17dB
- Video editig using http://fixounet.free.fr/avidemux/
- Slides done with Microsoft Powerpoint, recorded there. Then screencasted again with OBS to use its noise supression.
- Code with Mathworks Matlab 2016a
