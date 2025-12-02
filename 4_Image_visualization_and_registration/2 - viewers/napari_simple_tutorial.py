import napari
import nibabel as nib

path_nii_OG = r'..\phantom_for_tutorials\dice_test_phantom.nii.gz'
image = nib.load(path_nii_OG).get_fdata()

# Basic display of a 3D image
viewer = napari.Viewer()
viewer.add_image(image, name="Image")

# Create a mask and add as another layer
mask = image.copy()
mask[:] = 0.
mask[image > 2.5] = 1
mask = mask.astype(int)
viewer.add_labels(mask, name="Mask")
napari.run()
