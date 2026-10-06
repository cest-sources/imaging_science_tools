# NIfTI orientation: the right-hand rule phantom

A NIfTI of a hand posed as the usual drawing of the right-hand rule, with an interactive [NiiVue](https://github.com/niivue/niivue) viewer in [`right_hand_rule_phantom.ipynb`](right_hand_rule_phantom.ipynb).

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/cest-sources/imaging_science_tools/blob/main/4_nii_orientation/right_hand_rule_phantom.ipynb)

## Here you will learn

- **How a NIfTI knows where things are.** A NIfTI file holds a data array and a 4×4 affine that maps voxel indices to world millimetres. The world is RAS+: +X is the patient's Right, +Y Anterior, +Z Superior.
- **Why the axes are right-handed.** R × A = S, the same as thumb × index finger = middle finger of your right hand. The phantom is a hand in exactly that pose, so you can check any viewer against your own hand.
- **Why a hand is a good test object.** A right hand is chiral: its mirror image is a left hand. A brain is nearly symmetric, so a left-right flip can go unnoticed. With the hand, a flip is obvious.
- **The neurological and the radiological view, and why both exist.** This is the main point of the notebook (see below).
- **How to read a viewer safely.** Never trust which side of the screen something is on. Read the R/L letters at the slice edges.

## Neurological and radiological view

The file is the same in both views. The data and the affine do not change, and neither does the patient's right. Only the drawing on the screen is mirrored.

| | neurological | radiological |
| --- | --- | --- |
| patient's right appears | on the screen's right ("left is left") | on the screen's left |
| axial slice seen from | above (the head) | below (the feet) |
| coronal slice seen from | behind | the front |
| used by | neuroscience, brain atlases, SPM | clinical radiology, PACS, CT and MRI reports |

Why two conventions?

- **Radiological:** the doctor faces the patient, like shaking hands, so the patient's right is on the doctor's left. Axial CT and MRI slices are viewed from the foot of the bed, the same way the patient lies on the scanner table. Clinical images are always read this way.
- **Neurological:** the anatomist's view, as if you stand behind the patient and look down at the head. Left is left, which matches how brain atlases and anatomical drawings are made.

Why it matters: research software and clinical workstations often use different defaults. A wrong-side reading leads to wrong-side reports, and in the worst case to wrong-side surgery or a wrong lateralisation in a study. Sagittal slices look the same in both views.

With the phantom, the neurological axial view shows the palm, with the middle finger pointing at you. The radiological axial view shows the back of the hand, with the middle finger pointing away. In both, the thumb is next to the letter R, and it is still a right hand.

Do not confuse the display convention with "radiological storage order", which is the order of the voxels in the file. A correct affine makes the storage order irrelevant for a proper viewer.

## The phantom

The source is the healthy-adult hand template of Hegdé et al., the NIfTI behind [NIH 3D entry 17237](https://3d.nih.gov/entries/17237). The template file is `Hegde_etal_Healthy_adult_hand_template.nii` from [HegdeUSA/Hand_template](https://github.com/HegdeUSA/Hand_template) (MIT). NIH 3D hosts the surface made from that image; the notebook downloads the image.

The open hand is resampled to 1 mm and the digits are hinged:

| world axis | niivue label | digit |
| --- | --- | --- |
| +X | R (right) | thumb |
| +Y | A (anterior) | index finger |
| +Z | S (superior) | middle finger |

The palm, the index finger, and the thumb lie in the XY plane. +Z points out of the palm, the front of the hand. The middle finger follows +Z. The ring and little finger fold back and lie flat on the palm, and the little finger angles in toward the middle of the palm. Thumb × index = middle finger, the same right-handed triad as R × A = S.

[`right_hand_rule_phantom.nii.gz`](right_hand_rule_phantom.nii.gz) is that posed volume. The notebook rebuilds it, shows projections from three sides, compares the neurological and radiological views side by side, and opens the volume in NiiVue. Press the view button to switch between neurological and radiological and watch the axial and coronal panels flip. The panels stay blank on GitHub until the cells run.

```
pip install nibabel scipy matplotlib ipywidgets ipyniivue
jupyter lab right_hand_rule_phantom.ipynb
```
