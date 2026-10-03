# Right-hand rule phantom

A NIfTI of a hand posed as the usual drawing of the right-hand rule, with an interactive [NiiVue](https://github.com/niivue/niivue) viewer in [`right_hand_rule_phantom.ipynb`](right_hand_rule_phantom.ipynb).

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/cest-sources/imaging_science_tools/blob/main/6_right_hand_rule/right_hand_rule_phantom.ipynb)

The source is the healthy-adult hand template of Hegdé et al., the NIfTI behind [NIH 3D entry 17237](https://3d.nih.gov/entries/17237). The template file is `Hegde_etal_Healthy_adult_hand_template.nii` from [HegdeUSA/Hand_template](https://github.com/HegdeUSA/Hand_template) (MIT). NIH 3D hosts the surface made from that image; the notebook downloads the image.

The open hand is resampled to 1 mm and the digits are hinged:

| world axis | niivue label | digit |
| --- | --- | --- |
| +X | R (right) | thumb |
| +Y | A (anterior) | index finger |
| +Z | S (superior) | middle finger |

The palm, the index finger, and the thumb lie in the XY plane. +Z points out of the palm, the front of the hand. The middle finger follows +Z. The ring and little finger curl further toward the palm, and the little finger aims in toward the middle of the palm. Thumb × index = middle finger, the same right-handed triad as R × A = S.

[`right_hand_rule_phantom.nii.gz`](right_hand_rule_phantom.nii.gz) is that posed volume. The notebook rebuilds it and opens it in NiiVue: scroll a slice, drag the 3D panel. The panels stay blank on GitHub until the cells run.

```
pip install nibabel scipy matplotlib ipywidgets ipyniivue
jupyter lab right_hand_rule_phantom.ipynb
```
