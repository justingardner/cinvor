
#simCEintp

Toolbox to simulate channel encoding reconstruction of 8 stimulus orientations from 50 response instances of 2 voxels. This tests the case when stimulus orientations do not match the 
hypothetical channels orientation preferences. Channel responses are interpolated for 360 channels with orientation preferences 
ranging from 1:1:360 degs.

The code Interpolates channel responses on a space of 360 channel orientation preferences. This only works with orientation expressed on 360 deg space not on 180 deg space.

The code displays the clusters of voxels responses by orientations for the training set, the test set , the reconstructed orientations plotted against the true orientations and the test set average channel response normalized to the displayed orientation.

just run :

```matlab
  >> run slsimCEintp.m
```
