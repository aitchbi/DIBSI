# DIBSI
Domain-Informed B-Spline Interpolation (DIBSI).

An implementation of the method proposed in: 

Behjat, H., Dogan, Z., Van De Ville, D. and Sörnmo, L. Domain-informed spline interpolation. IEEE Transactions on Signal Processing, 67(15), pp. 3909-3921, 2019. doi.org/10.1109/TSP.2019.2922154 | preprint: arxiv.org/abs/1810.07502

An example application of this interpolation technique is for resampling brain functional Magnetic Resonance Imaging (fMRI) data in such way that the underlying inhomogenous brain structure of the individual is acconted for in the data interpolation phase. The following illustration is a replication of Fig. 8 from the above paper:  
![Fig. 8 from Behjat et al. 2019](figs/fig8_paper.jpg?raw=true)

(a) A slice of fMRI data of a subject, including a close-up of an ROI. (b) The subject’s brain anatomy at the same neurological coordinate as in (a). (c) Gray matter, (d) white matter and (e) cerebrospinal fluid segmented probability maps of the ROI shown in (b). (f) Description of the inhomogeneous domain along the marked line within the ROI in (b), from top to bottom. (g) DIBSI basis, of order three, associated to the domain shown in (f). (h) B-spline interpolation (BSI) and DIBSI of the functional samples along the marked line within the ROI in (a)
