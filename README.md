# Cell tracking for live-cell microscopy using an activity-prioritized assignment strategy

This repository contains the code to the publication:
“Cell tracking for live-cell microscopy using an activity-prioritized assignment strategy” Ruzaeva K. et al(in press, accepted at IEEE IPAS 2022)

and the post-processed *C.glutamicum* dataset from

[Noga Mosheiff et al., “Inheritance of cell-cycle duration in the presence of periodic forcing” Phys. Rev. X, vol. 8, pp. 021035, May 2018.](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.8.021035)


## Data
```
ACTIVITY-CELL-TRACKING
└───data  (contains preprocessed (segmented) datasets with different frame rates)
    └───N  (contains the dataset with N min frame rate)
        └───seg  (contains segmentation masks)
        └───im  (contains original data)
        └───res  (result will be stored here, the output is organized according to CTC naming convention)
    ...    
└───track.py (contains our tracking code)
└───utils.pt (contains some helper functions, i.e. Activity map and Gaussisan map generation)
```

