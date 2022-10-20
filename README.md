This repository contains the code to the publication:

and the post-processed C.glutamicum dataset from


##Data

ACTIVITY-CELL-TRACKING
└───data  (contains preprocessed (segmented) datasets with different frame rates)
    └───N  (contains the dataset with N min frame rate)
        └───seg  (contains segmentation masks)
        └───im  (contains original data)
        └───res  (result will be stored here, the output is organized according to CTC naming convention)
    ...    
└───track.py (contains our tracking code)
└───utils.pt (contains some helper functions, i.e. Activity map and Gaussisan map generation)


