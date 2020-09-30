# SpaceM -  method for single-cell metabolomics

This repository supplements the manuscript [(Rappez et al, 2020)](https://www.biorxiv.org/content/10.1101/510222v1) which presents SpaceM, a method for in situ single-cell metabolomics of cultured cells that integrates microscopy with MALDI-imaging mass spectrometry. The manuscript is currently under review. The SpaceM datasets presented in the manuscript are available on [MetaboLights](https://www.ebi.ac.uk/metabolights/reviewer417760fcbfbb6076b4ce5bd9a7e7c893).

For a detailed introduction into the method, presentation of the results of SpaceM investigations of various cells, and discussion of its advantages and capabilities, please read our [manuscript](https://www.biorxiv.org/content/10.1101/510222v1).

Here, we present the source code implementing the method. Moreover, we [present interactively using Google Collab](https://colab.research.google.com/drive/1CKdHDUkGIpAcBzrSfuCodMF_l2xbVAKT?usp=sharing) the downstream processing of the spatio-molecular matrices provided by SpaceM and replicate all main figures of the manuscript.


### Installation

We support `python3`. For the detailed requirements, see the file 'requirements.txt'. To install the dependencies, run:

`pip install -r requirements.txt`

Download and install [CellProfiler 3.0.0](https://cellprofiler.org/previous_releases/)
Download and install [Fiji release of December 22 2015](https://imagej.net/Fiji/Downloads)

In `paths.json` add the path of CellProfiler to `"CellProfiler path"` and the path of Fiji to `"Fiji path"`

The installation time is less than one hour.

### Data requirement

The main molder `MF` for SpaceM analysis should be created manually and organized as follows: 

```
MF
|
└─ Input
    |
    └─ MALDI
    |    *.RAW
    |    *.UDP
    |    *.imzML
    |    *.ibd
    |        
    └─ Microscopy
         |
         └─ preMALDI
         |    tile_1.tif
         |    tile_2.tif
         |    ...
         |    out.txt
         |   
         |
         └─ postMALDI
              tile_1.tif
              tile_2.tif
              ...
              out.txt
```


#### MALDI-imaging mass spectrometry

The `.RAW`, `.UDP`, `.imzML` and `.ibd` files are required. The data should be analyzed by using [METASPACE]( https://metaspace2020.eu/).
It is critical that ablation marks are visible and non-overlapping in the post-MALDI microscopy. For more details, see the Methods section in the manuscript.

#### Microscopy

For both pre- and post-MALDI images, a tiled acquisition of the cell culture area with black penmarks is required. 
The individual tiles should be stored with a text file `out.txt` containing the upper left pixel x-y coordinates in um of each tile. 
At the moment, SpaceM requires the stitching of the tiles and this particular code is optimized for the Nikon Ti-E microscope output format. However, this is a limitation of this implementation and not of the SpaceM method. The next version of the implementation will also be accepting pre-stitched images.

Once the processing is finished, add the path of the Main folder `MF` to the `"MF"` entry in the `path.json` file.

The CellProfiler (CP) software, CP pipeline for segmentation, and CP project files able to segment the cells are required. 
Both the path of the CP pipeline and project should be added to the `"CellProfiler pipeline path"` and 
`"CellProfiler project path"` in `paths.json`, respectively.

### Execute SpaceM

Run `python runAnalysis.py` and follow instructions when prompted. 
The produced results will be stored in the `Analysis` sub-folder. 
The final spatio-molecular matrix will be stored as `MORPHnMOL.csv` and can be found inside the `scAnalysis` sub-folder.

```
MF
|
└─ Analysis
    |
    └─ scAnalysis
    |    MORPHnMOL.csv
    
    ...
```




