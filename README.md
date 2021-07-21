
# Single-pixel based Data Fusion algorithm for spectral-temporal-spatial reconstruction

![fusion](https://user-images.githubusercontent.com/19323057/117669638-68276880-b1a7-11eb-975c-fa23b8beda13.png)


Codes for the data fusion procedure explained in: 

"Giga-voxel multidimensional fluorescence imaging combining single-pixel detection and data fusion"

F. SOLDEVILA, A. LENZ, A. GHEZZI, A. FARINA, C. D’ANDREA, AND E. TAJAHUERCE

*** To be Published ***

Codes by F.SOLDEVILA and A.LENZ

### Files explanation ###

data_fusion.m contains the main algorithm of data fusion, with some calls to the routines in \routines. It solves the gradient descent problem explained in the manuscript for a simulated 4D dataset.


datasets16-64-128-128.mat contains simulated measurements for the three detection systems in the manuscript (a single-pixel spectral camera, a single-pixel time-resolved camera, and a high-resolution pixelated detector) AND calibration measurements for the experimental devices (wavelength callibration).

Data.CCD: pixelated detector measurement

Data.L16: single-pixel spectral data

Data.PMT: single-pixel time-resolved data

Data.lmabda: wavelength vector for the 16 channels of the single-pixel spectral data

Data.time: time vector for the single-pixel time-resolved data


Code is provided as is. If you find it useful, please cite the paper
