<img src="./extra/logo.png" width="128">

# artscreenEEG

**Authors: Michael D. Nunez, Cort Horton, Siyi Deng, William Winter, Mark Dennison, Jacky Au, and Ramesh Srinivasan from the Human Neuroscience Lab at the University of California, Irvine**

artscreenEEG is a MATLAB function package to perform basic artifact correction on electroencephalographic (EEG) data. Please see demo/example_steps.m for an example of how to use these functions. 

This software is intended for users who wish to mitigate the muscle and electrical artifact found in all electroencephalographic (EEG) recordings, whether collected from the users' own labs or collected at other locations. Note that it is the view of the authors that all EEG recordings contain some amount of artifact and software tools can only mitigate this problem, despite claims by many EEG scientists.

## Getting Started

### Prerequisites

This software probably requires [MATLAB R2014a](http://www.mathworks.com/products/matlab/) or slighly earlier. Currently, the software was only extensively tested on R2014a. Later versions of MATLAB have not been extensively tested. On later versions some functions may be slow due to GUI drawing changes made by MathWorks.

This software should work on any operating system that runs the correct version of MATLAB.

### Downloading

The repository can be cloned with `git clone git@github.oit.uci.edu:mdnunez1/artscreenEEG.git`

The repository can also be may download via the _Download zip_ button above.

### Installation

After downloading/unzipping the repository, users will need to add these functions to the MATLAB path. In MATLAB, add the repository to the PATH with

```matlab
%Set 'artloc' to full directory path
artloc = 'C:\Users\MATLAB\artscreenEEG';
addpath(genpath(artloc));
```

### Usage

For a demonstration of useage, run the following. IMPORTANT NOTE: This script will ask your premission to download example EEG data from figshare.com

```matlab
%Open example script
edit example_steps.m;
%Run example script
example_steps;
```

### License

artscreenEEG is licensed under the GNU General Public License v3.0 and written by Michael D. Nunez, Cort Horton, Siyi Deng, William Winter, Mark Dennison, Jacky Au, and Ramesh Srinivasan from the Human Neuroscience Lab at the University of California, Irvine.

### Further Reading

Possible EEG artifact correction reference and citation:
Section 6 and Figure 4 of [Electroencephalography (EEG): neurophysics, experimental methods, and signal processing](https://www.researchgate.net/publication/290449135_Electroencephalography_EEG_neurophysics_experimental_methods_and_signal_processing)

### References

Makeig, S., Bell, A.J., Jung, T-P and Sejnowski, T.J.,
[Independent component analysis of electroencephalographic data.](https://www.researchgate.net/publication/2242002_Independent_Component_Analysis_of_Electroencephalographic_Data)
In: D. Touretzky, M. Mozer and M. Hasselmo (Eds). Advances in Neural 
Information Processing Systems 8:145-151, MIT Press, Cambridge, MA (1996).

Mognon A, Bruzzone L, Jovicich J, Buiatti M, 
[ADJUST: An Automatic EEG artifact Detector based on the Joint Use of Spatial and Temporal features.](https://www.researchgate.net/publication/45268818_ADJUST_An_automatic_EEG_artifact_detector_based_on_the_joint_use_of_spatial_and_temporal_features) 
Psychophysiology 48 (2), 229-240 (2011).