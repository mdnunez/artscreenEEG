<img src="./extra/logo.png" width="128">

# artscreenEEG
### Version 0.14.0

**Authors: Michael D. Nunez, Cort Horton, Siyi Deng, William Winter, and Ramesh Srinivasan from the Human Neuroscience Lab at the University of California, Irvine**

**Head model contributors: Mark Dennison, Jacky Au, Inez Falcon, and Derek C. Monroe at the University of California, Irvine**

artscreenEEG is a MATLAB package to perform basic artifact correction and analysis on electroencephalographic (EEG) data. Please see demo/example_steps.m for an example of how to use the basic functions. 

This software is intended for users who wish to mitigate the muscle and electrical artifact found in all electroencephalographic (EEG) recordings and perform basic EEG analysis, whether collected from the users' own labs or collected at other locations. 

While it is the view of the authors that all EEG recordings contain some amount of artifact and existing software tools can only mitigate this problem, specific analyses (such as band-limited EEG analyses) are more immune to this artifact. See Further Reading for analysis help.

## Getting Started

### Prerequisites

This software should work on any operating system that runs the correct version of MATLAB.

This software works best with [MATLAB R2014a](http://www.mathworks.com/products/matlab/) or slighly older MATLAB versions. This software works with newer versions of MATLAB however some graphical functions may be slower due to GUI drawing changes made by MathWorks.

### Downloading

The repository can be cloned with `git clone https://github.com/mdnunez/artscreenEEG.git`

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

### Further Reading

EEG analysis and artifact correction (see Section 6 for a discussion on EEG artifact correction): [Electroencephalography (EEG): neurophysics, experimental methods, and signal processing](https://www.researchgate.net/publication/290449135_Electroencephalography_EEG_neurophysics_experimental_methods_and_signal_processing)

### References

Oja, E., & Hyvarinen, A., 
A fast fixed-point algorithm for independent component analysis. 
Neural computation, 9(7), 1483-1492. (1997).

Mognon A, Bruzzone L, Jovicich J, Buiatti M, 
[ADJUST: An Automatic EEG artifact Detector based on the Joint Use of Spatial and Temporal features.](https://www.researchgate.net/publication/45268818_ADJUST_An_automatic_EEG_artifact_detector_based_on_the_joint_use_of_spatial_and_temporal_features) 
Psychophysiology 48 (2), 229-240 (2011).

Makeig, S., Bell, A.J., Jung, T-P and Sejnowski, T.J.,
[Independent component analysis of electroencephalographic data.](https://www.researchgate.net/publication/2242002_Independent_Component_Analysis_of_Electroencephalographic_Data)
In: D. Touretzky, M. Mozer and M. Hasselmo (Eds). Advances in Neural 
Information Processing Systems 8:145-151, MIT Press, Cambridge, MA (1996).

Nunez, M. D., Vandekerckhove, J., & Srinivasan, R. 
[How attention influences perceptual decision making: Single-trial EEG correlates of drift-diffusion model parameters.](https://www.researchgate.net/publication/298275031_How_attention_influences_perceptual_decision_making_Single-trial_EEG_correlates_of_drift-diffusion_model_parameters) 
Journal of Mathematical Psychology, 76, 117-130. (2017).

## Contributing

### Authorship

Please submit issues and associated pull requests. All .m code contributors will be listed as authors. All other contributors will be mentioned in the documentation.

### License

artscreenEEG is licensed under the GNU General Public License v3.0 and written by Michael D. Nunez, Cort Horton, Siyi Deng, William Winter, and Ramesh Srinivasan from the Human Neuroscience Lab at the University of California, Irvine.

### Version
The newest version of artscreenEEG is 0.13.1 indicating that artscreenEEG is still being developed (version 0)
