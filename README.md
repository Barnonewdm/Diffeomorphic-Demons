# Requirement
Please compile ITK4.13, set ITKV3_COMPATIBILITY and Module_ITKReview on.
* Note there may be some issues when compiling ITK, please pay attention to the ITK version and gcc version.

# Installation
`git clone https://github.com/Barnonewdm/Diffeomorphic-Demons.git`

`cd Diffeomorphic-Demons`

`mkdir build`

`cmake ..`

`make`


# Usage
DemonsRegistration is used to perform predict the deformation field.
DemonsWarp is used to warp the moving image with the predicted deformation field.

# Acknowledgements
Thanks to Tom Vercauteren for his initial work at 2007. Any issue is welcomed to report.
