# TeamsProject
Codes for analysis and figure making of teams paper
## Folder structure:

The data required to generate the figures (includes raw simulation files as well as compilations of various analyses) is available for download at the link:
Download this repo, then download the rar file at the above and extract it to the same folder as the repo. The data should all be present in a folder named "Final_Results". Once the folder structure is ready, the codes can be run to generate the figures. 

## Code structure
The topo files for all EMP as well as random networks are available in the folder "Raw_data" in the "Final_Results" folder. The networks were simulated using the julia codes at :https://github.com/csbBSSE/bmodel_julia. 
The codes are provided in two folders: Figures and Analysis. The analysis folder has functions and scripts used to perform various analysis on the data including calculation of stability metrics for all networks. The Figures folder contains all the figure related codes as well as the outputs. To run the figure codes, please change the variable "mainFolder" to the directory in which the Final_Results folder resides in all figure related scripts. 
