# rgbifanimals
These 4 R scripts will use BOLDigger output (or any other formatted sequence data containing taxonomic information), and determine the presence and diversity of alien species using GBIF occurrence data.

## User guide
1. Copy the four R scripts in a new folder
2. Copy the folder _Inputs_ in the just created folder
3. The _Inputs_ folder contains a tr.rdata file, which is used to calculate distances in 3_Main_script.R. Using this file saves some computing time.
4. Necessary input files have to be placed in the _Inputs_ filder and consist of
  - The _Boldigger_output.csv_, an example of which can be found in the repository
  - A _MetaData.csv_, of which the structure can be found in the example file
  - A _Synonyms.csv_ file, which contains taxonomic synonyms. This should be redundant as 3_Functions.R contains a function which checks the official name on WORMS, but making determined list can speed up processing.
5. Run the R scripts in order.
  - 1_Preperation.R will prepare the data for use in 3_Main_script.R and 4_Visualisation.R by creating several additional dataframes from the BOLDigger output
    All necessary files from 1_Preperation.R will be saved into a new folder _Output_.
  - 2_Functions.R contains all the functions used in 3_Main_script.R.
    - There is no need to actually run this script, everything will be run in the next script
  - 3_Main_script.R will calculate the shortest for each species found at each sampling location, to the nearest GBIF occurrence. 
    - Both the distance as the crow flies, as well as the distance over sea are calculated.
    - Depending on the amount if input data, the latter can take several hours of computing time.
  - 4_Visualisation.R will use the dataframes from 1_Preparation.R and the distances calculated in 3_Main_script.R to create several plots regarding native and alien species distribution, and diversity indices.
  - 5_Results.R is not necessary to run, but will provide insight in the ShortestPath versus DistanceOverSea differences
