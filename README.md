# rgbifanimals_shiny

## Origin of this project
This project is a fork of the rgbif project of [DvBMolEc](https://github.com/DvBMolEc/rgbifanimals). That project focussed on data analysis, where this project makes a shiny dashboard for the data that has been analysed in the main project. Information in this repository originates from the original project, [GBIF database](https://www.gbif.org/), and [NCBI taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy)

## Run locally
```bash
Rscript shiny.R
```

## Run in a docker container
```bash
docker build -f dockerfile -t timbosman/shiny .
docker run -id  -p 80:80 --name shiny timbosman/shiny
```


