# rgbifanimals_shiny

## Run locally
```bash
Rscript shiny.R
```

## Run in a docker container
```bash
docker build -f dockerfile -t timbosman/shiny .
docker run -id  -p 80:80 --name shiny timbosman/shiny
```
