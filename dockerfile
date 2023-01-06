FROM r-base
COPY . /usr/local/src/shinyapp
WORKDIR /usr/local/src/shinyapp
RUN apt-get update && apt-get install -y libgdal-dev libxml2-dev
RUN Rscript -e 'install.packages("requiRements")'
RUN Rscript -e 'requiRements::install(path_to_requirements = "/usr/local/src/shinyapp/requirements.txt")'
CMD ["Rscript", "shiny.R"]
