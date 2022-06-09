library(testthat)
source("functions.R")
############## create_rastered_world ##############
############## get_occurrence_data ##############
species1 = "Patella aspera"
species2 = "Nomen nescio"

df1 <- get_occurrence_data(species1)
expect_equal(nrow(df1), 217)
expect_error(get_occurrence_data(species2))

rm(species1, species2, df1)
############## check_occurrence_data ##############
species1 = "Patella aspera"
species2 = "Nomen nescio"

df1 <- get_occurrence_data(species1)
expect_equal(nrow(df1), 217)
expect_error(get_occurrence_data(species2))

rm(species1, species2, df1)
############## plot_distribution ##############
############## find_closest_registered_place ##############
############## filter_n_closest_coordinate_ceiling ##############
############## sp_format ##############
############## find_shortest_route_in_sea ##############
############## plot_shortest_path ##############
############## filter_on_distance ##############