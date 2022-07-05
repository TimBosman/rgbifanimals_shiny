#Input data
straightLine <- read.csv("Output/ShortestPath.csv")
throughSea <- read.csv("Output/DistanceOverSea.csv")

test <- merge(throughSea, straightLine, by.y = "Speciesname", by.x = "Specieslist")
for(item in levels(test$name)){
  test$straight[test$name == item] <- test[test$name == item, item] 
}

#remove unused sampling locations
test <- test[,!names(test) %in% c(levels(test$name), "Galway", "Gdynia")]

#plot
ggplot(test)+
  geom_point(aes(x = straight , y = distance), col = "green") +
  geom_abline()

ggplot(test) + 
  geom_density(aes(x = log10(straight), fill = is.na(distance)), alpha = 0.5)
