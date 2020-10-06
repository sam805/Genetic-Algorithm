# Genetic Algorithm
This code searches and finds near-optimal sets of auxiliary variables (called closure relationship) to increase the prediction accuracy of a mechanistic model (the unified model).
The GA function accepts 5 input parameters: filePath, pop.size, iteration_num, exePath, and resultPath.
filePath refers to the csv file address, which is considered a raw input values for the mechanist model (unified model). There is a method that converts this file into an acceptable format for a unified model. Unified model (stored in exePath address is in the form of a .exe file) executed, and it generated an output file. There is a function that parses the output file to export the required information. This information is used to calculate the fitness function.
Mutation probability is applied with the 0.03 of probability.
A two-point-crossover is used in this genetic algorithm.
pop.size refers to the number of chromosomes in the population.
iteration_num refers to the number of epochs to run the genetic algorithm.
Results path is the address that includes the sorted chromosomes for each epoch together with the corresponding fitness values and the source of chromosomes (if they are results of direct transfer,  crossover, or mutation).
