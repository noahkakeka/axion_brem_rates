# axion_brem_rates
This is the c++ code used in my master's thesis. It calculates the nucleon-nucleon-axion bremsstrahlung by inputting a density and mass to the code and outputting the rate.

To run, first download the gsl library. 

Create the exicutable by  g++ -o exicutable_name axion_brem_rates.cpp -lgsl -lgslcblas -lm

Then input a density and mass when running the executable and the location where you want the results to be placed 
./exicutable_name density mass /file_location/
