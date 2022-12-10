# EMEC436-Final-Project

## dat_to_sdf.jl
Contains functions for converting dat airfoil files to an SDF on an arbitrary grid. The most convinient for testing is testrun(;dir, n, m)

Can convert dat files to Star-CCM readable CSV with Airfoil(address_of_dat_file)

Can interpolate CSV data tables from Star-CCM to a Julia array with get_data_arrays()


## misc_functions.jl

Contains several graphing and interpolation functions.

Can convert data into a format fit for training with convert_to_trainable() and convert_to_small()


## neural_model.jl

Contains the code to train the neural network

Working models can be saved with:

@save "training_data/working_model_V6.jld2" model test_losses train_losses

And loaded with:

@load "training_data/working_model_V5.jld2" model test_losses

And run on correctly formatted data, (x), with:

y = model(x)

view_progress(test_x, test_y, best_model; n=3) plots a nice grid of data, n determines which airfoil is run

view_final(x,y,u;n=1) saves plots to drive
