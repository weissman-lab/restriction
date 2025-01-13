# Development and validation are csv files in which each row represents the 
# result of a single pulmonary function test, containing the following columns:
# site == pulmonary diagnostic lab at which PFT was performed
# age == patient age in years at the time of pulmonary function testing
# height = patient height (cm)
# weight = patient weight (kg)
# sex == 1 if male, 2 if female
# race == 1 if Asian, 2 if Hispanic, 3 if non-Hispanic Black, 4 if non-Hispanic
#   White, 5 if Other
# fev05 = forced expiratory volume in 0.5 seconds (L)
# fev1 = forced expiratory volume in 1 second (L)
# fev3 = forced expiratory volume in 3 seconds (L)
# fev6 = forced expiratory volume in 6 seconds (L)
# fvc = forced vital capcity (L)
# fev1_fvc = fev1 / fvc
# fef2575 = forced expiratory flow between 25% and 75% of vital capacity (L/s)
# fefmax = maximum forced expiratory flow (L/s)
# expiratory_time = expiratory time (s)
# tlc = total lung capacity (L)

development <- read_csv ("../data/development.csv")

validation <- read_csv ("../data/validation.csv")

