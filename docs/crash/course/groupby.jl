using DimensionalData
using Dates

# # DateTime operations

# Make a demo DimArray
tempo = Ti(range(DateTime(2000), step=Hour(6), length=365*5))
ds = rand(X(1:0.5:2), tempo)

# The following are `functions` that can be applied on a DateTime array.

# ## hour
hour.(tempo)

# ## dayofweek
dayofweek.(tempo)

# ## month
month.(tempo)

# ## dayofyear
dayofyear.(tempo)

# ## season
## TODOS: We will need several functions.

## TODO, we need a new function that can return 'DJF', 'DJF', 'DJF', ..., 'DJF', 'DJF', 'DJF'.

# ## select by month, days, years and seasons
## TODO, how do we select month 1 or 2, and even a group of them, i.e. [1,3,5]? Same for days, years and seasons.

# # groupby operations
# Here we use the same time functions from above

mean.(groupby(ds, Ti=month)) # is combining month from different years
## TODO. How do we aggregate by arbitrary months, let' say I want to aggregate every 3 months?

#
mean.(groupby(ds, Ti=year))

#
mean.(groupby(ds, Ti=yearmonth))

#
mean.(groupby(ds, Ti=hour))

#
mean.(groupby(ds, Ti=Dates.hour12))

## TODO. How do could we incorporate resample? Let's say if we have hour resolution I want to resample every 3,6,12.. hours?

mean.(groupby(ds, Ti=dayofyear)) # it will combine the same day from different year.

mean.(groupby(ds, Ti=yearmonthday)) # this does the a daily mean aggregation.

mean.(groupby(ds, Ti=yearmonth)) # this does a monthly mean aggregation

## TODO. Similar to the hourly resample, how do we do it for more than 1 day, let's say 8daily?

## TODO: Group by Dims. This should include the rasters input sampling.

# ## Apply custom function (i.e. normalization) to grouped output.
