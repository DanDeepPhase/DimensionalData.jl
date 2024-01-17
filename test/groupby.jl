using DimensionalData, Test, Dates
using DimensionalData: Dimensions
using DimensionalData: LookupArrays
using Statistics
const DD = DimensionalData

# Make a demo DimArray
A = rand(X(1:0.1:2), Y(1:20), Ti(DateTime(2000):Day(10):DateTime(2003)))
# Group by month and even/odd Y axis values
groupmeans = groupby(A, Ti=month, Y=isodd)
groupmeans = mean.(groupby(A, Ti=month, Y=isodd))
