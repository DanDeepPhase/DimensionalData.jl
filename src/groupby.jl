"""
    groupby(A::AbstractDimArray, dim::Dimension)


# Example

```julia
using DimensionalData, Dates
A = rand(X(1:0.1:20), Y(1:20), Ti(DateTime(2000):Day(3):DateTime(2003)))
# Group by month and even/odd Y axis values
groupmeans = mean.(groupby(A, Ti=month))
# or just means along the time dimension within groups
groupmeans = mean.(groupby(A, Ti=month); dims=Y)
# Or do something else with Y
groupmeans = mean.(groupby(A, Ti=month, Y=isodd))
```
"""
groupby(A::AbstractDimArray; kw...) = groupby(A, DD.kwdims(kw))
groupby(A::AbstractDimArray, dimfunc::Dimension...) = groupby(A, (dimfunc,))
function groupby(A::AbstractDimArray, dimfuncs::DimTuple)
    dims_indices = map(dimfuncs) do d
        group_lookup, group_indices = _group_indices(lookup(A, d), DD.val(d)) 
        rebuild(d, group_lookup), rebuild(d, group_indices)
    end
    dims = map(first, dims_indices)
    indices = map(last, dims_indices)
    # Hack DimPoints to get tuples of all view indices combinations
    views = map(D -> view(A, map(rebuild, dims, D)...), DimPoints(indices))
    return rebuild(views; dims=format(dims, views))
end

function _group_indices(lookup, f::Function)
    indices_dict = Dict()
    for (i, x) in enumerate(lookup)
         k = f(x)
         inds = get!(() -> Int[], indices_dict, k)
         push!(inds, i)
    end
    sorted = sort!(collect(pairs(indices_dict)))
    return first.(sorted), last.(sorted)
end
