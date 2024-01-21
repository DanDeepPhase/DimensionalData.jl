"""
    DimGroupByArray <: AbstractDimArray

`DimGroupByArray` is essentially a `DimArray` but holding
the results of a `groupby` operation.

This wrapper allows for specialisations on later broadcast or
reducing operations, e.g. for chunk reading with DiskArrays.jl,
because we know the data originates from a single array.
"""
struct DimGroupByArray{T,N,D<:Tuple,R<:Tuple,A<:AbstractArray{T,N},Na,Me} <: AbstractDimArray{T,N,D,A}
    data::A
    dims::D
    refdims::R
    name::Na
    metadata::Me
    function DimGroupByArray(
        data::A, dims::D, refdims::R, name::Na, metadata::Me
    ) where {D<:Tuple,R<:Tuple,A<:AbstractArray{T,N},Na,Me} where {T,N}
        checkdims(data, dims)
        new{T,N,D,R,A,Na,Me}(data, dims, refdims, name, metadata)
    end
end
function DimGroupByArray(data::AbstractArray, dims::Union{Tuple,NamedTuple};
    refdims=(), name=NoName(), metadata=NoMetadata()
)
    DimGroupByArray(data, format(dims, data), refdims, name, metadata)
end
@inline function rebuild(
    A::DimGroupByArray, data::AbstractArray, dims::Tuple, refdims::Tuple, name, metadata
)
    DimArray(data, dims, refdims, name, metadata) # Rebuild as a reguilar DimArray
end

function Base.summary(io::IO, A::DimGroupByArray{T,N}) where {T,N}
    print(io, _array_dim_description(A))
    print(io, string(nameof(typeof(A)), "{$(nameof(T)),$N}"))
end

function show_after(io::IO, mime, A::DimGroupByArray)
    A1 = DimArray(map(x -> map(rebuild, basedims(x), size(x)), A))
    printstyled(io, "  groups dimensions and sizes:"; color=:light_black)
    println()
    show_after(IOContext(io, :compact=>true, :typeinfo=>eltype(A1)), mime, A1)
end

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

Currently this returns a `DimArray` or `DimArray` of the original `AbstractDimArray`
or `AbstractDimStack` types. This may change to a
"""
groupby(A::AbstractDimArray, dimfunc::Dimension...; kw...) = groupby(A, (dimfunc,); kw...)
function groupby(A::AbstractDimArray, pairs::Pair...; kw...)
    dims = map(pairs) do (d, v)
        rebuild(basedims(d), v)
    end
    groupby(A, dims; kw...)
end
function groupby(A::AbstractDimArray, dimfuncs::DimTuple; kw...)
    dims_indices = map(dimfuncs) do d
        group_lookup, group_indices = _group_indices(lookup(A, d), DD.val(d))
        rebuild(d, group_lookup), rebuild(d, group_indices)
    end
    dims = map(first, dims_indices)
    indices = map(last, dims_indices)
    # Hack DimPoints to get tuples of all view indices combinations
    views = parent(map(D -> view(A, map(rebuild, dims, D)...), DimPoints(indices)))
    metadata = Metadata(:groupby => dimfuncs)
    return DimGroupByArray(views, format(dims, views), (), :groupby, metadata)
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
