
const DimOrDimTuple = Union{Dimension,Tuple{Vararg{<:Dimension}}}
const IntOrIntTuple = Union{Int,Tuple{Vararg{<:Int}}}

struct Ag end
struct DisAg end

struct Every{V,S}
    val::V
    start::S
end
Every(val; start=nothing) = Every(val, start)


const SCALE = """
- `scale`: the aggregation factor, which can be an `Integer`, a tuple of `Integer`
    for each dimension, or a `Dimension`
"""

const METHOD = """
- `method`: a function such as `mean` or `sum` that can combine the
    value of multiple cells to generate the aggregated cell, or a `Symbol` 
    `:start`, `:center`, or `:end` or `Locus` like `Start()` `Center()` or `End`
    to simply take a use single value form that position in the aggregation area.
"""

"""
    aggregate(method, object, scale; skipmissing)

Aggregate an `AbstractDimArray` or `AbstractDimStack` by `scale` using `method`.

# Arguments

$METHOD
- `object`: Object to aggregate, like `AbstractDimStack`, `AbstractDimArray` or `Dimension`.
$SCALE

When the aggregation `scale` of is larger than the array axis, the length of the axis is used.

# Keywords

- `skipmissing`: if `true`, any `missing` will be skipped during aggregation, so that
    only areas of all missing values will be aggregated to `missing`. If `false`, any
    aggegrated area containing a `missing` will be assigned `missing`.
"""
function aggregate end
aggregate(method, x, scale::Union{Integer,Every}; kw...) =
    aggregate(method, x, map(d -> rebuild(d, scale), dims(x)); kw...) 
aggregate(method, x, scale::Dimension; kw...) = aggregate(method, x, map(_ -> scale, size(x)); kw...) 
function aggregate(method, stack::AbstractDimStack, scale; kw...)
    f(A) = aggregate(method, A, scale; kw...)
    layers = map(f, values(stack))
    return DD.rebuild_from_arrays(stack, layers)
end
function aggregate(method, src::AbstractDimArray, scale::Tuple; kw...)
    dst = alloc_ag(method, src, scale)
    aggregate!(method, dst, src, scale; kw...)
end
aggregate(method, d::Dimension, scale::Integer) = rebuild(d, aggregate(method, lookup(d), scale))
aggregate(method, d::Dimension, scale::Every) = rebuild(d, aggregate(method, lookup(d), scale))
function aggregate(method, lookup::LookupArray, scale::Every)
    start = isnothing(scale.start) ? first(bounds(lookup)) : scale.start
    lookup = 
end
function aggregate(method, lookup::LookupArray, scale::Integer)
    scale == 1 && return lookup
    start, stop = _endpoints(method, lookup, scale)
    newlookup = lookup[start:scale:stop]
    if lookup isa AbstractSampled 
        sp = _aggregate(method, span(lookup), lookup, scale)
        return rebuild(newlookup; span=sp)
    end
    return newlookup
end

function _aggregate(method, span::Explicit, lookup, scale::Integer)
    m = val(span)
    range = firstindex(lookup):scale:lastindex(lookup)
    newm = similar(m, (2, length(range)))
    map(axies(newm, 2)) do i
        newm[1, i] = m[1, (i - 1) * scale + 1]
        newm[2, i] = m[2, i * scale]
    end
end
_aggregate(method, span::Regular, lookup, scale::Integer) = Regular(val(span) * scale)
function _aggregate(method, span::Span, lookup, scale::Irregular)
    range = firstindex(lookup):scale:lastindex(lookup)
    slicesspan(lookup, range)
end


"""
    aggregate!(method, dst::AbstractDimArray, src::AbstractDimArray, scale; skipmissingval=false)

Aggregate array `src` to array `dst` by `scale`, using `method`.

# Arguments

$METHOD
- `scale`: the aggregation factor, which can be an integer, a tuple of integers
  for each dimension, or any `Dimension`, `Selector` or `Int` combination you can
  usually use in `getindex`. Using a `Selector` will determine the scale by the
  distance from the start of the index in the `src` array.

When the aggregation `scale` of is larger than the array axis, the length of the axis is used.

# Keywords

- `progress`: show a progress bar.
- `skipmissingval`: if `true`, any `missingval` will be skipped during aggregation, so that
    only areas of all missing values will be aggregated to `missingval`. If `false`, any
    aggegrated area containing a `missingval` will be assigned `missingval`.

Note: currently it is _much_ faster to aggregate over memory-backed arrays.
Use [`read`](@ref) on `src` before use where required.
"""
aggregate!(method, dst::AbstractDimArray, src, scale::Pair...; kw...) =
    aggregate!(method, dst, src, map(rebuild, scale))
aggregate!(method, dst::AbstractDimArray, src, scale::Integer; kw...) =
    aggregate!(method, dst, src, map(d -> rebuild(d, scale), dims(src)))
aggregate!(method, dst::AbstractDimArray, src, scale::Tuple; kw...) =
    aggregate!(method, dst, src, map(rebuild, dims(src), scale))
aggregate!(method::Locus, dst::AbstractDimArray, src, scale::DimTuple; kw...) =
    aggregate!((method,), dst, src, scale)
# function aggregate!(loci::Tuple{Locus,Vararg}, dst::AbstractDimArray, src, scale::DimTuple; kw...)
#     intscale = _scale2int(Ag(), dims(src, dims(dst)), scale)
#     offsets = _agoffset.(loci, intscale)
#     indices = map(dims(dst), axes(dst), offsets, intscale) do d, ax, o, isc
#         a = upsample(first(ax), isc) + o
#         b = a + intscale - 1
#         rebuild(d, map(:, a, isc, b))
#     end
#     dst .= src[indices...]
# end
# Function/functor methods
function aggregate!(f, dst::AbstractDimArray, src, scale::DimTuple; skipmissing=false)
    aggregator = _dim_aggregator(dst, src, scale)
    if skipmissing
        broadcast!(dst, aggegator) do D 
            f(skipmissing(view(src, D...)))
        end
    else
        broadcast!(dst, aggregator) do D 
            f(view(src, D...))
        end
    end
end

function _dim_aggregator(dst, src, scale::DimTuple)
    map(dims(dst), dims(src), sortdims(scale, dims(dst))) do d, s, sc
        rebuild(d, _agg_indices(d, s, sc))
    end |> DimPoints
end

_agg_indices(dst, src, scale::Nothing) = rebuild(dst, :)
_agg_indices(dst, src, scale::Dimension) = _agg_indices(dst, src, val(scale))
function _agg_indices(dst::Dimension, src::Dimension, scale::Integer)
    map(eachindex(dst)) do i
        a = upsample(i, scale)
        b = a + scale - 1
        rebuild(dst, map(:, a, b))
    end
end
function _agg_indices(src::Dimension, dst::Dimension, scale::Every)
    @assert step(src) == val(scale) "lookup step: $(step(src)) $does not match scale: $(val(scale))"
    map(a) do i
        range = dims2indices(src, Interval(intervalbounds(d, i)...)) 
        rebuild(d, range) 
    end
end

"""
    disaggregate([method,] object, scale; filename, progress, keys)

Disaggregate array, or all arrays in a stack or series, by some scale.

# Arguments

- `object`: Object to aggregate, like `AbstractStack`, `AbstractDimArray` or a `Dimension`.
$SCALE
- `method`: for `disaggregate` this is used only to set the `Locus` of each `Dimension`.
    by default it will use the existing locii of the array.
"""
function disaggregate end
disaggregate(x, scale; kw...) = disaggregate(map(locus, dims(x)), x, scale; kw...)

function disaggregate(locus, stack::AbstractDimStack, scale)
    f(A) = disaggregate(locus, A, scale)
    layers = map(f, values(stack), Tuple(suffix))
    return DD.rebuild_from_arrays(stack, layers)
end
function disaggregate(locus, src::AbstractDimArray, scale; kw...)
    dst = alloc_disag(locus, src, scale; filename, suffix, kw...)
    disaggregate!(locus, dst, src, scale)
end
function disaggregate(locus::Locus, dim::Dimension, scale)
    rebuild(dim, disaggregate(locus, lookup(dim), scale))
end
function disaggregate(locus, lookup::LookupArray, scale)
    intscale = _scale2int(DisAg(), lookup, scale)
    intscale == 1 && return lookup

    len = length(lookup) * intscale
    step_ = step(lookup) / intscale
    start = lookup[1] - _agoffset(locus, intscale) * step_
    stop = start + (len - 1)  * step_
    index = LinRange(start, stop, len)
    if lookup isa AbstractSampled
        sp = disaggregate(locus, span(lookup), intscale)
        return rebuild(lookup; data=index, span=sp)
    else
        return rebuild(lookup; data=index)
    end
end

disaggregate(locus, span::Span, scale) = span
disaggregate(locus, span::Regular, scale) = Regular(val(span) / scale)

"""
    disaggregate!(dst, src, scale)

Disaggregate array `src` to array `dst` by some `scale`.

$SCALE
"""
function disaggregate!(dst::AbstractDimArray, src, scale)
    disaggregate!(dst, src, scale)
end
function disaggregate!(dst::AbstractDimArray, src, scale)
    intscale = _scale2int(DisAg(), dims(src), scale)
    broadcast!(dst, CartesianIndices(dst)) do I
        src[downsample.(Tuple(I), intscale)...]
    end
end

# Symbol locus to skip imports
for f in (:aggregate, :aggregate!, :disaggregate, :disaggregate!)
    @eval begin
        $f(::Type{L}, A::AbstractDimArray, args...) where L<:Locus = $f(L(), A, args...)
        function $f(locus::Symbol, A, args...)
            if locus in (:start, :Start)
                $f(Start(), A, args...)
            elseif locus in (:center, :Center)
                $f(Center(), A, args...)
            elseif locus in (:end, :End)
                $f(End(), A, args...)
            else
                throw(ArgumentError("$locus is not an accepted Symbol input. Use :start, :center or :end"))
            end
        end
    end
end

# Utils

# Allocate an array of the correct size to aggregate `A` by `scale`
function alloc_ag(method, A::AbstractDimArray, scale::Tuple)
    # Aggregate the dimensions
    aggdims = if method isa Tuple
        map(aggregate, method, dims(A), scale)
    else
        map(dims(A), scale) do d, sc
            aggregate(method, d, sc)
        end
    end
    # Dim aggregation determines the array size
    aggsize = map(length, aggdims)
    T = ag_eltype(method, A)
    return rebuild(A; data=similar(parent(A), T, aggsize), dims=aggdims)
end

# Allocate an array of the correct size to disaggregate `A` by `scale`
function alloc_disag(method::Tuple, A::AbstractDimArray, scale)
    intscale = _scale2int(DisAg(), dims(A), scale)
    dims_ = disaggregate.(method, dims(A), intscale)
    # Dim aggregation determines the array size
    aggsize = map(length, dims_)
    T = ag_eltype(method, A)
    return 
end

# Handle how methods like `mean` can change the type
ag_eltype(method::Tuple{<:Locus,Vararg}, A) = eltype(A)
ag_eltype(method::Locus, A) = eltype(A)
function ag_eltype(method, A)
    method_returntype = typeof(method(zero(eltype(A))))
    promote_type(eltype(A), method_returntype)
end

# Convert indicies from the aggregated array to the larger original array.
upsample(index::Int, scale::Int) = (index - 1) * scale + 1
upsample(index::Int, scale::Colon) = index

# Convert indicies from the original array to the aggregated array.
downsample(index::Int, scale::Int) = (index - 1) ÷ scale + 1
downsample(index::Int, scale::Colon) = index

# Convert scale or tuple of scale to integer using dims2indices
function _scale2int(x, dims::DimTuple, scale::Tuple)
    map((d, s) -> _scale2int(x, d, s), dims, DD.dims2indices(dims, scale))
end
_scale2int(x, dims, scale::Dimension) = _scale2int(x, dims, val(scale))
_scale2int(x, dims::DimTuple, scale::Int) = map(d -> _scale2int(x, d, scale), dims)
_scale2int(x, dims::DimTuple, scale::Colon) = map(d -> _scale2int(x, d, scale), dims)
_scale2int(x, dim::Dimension, scale::Int) = _scale2int(x, lookup(dim), scale)
_scale2int(::Ag, l::LookupArray, scale::Int) = scale > length(l) ? length(l) : scale
_scale2int(::DisAg, l::LookupArray, scale::Int) = scale
_scale2int(x, dim::Dimension, scale::Colon) = 1

_agoffset(locus::Locus, l::LookupArray, scale) = _agoffset(locus, scale)
_agoffset(method, l::LookupArray, scale) = _agoffset(locus(l), scale)
_agoffset(locus::Start, scale) = 0
_agoffset(locus::End, scale) = scale - 1
_agoffset(locus::Center, scale) = scale ÷ 2

function _endpoints(method, l::LookupArray, scale)
    start = firstindex(l) + _agoffset(method, l, scale)
    stop = (length(l) ÷ scale) * scale
    return start, stop
end
