function compress_intarray(arr::AbstractArray{<:Union{Integer, Missing}})
    if all(ismissing, arr)
        return arr
    end
    min_val, max_val = extrema(skipmissing(arr))
    # Determine the smallest integer type that fits the range of values
    target_type = if min_val >= typemin(Int8) && max_val <= typemax(Int8)
        Int8
    elseif min_val >= typemin(Int16) && max_val <= typemax(Int16)
        Int16
    elseif min_val >= typemin(Int32) && max_val <= typemax(Int32)
        Int32
    else
        Int64
    end
    # Convert non-missing elements to target_type, leaving missing values as-is
    return map(x -> ismissing(x) ? missing : convert(target_type, x), arr)
end


function compress_scalarint(x:: Union{Int64, Int32, Int16, Int8, Missing})
    ismissing(x) ? missing :
    x >= typemin(Int8) && x <= typemax(Int8) ? Int8(x) :
    x >= typemin(Int16) && x <= typemax(Int16) ? Int16(x) :
    x >= typemin(Int32) && x <= typemax(Int32) ? Int32(x) : Int64(x)
end

# FAILS AND I DON'T KNOW WHY
# function compress_tupleintarray(arr::AbstractArray{<:Union{Tuple, Missing}})
#     if all(x -> x isa Int, Iterators.flatten(skipmissing(arr)))==false
#         throw(ArgumentError("All data in tuples must be one of Int64, Int32, Int16, or Int8"))
#     end
#     min_val, max_val = extrema(Iterators.flatten(skipmissing(arr)))
#     target_type = if min_val >= typemin(Int8) && max_val <= typemax(Int8)
#         Int8
#     elseif min_val >= typemin(Int16) && max_val <= typemax(Int16)
#         Int16
#     elseif min_val >= typemin(Int32) && max_val <= typemax(Int32)
#         Int32
#     else
#         Int64
#     end
#     # Apply the type minimization across each element, preserving array structure
#     return map(t -> ismissing(t) ? missing : Tuple(target_type(x) for x in t), arr)
# end

function compress_tuplefloatarray(arr::AbstractArray{<:Union{Tuple, Missing}}, target_type::Type{<:AbstractFloat})
    if all(ismissing, arr)
        return arr
    end
    if all(x -> x isa AbstractFloat, Iterators.flatten(skipmissing(arr)))==false
        throw(ArgumentError("All data in tuples must be one of BigFloat, Float64, Float32, or Float16"))
    end
    min_value, max_value = extrema(Iterators.flatten(skipmissing(arr)))
    if min_value<=nextfloat(typemin(target_type)) || max_value>=prevfloat(typemax(target_type)) 
        throw(ArgumentError("Input data outside of target data type range"))
    end
    # Apply the type minimization across each element, preserving array structure
    return map(t -> ismissing(t) ? missing : Tuple(target_type(x) for x in t), arr)
end

function compress_floatarray(arr:: AbstractArray{<:Union{AbstractFloat, Missing}}, target_type::Type{<:AbstractFloat}) # This doesn't work on arrays that have no missing values, for some reason
    if all(ismissing, arr)
        return arr
    end
    # Check if target_type is one of the allowed types
    if !(target_type in (BigFloat, Float64, Float32, Float16))
        throw(ArgumentError("Target type must be one of BigFloat, Float64, Float32, or Float16"))
    end
    min_value, max_value = extrema(skipmissing(arr))
    if min_value<=nextfloat(typemin(target_type)) || max_value>=prevfloat(typemax(target_type)) 
        throw(ArgumentError("Input data outside of target data type range"))
    end
    return map(x -> ismissing(x) ? missing : convert(target_type, x), arr)
end

function compress_structure(s, target_type::Type{<:AbstractFloat}=Float64)
    if target_type==Float64
        return s
    end
    # Iterate over each field in the structure
    for name in fieldnames(typeof(s))
        # println(name)
        field = getfield(s, name)
        # Check the type and apply the appropriate function
        if isa(field, AbstractFloat)  # Scalar float
            converted_field=convert(target_type, field)
            # println("Converting field $(name) to type $(typeof(converted_field))")
            setfield!(s, name, converted_field)
        elseif isa(field, AbstractArray{<:Union{AbstractFloat, Missing}})  # Array of floats
            converted_field = compress_floatarray(field, target_type)
            # println("Converting field $(name) to type $(typeof(converted_field))")
            setfield!(s, name, converted_field)
        elseif isa(field, AbstractArray{<:Union{Tuple, Missing}})  # Array of tuples
            eltype_t = eltype(field)  # Type inside tuples
            if all(isa(t, Tuple{Vararg{<:Integer}}) for t in field)
                # FAILS AND I DON'T KNOW WHY
                # setfield!(s, name, compress_tupleintarray(field)) # Tuples of integers
            elseif all(isa(t, Tuple{Vararg{<:AbstractFloat}}) for t in field)
                # setfield!(s, name, compress_tuplefloatarray(field, target_type)) # Tuples of floats
            end
        end
    end
    return s
end
