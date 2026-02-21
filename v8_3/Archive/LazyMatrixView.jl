#######################################################################################
# CREATE ARRAY LAZY VIEWS TO SAVE ON MEMORY
#######################################################################################
# I need to create a massive matrix of cartesian indexes to create the gradient of the probability of staying single
# Instead of creating the array right away, I build a view.
struct LazyMatrixView{T, A <: Tuple, B <: Tuple} <: AbstractMatrix{T}
    A :: Vector{A}
    B :: Vector{B}
    dims :: Tuple{Vararg{Int}}  # Tuple for dimensions
end
# Define the size of the matrix
Base.size(lmv::LazyMatrixView) = (length(lmv.B), length(lmv.A))
# Define the axes for the matrix
Base.axes(lmv::LazyMatrixView) = (1:length(lmv.B), 1:length(lmv.A))
# Efficient single index getindex using LinearIndices without intermediate array creation
Base.getindex(lmv::LazyMatrixView, i::Int, j::Int) = LinearIndices(lmv.dims)[CartesianIndex(
    lmv.A[j]...,    # Spread the elements of the tuple from A
    lmv.B[i][4:end]...  # Spread the elements from the fourth to the end of the tuple from B
)]
# Overloaded getindex for slices, returning a matrix of linear indices
Base.getindex(lmv::LazyMatrixView, I::Union{Int, AbstractRange}, J::Union{Int, AbstractRange}) = [
    LinearIndices(lmv.dims)[CartesianIndex(lmv.A[j]..., lmv.B[i][4:end]...)] for i in I, j in J
]
# Implement checkbounds to ensure indices are valid
Base.checkbounds(lmv::LazyMatrixView, i, j) = 
    checkbounds(Bool, axes(lmv, 1), i) && checkbounds(Bool, axes(lmv, 2), j)


    
    
    