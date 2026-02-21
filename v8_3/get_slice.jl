# Use views to avoid copying slices
function get_slice(data, start, len)
    return @view data[start:(start+len-1)]
end