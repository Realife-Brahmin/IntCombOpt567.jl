module helperFunctions

export 
    argmax_smallestkey,
    myprintln

"""
    argmax_smallestkey(dict::Dict)

Returns the key with the maximum value. If multiple keys have the same
maximum value, returns the smallest key numerically among them.
"""
function argmax_smallestkey(dict::Dict)
    max_val = maximum(values(dict))
    keys_with_max = [k for (k, v) in dict if v == max_val]
    return minimum(keys_with_max)
end


#region myprintln
"""
    myprintln(verbose::Bool, msg::String)

Prints a message to the console if the `verbose` flag is set to `true`.

# Arguments
- `verbose::Bool`: A flag indicating whether to print the message.
- `msg::String`: The message to be printed.
"""
function myprintln(verbose::Bool, msg)
    if verbose
        println(msg)
    end
end
#endregion

end # module helperFunctions