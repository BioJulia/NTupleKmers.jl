###
### Indexing
###

@inline function inbounds_getindex2(x::Kmer{A,K,N}, i::Integer) where {A,K,N}
    # Emulation of BitIndex type
    i′ = i + div(64N - 2K, 2)
    val = (i′ - 1) << 1
    idx = (val >> 6) + 1
    off = 62 - (val & (UInt8(64) - 0x01))
    @inbounds begin
        chunk = x.data[idx]
    end
    bits = (chunk >> off) & UInt64(3)
    return reinterpret(eltype(x), 0x01 << bits)
end


"""
Get the Vector or NTuple of unsigned integers that the encoded  
"""
function packed_data end

@inline packed_data(kmer::Kmer) = kmer.data













@inline function extract_encoded_element(bidx::BitIndex{2,UInt64}, x::NTuple{N,UInt64}) where {N}
    @inbounds chunk = x[index(bidx)]
    offchunk = chunk >> (62 - offset(bidx))
    return offchunk & bitmask(bidx)
end

"Extract the element stored in a packed bitarray referred to by bidx."
@inline function extract_encoded_element(bidx::BitIndex{N,W}, data::AbstractArray{W}) where {N,W}
    @inbounds chunk = data[index(bidx)]
    offchunk = chunk >> offset(bidx)
    return offchunk & bitmask(bidx)
end