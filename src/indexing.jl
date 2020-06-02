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

@inline packed_data(seq::Kmer) = seq.data


@inline function inbounds_getindex(seq::Kmer, i::Integer)
    bidx = bitindex(seq, i)
    elem = extract_encoded_element(bidx, packed_data(seq))
    return decode(Alphabet(seq), elem)
end

