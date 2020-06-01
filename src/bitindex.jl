# BitIndex
# --------
#
# Utils for indexing bits in a vector of unsigned integers (internal use only).
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md
#
#         index(i)-1        index(i)        index(i)+1
# ....|................|..X.............|................|....
#                          |<-offset(i)-|
#                      |<--- 64 bits -->|

"""
A type trait of any bitvector-like type, that specifies how many bits are
required per element.

In this case out bitvector-like types are BioSequences, as they store sequences
of biological symbols in a succinct, packed form.

This trait should be defined and overloaded for your own bitvector-like types
or own concrete BioSequence subtypes.

!!! tip
    In BioSequences.jl, this trait is defined as a trait of the sequence type.
    
    But as a generic fallback, a default `BitsPerElem` trait for a sequence type
    is determined by the sequence type's Alphabet type trait.
    
    For example if you made a new sequence type, with an `Alphabet` trait of
    `DNAAlphabet{4}` and did not overload a `BitsPerElem` method specifically
    for the type, calling `BitsPerElem` on a sequence of that type results in
    a fallback result of `BitsPerElem{4}`, as dictated by the `DNAAlphabet{4}`.
    
    This needs repeating: BitsPerElem is a trait of sequence types, not of
    alphabet types. The alphabet types determine default fallbacks which are
    typically used because they are reasonable, but they may be overruled by a
    developer writing a `BitsPerElem` overload specifically for their sequence
    type. 
"""
struct BitsPerElem{N} end

BitsPerElem(::A) where {A<:NucleicAcidAlphabet{2}} = BitsPerElem{2}()
BitsPerElem(::A) where {A<:NucleicAcidAlphabet{4}} = BitsPerElem{4}()
BitsPerElem(s::BioSequence) = BitsPerElem(Alphabet(s))

"""
A type trait of any bitvector-like type, that specifies whether packed elements
are stored in an unsigned integer in order from left to right, or from right to
left.

```
Left to right element order:

             1----------->N    
MSB of UInt |.............X..| LSB of UInt

Right to left element order:

                  N<--------1
MSB of UInt |.....X..........| LSB of UInt
```                 
"""
abstract type ElementOrder end
struct LeftToRight <: ElementOrder end
struct RightToLeft <: ElementOrder end

ElementOrder(::Kmer{A,K,N}) where {A,K,N} = LeftToRight()

function registertype end

registertype(::Kmer{A,K,N}) where {A,K,N} = UInt64

struct BitIndex{N,W<:Unsigned,O<:ElementOrder}
    val::Int64
end

@inline function bitindex(::BitsPerSymbol{N}, ::Type{W}, ::Type{O}, i) where {N,W,O<:ElementOrder}
    return BitIndex{N,W,O}((i - 1) << trailing_zeros(N))
end

@inline function bitindex(x, i)
    return bitindex(BitsPerSymbol(x), registertype(x), ElementOrder(x), i)
end

index_shift(i::BitIndex{N,UInt64,O}) where {N,O}   = 6
index_shift(i::BitIndex{N,UInt32,O}) where {N,O}   = 5
index_shift(i::BitIndex{N,UInt16,O}) where {N,O}   = 4
index_shift(i::BitIndex{N,UInt8,O})  where {N,O}   = 3
offset_mask(i::BitIndex{N,W,O})      where {N,W,O} = UInt8(8 * sizeof(W)) - 0x01

index(i::BitIndex)  = (i.val >> index_shift(i)) + 1
offset(i::BitIndex) = i.val & offset_mask(i)

Base.:+(i::BitIndex, n::Int) = typeof(i)(i.val + n)
Base.:-(i::BitIndex, n::Int) = typeof(i)(i.val - n)
Base.:-(i1::T, i2::T)     where {T<:BitIndex} = i1.val - i2.val
Base.:(==)(i1::T, i2::T)  where {T<:BitIndex} = i1.val == i2.val
Base.isless(i1::T, i2::T) where {T<:BitIndex} = isless(i1.val, i2.val)
Base.cmp(i1::T, i2::T)    where {T<:BitIndex} = cmp(i1.val, i2.val)

@inline function nextposition(i::BitIndex{N,W,O}) where {N,W,O}
    return i + N
end

@inline function prevposition(i::BitIndex{N,W,O}) where {N,W,O}
    return i - N
end

function Base.iterate(i::BitIndex, s = 1)
    if s == 1
        return (index(i), 2)
    elseif s == 2
        return (offset(i), 3)
    else
        return nothing
    end
end

Base.show(io::IO, i::BitIndex) = print(io, '(', index(i), ", ", offset(i), ')')

"Extract the element stored in a packed bitarray referred to by bidx."
@inline function extract_encoded_element(bidx::BitIndex{B,W,RightToLeft}, data::T) where {B,W,N,T<:Union{Vector{W},NTuple{N,W}}}
    @inbounds chunk = data[index(bidx)]
    offchunk = chunk >> offset(bidx)
    return offchunk & bitmask(bidx)
end

"Extract the element stored in a packed bitarray referred to by bidx."
@inline function extract_encoded_element(bidx::BitIndex{B,W,LeftToRight}, data::T) where {B,W,N,T<:Union{Vector{W},NTuple{N,W}}}
    @inbounds chunk = data[index(bidx)]
    offchunk = chunk >> offset(bidx)
    return offchunk & bitmask(bidx)
end

# Create a bit mask that fills least significant `n` bits (`n` must be a
# non-negative integer).
"Create a bit mask covering the least significant `n` bits."
bitmask(::Type{T}, n::Integer) where {T} = (one(T) << n) - one(T)

# Create a bit mask filling least significant N bits.
# This is used in the extract_encoded_element function.
bitmask(bidx::BitIndex{N,W}) where {N, W} = bitmask(W, N)
bitmask(n::Integer) = bitmask(UInt64, n)
bitmask(::Type{T}, ::Val{N}) where {T, N} = (one(T) << N) - one(T)


# TODO: Work out places this is used and see if it is really nessecery given the
# bitmask methods above.
# TODO: Resolve this use of bits_per_symbol and A().
bitmask(::A) where {A<:Alphabet} = bitmask(bits_per_symbol(A()))
