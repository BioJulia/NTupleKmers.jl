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
BitsPerElem(k::Kmer) = BitsPerElem(Alphabet(k))

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

@inline function bitindex(::BitsPerElem{N}, ::Type{W}, ::O, i::Integer) where {N,W,O<:ElementOrder}
    return BitIndex{N,W,O}((i - 1) << trailing_zeros(N))
end

@inline function bitindex(x, i::Integer)
    return bitindex(BitsPerElem(x), registertype(x), ElementOrder(x), i)
end

"Compute the bitindex for the `i'th` base in `seq`'s bitvector-like datastore."
@inline function bitindex(seq::Kmer, i::Integer)
    return bitindex(BitsPerElem(seq), registertype(seq), ElementOrder(seq), i + n_unused(seq))
end

# TODO: Decide whether to define these in terms of firstindex and lastindex,
# instead of 1 and length(seq).
@inline firstbitindex(seq::Kmer) = bitindex(seq, 1 + n_unused(seq))
@inline lastbitindex(seq::Kmer)  = bitindex(seq, length(seq) + n_unused(seq))

@inline function eachbitindex(from::BitIndex{N,W,O}, to::BitIndex{N,W,O}) where {N,W,O}
    return from:N:to
end

@inline function eachbitindex(seq::Kmer, from = 1, to = length(seq))
    return eachbitindex(from, to)
end

index_shift(i::BitIndex{N,UInt64,O}) where {N,O} = 6
index_shift(i::BitIndex{N,UInt32,O}) where {N,O} = 5
index_shift(i::BitIndex{N,UInt16,O}) where {N,O} = 4
index_shift(i::BitIndex{N,UInt8,O})  where {N,O} = 3

offset_mask(i::BitIndex{N,W,O}) where {N,W,O} = UInt8(8 * sizeof(W)) - 0x01

"Get the index of the register the element pointed to by BitIndex `i`"
register(i::BitIndex) = (i.val >> index_shift(i)) + 1

"""
    offset(i::BitIndex{N,W,RightToLeft}) where {N,W}

Get the size of the right-shift of the register required to bring the element
pointed to by BitIndex `i`, to the `N` least significant digits.

```
               Before:
MSB of UInt |.....X..........| LSB of UInt
                  |--------->| Size of the right-shift required = offset
MSB of UInt |...............X| LSB of UInt
               After:
```
"""
@inline offset(i::BitIndex{N,W,RightToLeft}) where {N,W} = i.val & offset_mask(i)

"""
    offset(i::BitIndex{N,W,LeftToRight}) where {N,W}

Get the size of the right-shift of the register required to bring the element
pointed to by BitIndex `i`, to the `N` least significant digits.

```
               Before:
MSB of UInt |.....X..........| LSB of UInt
                  |--------->| Size of the right-shift required = offset
MSB of UInt |...............X| LSB of UInt
               After:
```
"""
@inline offset(i::BitIndex{N,W,LeftToRight}) where {N,W} = (8 * sizeof(W) - N) - (i.val & offset_mask(i))

"Create a bit mask covering the least significant `n` bits."
bitmask(bidx::BitIndex{N,W,O}) where {N,W,O} = (one(W) << N) - one(W)

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
        return (register(i), 2)
    elseif s == 2
        return (offset(i), 3)
    else
        return nothing
    end
end

Base.show(io::IO, i::BitIndex) = print(io, '(', register(i), ", ", offset(i), ')')

"Extract the element stored in a packed bitarray referred to by bidx."
@inline function extract_encoded_element(bidx::B, data) where {B<:BitIndex}
    @inbounds chunk = data[register(bidx)]
    offchunk = chunk >> offset(bidx)
    return offchunk & bitmask(bidx)
end






