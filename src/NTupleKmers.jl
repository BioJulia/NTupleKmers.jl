module NTupleKmers

export
    Kmer,
    DNAKmer,
    DNAKmer27,
    DNAKmer31,
    DNAKmer63,
    kmertype,
    capacity,
    n_unused
    

using BioSequences

import BioSequences: twobitnucs, inbounds_getindex

###
### Type definition
###

# In BioSequences.jl kmers are called Mer, not Kmer so there is no name clash.
"""
    Kmer{A<:NucleicAcidAlphabet{2},K,N} <: BioSequence{A}

A parametric, immutable, bitstype for representing Kmers - short sequences.

Given the number of Kmers generated from raw sequencing reads, avoiding
repetetive memory allocation and triggering of garbage collection is important,
as is the ability to effectively pack Kmers into arrays and similar collections.

In julia this means an immutable bitstype must represent such shorter Kmer
sequences. Thankfully this is not much of a limitation - kmers are rarely
manipulated and so by and large don't have to be mutable like `LongSequence`s.

Excepting their immutability, they fulfill the rest of the API and behaviours
expected from a concrete `BioSequence` type.

!!! warning
    Given their immutability, `setindex` and mutating sequence transformations
    are not implemented for kmers e.g. `reverse_complement!`. 

!!! tip
    Note that some sequence transformations that are not mutating are
    available, since they can return a new kmer value as a result e.g.
    `reverse_complement`. 
"""
struct Kmer{A<:NucleicAcidAlphabet{2},K,N} <: BioSequence{A}
    data::NTuple{N,UInt64}
end


###
### Shortcuts
###

"Shortcut for the type `Kmer{DNAAlphabet{2},K,N}`"
const DNAKmer{K,N} = Kmer{DNAAlphabet{2},K,N}

"Shortcut for the type `DNAKmer{27,1}`"
const DNAKmer27 = DNAKmer{27,1}

"Shortcut for the type `DNAKmer{31,1}`"
const DNAKmer31 = DNAKmer{31,1}

"Shortcut for the type `DNAKmer{63,2}`"
const DNAKmer63 = DNAKmer{63,2}


###
### Base Functions
###

"""
    kmertype(::Type{Kmer{A,K}}) where {A,K}

Fully resolve and incomplete kmer typing, computing the N parameter of
`Kmer{A,K,N}`, given only `Kmer{A,K}`.

## Example

```julia
julia> DNAKmer{63}
Kmer{DNAAlphabet{2},63,N} where N

julia> kmertype(DNAKmer{63})
Kmer{DNAAlphabet{2},63,2}

```
"""
@inline function kmertype(::Type{Kmer{A,K}}) where {A,K}
    return Kmer{A,K,ifelse(rem(2K, 64) != 0, div(2K, 64) + 1, div(2K, 64))}
end

@inline capacity(::Type{Kmer{A,K,N}}) where {A,K,N} = div(64N, 2)
@inline n_unused(::Type{Kmer{A,K,N}}) where {A,K,N} = capacity(Kmer{A,K,N}) - K

@inline function checkmer(::Type{Kmer{A,K,N}}) where {A,K,N}
    c = capacity(Kmer{A,K,N})
    if !(1 ≤ K ≤ c)
        throw(ArgumentError("K must be within 1..$c"))
    end
end

@inline Base.eltype(::Type{Kmer{A,K,N}}) where {A,K,N} = eltype(A)
@inline Base.length(x::Kmer{A,K,N}) where {A,K,N} = K
@inline Base.summary(x::Kmer{DNAAlphabet{2},K,N}) where {K,N} = string("DNA ", K, "-mer")
@inline Base.summary(x::Kmer{RNAAlphabet{2},K,N}) where {K,N} = string("RNA ", K, "-mer")

function Base.typemin(::Type{Kmer{A,K,N}}) where {A,K,N}
    checkmer(T)
    return T(ntuple(i -> zero(UInt64), N))
end

function Base.typemax(::Type{Kmer{A,K,N}}) where {A,K,N}
    checkmer(Kmer{A,K,N})
    return Kmer{A,K,N}((typemax(UInt64) >> (64N - 2K), ntuple(i -> typemax(UInt64),N - 1)...))
end

function Base.rand(::Type{Kmer{A,K,N}}) where {A,K,N}
    return Kmer{A,K,N}(ntuple(i -> rand(UInt64), N))
end


###
### Constructors
###

# Create a Mer from a sequence.
function (::Type{Kmer{A,K,N}})(seq::LongSequence{A}) where {A<:NucleicAcidAlphabet{2},K,N}
    seqlen = length(seq)
    if seqlen != K
        throw(ArgumentError("seq does not contain the correct number of nucleotides ($seqlen ≠ $K)"))
    end
    if seqlen > capacity(Kmer{A,K,N})
        throw(ArgumentError("Cannot build a mer longer than $(capacity(Kmer{A,K,N}))bp long"))
    end

    # Construct the head.
    bases_in_head = div(64 - (64N - 2K), 2)
    head = zero(UInt64)
    @inbounds for i in 1:bases_in_head
        nt = convert(eltype(typeof(seq)), seq[i])
        head = (head << 2) | UInt64(twobitnucs[reinterpret(UInt8, nt) + 0x01])
    end
    
    # And the rest of the sequence
    idx = Ref(bases_in_head + 1)
    
    tail = ntuple(Val{N - 1}()) do i
        Base.@_inline_meta
        body = zero(UInt64)
        @inbounds for i in 1:32
            nt = convert(eltype(typeof(seq)), seq[idx[]])
            body = (body << 2) | UInt64(twobitnucs[reinterpret(UInt8, nt) + 0x01])
            idx[] += 1
        end
        return body
    end

    return Kmer{A,K,N}((head, tail...))
end

function (::Type{BigMer{A,K}})(seq::LongSequence{A}) where {A<:NucleicAcidAlphabet{2},K}
    seqlen = length(seq)
    if seqlen != K
        throw(ArgumentError("seq does not contain the correct number of nucleotides ($seqlen ≠ $K)"))
    end
    if seqlen > BioSequences.capacity(BigMer{A,K})
        throw(ArgumentError("Cannot build a mer longer than $(BioSequences.capacity(BigMer{A,K}))bp long"))
    end

    x = zero(BioSequences.encoded_data_type(BigMer{A,K}))
    for c in seq
        nt = convert(eltype(BigMer{A,K}), c)
        x = (x << 2) | BioSequences.encoded_data_type(BigMer{A,K})(twobitnucs[reinterpret(UInt8, nt) + 0x01])
    end

    return BigMer{A,K}(x)
end


###
### Indexing
###

@inline function inbounds_getindex(x::Kmer{A,K,N}, i::Integer) where {A,K,N}
    # Emulation of BitIndex type
    i′ = i + div(64N - 2K, 2)
    val = (i′ - 1) << 1
    index = (val >> 6) + 1
    offset = 62 - (val & (UInt8(64) - 0x01))
    @inbounds begin
        chunk = x.data[index]
    end
    bits = (chunk >> offset) & UInt64(3)
    return reinterpret(eltype(x), 0x01 << bits)
end




@inline function choptail(x::NTuple{N,UInt64}) where {N}
    ntuple(Val{N - 1}()) do i
        Base.@_inline_meta
        return @inbounds x[i]
    end
end

@inline function setlast(x::NTuple{N,UInt64}, nt::DNA) where {N}
    @inbounds begin
        bits = UInt64(twobitnucs[reinterpret(UInt8, nt) + 0x01])
        tail = (x[N] & (typemax(UInt64) - UInt64(3))) | bits
    end
    return (choptail(x)..., tail)
end

"""
    shiftright

It is important to be able to efficiently shift the all the nucleotides in a kmer
one space to the left or right, as it is a key operation in iterating through
de bruijn graph neighbours or in building kmers a nucleotide at a time.
"""
function shiftright end

@inline shiftright(x::BigDNAMer{K}) where {K} = BigDNAMer{K}(reinterpret(UInt128, x) >> 2)

@inline function shiftright(x::Kmer{A,K,N}) where {A,K,N}
    return Kmer{A,K,N}(_shiftright(zero(UInt64), x.data...))
end

@inline function _shiftright(carry::UInt64, head::UInt64, tail...)
    return ((head >> 2) | carry, _shiftright((head & UInt64(3)) << 62, tail...)...)
end
@inline _shiftright(carry::UInt64) = ()


"""
    shiftleft

It is important to be able to efficiently shift the all the nucleotides in a kmer
one space to the left or right, as it is a key operation in iterating through
de bruijn graph neighbours or in building kmers a nucleotide at a time.
"""
function shiftleft end

@inline shiftleft(x::BigDNAMer{K}) where {K} = BigDNAMer{K}(reinterpret(UInt128, x) << 2)

@inline function shiftleft(x::Kmer{A,K,N}) where {A,K,N}
    _, newbits = _shiftleft(x.data...)
    # TODO: The line below is a workaround for julia issues #29114 and #3608
    newbits′ = newbits isa UInt64 ? (newbits,) : newbits
    return Kmer{A,K,N}(_cliphead(64N - 2K, newbits′...))
end

@inline function shiftleft(x::NTuple{N,UInt64}) where {N}
    _, newbits = _shiftleft(x...)
    # TODO: The line below is a workaround for julia issues #29114 and #3608
    return newbits isa UInt64 ? (newbits,) : newbits
end

@inline function _cliphead(by::Integer, head::UInt64, tail...)
    return (head & (typemax(UInt64) >> by), tail...)
end

@inline function _shiftleft(head::UInt64, tail...)
    carry, newtail = _shiftleft(tail...)
    # TODO: The line below is a workaround for julia issues #29114 and #36087
    newtail′ = newtail isa UInt64 ? (newtail,) : newtail 
    return head >> 62, ((head << 2) | carry, newtail′...)
end

@inline _shiftleft(head::UInt64) = (head & 0xC000000000000000) >> 62, head << 2



#=
@inline function shiftleft2(x::Kmer{A,K,N}) where {A,K,N}
    return _cliphead(64N - 2K, _shiftleft2(x.data...)...)
end

@inline function _shiftleft2(head::UInt64, next::UInt64, tail...)
    return (head << 2 | (next >> 62), _shiftleft2(next, tail...)...)
end
@inline _shiftleft2(head::UInt64) = head << 2
=#





end # module