module BenchNaiveKmers

using BenchmarkTools
using BioSequences
using BioSymbols

suite = BenchmarkGroup()

# Number of sequences.
N = 10_000

# Length of sequences.
L = 200

# Length of Kmer.
K = 63

"""
High level bit packing - lots to do here.
Just following TwoBit for the moment.
T: 00, C: 01, A: 10, G: 11
"""
function Base.convert(::Type{BitVector}, nt::DNA) #TODO: DNAAlphabet{2}

    if DNA_T == nt
        return BitVector((false,false))
    end

    if DNA_C == nt
        return BitVector((false,true))
    end

    if DNA_A == nt
        return BitVector((true,false))
    end

    if DNA_G == nt
        return BitVector((true,true))
    end

    error("Unknown.")
end

"High level bit packaging - lots to do here."
function Base.convert(::Type{BitVector}, nts::LongSequence) #TODO: DNAAlphabet{2}

    n = length(nts) * 2 #TODO: use known letter/symbol size instead of bit vector lengths.

    bs = BitVector(undef,n)
    i = 1
    for nt in nts
        bs[i:i+1] = convert(BitVector, nt)
        i = i + 2
    end
    return bs
end

"Slinding window iterator."
struct Window{I,K,S}
    n::Int
    m::I # Base.ReshapedArray{Bool,S,BitArray{1},Tuple{}}

    function Window{I,K,S}(m::I) where {I,K,S}
        n = size(m,2) - K + 1 # n kmers (works with vectors)
        return new{I,K,S}(n, m)
    end
end

function Window{K,S}(xs::I) where {K,S,I<:BitArray{1}} #TODO: The letter size S may be contained in I.
    m = reshape(xs, S, :)
    return Window{typeof(m),K,S}(m)
end

# Window(xs::I, n::Int) where {I} = Window{I,n}(xs)
Window{K}(xs) where {K} = Window{K,1}(xs)
Window{K,S}(xs::I) where {I,K,S} = Window{I,K,S}(xs)
Window{I,K,S}(xs) where {I,K,S} = Window{K,S}(convert(I, xs)) # Note converting early for correct n calculation.

Base.eltype(::Type{Window{I,K,S}}) where {I,K,S} = SubArray{eltype(I), 1, I, Tuple{UnitRange{Int64}}, true} # Note: applies to types like Vector{Char}.
Base.eltype(::Type{Window{I,K,S}}) where {I<:LongSequence,K,S} = I
Base.length(it::Window) = it.n


Base.size(seq::LongSequence, d) = length(seq) #TODO: this is a quick fix.

# Note: ReshapedArrayIterator.
function Base.iterate(it::Window{I,K,S}, state=0) where {I<:Base.ReshapedArray,K,S} #Note: zero based state (state behaviour matches eachcol iterator).

    if state < it.n
        next_state = state+1
        return (view(it.m, :, next_state:state + K), next_state)
    end

    return nothing
end

function Base.iterate(it::Window{I,K,S}, state=0) where {I,K,S} #Note: zero based state (state behaviour matches eachcol iterator).

    if state < it.n
        next_state = state+1
        return (view(it.m, next_state:state + K), next_state) #Note: assuming vector.
    end

    return nothing
end

"Frequency accumulator."
function frequencies(xs, fs = Dict{eltype(xs),Int}())
  for x in xs
    fs[x] = get(fs, x, 0) + 1
  end
  return fs
end

"Frequency accumulator with key change."
function frequencies(f::Function, xs, fs = Dict{eltype(xs),Int}())
  for x in xs
    fs[f(x)] = get(fs, x, 0) + 1
  end
  return fs
end

# Setup sequence generator/loader (keep isolated from benchmarks).
seqgen = Base.Generator(1:N) do _
    return randdnaseq(L)
end

let s = suite["LongSequence"] = BenchmarkGroup()

    seq = randdnaseq(L)

    s["collect"] = @benchmarkable collect(Window{K}($seq))

    # Load compatible seqs into memory for test.
    seqs = collect(seqgen)

    # Setup an iterator to generate window iterators.
    # iter = Base.Generator(seqs) do seq #TODO: write custom generator. Type is lost here.
    #     return Window{K}(seq)
    # end

    # Count kmers.
    s["canonical"] = @benchmarkable frequencies(canonical, Base.Iterators.flatten(iter)) setup=(iter = Base.Generator(Window{K}, $seqs))
    s["non-canonical"] = @benchmarkable frequencies(Base.Iterators.flatten(iter)) setup=(iter = Base.Generator(Window{K}, $seqs))

end

let s = suite["BitVector"] = BenchmarkGroup()

    seq = randdnaseq(L)

    s["collect"] = @benchmarkable collect(Window{K,2}($seq))

    # Load compatible seqs into memory for test.
    seqs = collect(BitVector, seqgen)

    # Setup an iterator to generate window iterators.
    # iter = Base.Generator(seqs) do seq #TODO: write custom generator. Type is lost here.
    #     return Window{K,2}(seq) #TODO: get letter size from alphabet type.
    # end

    # Count kmers.
    s["non-canonical"] = @benchmarkable frequencies(Base.Iterators.flatten(iter)) setup=(iter = Base.Generator(Window{K,2}, $seqs))

end

end # module BenchNaiveKmers
BenchNaiveKmers.suite

# # function loader() :: BitVector
# #     return randdnaseq(L)
# # end
#
# seq = randdnaseq(10)
#
# # k = 3
# # 1: 123-------
# # 2: -234------
# # 3: --345-----
# # 4: ---456----
# # 5: ----567---
# # 6: -----678--
# # 7: ------789-
# # 8: -------890
#
# Window{3}(seq)
# Window{BitVector,3,2}(seq)
#
#
# seqs = fill(randdnaseq(20), 1000)
#
# iter1 = Base.Generator(seqs) do seq
#     return Window{3}(seq)
# end
#
# fs1 = frequencies(Base.Iterators.flatten(iter1))
#
# iter2 = Base.Generator(seqs) do seq
#     return Window{BitVector,3,2}(seq)
# end
#
# fs2 = frequencies(Base.Iterators.flatten(iter2))
#
# sort(collect(BitVector, keys(fs1))) == sort(collect(keys(fs2)))
