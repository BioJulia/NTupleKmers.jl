module BenchNaiveKmers

using BenchmarkTools
using BioSequences

suite = BenchmarkGroup()

# Number of sequences.
N = 10_000

# Length of sequences.
L = 200

# Length of Kmer.
K = 63

"Slinding window iterator."
struct Window{I,K}
    ln::Int
    xs::I

    function Window{I,K}(xs) where {I,K}

        n = length(xs)

        ln = n - K + 1

        if ln < 1
            error("K is too large.")
        end

        return new(ln, xs)
    end
end

Window(xs::I, n::Int) where {I} = Window{I,n}(xs)
Window{K}(xs::I) where {I,K} = Window{I,K}(xs)

Base.eltype(::Type{Window{I, K}}) where {I,K} = SubArray{eltype(I), 1, I, Tuple{UnitRange{Int64}}, true} # Note: applies to types like Vector{Char}.
Base.eltype(::Type{Window{I, K}}) where {I<:LongSequence, K} = I
Base.length(it::Window) = it.ln

function Base.iterate(it::Window{I,K}, state=1) where {I, K}

    if state <= it.ln
        return (view(it.xs, state:state + K - 1), state+1)
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
seqs = Base.Generator(1:N) do _
    return randdnaseq(L)
end

# # Setup an iterator to generate window iterators.
# iter = Base.Generator(seqs) do seq
#     return Window{K}(seq)
# end

s = suite["LongSequence"] = BenchmarkGroup()

seq = randdnaseq(L)

s["collect"] = @benchmarkable collect(Window{K}($seq))

# Count kmers.
s["canonical"] = @benchmarkable frequencies(canonical, Base.Iterators.flatten(iter)) setup = (iter = Base.Generator(Window{K}, $seqs))
s["non-canonical"] = @benchmarkable frequencies(Base.Iterators.flatten(iter)) setup = (iter = Base.Generator(Window{K}, $seqs))

end # module BenchNaiveKmers
BenchNaiveKmers.suite
