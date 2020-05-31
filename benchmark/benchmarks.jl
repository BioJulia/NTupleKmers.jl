using BenchmarkTools, NTupleKmers, BioSequences

SUITE = BenchmarkGroup()

SUITE["creation"] = BenchmarkGroup(["constructors"])

# Add some benchmarks to the "creation" group
dnaseq = randdnaseq(63)
SUITE["creation"]["63mer"] = @benchmarkable DNAKmer{63,2}($dnaseq)










dnaseq = LongSequence{DNAAlphabet{2}}(randdnaseq(63))
m = DNAKmer{63, 2}(dnaseq)
oldm = BigDNAMer{63}(dnaseq)

@benchmark shiftright($m)
@benchmark shiftright($oldm)