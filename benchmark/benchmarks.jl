using BenchmarkTools, NTupleKmers, BioSequences

SUITE = BenchmarkGroup()

SUITE["creation"] = BenchmarkGroup(["constructors"])

# Add some benchmarks to the "creation" group
dnaseq = randdnaseq(63)
SUITE["creation"]["63mer"] = @benchmarkable DNAKmer{63,2}($dnaseq)











onem   = NTupleKmers.DNAKmer{31, 1}(LongSequence{DNAAlphabet{2}}(randdnaseq(31)))
twom   = NTupleKmers.DNAKmer{63, 2}(LongSequence{DNAAlphabet{2}}(randdnaseq(63)))
threem = NTupleKmers.DNAKmer{95, 3}(LongSequence{DNAAlphabet{2}}(randdnaseq(95)))
fourm  = NTupleKmers.DNAKmer{127, 4}(LongSequence{DNAAlphabet{2}}(randdnaseq(127)))


m = NTupleKmers.DNAKmer{63, 2}(dnaseq)
oldm = BigDNAMer{63}(dnaseq)

@benchmark shiftright($m)
@benchmark shiftright($oldm)