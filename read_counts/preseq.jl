#= 
    preseq.jl works in concert with the preseq toolbox to compute library complexity 
    as a function of sequence read count. It parallelizes the computation of read counts,
    which is the most expensive part of the complexity curve calculation in the original
    preseq toolbox
=#

using ArgParse

# Path to preseq executable
preseq_bin  = "preseq/preseq"
# Arguments to preseq (perhaps these should eventually be args to preseq.jl)
preseq_func = "c_curve"
preseq_s    = 10000

# Parse command line arguments
function parse_arguments()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--bam"
            help = "Input BAM file"
        "--output"
            help = "Output complexity file"
    end

    return parse_args(s)
end

# Essentially a Julia implementation of `uniq -c`
function unique_occurences(values)
    counts = []
    for value in values
        if size(counts, 1) == 0
            counts = [value 1]
        elseif counts[end, 1] != value
            counts = vcat(counts, [value 1])
        else
            counts[end, 2] += 1
        end
    end
    return counts
end

#function merge_counts(count_arrays)
#    # For a given count_arrays[i], count_arrays[i + 1]
#    # If last item in count_arrays[i] is identical to first item in count_arrays[i + 1]
#    # then sum up the two counts, update count_arrays[i + 1]'s count, and delete last row
#    # of count_arrays[i]
#    for count_array in count_arrays
#        println(count_array)
#    end
#end
#
#test_data = [1 1 2 2 2 2 2 3 3 4 5 6 7 7 8 8 8 8 9 9 9 10 11 11 11 12 13 14 14]
#
#c = Any[]
#push!(c, unique_occurences(test_data[1:10]))
#push!(c, unique_occurences(test_data[11:20]))
#push!(c, unique_occurences(test_data[21:end]))
#merge_counts(c)
#println(unique_occurences(test_data))
#exit()

# Combination of bitwise flags are used for determining properties of a read
is_primary(flag)                            = ~(flag & 0x100 != 0)
is_mapped(flag)                             = ~(flag & 0x4 != 0)
is_pairend(flag)                            = (flag & 0x1 != 0)
is_Trich(flag)                              = (is_pairend(flag) ? (flag & 0x40 != 0) : true)
is_mapping_paired(flag, mate_name, chrom)   = (((flag & 0x2 != 0) && ~((mate_name == "=") || (mate_name == chrom))) ? false : (flag & 0x2 != 0))

# Process command line arguments
bam_filename, bai_filename, output_filename = "", "", ""
parsed_args = parse_arguments()
for pa in parsed_args
    if pa[1] == "bam"
        bam_filename = pa[2]
        bai_filename = string(pa[2], ".bai")
    elseif pa[1] == "output"
        output_filename = pa[2]
    end
end

# Check if given BAM file and a corresponding index (BAI) file exist
if ~isfile(bam_filename)
    error("BAM file does not exist. Check file path and try again")
end

if ~isfile(bai_filename)
    error("BAI file does not exist. Make sure an index exists for the BAM file, or generate one by using: samtools index <filename.bam>")
end

tic()

# Determine number of entries in BAM file (using its index file)
idxstats = split(readall(`samtools idxstats $bam_filename`), '\n')
num_mapped, num_unmapped = 0, 0
for entry in idxstats
    fields          = split(entry)
    if size(fields, 1) > 0
        num_mapped      += parse(Int, fields[3])
        num_unmapped    += parse(Int, fields[4])
    end
end

# Read BAM file using samtools
starts = []
lines = split(readall(`samtools view $bam_filename`), '\n')
for line in lines
    # Read relevant fields from SAM file entry 
    fields      = split(line)
    if size(fields, 1) > 0
        flag        = parse(Int, fields[2])
        chrom       = fields[3]
        pos         = parse(Int, fields[4])
        mate_name   = fields[7]

        # Filter out secondary and unmapped reads
        if (is_primary(flag) && is_mapped(flag))
            # Additional filters
            if (~(is_mapping_paired(flag, mate_name, chrom)) || (is_mapping_paired(flag, mate_name, chrom) && is_Trich(flag)))
                push!(starts, pos + 1)
            end
        end
    end
end

toc()

# Write read counts to temporary file
count_filename = string(bam_filename, ".counts")
writedlm(count_filename, unique_occurences(starts)[1:end-1, 2])

# Use preseq to compute complexity from read counts
run(`./$preseq_bin $preseq_func -o $output_filename -V $count_filename -s $preseq_s`)

# Remove read counts file
rm(count_filename)

