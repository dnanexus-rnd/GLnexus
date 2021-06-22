
import argparse
import sys
import vcf


STRICT = False
QUIET = False

def ERROR(msg):
    sys.stderr.write("ERROR " + msg + "\n")
    exit(-1)

def WARNING(msg):
    if not QUIET:
        sys.stderr.write("WARNING " + msg + "\n")
    if STRICT:
        exit(-1)

def INFO(msg):
    if not QUIET:
        sys.stderr.write("INFO " + msg + "\n")

class range:

    def __init__(self, chrom, start, end):
        self.chrom = chrom
        self.start = start
        self.end = end

    def __init__(self, vcf_record):
        self.chrom = vcf_record.CHROM
        self.start = vcf_record.POS
        self.end = self.start + len(vcf_record.REF)

    def __eq__(self, other):
        return self.chrom == other.chrom and self.start == other.start and self.end == other.end

    def __ne__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):
        if (self.chrom < other.chrom):
            return True
        if (self.start < other.start or
            (self.start == other.start and self.end < other.end)):
            return True

        return False

    def __le__(self, other):
        return (self.__eq__(other) or self.__lt__(other))

    def __gt__(self, other):
        return not self.__le__(other)

    def __ge__(self, other):
        return not self.__lt__(other)

    def __hash__(self):
        return hash((self.chrom, self.start, self.end))

    def __str__(self):
        return "<" + self.chrom + "> " + self.start + ".." + self.end

# For a given set of (INFO or FORMAT) fields, returns its subset that is defined in a list of files.
# Params:
# @fields: a list of field names, or the string "<ALL>"
# @fields_from_files: a list of sets, where each element is a set of format field names defined in the file

def vcf_fields_defined(fields, fields_from_files):

    if isinstance(fields, str):
        if fields == "<ALL>":
            # Special case to include all possible fields found in intersection of all files
            ans = set.intersection(*fields_from_files)
        else:
            # A string that we don't understand, so we return an empty set
            return set()

    # If fields is not a string, expect it to be a list of strings
    else:
        fields = set(fields)
        ans = set.intersection(fields, *fields_from_files)

        ignored_fields = set(fields) ^ ans

        if ignored_fields:
            WARNING("The following fields were not found in all files and are ignored: {0}".format(
                ", ".join(ignored_fields)))
    return ans

# Wrapper around vcf_fields_defined to find
# the info and format fields defined for comparison
# Returns a tuple of (info_fields, format_fields) where each is
# a set of strings corresponding to the INFO and FORMAT fields to be used
# for comparison
def find_cmp_fields(input_reader, truth_reader, info_f, format_f):
    input_formats = set(input_reader.formats.keys()) # reader.formats returns a dict of name: format class
    input_infos = set(input_reader.infos.keys())

    truth_formats = set(truth_reader.formats.keys())
    truth_infos = set(truth_reader.infos.keys())

    info_fields = vcf_fields_defined(info_f, [input_infos, truth_infos])
    format_fields = vcf_fields_defined(format_f, [input_formats, truth_formats])

    return (info_fields, format_fields)

# Compare the samples defined in 2 vcf readers to make sure
# they match; raising an ERROR otherwise
# Params
# @input_reader: vcf_reader opened for the input file
# @truth_reader: vcf_reader opened for the truth file
# This function is not expected to return any values
def compare_vcf_samples(input_reader, truth_reader):
    if not (set(input_reader.samples) == set(truth_reader.samples)):
        ERROR("Input and truth files have different samples! Input samples: {0}; truth samples: {1}".format(
                input_reader.samples, truth_reader.samples))

# Compare a single vcf row from input and truth files, validating the INFO
# and FORMAT field specified in info_fs and format_fs; in addition to REF
# and ALT fields.
# Params
# @i_row: a vcf row (entry) from the input file to be compared
# @t_row: a vcf row (entry) from the truth file to be compared
# @info_fs: a list of strings specifying the INFO fields to be compared
# @format_fs: a list of strings specifying the FORMAT fields to be compared
# This function assumes taht i_row and t_row contain the same range and samples
# Returns:
# None if no differences were found based on search criteria
# Tuple of (i_row, t_row) if differences were found
def compare_vcf_row(i_row, t_row, info_fs, format_fs):

    # Expect this function to be called only for rows with identical ranges
    assert(range(i_row) == range(t_row))

    # Expect comparison for identical samples
    assert [x.sample for x in i_row.samples] == [x.sample for x in t_row.samples], (str(i_row.samples) + ' / ' + str(t_row.samples))
    if (i_row.REF != t_row.REF):
        WARNING("Found differing ref! {0} {1}".format(i_row.REF, t_row.REF))
        return (i_row, t_row)

    if (i_row.ALT != t_row.ALT):
        WARNING("Found differing alts! {0} {1}".format(i_row.ALT, t_row.ALT))
        return (i_row, t_row)

    if (i_row.QUAL != t_row.QUAL):
        WARNING("Found differing qualities! {0} {1}".format(i_row.QUAL, t_row.QUAL))
        return (i_row, t_row)

    if (i_row.FILTER != t_row.FILTER):
        WARNING("Found differing filters! {0} {1}".format(i_row.FILTER, t_row.FILTER))
        return (i_row, t_row)

    for info_f in info_fs:
        n_present = sum([1 for row in (i_row, t_row) if info_f in row.INFO])

        if not n_present:
            # INFO field is missing in both truth and output
            # We consider this as agreement
            continue
        elif n_present == 1:
            # INFO field is missing in one of the files, we consider
            # this as a mismatch
            WARNING("Missing INFO field: {0} in matched rows: {1}; {2}".format(info_f, i_row, t_row))
            return (i_row, t_row)
        else:
            # INFO field is in both files, we compare the values for
            # mismatch
            if i_row.INFO[info_f] != t_row.INFO[info_f]:
                return (i_row, t_row)

    for format_f in format_fs:

        n_present = sum([1 for row in (i_row, t_row) if format_f in row.FORMAT])

        for (i_sample, t_sample) in zip(i_row.samples, t_row.samples):

            if not n_present:
                # FORMAT field is missing in both truth and output
                # We consider this as agreement
                continue
            elif n_present == 1:
                # INFO field is missing in one of the files, we consider
                # this as a mismatch
                WARNING("Missing FORMAT field: {0} in matched rows: {1}; {2}".format(format_f, i_row, t_row))
                return (i_row, t_row)
            else:
                # FORMAT field is in both files, we compare the values for
                # mismatch
                if i_sample[format_f] != t_sample[format_f]:
                    WARNING("Found differing format value for {0}: {1}, {2}".format(format_f, i_sample[format_f], t_sample[format_f]))
                    return (i_row, t_row)

    # Finished comparison without any differences noted
    return None

# Compare the vcf rows (entries) of input and truth files
# Returns a tuple of (unmatched_i_rows, unmatched_t_rows, differing rows)
# where each element is a list of rows
def compare_vcf_rows(input_reader, truth_reader, info_fields, format_fields):
    t_entries = {}

    for row in truth_reader:
        t_entries[range(row)] = row

    unmatched_i_rows = []
    unmatched_t_rows = []
    differing_rows = []

    for i_row in input_reader:
        i_range = range(i_row)
        if i_range in t_entries:
            diff = compare_vcf_row(i_row, t_entries[i_range],
                                    info_fields, format_fields)
            if diff:
                differing_rows.append(t_entries[i_range])
                differing_rows.append(i_row)

            del t_entries[i_range]

        else:
            unmatched_i_rows.append(i_row)

    for (t_range, t_row) in t_entries.items():
        unmatched_t_rows.append(t_row)

    return (unmatched_i_rows, unmatched_t_rows, differing_rows)

# Write to stdout a vcf file containing all
# unmatched and differing vcf rows between the
# input and truth file
# Params
# @diff_entries: a 3-tuple of differing rows
# @template: a vcf_reader template neede by vcf_writer for
#            writing header lines

def write_diff(diff_entries, template):
    (unmatched_i_rows, unmatched_t_rows, differing_rows) = diff_entries

    writer = vcf.Writer(sys.stdout, template)

    if unmatched_i_rows:
        print("## Unmatched vcf row(s) in Input File")
        for row in unmatched_i_rows:
            writer.write_record(row)

    if unmatched_t_rows:
        print("## Unmatched vcf row(s) in Truth File")
        for row in unmatched_t_rows:
            writer.write_record(row)

    if differing_rows:
        print("## Matched row(s) which differ in Input and Truth Files")
        for row in differing_rows:
            writer.write_record(row)

# Overall routine for comparing input and truth vcf files
def compare_vcfs(input_file, truth_file, info_f, format_f):
    input_reader = vcf.Reader(open(input_file, 'r'))
    truth_reader = vcf.Reader(open(truth_file, 'r'))

    compare_vcf_samples(input_reader, truth_reader)
    info_fields, format_fields = find_cmp_fields(input_reader, truth_reader, info_f, format_f)

    if not info_fields:
        INFO("No INFO fields will be compared!")
    if not format_fields:
        INFO("No FORMAT fields will be compared!")

    if not QUIET:
        sys.stderr.write("Performing comparison on INFO fields: {0} and FORMAT fields: {1}\n".format(
        ','.join(info_fields) if info_fields else "-None-",
        ','.join(format_fields) if format_fields else "-None-"))

    # diff_entries = unmatched_i_rows, unmatched_t_rows, differing_rows
    diff_entries = compare_vcf_rows(input_reader, truth_reader, info_fields, format_fields)

    if any(diff_entries):
        write_diff(diff_entries, input_reader)
        exit(1)
    else:
        exit(None)

def parse_args(cmd=None):
    parser = argparse.ArgumentParser(description="Semantic comparison of 2 vcf files using user-specified criteria")
    parser.add_argument("--input", "-a", help="input file to be used for comparison")
    parser.add_argument("--truth", "-b", help="truth file to be used for comparison")
    parser.add_argument("--formats", "-f", nargs='+', default="<ALL>",
                        help="list of FORMAT field(s) that are to be examined for match, use <ALL> to compare all fields, default=<ALL>")
    parser.add_argument("--infos", "-i", nargs='+', default="<ALL>",
                        help="list of INFO field(s) that are to be examined for match, use <ALL> to compare all fields, default=<ALL>")
    parser.add_argument("--strict", default=False, action="store_true",
                       help="when strict mode is activated, the program will fail on all warnings in addition to errors")
    parser.add_argument("--quiet", default=False, action="store_true",
                       help="Do not print helpful debugging messages")
    if (cmd):
        return parser.parse_args(cmd)
    else:
        return parser.parse_args()

args = parse_args()
STRICT = args.strict
QUIET = args.quiet
compare_vcfs(args.input, args.truth, args.infos, args.formats)
