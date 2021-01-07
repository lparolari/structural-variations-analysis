FLAG_UNSET = 0
FLAG_SEGMENT_UNMAPPED = 4
FLAG_NEXT_SEGMENT_UNMAPPED = 8


def read_mates(input_file: str, raw_fields: bool = False, keep_comments: bool = False):
    """
    Return a list of dict where every entry is representing a field from the
    sam file.

    :param raw_fields: Allow to do not unpack fields and get the whole text line
    :param keep_comments: Allow to keep comments in sam file, i.e, line with `@` at the beginning
    """
    with open(input_file) as f:
        lines = f.readlines()

        if not keep_comments:
            lines = list(filter(lambda line: line[0] != "@", lines))
        
        if not raw_fields:
            def split_fields(line):
                values = line.split("\t")

                fields = {
                    "qname": values[0],
                    "flag": int(values[1]),
                    "rname": values[2],
                    "pos": int(values[3]),
                    "mapq": int(values[4]),
                    "cigar": values[5],
                    "rnext": values[6],
                    "pnext": int(values[7]),
                    "tlen": int(values[8]),
                    "seq": values[9],
                    "qual": values[10]
                }

                return fields

            lines = list(map(split_fields, lines))
        
        return lines



def to_wig(ls, step_type: str = "fixedStep", chrom: str = "genome", start: int = 1, step: int = 1, span: int = 1):
    """
    Prints to the standard output a list in wig format.
    """
    print(f"{step_type} chrom={chrom} start={start} step={step} span={span}")

    for item in ls:
        print(item)