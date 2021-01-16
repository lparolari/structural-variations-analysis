FLAG_UNSET = 0
FLAG_SEGMENT_UNMAPPED = 4
FLAG_NEXT_SEGMENT_UNMAPPED = 8


MATE_LENGTH = 100


def read_mates(
    input_file: str, 
    raw_fields: bool = False, 
    keep_comments: bool = False, 
    keep_fields = ["qname", "flag", "rname", "pos", "mapq", "cigar", "rnext", "pnext", "tlen", "seq", "qual", ]):
    
    """
    Return a list of dict where every entry is representing a field from the
    sam file.

    :param raw_fields: Allow to do not unpack fields and get the whole text line
    :param keep_comments: Allow to keep comments in sam file, i.e, line with `@` at the beginning
    :param keep_fields: Allot to keep only a subset of fields. By default it keeps all fields.
    """
    
    with open(input_file) as f:
        lines = f.readlines()

        if not keep_comments:
            lines = list(filter(lambda line: line[0] != "@", lines))
        
        if not raw_fields:
            def split_fields(line):
                values = line.split("\t")

                fields = {
                    **( { "qname":      values[0]  }  if "qname" in keep_fields else {}),
                    **( { "flag":   int(values[1]) }  if "flag"  in keep_fields else {}),
                    **( { "rname":      values[2]  }  if "rname" in keep_fields else {}),
                    **( { "pos":    int(values[3]) }  if "pos"   in keep_fields else {}),
                    **( { "mapq":   int(values[4]) }  if "mapq"  in keep_fields else {}),
                    **( { "cigar":      values[5]  }  if "cigar" in keep_fields else {}),
                    **( { "rnext":      values[6]  }  if "rnext" in keep_fields else {}),
                    **( { "pnext":  int(values[7]) }  if "pnext" in keep_fields else {}),
                    **( { "tlen":   int(values[8]) }  if "tlen"  in keep_fields else {}),
                    **( { "seq":        values[9]  }  if "seq"   in keep_fields else {}),
                    **( { "qual":       values[10] }  if "qual"  in keep_fields else {}),
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


def is_plausible_tlen(min_tlen = 0, max_tlen = 20000):
    """
    Returns true whether the tlen of mates is geq `min_tlen` and
    leq `max_tlen`, false otherwise.
    """

    def is_plausible_tlen_fixed(mate):
        return (mate["tlen"] >= min_tlen) and (mate["tlen"] <= max_tlen)

    return is_plausible_tlen_fixed
    

def is_first_read_exlusively_mapped(mate):
    """
    Returns true whether the first read maps exclusively, 
    false otherwise.
    """
    return (mate["flag"] & (FLAG_SEGMENT_UNMAPPED | FLAG_NEXT_SEGMENT_UNMAPPED)) == FLAG_SEGMENT_UNMAPPED


def is_second_read_exlusively_mapped(mate): 
    """
    Returns true whether the second read maps exclusively, 
    false otherwise.
    """
    return (mate["flag"] & (FLAG_SEGMENT_UNMAPPED | FLAG_NEXT_SEGMENT_UNMAPPED)) == FLAG_NEXT_SEGMENT_UNMAPPED
    

def is_first_and_second_read_mapped(mate): 
    """
    Returns true whether first and second read are mapped, 
    false otherwise.
    """
    return (mate["flag"] & (FLAG_SEGMENT_UNMAPPED | FLAG_NEXT_SEGMENT_UNMAPPED)) == FLAG_UNSET


def filter_out_invalid_mates(mates):
    """
    Remove mates with implausible tlen.
    """
    return list(filter(is_plausible_tlen(), mates))