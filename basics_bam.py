import sys
import pysam
import matplotlib.pyplot as plt
from collections import Counter


def get_bam(bam_file):
    return pysam.AlignmentFile(bam_file, "rb", check_sq=False)


def get_mapping_num_per_reads(pysam_ali):
    map_count = Counter([(x.qname, x.is_read1)
                         for x in pysam_ali
                         if not x.is_unmapped])
    return map_count


if __name__ == "__main__":

    my_bam = get_bam(sys.argv[1])
    counts = get_mapping_num_per_reads(my_bam)
    plt.hist(counts.values(), 50)

    plt.show()
