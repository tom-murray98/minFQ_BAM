import pysam
import logging


class ReadBam:
    """
        :param bam_file: A bam format alignment file.
        :return:
            name: sequence name
            seq: sequence
            qual: quality of reads
            start_time_pr: start time of run (per read)
            read_basecall_id: read ID and basecall ID
    `       channel: channel ID
            read_number: read number
            seq_to_signal :sequence to signal move table
            rg_id : read group identifier
            dt: time of run
            ds: run ID in basecall model
            bm: basecall model
            ri: run ID
            lb: sample ID
            pl: read platform (e.g ONT, PacBio)
            pm: device position
            pu: flow cell ID
            al: unknown (outputs unclassified)
    """

    def __init__(self, bam_file):
        """
        constructor: Initializes class with the path to the BAM
        """
        # self.sam_file = None
        self.bam_file = pysam.AlignmentFile(bam_file, "rb", check_sq=False)
        self.tags = self._get_rg_tags(self.bam_file)
        self.reads = {}
        ### This is going to be very slow for a big file.
        for read in self.process_reads():
            self.reads[read['read_id']] = read

    def bam_read_count(self):
        """Counts the number of reads in a BAM file."""
        return len(self.reads)

    def _get_rg_tags(self, sam_file):
        """Extracts relevant information from read group (RG) tags."""
        rg_tags = sam_file.header.get("RG", [])
        # Exits script if no RG field (for Graeme)
        if not rg_tags:
            print("This file does not contain an @RG field")
            exit()
        rg_tags_dict = {}
        for rg_tag in rg_tags:
            # print(rg_tag)
            rg_tags_dict['rg_id'] = rg_tag.get("ID", None)
            rg_tags_dict['dt'] = rg_tag.get("DT", None)
            ds = rg_tag.get("DS", None)
            # DS contains bm and ri tags
            ds_tags = ds.split(" ")
            # Splits into separate variables, removes proceeding label
            rg_tags_dict['bm'] = ds_tags[0].replace("basecall_model=", "")
            rg_tags_dict['ri'] = ds_tags[1].replace("runid=", "")
            rg_tags_dict['lb'] = rg_tag.get("LB", None)
            rg_tags_dict['pl'] = rg_tag.get("PL", None)
            rg_tags_dict['pm'] = rg_tag.get("PM", None)
            rg_tags_dict['pu'] = rg_tag.get("PU", None)
            rg_tags_dict['al'] = rg_tag.get("al", None)

            # return rg_id, dt, bm, ri, lb, pl, pm, pu, al
            return rg_tags_dict

    def read_bam(self):
        """Reads BAM file, iterates through information from reads and yields it."""

        # sam_file = pysam.AlignmentFile(self.bam_file, "rb", check_sq=False)
        sam_file = self.bam_file

        """
        # Check if the file is aligned or unaligned, redundant??

        if sam_file.nreferences != 0:
            print("This file is aligned")
        else:
            print("This file is unaligned")
        """
        # Get info from RG tags
        # rg_id, dt, bm, ri, lb, pl, pm, pu, al = self.get_rg_tags(sam_file)
        # rg_id, dt, bm, ri, lb, pl, pm, pu, al = self.tags
        # Iterates through reads
        for read in sam_file:
            name = read.query_name
            seq = read.seq
            qual = read.qual
            start_time_pr = read.get_tag("st")
            channel = read.get_tag("ch")
            read_number = read.get_tag("rn")
            read_basecall_id = read.get_tag("RG")
            # Yields a tuple containing RG and read data
            yield (
                name,
                seq,
                qual,
                start_time_pr,
                channel,
                read_number,
                read_basecall_id,
                self.tags['rg_id'],
                self.tags['dt'],
                self.tags['bm'],
                self.tags['ri'],
                self.tags['lb'],
                self.tags['pl'],
                self.tags['pu'],
                self.tags['pm'],
                self.tags['al'],
            )

    def process_reads(self):
        """Makes dictionary containing current read info"""
        for (
                name,
                seq,
                qual,
                start_time_pr,
                channel,
                read_number,
                read_basecall_id,
                rg_id,
                dt,
                bm,
                ri,
                lb,
                pl,
                pu,
                pm,
                al,
        ) in self.read_bam():
            bam_read = {
                "read_group_id": rg_id,
                "time_of_run": dt,
                "basecall_model": bm,
                "run_id": ri,
                "sample_id": lb,
                "platform": pl,
                "flowcell_id": pu,
                "device_position": pm,
                "al": al,
                "read_id": name,
                "sequence": seq,
                "sequence_length": len(str(seq)),
                "quality": qual,
                "start_time_per_run": start_time_pr,
                "channel": channel,
                "read_number": read_number,
                "read_basecall_id": read_basecall_id,
            }
            yield bam_read



if __name__ == "__main__":
    # Uses ReadBam to read ban file in and yield results
    read_bam_try = ReadBam(bam_file="test_data/230498_pass_68a2633f_385676fd_0.bam")
    # Uses .process_reads to move into library
    # read_bam_try.process_reads()
    # for read in read_bam_try.read_bam():
    #    print(read)
    # print (len(read_bam_try.reads))
    # print (read_bam_try.reads.keys())
    # print (len(read_bam_try.reads[2]['sequence']))
    # print (len(read_bam_try.reads))
    # print (read_bam_try.bam_read_count())

    """
    PSEUDOCODE:
    We want to replace: for desc, name, seq, qual in readfq(fh):
    with our BAM file handler.

    Ideally, we want to be able to do something like:

    for desc, name, seq, qual in ReadBam(bam_file=fh).process_reads():
        description_dict = parse_fastq_description(desc

    where desc is the fastq style formatted read description as a string with style RunID=XXXX ReadID=XXXX etc.

    or you could do:

    for desc, name, seq, qual in ReadBam(bam_file=fh).process_reads():
        description_dict = desc

    i.e the desc object is already a dictionary.


    """
