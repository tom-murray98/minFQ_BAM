import pysam


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

    def __init__(self, bam_file=None):
        """constructor: Initializes class with the path to the BAM"""
        self.bam_file = bam_file

    def get_rg_tags(self):
        """Extracts relevant information from read group (RG) tags."""
        rg_tags = self.sam_file.header.get("RG", [])
        for rg_tag in rg_tags:
            rg_id = rg_tag.get("ID", None)
            dt = rg_tag.get("DT", None)
            ds = rg_tag.get("DS", None)
            # DS contains bm and ri tags
            ds_tags = ds.split(" ")
            # Splits into separate variables, removes proceeding label
            bm = ds_tags[0].replace("basecall_model=", "")
            ri = ds_tags[1].replace("runid=", "")
            lb = rg_tag.get("LB", None)
            pl = rg_tag.get("PL", None)
            pm = rg_tag.get("PM", None)
            pu = rg_tag.get("PU", None)
            al = rg_tag.get("al", None)

            return rg_id, dt, bm, ri, lb, pl, pm, pu, al

    def read_bam(self):
        """Reads BAM file, iterates through information from reads and yields it."""
        if self.bam_file:
            self.sam_file = pysam.AlignmentFile(self.bam_file, "rb", check_sq=False)
            # Check if the file is aligned or unaligned, redundant??
            if self.sam_file.nreferences != 0:
                print("This file is aligned")
            else:
                print("This file is unaligned")

            # Get info from RG tags
            rg_id, dt, bm, ri, lb, pl, pm, pu, al = self.get_rg_tags()
            # Iterates through reads
            for read in self.sam_file:
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
                    rg_id,
                    dt,
                    bm,
                    ri,
                    lb,
                    pl,
                    pu,
                    pm,
                    al,
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

            print(bam_read)
            break


if __name__ == "__main__":
    # Uses ReadBam to read ban file in and yield results
    read_bam_try = ReadBam(bam_file="/home/p2solo/PycharmProjects/pythonProject/sort_ds1263_NUH7_M1.sup.meth.hg38.bam")
    # Uses .process_reads to move into library
    read_bam_try.process_reads()
