import pysam


class ReadBam:
    """
    :param bam_file: A bam format alignment file.
    :return:
        id_tag : ID
        dt_tag: time of run
        basecall_model_tag: Basecall Model ID
        runid_tag: run ID
        lb: sample ID
        pl: outputs 'ONT', device organisation?
        pm: device position
        pu: flow cell ID
        al: unknown (outputs unclassified when present)
    """

    def __init__(self, bam_file=None):
        self.sam_file = None
        self.bam_file = bam_file

    def get_rg_tags(self):
        rg_tags = self.sam_file.header.get("RG", [])
        if not rg_tags:
            print("This BAM file does not contain an @RG field")
            exit()
        else:
            for rg_tag in rg_tags:
                id_tag = rg_tag.get("ID", None)
                dt_tag = rg_tag.get("DT", None)
                ds_tag = rg_tag.get("DS", None)
                ds_tags = ds_tag.split(" ")
                basecall_model_tag = ds_tags[0].replace("basecall_model=", "")
                runid_tag = ds_tags[1].replace("runid=", "")
                lb_tag = rg_tag.get("LB", None)
                pl_tag = rg_tag.get("PL", None)
                pm_tag = rg_tag.get("PM", None)
                pu_tag = rg_tag.get("PU", None)
                al_tag = rg_tag.get("al", None)

                return (
                    id_tag,
                    dt_tag,
                    basecall_model_tag,
                    runid_tag,
                    lb_tag,
                    pl_tag,
                    pm_tag,
                    pu_tag,
                    al_tag,
                )

    def read_bam(self):
        if self.bam_file:
            self.sam_file = pysam.AlignmentFile(self.bam_file, "rb", check_sq=False)

            (
                id_tag,
                dt_tag,
                basecall_model_tag,
                runid_tag,
                lb_tag,
                pl_tag,
                pm_tag,
                pu_tag,
                al_tag,
            ) = self.get_rg_tags()

            # for read in self.sam_file:
            #     name = read.query_name
            #     seq = read.seq
            #     qual = read.qual
            #     start_time_pr = read.get_tag("st")
            #     channel = read.get_tag("ch")
            #     read_number = read.get_tag("rn")
            #     read_basecall_id = read.get_tag("RG")

            yield (
                # name,
                # seq,
                # qual,
                # start_time_pr,
                # channel,
                # read_number,
                # read_basecall_id,
                id_tag,
                dt_tag,
                basecall_model_tag,
                runid_tag,
                lb_tag,
                pl_tag,
                pm_tag,
                pu_tag,
                al_tag,
            )

    def process_reads(self):
        for (
            # name,
            # seq,
            # qual,
            # start_time_pr,
            # channel,
            # read_number,
            # read_basecall_id,
            id_tag,
            dt_tag,
            basecall_model_tag,
            runid_tag,
            lb_tag,
            pl_tag,
            pm_tag,
            pu_tag,
            al_tag,
        ) in self.read_bam():
            bam_read = {
                "ID": id_tag,
                "time_of_run": dt_tag,
                "sample_id": lb_tag,
                "basecall_model": basecall_model_tag,
                "runid": runid_tag,
                "platform": pl_tag,
                "flow_cell_id": pu_tag,
                "device_position": pm_tag,
                "al": al_tag,
                # "read_id": name,
                # "sequence": seq,
                # "sequence_length": len(str(seq)),
                # "quality": qual,
                # "start_time_per_run": start_time_pr,
                # "channel": channel,
                # "read_number": read_number,
                # "read_basecall_id": read_basecall_id,
            }

            print(bam_read)


if __name__ == "__main__":
    # Uses ReadBam to read ban file in and yield resultsds1263_single_OLD.bam
    read_bam_try = ReadBam(bam_file="230498_pass_68a2633f_385676fd_0.bam")
    # Uses .process_reads to move into library and print
    read_bam_try.process_reads()
