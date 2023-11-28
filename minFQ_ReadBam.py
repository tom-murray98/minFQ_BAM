import pysam
import click


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
        """ Initializes bam_file and tags as instances """
        # Defines self.bam_file as an AlignmentFile object for reading input BAM file
        self.bam_file = pysam.AlignmentFile(bam_file, 'rb', check_sq=False)
        self.tags = self.get_rg_tags(self.bam_file)

    def get_rg_tags(self, bam_file):
        """ Detects and retrieves @RG header tags from a BAM file and assigns each one to a dictionary """
        # Gets @RG from header
        rg_tags = bam_file.header.get('RG', [])
        # Exits script if no @RG field present
        if not rg_tags:
            print('This file does not contain an @RG field')
            exit()
        # Only return the first RG tag for ease
        rg_tag = rg_tags[0]
        # Creates dictionary of @RG header tags and their values
        rg_tags_dict = {
            'rg_id': rg_tag.get('ID', None),
            'start_time': rg_tag.get('DT', None),
            'basecall_model_version_id': rg_tag.get('DS', '').split(' ')[1].replace('basecall_model=', ''),
            'run_id': rg_tag.get('DS', '').split(' ')[0].replace('runid=', ''),
            'sample_id': rg_tag.get('LB', None),
            'platform': rg_tag.get('PL', None),
            'position': rg_tag.get('PM', None),
            'flow_cell_id': rg_tag.get('PU', None),
            'al': rg_tag.get('al', None)
        }
        return rg_tags_dict

    def read_bam(self):
        """ Extracts read data from a BAM file iteratively and yields in a FASTQ friendly format"""
        for read in self.bam_file:
            name = read.query_name
            seq = read.seq
            qual = read.qual
            start_time_pr = read.get_tag('st')
            channel = read.get_tag('ch')
            read_number = read.get_tag('rn')
            read_basecall_id = read.get_tag('RG')

            # Combines read and @RG metrics to create description variable, using the same format as FASTQ files
            desc = '@' + str(name) + ' runid=' + str(self.tags['run_id']) + ' read=' + str(read_number) + ' ch=' + str(
                channel) + ' start_time=' + str(start_time_pr) + ' flow_cell_id=' + str(
                self.tags['flow_cell_id']) + ' protocol_group_id=' + 'UNKNOWN' + ' sample_id=' + str(
                self.tags['sample_id']) + ' parent_read_id=' + str(name) + ' basecall_model_version_id=' + str(
                self.tags['basecall_model_version_id'])
            # Yields FASTQ metrics
            yield desc, name, seq, qual


# This function was taken from https://github.com/LooseLab/minotourcli/blob/master/minFQ/fastq_handler_utils.py#L137
# Not intended for permanent use, testing purpose only
def parse_fastq_description(description):
    description_dict = {}
    descriptors = description.split(" ")
    for item in descriptors:
        if "=" in item:
            bits = item.split("=")
            description_dict[bits[0]] = bits[1]
    return description_dict


""" Click decorators for creation of file handle input in the command line, as well as an @RG tag filter """


@click.command()
@click.argument('fh', type=click.Path(exists=True))
@click.option('--rg', is_flag=True, show_default=True, default=False,
              help="Prints @RG tags only found in the BAM file header.")
def main(fh, rg):
    for desc, name, seq, qual in ReadBam(bam_file=fh).read_bam():
        if rg:
            var = ReadBam(bam_file=fh).tags
            click.echo(var)
            break

        else:
            description_dict = parse_fastq_description(desc)
            click.echo(description_dict)


if __name__ == '__main__':
    main()

