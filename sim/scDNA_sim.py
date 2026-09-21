import argparse
from Bio import SeqIO
import argparse
import re
import numpy as np
import pandas as pd
import os
import random
import subprocess
from pathlib import Path
from bisect import bisect_right


"""
Inputs:

1. A file containing phylogeny and copy number events (with paths to ecDNA CN vectors if needed)
2. Reference genome files (paternal and maternal)


"""

parser = argparse.ArgumentParser()


parser.add_argument(
    "--mat",
    type = str,
    help = "Path to maternal genome"
)

parser.add_argument(
    "--pat",
    type = str,
    help = "Path to paternal genome (doesn't really matter which is paternal or maternal)"
)

parser.add_argument(
    "--cnv",
    type = str,
    help = "Path to file containing CNV information phylogenetic tree"
)

parser.add_argument(
    "--ado-freq",
    type = float,
    default = 0.00001,
    help = "Number of ADO events expected per base"
)


parser.add_argument(
    "--ado-mean",
    type = float,
    default = 3,
    help = "Average ADO size (log 10)"
)


parser.add_argument(
    "--ado-var",
    type = float,
    default = 1,
    help = "ADO Variance (log 10)"
)


parser.add_argument(
    "--coverage-mean",
    type = float,
    default = 0.1,
    help = "Mean of coverage per cell"
)

parser.add_argument(
    "--coverage-var",
    type = float,
    default = 0.5,
    help = "Variance of coverage per cell"
)

parser.add_argument(
    "--out",
    type = str,
    default = './out',
    help = "Path to output folder"
)


args = parser.parse_args()

########################################################

def run_sim(args) :
    os.makedirs(args.out, exist_ok = True)
    os.makedirs(f"{args.out}/dwgsim/", exist_ok = True)
    os.makedirs(f"{args.out}/cn_mat/", exist_ok = True)

    maternal, paternal = read_fasta(args.mat, args.pat)

    # First read: get the tree
    with open(args.cnv, 'r') as file :
        tree_str = file.readline()

        print("Parsed tree:")
        tree_dict_list = parse_subnodes(tree_str)
        print(tree_dict_list)

    # Get the list of CNVs
    cnv_df = pd.read_csv(args.cnv, skiprows=1, sep = '\t')
    # print(cnv_df)

    final_mat_gt_cnvs = []
    final_pat_gt_cnvs = []

    # Don't have to open each time
    ecDNA_dict = {}

    # Tracks current cell number
    cell_count = 0

    for cnv_set in tree_dict_list :
        mat_genome, pat_genome = generate_genomes(cnv_set, cnv_df, maternal, paternal)

        print("Generated Genomes")
        # print(mat_genome.sequences())
        # print(mat_genome.copy_number_matrix())

        no_cells = cnv_set['leaves']

        mat = mat_genome.sequences()
        pat = pat_genome.sequences()
        true_cn_matrix_mat = mat_genome.copy_number_matrix()
        true_cn_matrix_pat = pat_genome.copy_number_matrix()

        # Create a fasta for EACH cell
        for i in range(no_cells) :
            true_cn_matrix_mat['copies_final'] = true_cn_matrix_mat['copies']
            true_cn_matrix_pat['copies_final'] = true_cn_matrix_pat['copies']

            fname = f"{args.out}/dwgsim/TEMP.fa"
            with open(fname, "w") as f:

                for hap, sequence, matrix, final_list in [('mat', mat, true_cn_matrix_mat, final_mat_gt_cnvs), 
                                            ('pat', pat, true_cn_matrix_pat, final_pat_gt_cnvs)] :
                    for chrom, seq in sequence.items():
                        # Sample from ecDNA distribution
                        if 'ecDNA' in chrom :
                            path = chrom.split('____')[-1]

                            if path not in ecDNA_dict :

                                ecDNA_df = pd.read_csv(path, sep = '\t')
                                ecDNA_dict[path] = list(ecDNA_df['ecDNA_0'])
                            ecDNA_count = random.choice(ecDNA_dict[path])
                            if ecDNA_count > 0 :
                                full_seq = ""
                                for i in range(ecDNA_count) :
                                    full_seq += seq[0]


                                seq_after_ado = sim_ado(full_seq, args)
                                for seq_no, seq in enumerate(seq_after_ado) :
                                    f.write(f">{hap}_{chrom.split('____')[0]}_{seq_no}\n")
                                    f.write(f"{seq}\n")
                            # Add counts to copy number matrix
                            matrix['copies_final'] = matrix['copies_final'] + matrix[chrom] * ecDNA_count

                        else :
                            
                            seq_after_ado = sim_ado(seq[0], args)
                            for seq_no, seq in enumerate(seq_after_ado) :
                                f.write(f">{hap}_{chrom}_{seq_no}\n")
                                f.write(f"{seq}\n")

                    # Add to the gt matrix
                    cell_cn = {
                        f"{row.chrom}:{row.start}-{row.end}": row.copies_final
                        for row in matrix.itertuples()
                    }
                    cell_cn['cell'] = f"cell{cell_count}"
                    final_list.append(cell_cn)

                # Call the DWGSIM simulator 
                # TODO: format for server
                fname_wsl = windows_to_wsl(fname)
                outname_wsl = windows_to_wsl(f"{args.out}/dwgsim/cell{cell_count}_sim")
                coverage = np.random.lognormal(np.log(args.coverage_mean), args.coverage_var)
                coverage = max(coverage, 0.001)



                cmd = [
                    "wsl",
                    "/home/kmyut/miniconda3/envs/dwgsim/bin/dwgsim",
                    "-H",
                    "-C", str(coverage),
                    "-1", "100",
                    "-2", "100",
                    "-e", "0",
                    "-E", "0",
                    "-r", "0",
                    fname_wsl,
                    outname_wsl,
                ]
                print("Running:", " ".join(cmd))

                subprocess.run(cmd, check=True)

                # os.remove(fname)



            cell_count += 1
            # print(mat)
            # print(true_cn_matrix_mat)
    final_mat_df = pd.DataFrame(final_mat_gt_cnvs)
    final_pat_df = pd.DataFrame(final_pat_gt_cnvs)

    final_mat_df.to_csv(f"{args.out}/cn_mat/mat.tsv", sep = '\t', index = False)
    final_pat_df.to_csv(f"{args.out}/cn_mat/pat.tsv", sep = '\t', index = False)


def windows_to_wsl(path):
    path = Path(path).resolve()

    # C:\Users\... -> /mnt/c/Users/...
    drive = path.drive[0].lower()
    rest = path.as_posix()[2:]  # remove "C:"
    return f"/mnt/{drive}{rest}"


# Simulate ADO events
def sim_ado(seq, args) :
    seq_len = len(seq)

    # Number of ADO events
    no_events = np.random.poisson(args.ado_freq * seq_len)

    if no_events == 0:
        return [seq]

    # Sample ADO lengths in bp.
    # log10(L) ~ N(mean, variance)
    ado_lengths = 10 ** np.random.normal(
        loc=args.ado_mean,
        scale=np.sqrt(args.ado_var),
        size=no_events
    )

    # Generate random ADO intervals
    events = []

    for length in ado_lengths:
        length = int(round(length))

        # Don't allow an event longer than the sequence
        length = min(length, seq_len)

        if length == 0:
            continue

        start = np.random.randint(0, seq_len - length + 1)
        end = start + length

        events.append((start, end))

    if not events:
        return [seq]

    # Sort events by start position
    events.sort()

    # Merge overlapping ADO events
    merged = []
    start, end = events[0]

    for next_start, next_end in events[1:]:
        if next_start <= end:
            end = max(end, next_end)
        else:
            merged.append((start, end))
            start, end = next_start, next_end

    merged.append((start, end))

    # Cut out ADO regions
    fragments = []
    prev_end = 0

    for start, end in merged:
        if start > prev_end:
            fragments.append(seq[prev_end:start])

        prev_end = end

    # Sequence after final ADO event
    if prev_end < seq_len:
        fragments.append(seq[prev_end:])

    return fragments



def generate_genomes(cnv_set, cnv_df, maternal, paternal):

    results = []

    # --------------------------------------------
    # Get CNVs for this particular genome
    # --------------------------------------------

    cnv_ids = cnv_set["cnvs"]

    selected = cnv_df[
        cnv_df["id"].isin(cnv_ids)
    ]

    # --------------------------------------------
    # Fresh genomes
    # --------------------------------------------

    mat = CNVGenome(maternal)
    pat = CNVGenome(paternal)

    # Build nodes using ALL CNVs that will potentially
    # affect each haplotype.

    print("Building genomes")
    mat.build(cnv_df)
    pat.build(cnv_df)
    print("Finished building genomes")

    # --------------------------------------------
    # Apply EXACTLY in cnv_set["cnvs"] order
    # --------------------------------------------

    for cnv_id in cnv_ids:
        print(f"Applying {cnv_id}")

        cnv = selected[
            selected["id"] == cnv_id
        ].iloc[0]

        if cnv["hap"] == "M":
            mat.apply(cnv)

        elif cnv["hap"] == "P":
            pat.apply(cnv)

        elif cnv["hap"] == "D" :
            mat.apply(cnv)
            pat.apply(cnv)

        else:
            raise ValueError(
                f"Unknown haplotype: {cnv['hap']}"
            )



    return mat, pat




# Return a dict of {chr : seq} for maternal and paternal
def read_fasta(mat_file, pat_file) :
    mat_dict = {}

    for record in SeqIO.parse(mat_file, "fasta"):
        mat_dict[record.id] = str(record.seq)

    pat_dict = {}

    for record in SeqIO.parse(pat_file, "fasta"):
        pat_dict[record.id] = str(record.seq)

    return mat_dict, pat_dict


# Make list of dicts [{'leaves' : LEAF_NO, 'cnvs' : [LIST OF CNVs]}]
# CNV list order is based on the order that they happen
def parse_subnodes(s):
    result = []

    def parse(i):
        # Leaf
        if s[i] != "(":
            m = re.match(r"\d+", s[i:])
            if not m:
                raise ValueError(f"Expected leaf at position {i}")

            leaf = int(m.group())
            i += len(m.group())

            cnvs = []
            if i < len(s) and s[i] == "[":
                end = s.index("]", i)
                cnvs = [int(x) for x in s[i + 1:end].split(",")]
                i = end + 1

            result.append({
                "leaves": leaf,
                "cnvs": cnvs,
            })

            return i

        # Start of a subtree
        i += 1
        start = len(result)

        # Parse children
        while True:
            i = parse(i)

            if s[i] == ",":
                i += 1
            elif s[i] == ")":
                i += 1
                break

        # Check for annotation on this subtree
        cnvs = []
        if i < len(s) and s[i] == "[":
            end = s.index("]", i)
            cnvs = [int(x) for x in s[i + 1:end].split(",")]
            i = end + 1

            # Apply subtree CNVs to every leaf produced by this subtree.
            for j in range(start, len(result)):
                result[j]["cnvs"] = cnvs + result[j]["cnvs"]


        return i

    parse(0)
    return result



from dataclasses import dataclass
from typing import Optional


@dataclass
class Node:
    """
    A piece of sequence.

    start/end are ALWAYS reference coordinates.
    They never change when sequence is inserted/deleted.
    """

    chrom: str
    start: int
    end: int
    seq: str

    # Where this sequence came from in the reference.
    # For an original node, this is start/end.
    # For an amplified node, it still points to the original
    # reference coordinates.
    source_start: int
    source_end: int

    prev: Optional["Node"] = None
    next: Optional["Node"] = None

    # Original reference node vs amplified copy
    original: bool = True

from dataclasses import dataclass
from typing import Optional


@dataclass
class Node:
    chrom: str
    start: int
    end: int
    seq: str

    original: bool = True

    # Which original nodes are represented by this physical node?
    source_nodes: list = None

    prev: Optional["Node"] = None
    next: Optional["Node"] = None

    def __post_init__(self):
        if self.source_nodes is None:
            self.source_nodes = []



class CNVGenome:

    def __init__(self, genome):
        """
        genome:
            Your existing genome dictionary, e.g.

            {
                "chr1": ["AAATTTCCCGGG"],
                "chr2": ["TTTCCCAAA"],
                "chr3": ["GGGGTTTTCCCC"]
            }

        The input genome is NOT modified.
        """

        # Make our own copy.
        self.genome = {
            chrom: list(seq)
            for chrom, seq in genome.items()
        }

        # Linked-list heads/tails for the normal chromosomes.
        self.head = {}
        self.tail = {}

        # Maps an ORIGINAL coordinate to its original node.
        #
        # (chr1, 5) -> original node containing coordinate 5
        self.original_nodes = {}
        self.original_node_list = {}
        self.original_node_starts = {}


        self.ecDNA_sources = {}


        # Number of ecDNA events for THIS genome.
        self.ecDNA_count = 0

    # ========================================================
    # PARSE LOC
    # ========================================================
    @staticmethod
    def parse_locs(loc):
        """
        Parse one or more genomic intervals.

        Examples:
            chr1:5-10
            chr1:5-6,chr1:1-3

        Returns:
            [
                ("chr1", 5, 10),
                ("chr1", 1, 3),
            ]
        """

        fragments = []

        for part in loc.split(","):

            part = part.strip()

            chrom, coords = part.split(":")
            start, end = coords.split("-")

            fragments.append(
                (chrom, int(start), int(end))
            )

        return fragments


    # ========================================================
    # BUILD INITIAL LINKED LIST
    # ========================================================


    def build(self, cnvs):
        """
        Build the original genome nodes.

        cnvs is a pandas DataFrame.

        Every CNV boundary becomes a node boundary.

        Example:

            AMP chr1:56-90
            AMP chr1:70-80

        becomes:

            1-55
            56-69
            70-80
            81-90
            91-end

        Original nodes are stored as intervals rather than
        creating one dictionary entry per base pair.
        """

        boundaries = {}

        # --------------------------------------------
        # Start/end of every chromosome
        # --------------------------------------------

        for chrom, seq_list in self.genome.items():

            sequence = "".join(seq_list)

            boundaries[chrom] = {
                1,
                len(sequence) + 1,
            }

        # --------------------------------------------
        # Add every CNV boundary
        # --------------------------------------------

        for cnv in cnvs.itertuples(index=False):

            if cnv.type != "WGD":

                fragments = self.parse_locs(cnv.loc)

                for chrom, start, end in fragments:

                    if chrom not in boundaries:
                        continue

                    boundaries[chrom].add(start)
                    boundaries[chrom].add(end + 1)

        # --------------------------------------------
        # Create linked lists
        # --------------------------------------------

        for chrom, points in boundaries.items():

            points = sorted(points)

            sequence = "".join(self.genome[chrom])

            previous = None

            # Initialize permanent original-node structures
            self.original_node_list[chrom] = []
            self.original_node_starts[chrom] = []

            for i in range(len(points) - 1):

                start = points[i]
                end = points[i + 1] - 1

                seq = sequence[start - 1:end]

                node = Node(
                    chrom=chrom,
                    start=start,
                    end=end,
                    seq=seq,
                    original=True,
                )

                # This original node represents itself.
                node.source_nodes = [node]

                # ----------------------------------------
                # Linked list
                # ----------------------------------------

                if previous is None:
                    self.head[chrom] = node

                else:
                    previous.next = node
                    node.prev = previous

                previous = node

                # ----------------------------------------
                # Permanent list of original nodes
                # ----------------------------------------

                self.original_node_list[chrom].append(node)
                self.original_node_starts[chrom].append(start)

            # --------------------------------------------
            # Chromosome tail
            # --------------------------------------------

            self.tail[chrom] = previous


    def get_original_node(self, chrom, pos):
        """
        Return the original node containing genomic coordinate `pos`.

        Parameters
        ----------
        chrom : str
            Chromosome name.

        pos : int
            1-based genomic coordinate.

        Returns
        -------
        Node or None
            Original node containing this position.
        """

        nodes = self.original_node_list.get(chrom)

        if nodes is None:
            return None

        starts = self.original_node_starts[chrom]

        # Find the last node whose start <= pos
        i = bisect_right(starts, pos) - 1

        if i < 0:
            return None

        node = nodes[i]

        # Make sure the coordinate actually falls in this node.
        if node.start <= pos <= node.end:
            return node

        return None

    
    # ========================================================
    # GET CURRENT PATH
    # ========================================================

    def get_current_path(self, chrom, start, end):
        """
        THE CENTRAL OPERATION.

        Find the original node containing `start`.
        Find the original node containing `end`.

        Then walk the CURRENT linked list from start_node
        to end_node.

        Anything physically inserted between those anchors
        is included.

        This is the behavior you described.
        """

        start_node = self.get_original_node(
            chrom,
            start,
        )

        end_node = self.get_original_node(
            chrom,
            end,
        )

        # ----------------------------------------------------
        # Make sure the current topology still contains both
        # anchors.
        #
        # If one was deleted, there is no current path.
        # ----------------------------------------------------

        if not self.is_in_current_list(start_node):

            return []

        if not self.is_in_current_list(end_node):

            return []

        # ----------------------------------------------------
        # Walk current topology.
        # ----------------------------------------------------

        path = []

        node = start_node

        while node is not None:

            path.append(node)

            if node is end_node:
                return path

            node = node.next

        # We reached the end without finding end_node.
        #
        # This can happen if the topology was altered such
        # that the end anchor is no longer downstream.
        return []

    # ========================================================
    # CHECK WHETHER NODE STILL EXISTS
    # ========================================================

    def is_in_current_list(self, target):

        node = self.head[target.chrom]

        while node is not None:

            if node is target:
                return True

            node = node.next

        return False

    # ========================================================
    # GET SEQUENCE FOR CURRENT PATH
    # ========================================================

    def path_sequence(self, path):
        """
        Concatenate all sequence in a path.
        """

        return "".join(
            node.seq
            for node in path
        )

    # ========================================================
    # INSERT BEFORE
    # ========================================================

    def insert_before(self, node, new_node):

        previous = node.prev

        new_node.prev = previous
        new_node.next = node

        node.prev = new_node

        if previous is not None:

            previous.next = new_node

        else:

            self.head[node.chrom] = new_node

    # ========================================================
    # INSERT AFTER
    # ========================================================

    def insert_after(self, node, new_node):

        following = node.next

        new_node.prev = node
        new_node.next = following

        node.next = new_node

        if following is not None:

            following.prev = new_node

        else:

            self.tail[node.chrom] = new_node

    # ========================================================
    # REMOVE NODE
    # ========================================================

    def remove_node(self, node):

        previous = node.prev
        following = node.next

        if previous is not None:

            previous.next = following

        else:

            self.head[node.chrom] = following

        if following is not None:

            following.prev = previous

        else:

            self.tail[node.chrom] = previous

        node.prev = None
        node.next = None

    # ========================================================
    # AMP
    # ========================================================

    def amp(self, chrom, start, end):
        """
        Amplify the CURRENT path between original coordinates
        start and end.

        The copied sequence includes any material inserted by
        previous amplifications.
        """

        path = self.get_current_path(
            chrom,
            start,
            end,
        )

        if not path:
            return

        # --------------------------------------------
        # Snapshot the sequence BEFORE modifying list.
        # --------------------------------------------

        copies = []

        for node in path:

            copy = Node(
                chrom=node.chrom,
                start=node.start,
                end=node.end,
                seq=node.seq,
                original=False,
            )

            copy.source_nodes = list(node.source_nodes)
            copies.append(copy)

        # --------------------------------------------
        # Insert the copies immediately after the
        # current path.
        #
        # Since the original end_node is the final node
        # in the path, insert after it.
        # --------------------------------------------

        insertion_point = path[-1]

        for copy in copies:

            self.insert_after(
                insertion_point,
                copy,
            )

            insertion_point = copy

    # ========================================================
    # DELETE
    # ========================================================

    def delete(self, chrom, start, end):
        """
        Delete the CURRENT path between original coordinates.

        Everything physically between the original start/end
        anchors is removed.

        This includes amplified material.
        """

        path = self.get_current_path(
            chrom,
            start,
            end,
        )

        if not path:
            return

        for node in path:

            self.remove_node(node)

    # ========================================================
    # ecDNA
    # ========================================================

    def make_ecDNA(self, fragments, vec):

        pieces = []
        source_nodes = []

        for chrom, start, end in fragments:

            path = self.get_current_path(
                chrom,
                start,
                end,
            )

            if not path:
                continue

            pieces.append(
                self.path_sequence(path)
            )

            for node in path:
                source_nodes.extend(
                    node.source_nodes
                )

        sequence = "".join(pieces)

        if not sequence:
            return None

        self.ecDNA_count += 1

        name = f"ecDNA{self.ecDNA_count}____{vec}"

        self.genome[name] = [sequence]

        # Store provenance for truth matrix
        self.ecDNA_sources[name] = source_nodes

        return name



    # ========================================================
    # APPLY ONE CNV
    # ========================================================

    def apply(self, cnv):

        if cnv.type != 'WGD' :
            fragments = self.parse_locs(cnv["loc"])

        cnv_type = cnv["type"]

        # ============================================
        # AMP
        #
        # Each fragment is treated as its own AMP.
        # ============================================

        if cnv_type == "AMP":

            for chrom, start, end in fragments:

                self.amp(
                    chrom,
                    start,
                    end,
                )

        # ============================================
        # DEL
        #
        # Each fragment is treated as its own DEL.
        # ============================================

        elif cnv_type == "DEL":

            for chrom, start, end in fragments:

                self.delete(
                    chrom,
                    start,
                    end,
                )

        # ============================================
        # ecDNA
        #
        # All fragments are collected and concatenated
        # into ONE new chromosome.
        # ============================================

        elif cnv_type == "ecDNA":

            self.make_ecDNA(
                fragments,
                cnv["vec"],
            )


        
        elif cnv_type == "WGD":

            for chrom in list(self.head.keys()):

                # Find the current first and last nodes.
                first = self.head.get(chrom)
                last = self.tail.get(chrom)

                if first is None or last is None:
                    continue

                # Convert the current chromosome into an
                # AMP using the original-coordinate anchors.
                #
                # Find an original coordinate represented by
                # the first/last nodes.
                start = first.source_nodes[0].start
                end = last.source_nodes[-1].end

                self.amp(
                    chrom,
                    start,
                    end,
                )


        else:

            raise ValueError(
                f"Unknown CNV type: {cnv_type}"
            )

    # ========================================================
    # GET ALL CURRENT CHROMOSOME SEQUENCES
    # ========================================================

    def sequences(self):
        """
        Return the entire current genome.

        Example:

        {
            "chr1": ["AAAA..."],
            "chr2": ["TTTT..."],
            "ecDNA1_path": ["CCCC..."]
        }
        """

        result = {}

        # Normal chromosomes
        for chrom in self.head:

            result[chrom] = [
                self.get_chromosome_sequence(chrom)
            ]

        # ecDNA chromosomes
        for chrom, seq in self.genome.items():

            if chrom.startswith("ecDNA"):

                result[chrom] = seq

        return result

    # ========================================================
    # CURRENT CHROMOSOME SEQUENCE
    # ========================================================

    def get_chromosome_sequence(self, chrom):

        pieces = []

        node = self.head[chrom]

        while node is not None:

            pieces.append(node.seq)

            node = node.next

        return "".join(pieces)


    def copy_number_matrix(self):


        rows = []

        # ------------------------------------------------
        # Get every original section.
        # ------------------------------------------------

        for chrom, originals in self.original_node_list.items():

            for original in originals:

                # ----------------------------------------
                # Normal chromosome copy count
                # ----------------------------------------

                copies = 0

                node = self.head.get(chrom)

                while node is not None:

                    if original in node.source_nodes:
                        copies += 1

                    node = node.next

                row = {
                    "chrom": chrom,
                    "start": original.start,
                    "end": original.end,
                    "copies": copies,
                }

                # ----------------------------------------
                # ecDNA counts
                # ----------------------------------------

                for ecDNA_name, sources in self.ecDNA_sources.items():

                    row[ecDNA_name] = sources.count(
                        original
                    )

                rows.append(row)

        return pd.DataFrame(rows)


    # ========================================================
    # DEBUG
    # ========================================================

    def dump(self, chrom):

        node = self.head[chrom]

        while node is not None:

            print(
                f"{node.start}-{node.end}: "
                f"{node.seq} "
                f"{'ORIGINAL' if node.original else 'COPY'}"
            )

            node = node.next


run_sim(args)