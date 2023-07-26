import subprocess
import dataclasses
import re
import time
import argparse
import pathlib

from Bio import SeqIO

MATCH_STRUCTURE_INFO_PATTERN = re.compile(
    r"\s+(\d+)\s+\d+\s+([\.\(\)]+)\s+(-?\d+\.\d+).*ID\s=\s(\d+)"
)


@dataclasses.dataclass
class CotranscriptionStep:
    id_: str
    structure: str
    energy: str
    step: str


parser = argparse.ArgumentParser(prog="run_drtransformer")
parser.add_argument("-w", "--window_size", type=int)
parser.add_argument("fasta_file", type=pathlib.Path)


def main():
    args = parser.parse_args()

    with open(args.fasta_file) as handle:
        sequence_fasta = next(SeqIO.parse(handle, "fasta"))

    sequence = str(sequence_fasta.seq)

    for i in range(0, len(sequence) - args.window_size + 1):
        start, stop = (i, i + args.window_size)
        drtransformer_process = subprocess.Popen(
            f"DrTransformer",
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            shell=True,
        )

        start_time = time.perf_counter()
        sequence_window = sequence[start:stop]
        assert len(sequence_window) == args.window_size

        outs, _ = drtransformer_process.communicate(sequence_window.encode("utf-8"))

        stop_time = time.perf_counter()

        elapsed_time = stop_time - start_time
        print("Elapsed time", elapsed_time)

        result = outs.decode("utf-8")
        result_lines = list(
            filter(
                lambda x: not x.startswith("#") and not len(x.strip()) == 0,
                result.split("\n"),
            )
        )  # remove headers
        id_structure_map: dict[int, CotranscriptionStep] = {}

        for result_line in result_lines:
            structure_info_match = MATCH_STRUCTURE_INFO_PATTERN.match(result_line)
            if not structure_info_match:
                print(result_line)
            transcription_step = structure_info_match[1]
            structure = structure_info_match[2]
            energy = structure_info_match[3]
            id_ = structure_info_match[4]

            if id_ in id_structure_map:
                if len(id_structure_map[id_].structure) < len(structure):
                    id_structure_map[int(id_)] = CotranscriptionStep(
                        id_, structure, energy, transcription_step
                    )
            else:
                id_structure_map[int(id_)] = CotranscriptionStep(
                    id_, structure, energy, transcription_step
                )

        with open(f"{args.fasta_file.parts[-1][:-3]}_{start}-{stop}.drf", "w") as fh:
            for id_ in sorted(id_structure_map.keys()):
                cotrans_struct = id_structure_map[id_]
                cotrans_struct_tuple = dataclasses.astuple(cotrans_struct)
                reconstructed_line = " ".join(cotrans_struct_tuple)
                fh.write(f"{reconstructed_line}\n")


if __name__ == "__main__":
    main()
