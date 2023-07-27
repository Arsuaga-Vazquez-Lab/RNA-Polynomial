# DrTransformer-Runner

DrTransformer-Runner is a script used to segment the sequences of the plasmids *[data/pfc8_snrpn_coding_strand.fa](https://github.com/Arsuaga-Vazquez-Lab/RNA-Polynomial/blob/main/data/pfc53_airn_coding_strand.fa)* and *[pfc53_airn_coding_strand.fa](https://github.com/Arsuaga-Vazquez-Lab/RNA-Polynomial/blob/main/data/pfc53_airn_coding_strand.fa)* and generate RNA secondary structures of the sequence segments stored in *[pfc8_snrpn.zip](https://github.com/Arsuaga-Vazquez-Lab/RNA-Polynomial/blob/main/data/pfc8_snrpn.zip)* and *[pfc53_airn.zip](https://github.com/Arsuaga-Vazquez-Lab/RNA-Polynomial/blob/main/data/pfc53_airn.zip)*


## Setup

# Dependencies
* Python 3.8+
* [BioPython](https://biopython.org/)
* [DrTransformer](https://github.com/ViennaRNA/drtransformer)

Follow instructions to install [DrTransformer](https://github.com/ViennaRNA/drtransformer) and it's dependencies.
> [!IMPORTANT]  
> There is a requirement than you have [Vienna RNA](https://github.com/ViennaRNA/ViennaRNA) installed as well with the Python3 Bindings. The recommended way to install is via the binary packages provided on their website.

## Usage

```bash
python3 drtransformer-runner/run_drtransformer.py -w 200 data/pfc53_airn_coding_strand.fa
```
* `-w` or `--width` specifies the width of the sliding window used.

This will generate the folders found in the zip files found in the data directory.
