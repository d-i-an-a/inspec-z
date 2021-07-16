# Inspec-z
Interactive inspection and validation of 1D spectra.
- [Docker support](#Docker)
- [Conda Installation](#Conda-Installation)
- [Run sequence](#Run-sequence)


# Docker 
Docker instructions for Mac OS X can be found [here](Docker/README.md).

# Conda installation
Clone package locally or download the zip file:

```bash
git clone https://github.com/d-i-an-a/inspec-z.git
cd inspec-z
```

Inspec-z requires: `python`, `pyyaml`, `numpy`, `astropy`, `scipy`, and `matplotlib`.

A conda environment is recommended. Either install the latest version of the required packages:

```bash
conda create --name inspec-z
conda activate inspec-z
conda install python pyyaml numpy astropy matplotlib scipy
```
or install the lastest tested version:

`conda create --name  inspec-z --file requirements.txt`

# Run sequence
An example configuration and data files are provided.

1. Edit or create the `yaml` configuration file in `config_files` folder.
- Set up the file paths and file names. To follow the example replace `/path/to/` with the path to the folder containing the `inspec-z` code.
- Change the selection of the sample to be inspected.

2. Prepare catalog: `python prepare_catalog.py config_files/example.yaml`
3. Run main program: `python inspecz.py config_files/example.yaml`
4. Merge results with initial data file and get report: `python finalize.py config_files/example.yaml`

<!-- ![inspecz_main_screenshot](https://user-images.githubusercontent.com/48242007/116094996-e95b0780-a69f-11eb-8010-a1cc7c2b9633.png) -->

![Isnpecz_v1](https://user-images.githubusercontent.com/48242007/125954299-a3b4f8ad-2991-4a89-8843-363f468642ff.png)
