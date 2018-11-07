# Dataset

## Download

You can download pre-built dataset here:

# TODO: upload and fill in links
md5 `05ad85d871958a05c02ab51a4fde8530` [training](url)  
md5 `e53db4bff7dc4784123ae6df72e3b1f0` [validation](url)  
md5 `677b757ccec4809febd83850b43e1616` [test](url)  


## Generation

To generate the training data, run 
```
python get_data.py -o guacamol/data --chembl
```
which will download and process chembl for you in your current folder.

This script will use the molecules from `holdout_set_gcm_v1.smiles` as a
holdout set, and will exclude molecules very similar to these.

### Docker
To be sure that you have the right dependencies you can build a Docker image, run from the top-level directory:
```
docker build -t guacamol-deps dockers/
```
Then you can run:
```
docker run --rm -it  -v `pwd`:/guacamol -w /guacamol guacamol-deps python -m guacamol.data.get_data -o guacamol/data --chembl
```
