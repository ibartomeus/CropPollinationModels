[![DOI](https://zenodo.org/badge/240485535.svg)](https://zenodo.org/badge/latestdoi/240485535)

# Crop Pollination Models

This repo creates a dynamic paper to understand the effect of pollinators on crop production. It is part of the paper "Wild insects and honey bees are equally important to crop yields in a global analysis" by James Reilly et al. [DOI TBA](), but it will be updated as more data is available.  

# The question:  

We aim to provide a permanent and updated answer to the question: **"Which is the contribution of insect pollinators to crop production?"**.
Using the CropPol dataset in combination with the best models identified in Reilly et al 2023, we provide a dynamic document that runs three types of models ... This dynamic document can be found here: [TBA]().  

# Why?  

We believe that science is iterative. Since the first case studies reporting a positive effect of pollinators on crop yield, more and more data has accumulated. This allowed us to conduct synthesis summarizing what we know [e.g. 1,2,3]. However, the answer is not easy, as we find a huge variability in crop types, pollinator communities, agricultural practices and environmental contexts. As the question is not settled, we aim to embrace this uncertainty and periodically report updates as our knowledge increases.  

# How? 

Every time CropPol database has a new release (i.e. new datasets are added) this repo automatically opens a new branch, in this brach it automatically re-run all models and figures using the new CropPol version and updates the report. Then, an email is triggered informing the repo maintainers that a new version is ready to review. Upon personal revision, the maintainers manually merge the new branch into the master branch and create a new version of CropPollinationModels. All versions are stored in Zenodo [DOI TBA]() so changes over time can be tracked. 

As all models are open, models can be updated if necessary and new models can be added into the report. If you want to contribute to the modelling efforts, let us know in an [issues]() or directly make a pull request. 

# Details  

- `Reilly_models` folder contains the script which reproduce the original analisys and its useful to understand some decisions on how final modls where selected.
- `dynamic_report` contains the Quarto webpage ... 

# Refrences 

- Reilly et al. 2023 Submitted.
- Allen-Perkins A, et al. (2022) CropPol: a dynamic, open and global database on crop pollination. Ecology 103 (3): e3614.  
- 1: Garibaldi L et al. 2013 Wild pollinators enhance fruit set of crops regardless of honey-bee abundance. Science 339, 1608–1611. (doi:10.1126/science.1230200). 
- 2: Rader R et al. 2016 Non-bee insects are important contributors to global crop pollination. PNAS 113, 146–151. (doi:10.1073/pnas.1517092112). 
- 3: Dainese M et al. 2019 A global synthesis reveals biodiversity-mediated benefits for crop production. Science Advances 5: eaax0221.  

# To Do:  
[ ] Add Links to download a PDF and social media to share.   
[x] Add htmls to /docs automatically to make it visible in github pages.   
[ ] Add a news.md for versioning.   
[ ] make github actions to work.   
[ ] Re-think data cleaning.   
[ ] Re-think final models to show.  
[ ] Use Renv to freeze package versions
[ ] On commit -> [skip ci]

