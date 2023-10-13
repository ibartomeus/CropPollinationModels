[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7818147.svg)](https://doi.org/10.5281/zenodo.7818147)


# Crop Pollination Models

This repo creates a dynamic paper to understand the effect of pollinators on crop production. It is part of the paper "Wild insects and honey bees are equally important to crop yields in a global analysis" by James Reilly et al. (2023), but it will be further updated as more data is available. This repo can be cited with DOI: [10.5281/zenodo.7481551](https://doi.org/10.5281/zenodo.7481551).

# The question:  

We aim to provide a permanent and updated answer to the question: **"What is the contribution of insect pollinators to crop production?"**.
Using the [CropPol](https://github.com/ibartomeus/OBservData) database in combination with the best model identified in Reilly et al 2023, we provide a dynamic document that runs the full model and sumarize the effect sizes of all parameters. This dynamic document can be found here: [https://ibartomeus.github.io/CropPollinationModels](https://ibartomeus.github.io/CropPollinationModels).  

# Why?  

We believe that science is iterative. Since the first case studies reporting a positive effect of pollinators on crop yield, more and more data has accumulated. This allowed us to conduct synthesis summarizing what we know (e.g. Garibaldi et al 2013, Rader et al 2016, Dainese et al 2019). However, the answer is not easy, as we find a huge variability in crop types, pollinator communities, agricultural practices and environmental contexts. As the question is not settled, we aim to embrace this uncertainty and periodically report updates as our knowledge increases.   

# How? 

Every time CropPol database has a new release (i.e. new datasets are added) this repo automatically re-run all models and figures using the new CropPol version and updates the report. Then, an email is triggered informing the repo maintainers that a new version is ready to review. Upon personal revision, the maintainers manually create a new release of CropPollinationModels. All versions are stored in Zenodo [10.5281/zenodo.7481551](https://doi.org/10.5281/zenodo.7481551) so changes over time can be tracked.  

As all models are open, models can be updated if necessary and new models can be added into the report. If you want to contribute to the modelling efforts, let us know in an [issues](https://github.com/ibartomeus/CropPollinationModels/issues) or directly make a pull request. You can check the `News.md` file to see main version changes.

# But how?  

This repo is hosted in github and uses github actions (`/.github/workflows/publish.yml`) to trigger the document update. To build the webpage we used Quarto (`.qmd` files + `_quarto.yml`). To ensure consistency in package versions we use Renv packge (using `init()` to set up the R dependency management and `snapshot()` to update it).

To set up the github actions we first created a gh-pages branch for deployment (`quarto publish gh-pages` in terminal; details [here](https://quarto.org/docs/publishing/github-pages.html)) and connect it to github (on github `Settings/Pages/Deploy from branch/gh-pages/(root)`). You also need to add a `.nojekyll` empty file in the root and gitignore `/.quarto/` and `/_site/`. We scheduled the github action to ran once a month using `cron` but it can be triggered manually from the github webpage. Note that you can also re-build the webpage manually (`quarto publish gh-pages` in terminal). Also, R scripts ran on a server, as specified in the quarto and github actions ylm´s thanks to `Renv` package and github computer time generosity. Everything is licenced with a MIT licence (`LICENCE`) and all this couldn't be possible without the help of Paco Rodriguez.    

The `scripts` folder contains the script which reproduce the original analysis for the recods and its useful to understand some decisions on how final models were selected. It also contains scripts to read the data and store past results.

# References 

- Reilly et al. 2023 Submitted.
- Allen-Perkins A, et al. (2022) CropPol: a dynamic, open and global database on crop pollination. Ecology 103 (3): e3614.  
- Garibaldi L et al. 2013 Wild pollinators enhance fruit set of crops regardless of honey-bee abundance. Science 339, 1608–1611. (doi:10.1126/science.1230200). 
- Rader R et al. 2016 Non-bee insects are important contributors to global crop pollination. PNAS 113, 146–151. (doi:10.1073/pnas.1517092112). 
- Dainese M et al. 2019 A global synthesis reveals biodiversity-mediated benefits for crop production. Science Advances 5: eaax0221.  

