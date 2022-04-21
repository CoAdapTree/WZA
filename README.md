
         ________ __                        __       __  ________   ______  
        /        /  |                      /  |  _  /  |/        | /      \
        $$$$$$$$/$$ |____    ______        $$ | / \ $$ |$$$$$$$$/ /$$$$$$  |
          $$ |  $$      \  /      \       $$ |/$  \$$ |    /$$/  $$ |__$$ |
          $$ |  $$$$$$$  |/$$$$$$  |      $$ /$$$  $$ |   /$$/   $$    $$ |
          $$ |  $$ |  $$ |$$    $$ |      $$ $$/$$ $$ |  /$$/    $$$$$$$$ |
          $$ |  $$ |  $$ |$$$$$$$$/       $$$$/  $$$$ | /$$/____ $$ |  $$ |
          $$ |  $$ |  $$ |$$       |      $$$/    $$$ |/$$      |$$ |  $$ |
          $$/   $$/   $$/  $$$$$$$/       $$/      $$/ $$$$$$$$/ $$/   $$/




# The WZA

A repository containing the analysis scripts and simulation configuration files to accompany our project on the Weighted-Z Analysis for GEA (WZA).

Briefly, the WZA is a method to aggregate data across closely linked sites in GEA studies. In other words, the WZA is a method for performing GEA using analysis windows rather than individual genetic markers (e.g. SNPs).

There are two main parts to this repository, the first is the simulation study [SimulationStudy/](SimulationStudy/)and the second is the analysis of real data (```dataAnalysis/```).

The code in the [dataAnalysis/](dataAnalysis/) directory contains the functions to aggregate data across sites in GEA studies. For the pre-print describing the WZA, we re-analyzed data from Yeaman et al (2016). For ongoing work in the CoAdapTree project, we used data in a different format, so I wrote code to work specifically with that kind of data.

Plots and code to generate them for the data analysis and simulation study are available within the appropriate directoy. Additionally, the write up and code to make a poster I have presented for this project are available here (for anyone who needs help to sleep).

## Notes

Going back and forth from servers to my local machine, there was a lot of tweaking paths. If you want to implement any of these scripts, you'll need to take a good look at paths to make sure everything will work on your machine.

During the course of this project I referred to the Gradient map as the "cline" map. Lots of files and scripts use this nomenclature. However, during the writing process we decided that "cline" is perhaps too loaded a term to use, so we went with "gradient" instead.

I typically use English (UK) spelling, but sometimes just let spell-checkers do their thing so there is a bit of a combination of UK and US spellings of things.
...