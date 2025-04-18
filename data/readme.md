# Data files

## Pome_competition_experiment_data.csv
<ul>
  <li> This file contains data for the competition experiments</li>
  <li> It contains the counts of the evolved strain and the ancestor</li>
  <li> The columns are:</li>
  <ul>
  <li> 1: Competition.assay: competition assay ID </li>
  <li> 2: Evolution.environment: The environment where the evolved strain evolved </li>
  <li> 3: Population.size: population size of the evovled strain, either 24 (24-well plate) or 96 (96-well plate)</li>
  <li> 4: Evolved.strain.ID:  ID of the evolved strain </li>
  <li> 5: Evolved.strain.mating.type:  mating type of the evolved strain </li>
  <li> 6: Evolved.strain.ade6.marker: ade6 allele of the evolved strain </li>
  <li> 7: Evolved.strain.genotype: combined genotype of mating and ade6 marker of the evolved strain </li>
  <li> 8: Ancestor.genotype: combined genotype of the competitor</li>
  <li> 9: Ancestor.mating.type: mating type of the competitor </li>
  <li> 10: Ancestor.ade6.marker: ade6 allele of the competitor </li>
  <li> 11: Competition.replicate:  compettion replicate </li>
  <li> 12: Round: competition round </li>
  <li> 13: Picture.ID.initial: ID of the photograph taken of the plate from the starting time point </li>
  <li> 14: Evolved.strain.initial.CFU: Colony count of the evolved strain at the starting time point </li>
  <li> 15: Ancestor.initial.CFU: Colony count of the competitor at the starting time point </li>
  <li> 16: Competition.environment: The competition environment </li>
  <li> 17: Picture.ID.final: ID of the photograph taken of the plate from the final time point </li>
  <li> 18: Evolved.strain.final.CFU: Colony count of the evolved strain at the final time point </li>
  <li> 19: Ancestor.final.CFU: Colony count of the competitor from the final time point </li>
  </ul>
</ul>

## Pombe_competition_experiment_models_original.RData

This in an R data file contains the saved brms models for the analysis of competition experiment data. This intermediate step is saved so that the posterior distributions do not need to ba calculated again eveytime, since sampling for those models is a bit slow. This file can be loaded into R.

## ODdata_filtered.csv

<ul>
  <li> This file contain data for clone growth curve measurements, maximal growth rates and carrying capacities have been extracted from the fitted growth curves </li>
  <li> The columns are: </li>
  <ul>
    <li> 1: Temperature: Temperature treatment where the growth curve was measured </li>
    <li> 2: Clone.plate: ID of the cryoplate clones are stored in </li>
    <li> 3: Replicate: Replicate ID </li>
    <li> 4: Pregrowth.OD: Optical density measured from the pregrowth plate, uased as a covariate in the model </li>
    <li> 5: Date: Date of the measurement </li>
    <li> 6: Well number on the bioscreen plate </li>
    <li> 7: Population.ID: ID of the population from where the clone was sampled from </li>
    <li> 8: CLone:ID : ID of the clone </li>
    <li> 9: Mating.type: mating type of the clone </li>
    <li> 10: ade6.allele: ade6 allele of the clone </li>
    <li> 11: Evolution.treatment: The environment that the population experienced during the evolution treatment </li>
    <li> 12: Population.size: The population size treatment that the population had during the evolution treatment (either 24 96 well-plate) </li>
    <li> 13: Experiment.ID: experiment ID </li>
    <li> 14: k: value for the carrying capacity parameter extracted from the growth curve fit </li>
    <li> 15: r: value for the maximum growth rate paremeter extracted from the growth curve fit </li>
    <li> 16: maxOD: maximum optical density value in the growth curve measurement </li>
    <li> 17: Observation: notes from the growth curve fits </li>
    <li> 18: kasvu: notes from the growth curve fits </li>
  </ul>
</ul>
