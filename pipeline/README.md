# Event driven analytics pipeline for ABMI core

## Goal

To automate the sh*t out of everything that we do.

## Execution

### Hosting

The steps presented here would require their own project in
the GitHub/GitLab organization. GitHub is great, but GitLab has better
tooling for continuois integration (CI/CD) pipelines
(if needed we can set up own own Kubernetes cluster to run jobs) and 
no restrictions around private projects.

Private projects might be important if we have sensitive info
(species data we can't show, API endpoints we can't expose etc)

Public projects can be mirrored on GitHub, or made public on GitLab.

### Versioning

Major.Minor.Bugfix

* Major versions: each year would represent a major version, meaning that things are expected **not** to be backwards compatible.
* Minor versions: this represents backwards compatible changes, which can be a case when adding a feature to a previous version.
* Bugfix releases: these address specific issues.

Version life cycle:

1. Bump up major version in development branch.
2. Once working, merge to master branch.
3. Add new features and bump minor version (feature branches 1st, merge to development branch)
4. Fix bugs (track through issues: hotfix branch, merge to master)

The event driven in the title indicates that whenever code base
in master and development branches changes, it triggers the
CI/CD pipeline to run package and unit tests and 
update all the project artefacts (data bundles,
R packages, Docker containers).

### Steps

#### 1. Data ingest

Set up DB connections (safely store credentials/secrets/tokens in private repos)
for geospatial (Eric) and species data sets (Oracle, BU/WildTrax, BOREAL).

Use scripts to

1. pull data from DBs
2. run checks to make sure data is in the expected format
3. fix issues: WHO IS FIXING IT? ARE WE MAINTAINING ALL MASTER LOOKUP TABLES?
4. CI/CD pipeline to create data bundles for modeling

We need the following sets of scripts:

* `species` to turn raw data into sample x species tables + metadata
* `geo` to turn geospatial summaries into tables ready for modeling

he output is a bundle that I call `sdtin` (standard input).

#### 2. Modeling

Take scripts we already have and standardize across taxa.

Keep it modular: define model output to be atomic (1 bootstrap run for 1 species in one modeling region). These model outputs than can be summarized by other functions. Avoid loops.

Out put of the modeling is `stdout` (standard output) that is described in the file `StandardizedSpeciesOutputMetadata.Rmd`.

#### 3. Derivatives

This is the realm of the **cure4insect** package which takes the `stdout` objects and turn that into predictions.

Here we need to deal with GEOSPATIAL STUFF FOR PREDICTION that comes from `geo` scripts. Workflow can be:

1. Standardize geospatial raw data and store in same DB.
2. have specialized scripts for prediction by project and use appropriate inputs.



