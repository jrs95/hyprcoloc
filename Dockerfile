FROM rocker/r-ver:4.3.3

# Necessary libraries
RUN Rscript -e "install.packages('devtools')"
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/Archive/RcppEigen/RcppEigen_0.3.3.9.3.tar.gz', repos=NULL, type = 'source')"
RUN Rscript -e "devtools::install_github('jrs95/hyprcoloc', build_opts = c('--resave-data', '--no-manual'), upgrade = 'never')"

LABEL version="1.0" maintainer="diegomscoelho@gmail.com"

CMD ["Rscript", "-e", "library('hyprcoloc')"]