# Rmarkdown_workshop

Minicurso "Rmarkdown para Comunicação Estatística" apresentado na Jornada dos 
60 anos do IME UFRGS.

## Para visualizar os slides

https://markus-stein.github.io/Rmarkdown_workshop/index.html


## Para reproduzir os slides

Os slides foram escritos usando o pacote `xaringan` do R. Para reproduzi-los 
verifique:

1. se você o tem instalado em sua máquina, `library(xaringan)`, senão, para 
instalar a versão disponível no CRAN

```
install.packages("xaringan")
library("xaringan")
```

Ou, para a versão em desenvolvimento no Github (instalar o pacote `devtools` 
também caso não esteja instalado)

```
# install.packages("devtools")
devtools::install_github("r-lib/devtools")
library("xaringan")
```

2. Outros pacotes também são necessários `gapminder`, `knitr` e `dplyr`.

3. Clonar o repositório https://github.com/markus-stein/Rmarkdown_workshop e abrir o projeto no Rstudio (arquivo `.Rproj`). 

4. Dentro do projeto abra e `knit` o arquivo da apresentação, `index.Rmd`.

![][id]

[id]: images/cropped-logo-60-transparente.png


<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative Commons Attribution-NonCommercial 4.0 International License</a>.